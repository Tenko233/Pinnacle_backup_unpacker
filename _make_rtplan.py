# -*- coding: utf-8 -*-
"""
Generate a DICOM RT Plan (RTPLAN) object from Pinnacle treatment-plan data.

The beam geometry (gantry/collimator/couch angles, jaw and MLC leaf positions,
cumulative meterset weights) is read from the ``plan.Trial`` text header and
mapped to the DICOM RT Beams module.

Note: jaw/MLC leaf sign conventions follow the IEC 61217 beam limiting-device
coordinate system; magnitudes are taken directly from Pinnacle, but the +/- X/Y
sign convention should be cross-checked against a reference TPS export before
clinical use (research use only).
"""

import os
import re
import sys
import time

import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from _make_rtdose import (
    fmt_ds,
    fmt_ds_list,
    parse_patient_position,
    parse_patient_registry,
    read_image_uids,
)

RTPLAN_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.5"
IMPLEMENTATION_CLASS_UID = "1.2.826.0.1.3680043.8.498.75006884747854523615841001"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _find_block(text, open_idx):
    depth = 0
    i = open_idx
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def _top_level_blocks(text, entry_re):
    res = []
    i = 0
    n = len(text)
    while i < n:
        m = entry_re.search(text, i)
        if not m:
            break
        o = text.find("{", m.start())
        if o < 0:
            break
        e = _find_block(text, o)
        if e < 0:
            break
        res.append((text[o + 1:e], m))
        i = e + 1
    return res


def _get_float(text, key):
    m = re.search(re.escape(key) + r"\s*=\s*([-0-9.eE+]+)", text)
    return float(m.group(1)) if m else None


def _get_int(text, key):
    m = re.search(re.escape(key) + r"\s*=\s*(\d+)", text)
    return int(m.group(1)) if m else None


def _get_str(text, key):
    m = re.search(re.escape(key) + r'\s*=\s*"([^"]*)"', text)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def _parse_beams(trial_text):
    """Return (beams, mlc_bounds_mm) parsed from the plan.Trial header."""
    beams = []
    bl = re.search(r"BeamList\s*=\s*\{", trial_text)
    beamlist = ""
    if bl:
        blo = trial_text.find("{", bl.start())
        ble = _find_block(trial_text, blo)
        beamlist = trial_text[blo + 1:ble]
    for m in re.finditer(r"Beam\s*=\s*\{", beamlist):
        o = beamlist.find("{", m.start())
        e = _find_block(beamlist, o)
        b = beamlist[o + 1:e]
        if _get_str(b, "Name") and _get_field_num(b) is not None:
            beams.append(b)

    # MLC leaf boundaries (Varian Millennium 120: 60 pairs, 1 cm outer / 0.5 cm inner)
    mlc_bounds_mm = None
    row = re.search(r"RowLabelList\s*=\s*\{", trial_text)
    if row:
        ro = trial_text.find("{", row.start())
        re_ = _find_block(trial_text, ro)
        ys = [float(x) for x in re.findall(r"Y =\s*([-0-9.]+)\s*cm", trial_text[ro + 1:re_])]
        if ys:
            bounds = [ys[0] - 0.5]
            for i in range(len(ys)):
                w = 1.0 if (i < 10 or i >= 50) else 0.5
                bounds.append(ys[i] + w / 2)
            mlc_bounds_mm = [round(b * 10, 4) for b in bounds]
    return beams, mlc_bounds_mm


def _get_field_num(b):
    return _get_int(b, "NumberOfControlPoints")


def _parse_prescription(trial_text):
    """Return (prescription_name, prescription_dose_cGy, n_fractions, total_mu)."""
    presc_name = None
    presc_dose = None
    n_frac = None
    for m in re.finditer(r"Prescription\s*=\s*\{", trial_text):
        o = trial_text.find("{", m.start())
        e = _find_block(trial_text, o)
        blk = trial_text[o + 1:e]
        if presc_name is None:
            presc_name = _get_str(blk, "Name")
            presc_dose = _get_float(blk, "PrescriptionDose")
            n_frac = _get_int(blk, "NumberOfFractions")
    total_mu = _get_float(trial_text, "RequestedMonitorUnitsPerFraction")
    return presc_name or "", presc_dose or 212.0, n_frac or 1, total_mu


def _parse_iso(trial_text, plan_dir):
    """Return the DICOM isocenter (mm) from the first beam's prescription point."""
    # The isocenter POI name referenced by the first beam
    iso_name = None
    bl = re.search(r"BeamList\s*=\s*\{", trial_text)
    if bl:
        blo = trial_text.find("{", bl.start())
        ble = _find_block(trial_text, blo)
        first_beam = re.search(r"Beam\s*=\s*\{", trial_text[blo + 1:ble])
        if first_beam:
            o = trial_text.find("{", blo + 1 + first_beam.start())
            e = _find_block(trial_text, o)
            iso_name = _get_str(trial_text[o + 1:e], "IsocenterName")
    # look up the POI
    pts_path = os.path.join(plan_dir, "plan.Points")
    if os.path.exists(pts_path):
        pts = open(pts_path, "rb").read().decode("latin-1")
        for m in re.finditer(r"Poi\s*=\s*\{", pts):
            o = pts.find("{", m.start())
            e = _find_block(pts, o)
            blk = pts[o + 1:e]
            name = _get_str(blk, "Name")
            if name == iso_name or (iso_name is None and name == "ISO"):
                x = _get_float(blk, "XCoord")
                y = _get_float(blk, "YCoord")
                z = _get_float(blk, "ZCoord")
                if x is not None:
                    # convert_point HFS: X, -Y, -Z (mm)
                    return [x * 10, -y * 10, -z * 10]
    return [0.0, 0.0, 0.0]


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------

def generate_rtplan(patient_dir, plan_dir_name, output_dir, dose_ref_uid=None):
    """Generate a DICOM RTPLAN for one plan directory.

    Returns ``(out_path, plan_sop_uid)``.
    """
    os.makedirs(output_dir, exist_ok=True)
    plan_dir = os.path.join(patient_dir, plan_dir_name)
    info, full_name, plans = parse_patient_registry(patient_dir)
    plan = next((p for p in plans if f"Plan_{p['plan_id']}" == plan_dir_name), None)
    plan_name = plan["plan_name"] if plan else plan_dir_name
    image_set_id = plan["primary_ct"] if plan else 0

    study_uid, frame_uid = read_image_uids(patient_dir, str(image_set_id))
    patient_position = parse_patient_position(plan_dir)
    trial_text = open(os.path.join(plan_dir, "plan.Trial"), "rb").read().decode("latin-1")

    beams, mlc_bounds_mm = _parse_beams(trial_text)
    presc_name, presc_dose, n_frac, total_mu = _parse_prescription(trial_text)
    iso_dicom = _parse_iso(trial_text, plan_dir)

    dob = (info["dob"] or "").replace("-", "").replace("/", "").replace(" ", "")
    sex = (info["gender"] or "O")[:1].upper()
    # DICOM PatientID should be the hospital MRN (fall back to Pinnacle internal id)
    patient_id = str(info["mrn"] or info["patient_id"] or "")

    # ---- file meta --------------------------------------------------------
    sop_uid = generate_uid()
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = RTPLAN_SOP_CLASS_UID
    file_meta.MediaStorageSOPInstanceUID = sop_uid
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    file_meta.ImplementationClassUID = IMPLEMENTATION_CLASS_UID

    ds = FileDataset(None, {}, file_meta=file_meta, preamble=b"\x00" * 128)
    ds.is_implicit_VR = False
    ds.is_little_endian = True

    # ---- SOP Common / General --------------------------------------------
    ds.SpecificCharacterSet = "ISO_IR 100"
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")
    ds.SOPClassUID = RTPLAN_SOP_CLASS_UID
    ds.SOPInstanceUID = sop_uid
    ds.Modality = "RTPLAN"
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = generate_uid()
    ds.FrameOfReferenceUID = frame_uid
    ds.SeriesNumber = "1"

    # ---- Patient ----------------------------------------------------------
    ds.PatientName = full_name or "UNKNOWN"
    ds.PatientID = patient_id or "UNKNOWN"
    if dob:
        ds.PatientBirthDate = dob
    if sex:
        ds.PatientSex = sex

    # ---- RT Plan ----------------------------------------------------------
    ds.RTPlanLabel = plan_name
    ds.RTPlanName = plan_name
    ds.RTPlanDate = time.strftime("%Y%m%d")
    ds.RTPlanTime = time.strftime("%H%M%S")
    ds.PlanIntent = "CURATIVE"
    ds.RTPlanGeometry = "PATIENT"

    # ---- Prescription -----------------------------------------------------
    target_prescription_gy = presc_dose * n_frac / 100.0
    dose_ref = Dataset()
    dose_ref.DoseReferenceNumber = "1"
    dose_ref.DoseReferenceDescription = presc_name
    dose_ref.DoseReferenceType = "TARGET"
    dose_ref.DoseReferenceStructureType = "SITE"
    dose_ref.TargetPrescriptionDose = fmt_ds(target_prescription_gy)
    dose_ref.DoseReferenceUID = dose_ref_uid or generate_uid()
    ds.DoseReferenceSequence = Sequence([dose_ref])

    # ---- Fraction group ---------------------------------------------------
    fg = Dataset()
    fg.FractionGroupNumber = "1"
    fg.NumberOfFractionsPlanned = str(n_frac)
    fg.NumberOfBeams = str(len(beams))
    fg.NumberOfBrachyApplicationSetups = "0"
    refs = []
    per_beam_mu = (total_mu / len(beams)) if beams else 0.0
    for bi in range(len(beams)):
        bds = Dataset()
        bds.ReferencedBeamNumber = str(bi + 1)
        bds.BeamDose = fmt_ds(presc_dose / 100.0)
        bds.BeamMeterset = fmt_ds(per_beam_mu)
        refs.append(bds)
    fg.ReferencedBeamSequence = Sequence(refs)
    ds.FractionGroupSequence = Sequence([fg])

    # ---- beams ------------------------------------------------------------
    beam_seq = []
    for bi, b in enumerate(beams):
        name = _get_str(b, "Name")
        machine = _get_str(b, "MachineNameAndVersion") or "UNKNOWN"
        btype = _get_str(b, "SetBeamType") or "STATIC"
        ncp = _get_int(b, "NumberOfControlPoints") or 0

        cp_m = re.search(r"ControlPointList\s*=\s*\{", b)
        cp_entries = []
        if cp_m:
            cpo = b.find("{", cp_m.start())
            cpe = _find_block(b, cpo)
            cp_entries = _top_level_blocks(b[cpo + 1:cpe], re.compile(r"#(\d+)\s*=\s*\{"))
            cp_entries.sort(key=lambda t: int(t[1].group(1)))

        beam = Dataset()
        beam.BeamNumber = str(bi + 1)
        beam.BeamName = name
        beam.BeamType = "DYNAMIC" if "Arc" in btype else "STATIC"
        beam.RadiationType = "PHOTON"
        beam.TreatmentMachineName = machine.split(":")[0].strip() or "UNKNOWN"
        beam.PrimaryDosimeterUnit = "MU"
        beam.SourceAxisDistance = "1000"

        dev_seq = Sequence()
        for dt, npairs in [("X", "1"), ("Y", "1"), ("MLCX", "60")]:
            dev = Dataset()
            dev.RTBeamLimitingDeviceType = dt
            dev.SourceToBeamLimitingDeviceDistance = "1000"
            dev.NumberOfLeafJawPairs = npairs
            dev_seq.append(dev)
        beam.BeamLimitingDeviceSequence = dev_seq

        cp_seq = Sequence()
        cum = 0.0
        for inner, m in cp_entries:
            cp = Dataset()
            cp.ControlPointIndex = str(m.group(1))
            cp.CumulativeMetersetWeight = fmt_ds(cum)
            cp.NominalBeamEnergy = "6"
            cp.DoseRateSet = "600"
            cp.IsocenterPosition = fmt_ds_list(iso_dicom)
            g = _get_float(inner, "Gantry")
            couch = _get_float(inner, "Couch")
            coll = _get_float(inner, "Collimator")
            cp.GantryAngle = fmt_ds(g if g is not None else 0.0)
            cp.BeamLimitingDeviceAngle = fmt_ds(coll if coll is not None else 0.0)
            cp.PatientSupportAngle = fmt_ds(couch if couch is not None else 0.0)
            cp.GantryRotationDirection = "NONE"
            cp.BeamLimitingDeviceRotationDirection = "NONE"
            cp.PatientSupportRotationDirection = "NONE"

            pos_seq = Sequence()
            L = _get_float(inner, "LeftJawPosition")
            R = _get_float(inner, "RightJawPosition")
            T = _get_float(inner, "TopJawPosition")
            Btm = _get_float(inner, "BottomJawPosition")
            for dt, vals in [("X", [-L * 10, R * 10]), ("Y", [-Btm * 10, T * 10])]:
                d = Dataset()
                d.RTBeamLimitingDeviceType = dt
                d.LeafJawPositions = fmt_ds_list([v if v is not None else 0.0 for v in vals])
                pos_seq.append(d)
            mm = re.search(r"Points\[\]\s*=\s*\{", inner)
            if mm and mlc_bounds_mm:
                mo = inner.find("{", mm.start())
                me = _find_block(inner, mo)
                pts = [float(x) for x in re.findall(r"[-0-9.eE+]+", inner[mo + 1:me])]
                A = pts[0::2]
                B = pts[1::2]
                leaf_jaw = [a * 10 for a in A] + [b * 10 for b in B]
                d = Dataset()
                d.RTBeamLimitingDeviceType = "MLCX"
                d.LeafJawPositions = fmt_ds_list(leaf_jaw)
                d.LeafPositionBoundaries = [str(v) for v in mlc_bounds_mm]
                pos_seq.append(d)
            cp.BeamLimitingDevicePositionSequence = pos_seq

            w = _get_float(inner, "Weight")
            cum += (w if w else 0.0)
            cp_seq.append(cp)
        beam.ControlPointSequence = cp_seq
        beam_seq.append(beam)
    ds.BeamSequence = Sequence(beam_seq)

    out_path = os.path.join(output_dir, f"RTPLAN_{plan_dir_name}.dcm")
    pydicom.dcmwrite(out_path, ds, write_like_original=False)
    return out_path, sop_uid


# ---------------------------------------------------------------------------
# main (standalone)
# ---------------------------------------------------------------------------

def main():
    import _make_rtdose  # noqa: F401  (only for the module constants when run standalone)
    # standalone defaults kept for backwards compatibility
    from _make_rtdose import BASE as _BASE, PLAN_NAME as _PLAN, OUTPUT_DIR as _OUT
    out_path, uid = generate_rtplan(_BASE, _PLAN, _OUT)
    print("WROTE:", out_path)
    print("Plan UID:", uid)


if __name__ == "__main__":
    sys.exit(main())
