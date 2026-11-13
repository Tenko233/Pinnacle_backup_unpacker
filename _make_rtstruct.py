# -*- coding: utf-8 -*-
"""
Generate a DICOM RT Structure Set (RTSTRUCT) from Pinnacle ROI/POI data.

Sources:
  * ``plan.Points``   -- POI markers (e.g. Ref, ISO)  -> POINT contours
  * ``plan.roi``      -- ROI contours (Body, organs, PTV/GTV/CTV, ...)
  * ``ImageSet_N.DICOM`` -- the CT images that contours are referenced against
"""

import os
import re
import time
import glob

import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from _make_rtdose import (
    convert_point,
    parse_patient_position,
    parse_patient_registry,
)

RTSTRUCT_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.3"
CT_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.2"
IMPLEMENTATION_CLASS_UID = "1.2.826.0.1.3680043.8.498.75006884747854523615841001"

# Pinnacle colour name -> [R, G, B] (extracted from a Pinnacle TPS export)
COLOR_MAP = {
    "aquamarine": [127, 255, 212],
    "blue": [0, 0, 255],
    "brown": [165, 80, 65],
    "forest": [34, 139, 34],
    "green": [0, 255, 0],
    "lavender": [200, 180, 255],
    "lightblue": [0, 128, 255],
    "olive": [165, 161, 55],
    "orange": [255, 150, 0],
    "purple": [255, 0, 255],
    "red": [255, 0, 0],
    "skyblue": [0, 255, 255],
    "slateblue": [128, 0, 255],
    "steelblue": [70, 130, 180],
    "teal": [0, 160, 160],
    "tomato": [255, 83, 76],
    "yellow": [255, 255, 0],
    "yellowgreen": [192, 255, 0],
}


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


def _get_float(text, key):
    m = re.search(re.escape(key) + r"\s*=\s*([-0-9.eE+]+)", text)
    return float(m.group(1)) if m else None


def _get_int(text, key):
    m = re.search(re.escape(key) + r"\s*=\s*(\d+)", text)
    return int(m.group(1)) if m else None


def _get_str(text, key):
    m = re.search(re.escape(key) + r'\s*=\s*"([^"]*)"', text)
    return m.group(1) if m else None


def parse_pois(points_path):
    """Parse plan.Points -> list of POI dicts (name, x, y, z, type, color)."""
    if not os.path.exists(points_path):
        return []
    text = open(points_path, "rb").read().decode("latin-1")
    pois = []
    for m in re.finditer(r"Poi\s*=\s*\{", text):
        o = text.find("{", m.start())
        e = _find_block(text, o)
        blk = text[o + 1:e]
        name = _get_str(blk, "Name")
        if not name:
            continue
        pois.append({
            "name": name,
            "x": _get_float(blk, "XCoord") or 0.0,
            "y": _get_float(blk, "YCoord") or 0.0,
            "z": _get_float(blk, "ZCoord") or 0.0,
            "type": _get_str(blk, "PoiInterpretedType") or "MARKER",
            "color": _get_str(blk, "Color") or "red",
        })
    return pois


def parse_rois(roi_path):
    """Parse plan.roi -> list of ROI dicts (name, type, color, curves).

    Each curve is a list of (x, y, z) tuples in Pinnacle cm.
    """
    if not os.path.exists(roi_path):
        return []
    text = open(roi_path, "rb").read().decode("latin-1")
    rois = []
    for m in re.finditer(r"roi\s*=\s*\{", text):
        o = text.find("{", m.start())
        e = _find_block(text, o)
        blk = text[o + 1:e]
        name = re.search(r"name:\s*(\S+)", blk)
        if not name:
            continue
        roi = {
            "name": name.group(1),
            "type": (re.search(r"roiinterpretedtype:\s*(\S+)", blk).group(1)
                     if re.search(r"roiinterpretedtype:\s*(\S+)", blk) else "ORGAN"),
            "color": (re.search(r"color:\s*(\S+)", blk).group(1)
                      if re.search(r"color:\s*(\S+)", blk) else "red"),
            "curves": [],
        }
        for cm in re.finditer(r"curve\s*=\s*\{", blk):
            co = blk.find("{", cm.start())
            ce = _find_block(blk, co)
            cblk = blk[co + 1:ce]
            pm = re.search(r"points\s*=\s*\{", cblk)
            if not pm:
                continue
            po = cblk.find("{", pm.start())
            pe = _find_block(cblk, po)
            pts_txt = cblk[po + 1:pe]
            vals = re.findall(r"[-0-9.eE+]+", pts_txt)
            curve = []
            for i in range(0, len(vals) - 2, 3):
                try:
                    curve.append((float(vals[i]), float(vals[i + 1]), float(vals[i + 2])))
                except ValueError:
                    continue
            # Pinnacle closed curves repeat the first point as a closing point;
            # DICOM CLOSED_PLANAR contours must not duplicate it.
            if len(curve) > 1:
                d = sum((curve[0][k] - curve[-1][k]) ** 2 for k in range(3)) ** 0.5
                if d < 1e-3:
                    curve = curve[:-1]
            if curve:
                roi["curves"].append(curve)
        rois.append(roi)
    return rois


def read_ct_slices(patient_dir, image_set_id):
    """Read the CT image set; return (slices, study_uid, series_uid, frame_uid).

    ``slices`` is a list of (z_mm, SOPInstanceUID) sorted ascending by z.
    """
    dcm_dir = os.path.join(patient_dir, f"ImageSet_{image_set_id}.DICOM")
    files = sorted(glob.glob(os.path.join(dcm_dir, "*.dcm")))
    if not files:
        return None, None, None, None
    slices = []
    study_uid = series_uid = frame_uid = None
    for f in files:
        try:
            d = pydicom.dcmread(f, stop_before_pixels=True)
        except Exception:
            continue
        z = float(d.ImagePositionPatient[2])
        slices.append((z, d.SOPInstanceUID))
        if study_uid is None:
            study_uid = d.StudyInstanceUID
            series_uid = d.SeriesInstanceUID
            frame_uid = d.FrameOfReferenceUID
    slices.sort(key=lambda t: t[0])
    return slices, study_uid, series_uid, frame_uid


def _nearest_slice_uid(z_mm, slices):
    """Return the SOPInstanceUID of the CT slice closest to ``z_mm``."""
    if not slices:
        return None
    best = min(slices, key=lambda t: abs(t[0] - z_mm))
    return best[1]


def generate_rtstruct(patient_dir, plan_dir_name, output_dir, image_set_id=None):
    """Generate a DICOM RTSTRUCT for one plan. Returns ``out_path``."""
    os.makedirs(output_dir, exist_ok=True)
    plan_dir = os.path.join(patient_dir, plan_dir_name)
    info, full_name, plans = parse_patient_registry(patient_dir)
    plan = next((p for p in plans if f"Plan_{p['plan_id']}" == plan_dir_name), None)
    plan_name = plan["plan_name"] if plan else plan_dir_name
    if image_set_id is None:
        image_set_id = plan["primary_ct"] if plan else 0

    patient_position = parse_patient_position(plan_dir)

    # ---- CT reference -----------------------------------------------------
    slices, study_uid, series_uid, frame_uid = read_ct_slices(patient_dir, image_set_id)

    # ---- ROIs & POIs ------------------------------------------------------
    pois = parse_pois(os.path.join(plan_dir, "plan.Points"))
    rois = parse_rois(os.path.join(plan_dir, "plan.roi"))

    dob = (info["dob"] or "").replace("-", "").replace("/", "").replace(" ", "")
    sex = (info["gender"] or "O")[:1].upper()
    patient_id = str(info["mrn"] or info["patient_id"] or "")

    # ---- file meta --------------------------------------------------------
    sop_uid = generate_uid()
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = RTSTRUCT_SOP_CLASS_UID
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
    ds.SOPClassUID = RTSTRUCT_SOP_CLASS_UID
    ds.SOPInstanceUID = sop_uid
    ds.Modality = "RTSTRUCT"
    ds.StudyInstanceUID = study_uid or generate_uid()
    ds.SeriesInstanceUID = generate_uid()
    ds.SeriesNumber = "1"
    ds.Manufacturer = "Philips"
    ds.ManufacturerModelName = "Pinnacle3"

    # ---- Patient ----------------------------------------------------------
    ds.PatientName = full_name or "UNKNOWN"
    ds.PatientID = patient_id or "UNKNOWN"
    if dob:
        ds.PatientBirthDate = dob
    if sex:
        ds.PatientSex = sex

    # ---- Structure Set ----------------------------------------------------
    ds.StructureSetLabel = plan_name
    ds.StructureSetName = "POIandROI"
    ds.StructureSetDate = time.strftime("%Y%m%d")
    ds.StructureSetTime = time.strftime("%H%M%S")
    ds.ApprovalStatus = "UNAPPROVED"

    # ---- assemble ROI list: POI markers first, then ROI contours ----------
    roi_list = []  # dicts: name, type, color, curves (each curve a list of (x,y,z) cm)
    for p in pois:
        roi_list.append({
            "name": p["name"],
            "type": p["type"] or "MARKER",
            "color": p["color"],
            "points": [{"x": p["x"], "y": p["y"], "z": p["z"]}],
            "is_marker": True,
        })
    for r in rois:
        roi_list.append({
            "name": r["name"],
            "type": r["type"],
            "color": r["color"],
            "curves": r["curves"],
            "is_marker": False,
        })

    ssroi_seq = []
    roic_seq = []
    obs_seq = []
    for idx, r in enumerate(roi_list, start=1):
        # StructureSetROI
        ss = Dataset()
        ss.ROINumber = str(idx)
        ss.ReferencedFrameOfReferenceUID = frame_uid or ""
        ss.ROIName = r["name"]
        ss.ROIGenerationAlgorithm = ""
        ssroi_seq.append(ss)

        # ROI Contour
        rc = Dataset()
        rc.ReferencedROINumber = str(idx)
        rc.ROIDisplayColor = [int(x) for x in COLOR_MAP.get(r["color"], [255, 255, 255])]
        cont_seq = []
        if r["is_marker"]:
            for pt in r["points"]:
                c = Dataset()
                dcm = convert_point(pt, patient_position)
                c.ContourGeometricType = "POINT"
                c.NumberOfContourPoints = "1"
                c.ContourData = [dcm[0], dcm[1], dcm[2]]
                uid = _nearest_slice_uid(dcm[2], slices)
                if uid:
                    img = Dataset()
                    img.ReferencedSOPClassUID = CT_SOP_CLASS_UID
                    img.ReferencedSOPInstanceUID = uid
                    c.ContourImageSequence = Sequence([img])
                cont_seq.append(c)
        else:
            for curve in r["curves"]:
                c = Dataset()
                c.ContourGeometricType = "CLOSED_PLANAR"
                data = []
                z_dcm = None
                for (x, y, z) in curve:
                    dcm = convert_point({"x": x, "y": y, "z": z}, patient_position)
                    data += [dcm[0], dcm[1], dcm[2]]
                    z_dcm = dcm[2]
                c.NumberOfContourPoints = str(len(curve))
                c.ContourData = data
                uid = _nearest_slice_uid(z_dcm, slices)
                if uid:
                    img = Dataset()
                    img.ReferencedSOPClassUID = CT_SOP_CLASS_UID
                    img.ReferencedSOPInstanceUID = uid
                    c.ContourImageSequence = Sequence([img])
                cont_seq.append(c)
        rc.ContourSequence = Sequence(cont_seq)
        roic_seq.append(rc)

        # RT ROI Observation
        obs = Dataset()
        obs.ObservationNumber = str(idx)
        obs.ReferencedROINumber = str(idx)
        obs.RTROIInterpretedType = r["type"]
        obs.ROIInterpreter = ""
        obs_seq.append(obs)

    ds.StructureSetROISequence = Sequence(ssroi_seq)
    ds.ROIContourSequence = Sequence(roic_seq)
    ds.RTROIObservationsSequence = Sequence(obs_seq)

    # ---- Referenced Frame of Reference ------------------------------------
    if frame_uid and study_uid and series_uid and slices:
        fr = Dataset()
        fr.FrameOfReferenceUID = frame_uid
        study = Dataset()
        study.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.2"  # Detached Study Management
        study.ReferencedSOPInstanceUID = study_uid
        series = Dataset()
        series.SeriesInstanceUID = series_uid
        img_seq = []
        for z, uid in slices:
            im = Dataset()
            im.ReferencedSOPClassUID = CT_SOP_CLASS_UID
            im.ReferencedSOPInstanceUID = uid
            img_seq.append(im)
        series.ContourImageSequence = Sequence(img_seq)
        study.RTReferencedSeriesSequence = Sequence([series])
        fr.RTReferencedStudySequence = Sequence([study])
        ds.ReferencedFrameOfReferenceSequence = Sequence([fr])

    out_path = os.path.join(output_dir, f"RTSTRUCT_{plan_dir_name}.dcm")
    pydicom.dcmwrite(out_path, ds, write_like_original=False)
    return out_path
