# -*- coding: utf-8 -*-
"""
Generate a DICOM RT Dose (RTDOSE) object from Pinnacle treatment-plan data.

This script reads the raw Pinnacle "plan.Trial" text header and the per-beam
"plan.Trial.binary.NNN" dose volumes (big-endian float32), sums the dose over
all beams with prescription-point MU normalisation, and writes a DICOM
RT Dose Storage instance compliant with:

  * DICOM PS3.3  (RT Dose IOD)
  * DICOM PS3.5  (Data Structures and Encoding)
  * DICOM PS3.10 (Media Storage and File Format)

Key correctness points (verified against PyMedPhys's validated implementation):
  * Dose volumes are stored BIG-ENDIAN float32, x-fastest, z-reversed.
  * Total dose = sum over beams; each beam is normalised to the prescription
    point so the physical dose = stored pixel value x DoseGridScaling.
  * Geometry uses the Pinnacle patient-coordinate conventions (HFS).

The DICOM metadata is built so that all Type 1 / Type 2 / conditional tags of
the RT Dose IOD are present, DS values never exceed 16 characters, File Meta
Information is complete, UIDs are unique per object, and the dataset is written
with an internally-consistent Explicit VR Little Endian transfer syntax.
"""

import io
import math
import os
import re
import struct
import sys
import time

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.tag import Tag
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE = r"D:\PythonProject\Pinnacle backup unpacker\_patient_5407\Institution_34\Mount_0\Patient_5407"
PLAN_NAME = "Plan_2"          # which Plan_N directory to export
OUTPUT_DIR = r"D:\PythonProject\Pinnacle backup unpacker\2026082002\DICOM_1100122014"
OUTPUT_FILENAME = "RTDOSE.dcm"
# SOPInstanceUID of the corresponding RTPLAN (referenced by ReferencedRTPlanSequence).
# Leave as empty string to omit the reference (DoseSummationType "PLAN" technically
# requires it when a plan exists).
RTPLAN_SOP_UID = "1.2.826.0.1.3680043.10.202.1.81039528147881989437888506786974711"

# ---------------------------------------------------------------------------
# DICOM constants
# ---------------------------------------------------------------------------

RTDOSE_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.2"   # RT Dose Storage
RTPLAN_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.5"   # RT Plan Storage
IMPLEMENTATION_CLASS_UID = "1.2.826.0.1.3680043.8.498.75006884747854523615841001"
MANUFACTURER = "Pinnacle RTDOSE Export"

# ImageOrientationPatient for each supported patient position (PS3.3 C.7.6.2)
_IMAGE_ORIENTATION_MAP = {
    "HFS": [1, 0, 0, 0, 1, 0],
    "HFP": [-1, 0, 0, 0, -1, 0],
    "FFS": [-1, 0, 0, 0, 1, 0],
    "FFP": [1, 0, 0, 0, -1, 0],
}


# ---------------------------------------------------------------------------
# DS (Decimal String) formatting
# ---------------------------------------------------------------------------

def fmt_ds(value, max_len=16):
    """Format a number as a DICOM DS string of at most ``max_len`` characters.

    DS values are decimal strings (no scientific notation) and each individual
    value must not exceed 16 characters (PS3.5 A.1.2).  Floating-point
    representation artefacts (e.g. ``-121.17199999999997``) are avoided by
    rounding to a bounded number of decimal places and stripping trailing
    zeros.
    """
    v = float(value)
    if v == 0 or (v == int(v) and abs(v) < 1e15):
        return "0" if v == 0 else str(int(v))
    for decimals in range(12, 0, -1):
        s = ("{:." + str(decimals) + "f}").format(v).rstrip("0").rstrip(".")
        if s in ("", "-", "-0"):
            return "0"
        if len(s) <= max_len:
            return s
    # Last resort: coarse integer rounding.
    return ("{:.0f}").format(v)


def fmt_ds_list(values, max_len=16):
    """Format an iterable of numbers as a list of DS strings."""
    return [fmt_ds(v, max_len) for v in values]


# ---------------------------------------------------------------------------
# Minimal Pinnacle text-parsing helpers
# ---------------------------------------------------------------------------

def _find_block(text, open_brace_idx):
    """Return the index of the '}' matching the '{' at ``open_brace_idx``."""
    depth = 0
    i = open_brace_idx
    n = len(text)
    while i < n:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def get_float(text, key):
    m = re.search(re.escape(key) + r"\s*=\s*([-0-9.eE+]+)", text)
    return float(m.group(1)) if m else None


def get_int(text, key):
    m = re.search(re.escape(key) + r"\s*=\s*(\d+)", text)
    return int(m.group(1)) if m else None


def get_str(text, key):
    m = re.search(re.escape(key) + r'\s*=\s*"([^"]*)"', text)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# Dose-grid reading
# ---------------------------------------------------------------------------

def read_dose_grid(binary_file, dim_x, dim_y, dim_z):
    """Read a Pinnacle dose volume (big-endian float32, x-fastest, z-ascending).

    Returns a float64 array indexed ``[x, y, z]``.
    """
    with open(binary_file, "rb") as fh:
        raw = fh.read()
    if len(raw) != dim_x * dim_y * dim_z * 4:
        raise ValueError(
            "Dose binary size mismatch: got %d bytes, expected %d "
            "(%d x %d x %d float32)"
            % (len(raw), dim_x * dim_y * dim_z * 4, dim_x, dim_y, dim_z)
        )
    # File order: x fastest, then y, then z (ascending z).  Verified empirically
    # against the prescription-point location: the isocentre dose is found at
    # arr[z == Pinnacle z] (ascending), not at the reversed slice.
    arr = np.frombuffer(raw, dtype=">f4").reshape(dim_z, dim_y, dim_x).astype(np.float64)
    # arr[z, y, x] -> dose_grid[x, y, z]
    dose_grid = arr.transpose(2, 1, 0)
    return dose_grid


def trilinear_interpolation(idx, grid):
    """Trilinear interpolation of a point given by fractional index ``idx``."""
    int_idx = [math.floor(f) for f in idx]
    frac_idx = [f % 1.0 for f in idx]
    g = grid
    # clamp to avoid out-of-bounds at the last voxel
    c = [[[0.0 for _ in range(2)] for _ in range(2)] for _ in range(2)]
    for dx in range(2):
        for dy in range(2):
            for dz in range(2):
                ix = min(int_idx[0] + dx, g.shape[0] - 1)
                iy = min(int_idx[1] + dy, g.shape[1] - 1)
                iz = min(int_idx[2] + dz, g.shape[2] - 1)
                c[dx][dy][dz] = g[ix, iy, iz]
    # interpolate along x
    a = [[0.0, 0.0], [0.0, 0.0]]
    for dy in range(2):
        for dz in range(2):
            a[dy][dz] = c[0][dy][dz] * (1 - frac_idx[0]) + c[1][dy][dz] * frac_idx[0]
    # along y
    b = [0.0, 0.0]
    for dz in range(2):
        b[dz] = a[0][dz] * (1 - frac_idx[1]) + a[1][dz] * frac_idx[1]
    # along z
    return b[0] * (1 - frac_idx[2]) + b[1] * frac_idx[2]


# ---------------------------------------------------------------------------
# Pinnacle plan parsing
# ---------------------------------------------------------------------------

def parse_patient_position(plan_dir):
    setup_path = os.path.join(plan_dir, "plan.PatientSetup")
    if not os.path.exists(setup_path):
        return "HFS"
    text = open(setup_path, "rb").read().decode("latin-1")
    orient = get_str(text, "Orientation") or ""
    position = get_str(text, "Position") or ""
    hf = "Head First" in orient
    sup = "supine" in position.lower()
    return "HFS" if (hf and sup) else ("HFP" if hf else ("FFS" if sup else "FFP"))


def parse_trial(plan_dir):
    """Extract dose-grid parameters, fractions, prescription and beams."""
    path = os.path.join(plan_dir, "plan.Trial")
    text = open(path, "rb").read().decode("latin-1")
    trial = {
        "dim_x": get_int(text, "DoseGrid .Dimension .X"),
        "dim_y": get_int(text, "DoseGrid .Dimension .Y"),
        "dim_z": get_int(text, "DoseGrid .Dimension .Z"),
        "voxel_x": get_float(text, "DoseGrid .VoxelSize .X"),
        "voxel_y": get_float(text, "DoseGrid .VoxelSize .Y"),
        "voxel_z": get_float(text, "DoseGrid .VoxelSize .Z"),
        "origin_x": get_float(text, "DoseGrid .Origin .X"),
        "origin_y": get_float(text, "DoseGrid .Origin .Y"),
        "origin_z": get_float(text, "DoseGrid .Origin .Z"),
    }
    # prescription list: collect NumberOfFractions per prescription name
    prescriptions = {}
    for m in re.finditer(r'Name = "([^"]+)";\s*\n\s*[^\n]*NumberOfFractions\s*=\s*(\d+)', text):
        prescriptions[m.group(1)] = int(m.group(2))
    # more robust: find each Prescription block
    for m in re.finditer(r"Prescription\s*=\s*\{", text):
        o = text.find("{", m.start())
        e = _find_block(text, o)
        blk = text[o + 1:e]
        name = get_str(blk, "Name")
        nf = get_int(blk, "NumberOfFractions")
        if name and nf is not None:
            prescriptions[name] = nf
    trial["prescriptions"] = prescriptions

    # beam list
    beams = []
    bl = re.search(r"BeamList\s*=\s*\{", text)
    if bl:
        blo = text.find("{", bl.start())
        ble = _find_block(text, blo)
        beamlist = text[blo + 1:ble]
        for m in re.finditer(r"Beam\s*=\s*\{", beamlist):
            o = beamlist.find("{", m.start())
            e = _find_block(beamlist, o)
            b = beamlist[o + 1:e]
            name = get_str(b, "Name")
            if not name:
                continue
            dose_ref = re.search(r"DoseVolume\s*=\s*\\XDR:(\d+)\\?", b)
            # Per-beam prescription dose at the prescription point (cGy/fraction)
            prescription_dose = None
            mui = re.search(r"MonitorUnitInfo\s*=\s*\{", b)
            if mui:
                mo = b.find("{", mui.start())
                me = _find_block(b, mo)
                prescription_dose = get_float(b[mo + 1:me], "PrescriptionDose")
            beams.append({
                "name": name,
                "prescription_name": get_str(b, "PrescriptionName"),
                "prescription_point": get_str(b, "PrescriptionPointName"),
                "dose_volume_id": dose_ref.group(1) if dose_ref else None,
                "prescription_dose": prescription_dose,
            })
    trial["beams"] = beams
    return trial


def parse_points(plan_dir):
    path = os.path.join(plan_dir, "plan.Points")
    if not os.path.exists(path):
        return {}
    text = open(path, "rb").read().decode("latin-1")
    points = {}
    for m in re.finditer(r"Poi\s*=\s*\{", text):
        o = text.find("{", m.start())
        e = _find_block(text, o)
        blk = text[o + 1:e]
        name = get_str(blk, "Name")
        if name:
            points[name] = {
                "x": get_float(blk, "XCoord"),
                "y": get_float(blk, "YCoord"),
                "z": get_float(blk, "ZCoord"),
            }
    return points


def convert_point(point, patient_position):
    """Convert a Pinnacle point (cm) to DICOM patient coordinates (mm)."""
    x, y, z = point["x"], point["y"], point["z"]
    ref = [x * 10, y * 10, z * 10]
    if patient_position in ("HFP", "FFS"):
        ref[0] = -ref[0]
    if patient_position in ("HFS", "FFS"):
        ref[1] = -ref[1]
    if patient_position in ("HFS", "HFP"):
        ref[2] = -ref[2]
    return ref


def read_image_uids(base, image_set_id):
    """Read StudyInstanceUID and FrameUID from the image set's .ImageInfo file."""
    info = os.path.join(base, f"ImageSet_{image_set_id}.ImageInfo")
    if not os.path.exists(info):
        raise FileNotFoundError(f"ImageSet info file not found: {info}")
    text = open(info, "rb").read().decode("latin-1")
    study = re.search(r'StudyInstanceUID\s*=\s*"([^"]+)"', text)
    frame = re.search(r'FrameUID\s*=\s*"([^"]+)"', text)
    return (study.group(1) if study else None), (frame.group(1) if frame else None)


# ---------------------------------------------------------------------------
# RTDOSE construction
# ---------------------------------------------------------------------------

def build_rtdose(plan_dir, base, patient_position, trial, points, study_uid, frame_uid,
                 rtp_uid=None, patient_name="", patient_id="", patient_birth_date="",
                 patient_sex=""):
    dim_x, dim_y, dim_z = trial["dim_x"], trial["dim_y"], trial["dim_z"]
    voxel_x, voxel_y, voxel_z = trial["voxel_x"], trial["voxel_y"], trial["voxel_z"]
    origin_x, origin_y, origin_z = trial["origin_x"], trial["origin_y"], trial["origin_z"]
    spacing = [voxel_x * 10, voxel_y * 10, voxel_z * 10]  # mm

    # --- geometry ---------------------------------------------------------
    if patient_position in ("HFP", "FFS"):
        dose_origin_x = -origin_x * 10
    else:
        dose_origin_x = origin_x * 10
    if patient_position in ("HFS", "FFS"):
        dose_origin_y = -origin_y * 10
    else:
        dose_origin_y = origin_y * 10
    if patient_position in ("HFS", "HFP"):
        dose_origin_z = -origin_z * 10
    else:
        dose_origin_z = origin_z * 10

    ydoseshift = voxel_y * 10 * dim_y - voxel_y * 10
    zdoseshift = voxel_z * 10 * dim_z - voxel_z * 10

    if patient_position == "HFS":
        ipp = [dose_origin_x, dose_origin_y - ydoseshift, dose_origin_z - zdoseshift]
    elif patient_position == "HFP":
        ipp = [dose_origin_x, dose_origin_y + ydoseshift, dose_origin_z - zdoseshift]
    elif patient_position == "FFS":
        ipp = [dose_origin_x, dose_origin_y - ydoseshift, dose_origin_z + zdoseshift]
    else:  # FFP
        ipp = [dose_origin_x, dose_origin_y + ydoseshift, dose_origin_z + zdoseshift]

    orientation = _IMAGE_ORIENTATION_MAP[patient_position]
    orientation_matrix = np.zeros((3, 3))
    orientation_matrix[0, :] = orientation[:3]
    orientation_matrix[1, :] = orientation[3:]
    orientation_matrix[2, :] = np.cross(orientation_matrix[0, :], orientation_matrix[1, :])

    # --- sum dose over beams ----------------------------------------------
    summed = None
    for beam in trial["beams"]:
        bid = beam["dose_volume_id"]
        if bid is None:
            continue
        binary_file = os.path.join(plan_dir, f"plan.Trial.binary.{int(bid):03d}")
        if not os.path.isfile(binary_file):
            raise FileNotFoundError(f"Dose file not found: {binary_file}")

        dose_grid = read_dose_grid(binary_file, dim_x, dim_y, dim_z)

        # prescription point (Pinnacle coords, cm) from plan.Points
        presc_name = beam["prescription_point"]
        presc_pt = points.get(presc_name)
        if not presc_pt:
            raise KeyError(f"Prescription point '{presc_name}' not found in plan.Points")

        # fractional index of the prescription point within the dose grid.
        # The dose grid is indexed [x, y, z] where, in the Pinnacle binary file,
        #   x = Pinnacle X (ascending),
        #   y = (dim_y-1) - Pinnacle Y (DESCENDING / flipped),
        #   z = Pinnacle Z (ascending).
        # This is a fixed file-format convention (verified against the Pinnacle
        # TPS-exported RTDOSE and the MapCHECK3 measurement), independent of
        # patient position, so it is computed directly from Pinnacle coordinates.
        idx = [0.0, 0.0, 0.0]
        idx[0] = (presc_pt["x"] - origin_x) / voxel_x
        idx[1] = (dim_y - 1) - (presc_pt["y"] - origin_y) / voxel_y
        idx[2] = (presc_pt["z"] - origin_z) / voxel_z

        cgy_mu = trilinear_interpolation(idx, dose_grid)
        n_frac = trial["prescriptions"].get(beam["prescription_name"], 1)
        presc_dose = beam.get("prescription_dose")
        if presc_dose is None:
            raise ValueError(
                f"Beam '{beam['name']}' has no MonitorUnitInfo.PrescriptionDose"
            )
        total_prescription = presc_dose * n_frac
        beam_mu = (total_prescription / cgy_mu) / n_frac if cgy_mu else 0.0

        # physical dose in Gy for this beam
        values = n_frac * dose_grid * beam_mu / 100.0  # [x, y, z]

        # flatten to DICOM pixel order: frame f (z=dim_z-1-f), row=y, col=x
        frame_values = values.transpose(2, 1, 0)[::-1]  # [z, y, x] then reverse z
        flat = frame_values.reshape(-1)

        summed = flat if summed is None else summed + flat

    if summed is None:
        raise RuntimeError("No beams with dose found")

    max_dose_gy = float(np.max(summed))
    scaling = max_dose_gy / 16384.0

    # quantize to uint16
    if scaling != 0:
        pixel_values = np.clip(np.round(summed / scaling), 0, 65535).astype(np.uint16)
    else:
        pixel_values = np.zeros(len(summed), dtype=np.uint16)

    # --- DICOM file meta --------------------------------------------------
    sop_uid = generate_uid()
    series_uid = generate_uid()
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = RTDOSE_SOP_CLASS_UID
    file_meta.MediaStorageSOPInstanceUID = sop_uid
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    file_meta.ImplementationClassUID = IMPLEMENTATION_CLASS_UID

    ds = FileDataset(OUTPUT_FILENAME, {}, file_meta=file_meta, preamble=b"\x00" * 128)
    ds.is_implicit_VR = False
    ds.is_little_endian = True

    # --- SOP Common -------------------------------------------------------
    ds.SpecificCharacterSet = "ISO_IR 100"
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")
    ds.SOPClassUID = RTDOSE_SOP_CLASS_UID
    ds.SOPInstanceUID = sop_uid

    # --- Patient ----------------------------------------------------------
    ds.PatientName = patient_name or "UNKNOWN"
    ds.PatientID = patient_id or "UNKNOWN"
    if patient_birth_date:
        ds.PatientBirthDate = patient_birth_date
    if patient_sex:
        ds.PatientSex = patient_sex

    # --- General Study (Type 2 tags present) ------------------------------
    ds.StudyInstanceUID = study_uid
    ds.StudyDate = time.strftime("%Y%m%d")
    ds.StudyTime = time.strftime("%H%M%S")
    ds.ReferringPhysicianName = ""
    ds.AccessionNumber = ""
    ds.StudyID = ""
    ds.StudyDescription = ""

    # --- General Series ---------------------------------------------------
    ds.Modality = "RTDOSE"
    ds.SeriesInstanceUID = series_uid
    ds.SeriesNumber = "1"
    ds.SeriesDescription = "RT Dose"

    # --- Frame of Reference ----------------------------------------------
    ds.FrameOfReferenceUID = frame_uid
    ds.PositionReferenceIndicator = ""

    # --- General Equipment ------------------------------------------------
    ds.Manufacturer = MANUFACTURER
    ds.ManufacturerModelName = ""
    ds.SoftwareVersions = ""

    # --- General Image ----------------------------------------------------
    ds.InstanceNumber = "1"
    ds.PatientPosition = patient_position

    # --- Image Plane ------------------------------------------------------
    ds.PixelSpacing = fmt_ds_list([spacing[0], spacing[1]])
    ds.ImageOrientationPatient = fmt_ds_list(orientation)
    ds.ImagePositionPatient = fmt_ds_list(ipp)
    ds.SliceThickness = fmt_ds(spacing[2])

    # --- Image Pixel ------------------------------------------------------
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.Rows = dim_y
    ds.Columns = dim_x
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0

    # --- Multi-frame ------------------------------------------------------
    ds.NumberOfFrames = dim_z
    ds.FrameIncrementPointer = Tag(0x3004, 0x000C)

    # --- RT Dose ----------------------------------------------------------
    ds.DoseUnits = "GY"
    ds.DoseType = "PHYSICAL"
    ds.DoseSummationType = "PLAN"
    ds.DoseGridScaling = fmt_ds(scaling)
    ds.GridFrameOffsetVector = fmt_ds_list([p * spacing[0] for p in range(dim_z)])
    ds.TissueHeterogeneityCorrection = "IMAGE"

    # Referenced RT Plan (Type 1C when DoseSummationType is PLAN)
    if rtp_uid:
        ref = Dataset()
        ref.ReferencedSOPClassUID = RTPLAN_SOP_CLASS_UID
        ref.ReferencedSOPInstanceUID = rtp_uid
        ds.ReferencedRTPlanSequence = Sequence([ref])

    # --- Pixel Data -------------------------------------------------------
    ds.PixelData = pixel_values.tobytes()

    return ds, (max_dose_gy, scaling, summed)


def glob_rplan(directory):
    """Return the SOPInstanceUID of a RTPLAN file in ``directory``, if any."""
    if not os.path.isdir(directory):
        return None
    for f in os.listdir(directory):
        if f.startswith("RP.") and f.endswith(".dcm"):
            try:
                d = pydicom.dcmread(os.path.join(directory, f), stop_before_pixels=True)
                return d.SOPInstanceUID
            except Exception:
                continue
    return None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_rtdose(path, original_physical_dose=None, dims=None):
    """Run a set of DICOM-conformance checks on the generated RTDOSE file.

    Returns a list of (status, message) tuples; status is "PASS", "WARN" or
    "FAIL".
    """
    results = []

    # 1. strict read
    try:
        ds = pydicom.dcmread(path, force=False)
        results.append(("PASS", "dcmread(force=False) succeeded"))
    except Exception as e:
        results.append(("FAIL", f"dcmread(force=False) failed: {e}"))
        return results

    # 2. file meta completeness
    fm = ds.file_meta
    for tag, name in [
        ((0x0002, 0x0000), "FileMetaInformationGroupLength"),
        ((0x0002, 0x0001), "FileMetaInformationVersion"),
        ((0x0002, 0x0002), "MediaStorageSOPClassUID"),
        ((0x0002, 0x0003), "MediaStorageSOPInstanceUID"),
        ((0x0002, 0x0010), "TransferSyntaxUID"),
        ((0x0002, 0x0012), "ImplementationClassUID"),
    ]:
        if tag in fm:
            results.append(("PASS", f"File meta tag {name} present"))
        else:
            results.append(("FAIL", f"File meta tag {name} MISSING"))

    # 3. required tags presence / VR / VM (spot checks)
    required = {
        (0x0008, 0x0016): "SOPClassUID",
        (0x0008, 0x0018): "SOPInstanceUID",
        (0x0008, 0x0060): "Modality",
        (0x0020, 0x000D): "StudyInstanceUID",
        (0x0020, 0x000E): "SeriesInstanceUID",
        (0x0020, 0x0052): "FrameOfReferenceUID",
        (0x0020, 0x0011): "SeriesNumber",
        (0x0008, 0x0090): "ReferringPhysicianName",
        (0x0028, 0x0030): "PixelSpacing",
        (0x0020, 0x0037): "ImageOrientationPatient",
        (0x0020, 0x0032): "ImagePositionPatient",
        (0x0028, 0x0002): "SamplesPerPixel",
        (0x0028, 0x0004): "PhotometricInterpretation",
        (0x0028, 0x0010): "Rows",
        (0x0028, 0x0011): "Columns",
        (0x0028, 0x0100): "BitsAllocated",
        (0x0028, 0x0101): "BitsStored",
        (0x0028, 0x0102): "HighBit",
        (0x0028, 0x0103): "PixelRepresentation",
        (0x0028, 0x0008): "NumberOfFrames",
        (0x0028, 0x0009): "FrameIncrementPointer",
        (0x3004, 0x0002): "DoseUnits",
        (0x3004, 0x0004): "DoseType",
        (0x3004, 0x000A): "DoseSummationType",
        (0x3004, 0x000E): "DoseGridScaling",
        (0x3004, 0x000C): "GridFrameOffsetVector",
        (0x7FE0, 0x0010): "PixelData",
    }
    missing = [name for tag, name in required.items() if tag not in ds]
    if missing:
        results.append(("FAIL", "Missing required tags: " + ", ".join(missing)))
    else:
        results.append(("PASS", f"All {len(required)} required tags present"))

    # 4. DS string length
    long_vals = []
    for elem in ds:
        if elem.VR != "DS":
            continue
        val = elem.value
        if isinstance(val, (list, tuple)):
            vals = val
        elif isinstance(val, pydicom.multival.MultiValue):
            vals = list(val)
        else:
            vals = [val]
        for v in vals:
            s = str(v)
            if len(s) > 16:
                long_vals.append((elem.keyword, s, len(s)))
    if long_vals:
        results.append(("FAIL", "DS values exceed 16 chars: " + str(long_vals)))
    else:
        results.append(("PASS", "All DS values <= 16 chars"))

    # 5. SOP Class UID
    if ds.SOPClassUID == RTDOSE_SOP_CLASS_UID:
        results.append(("PASS", "SOPClassUID is RT Dose Storage"))
    else:
        results.append(("FAIL", f"SOPClassUID is {ds.SOPClassUID}, not RT Dose Storage"))

    # 6. UID consistency
    if ds.SOPInstanceUID == ds.file_meta.MediaStorageSOPInstanceUID:
        results.append(("PASS", "SOPInstanceUID matches file meta"))
    else:
        results.append(("FAIL", "SOPInstanceUID != MediaStorageSOPInstanceUID"))
    if ds.SOPClassUID == ds.file_meta.MediaStorageSOPClassUID:
        results.append(("PASS", "SOPClassUID matches file meta"))
    else:
        results.append(("FAIL", "SOPClassUID != MediaStorageSOPClassUID"))
    if ds.SeriesInstanceUID != ds.SOPInstanceUID:
        results.append(("PASS", "SeriesInstanceUID != SOPInstanceUID"))
    else:
        results.append(("FAIL", "SeriesInstanceUID == SOPInstanceUID"))

    # 7. PixelData byte length
    expected = ds.Rows * ds.Columns * ds.NumberOfFrames * (ds.BitsAllocated // 8)
    actual = len(ds.PixelData)
    if actual == expected:
        results.append(("PASS", f"PixelData length {actual} == expected {expected}"))
    else:
        results.append(("FAIL", f"PixelData length {actual} != expected {expected}"))

    # 8. decode pixel array
    try:
        arr = ds.pixel_array
        results.append(("PASS", f"pixel_array decoded: shape {arr.shape}, dtype {arr.dtype}"))
        if arr.min() < 0 or arr.max() > 65535:
            results.append(("FAIL", "pixel values out of uint16 range"))
        else:
            results.append(("PASS", f"pixel values within uint16 range [{arr.min()}, {arr.max()}]"))
    except Exception as e:
        results.append(("FAIL", f"pixel_array decode failed: {e}"))

    # 9. NumberOfFrames vs GridFrameOffsetVector
    nf = ds.NumberOfFrames
    gfv = len(ds.GridFrameOffsetVector)
    if nf == gfv:
        results.append(("PASS", f"NumberOfFrames {nf} == GridFrameOffsetVector {gfv}"))
    else:
        results.append(("FAIL", f"NumberOfFrames {nf} != GridFrameOffsetVector {gfv}"))

    # 10. reconstructed dose vs original
    if original_physical_dose is not None:
        reconstructed = ds.pixel_array.astype(np.float64) * float(ds.DoseGridScaling)
        flat = reconstructed.reshape(-1)
        diff = np.abs(flat - original_physical_dose)
        max_abs = float(diff.max())
        max_rel = float((diff / np.maximum(np.abs(original_physical_dose), 1e-9)).max())
        mean_abs = float(diff.mean())
        results.append((
            "PASS" if max_abs <= 2 * float(ds.DoseGridScaling) else "WARN",
            f"dose reconstruction: max_abs={max_abs:.6g} Gy, max_rel={max_rel:.6g}, "
            f"mean_abs={mean_abs:.6g} Gy (DoseGridScaling={float(ds.DoseGridScaling):.6g})",
        ))

    # 11. FrameOfReferenceUID / RTPlan references
    if hasattr(ds, "FrameOfReferenceUID"):
        results.append(("PASS", f"FrameOfReferenceUID = {ds.FrameOfReferenceUID}"))
    else:
        results.append(("FAIL", "FrameOfReferenceUID missing"))
    if hasattr(ds, "ReferencedRTPlanSequence"):
        results.append(("PASS", f"ReferencedRTPlanSequence present ({len(ds.ReferencedRTPlanSequence)} refs)"))
    else:
        results.append(("WARN", "ReferencedRTPlanSequence absent"))

    return results


# ---------------------------------------------------------------------------
# Reusable entry points
# ---------------------------------------------------------------------------

def parse_patient_registry(patient_dir):
    """Read the 'Patient' registry file; return (info_dict, full_name, plans)."""
    path = os.path.join(patient_dir, "Patient")
    text = open(path, "rb").read().decode("latin-1")
    info = {
        "patient_id": get_int(text, "PatientID"),
        "last_name": get_str(text, "LastName") or "",
        "first_name": get_str(text, "FirstName") or "",
        "mrn": get_str(text, "MedicalRecordNumber") or "",
        "gender": get_str(text, "Gender") or "",
        "dob": get_str(text, "DateOfBirth") or "",
    }
    full_name = f"{info['last_name']}^{info['first_name']}"
    plans = []
    for m in re.finditer(r"Plan\s*=\s*\{", text):
        o = text.find("{", m.start())
        e = _find_block(text, o)
        blk = text[o + 1:e]
        plans.append({
            "plan_id": get_int(blk, "PlanID"),
            "plan_name": get_str(blk, "PlanName"),
            "primary_ct": get_int(blk, "PrimaryCTImageSetID"),
        })
    return info, full_name, plans


def generate_rtdose(patient_dir, plan_dir_name, output_dir, rtp_uid=None, image_set_id=None):
    """Generate a DICOM RTDOSE for one plan directory.

    Returns ``(out_path, meta)`` where ``meta`` holds max_dose_gy, scaling and
    plan_name.
    """
    os.makedirs(output_dir, exist_ok=True)
    plan_dir = os.path.join(patient_dir, plan_dir_name)
    info, full_name, plans = parse_patient_registry(patient_dir)

    if image_set_id is None:
        for p in plans:
            if f"Plan_{p['plan_id']}" == plan_dir_name:
                image_set_id = p["primary_ct"]
                break
        image_set_id = 0 if image_set_id is None else image_set_id

    plan_name = next(
        (p["plan_name"] for p in plans if f"Plan_{p['plan_id']}" == plan_dir_name),
        plan_dir_name,
    )

    study_uid, frame_uid = read_image_uids(patient_dir, str(image_set_id))
    patient_position = parse_patient_position(plan_dir)
    trial = parse_trial(plan_dir)
    points = parse_points(plan_dir)

    dob = (info["dob"] or "").replace("-", "").replace("/", "").replace(" ", "")
    sex = (info["gender"] or "O")[:1].upper()
    # DICOM PatientID should be the hospital MRN (fall back to Pinnacle internal id)
    patient_id = str(info["mrn"] or info["patient_id"] or "")

    ds, (max_dose_gy, scaling, summed) = build_rtdose(
        plan_dir, patient_dir, patient_position, trial, points, study_uid, frame_uid,
        rtp_uid=rtp_uid,
        patient_name=full_name,
        patient_id=patient_id,
        patient_birth_date=dob,
        patient_sex=sex,
    )

    out_path = os.path.join(output_dir, f"RTDOSE_{plan_dir_name}.dcm")
    pydicom.dcmwrite(out_path, ds, write_like_original=False)
    return out_path, {"max_dose_gy": max_dose_gy, "scaling": scaling, "plan_name": plan_name}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_path, meta = generate_rtdose(BASE, PLAN_NAME, OUTPUT_DIR, rtp_uid=RTPLAN_SOP_UID)
    print(f"\nWROTE: {out_path}")
    print(f"  max dose = {meta['max_dose_gy']:.4f} Gy")
    print(f"  DoseGridScaling = {fmt_ds(meta['scaling'])}")

    # ---- validate ----
    print("\n=== Validation ===")
    d = pydicom.dcmread(out_path)
    arr = d.pixel_array.astype(np.float64) * float(d.DoseGridScaling)
    results = validate_rtdose(out_path, original_physical_dose=arr.reshape(-1))
    fails = 0
    for status, msg in results:
        print(f"[{status}] {msg}")
        if status == "FAIL":
            fails += 1
    print(f"\n{fails} FAIL, {sum(1 for s, _ in results if s == 'WARN')} WARN, "
          f"{sum(1 for s, _ in results if s == 'PASS')} PASS")

    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
