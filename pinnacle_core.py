# -*- coding: utf-8 -*-
"""
Core back-end for the "Pinnacle backup -> DICOM export" tool.

Provides:
  * scan_tar()        -- list patients and their plans inside a Pinnacle .tar backup
  * extract_patient() -- extract a single patient directory from the tar
  * export_plan()     -- generate RTDOSE + RTPLAN for one plan
"""

import glob
import os
import re
import shutil
import tarfile


def parse_patient_registry_text(text):
    """Parse a 'Patient' registry file's text; return (info, plans)."""
    def _get_int(key):
        m = re.search(re.escape(key) + r"\s*=\s*(\d+)", text)
        return int(m.group(1)) if m else None

    def _get_str(key):
        m = re.search(re.escape(key) + r'\s*=\s*"([^"]*)"', text)
        return m.group(1) if m else ""

    def _find_block(start):
        depth = 0
        i = start
        while i < len(text):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
                if depth == 0:
                    return i
            i += 1
        return -1

    info = {
        "patient_id": _get_int("PatientID"),
        "last_name": _get_str("LastName"),
        "first_name": _get_str("FirstName"),
        "mrn": _get_str("MedicalRecordNumber"),
        "gender": _get_str("Gender"),
        "dob": _get_str("DateOfBirth"),
    }
    plans = []
    for m in re.finditer(r"Plan\s*=\s*\{", text):
        o = text.find("{", m.start())
        e = _find_block(o)
        blk = text[o + 1:e]
        pid = re.search(r"PlanID\s*=\s*(\d+)", blk)
        pname = re.search(r'PlanName\s*=\s*"([^"]*)"', blk)
        ct = re.search(r"PrimaryCTImageSetID\s*=\s*(\d+)", blk)
        plans.append({
            "plan_id": int(pid.group(1)) if pid else None,
            "plan_name": pname.group(1) if pname else "",
            "primary_ct": int(ct.group(1)) if ct else 0,
        })
    return info, plans


def scan_tar(tar_path):
    """Scan a Pinnacle backup tar; return a list of patient dicts (each with plans)."""
    patients = {}
    with tarfile.open(tar_path, "r") as tar:
        for m in tar.getmembers():
            if not m.isfile():
                continue
            parts = m.name.rstrip("/").split("/")
            if len(parts) < 4 or parts[-1] != "Patient" or not parts[-2].startswith("Patient_"):
                continue
            patient_dir_name = parts[-2]
            if patient_dir_name in patients:
                continue
            f = tar.extractfile(m)
            text = f.read().decode("latin-1", errors="replace")
            f.close()
            info, plans = parse_patient_registry_text(text)
            patients[patient_dir_name] = {
                "dir_name": patient_dir_name,
                "institution": parts[-4],
                "mount": parts[-3],
                "patient_id": info["patient_id"],
                "name": f"{info['last_name']}^{info['first_name']}",
                "mrn": info["mrn"],
                "gender": info["gender"],
                "dob": info["dob"],
                "plans": plans,
                "tar_path": tar_path,
            }
    return list(patients.values())


def extract_patient(tar_path, patient_dir_name, dest_root):
    """Extract a patient directory from the tar into ``dest_root``.

    Only the files needed for RTDOSE/RTPLAN export are extracted (raw ``.img``
    image volumes and timestamped log files are skipped).  Colons in filenames
    are replaced with underscores for Windows compatibility.

    Returns the absolute path of the extracted patient directory.
    """
    os.makedirs(dest_root, exist_ok=True)
    marker = "/" + patient_dir_name + "/"
    with tarfile.open(tar_path, "r") as tar:
        members = []
        for m in tar.getmembers():
            if marker not in m.name:
                continue
            base = m.name.rstrip("/").split("/")[-1]
            if ":" in base:            # timestamped ErrorLog/Transcript (invalid on Windows)
                continue
            if base.endswith(".img"):  # large raw image volume, not needed
                continue
            members.append(m)
        if not members:
            raise ValueError(f"Patient directory '{patient_dir_name}' not found in {tar_path}")
        for m in members:
            target = os.path.join(dest_root, m.name.replace(":", "_"))
            if m.isdir():
                os.makedirs(target, exist_ok=True)
            else:
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with tar.extractfile(m) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)

    for root, _dirs, _files in os.walk(dest_root):
        if os.path.basename(root) == patient_dir_name:
            return root
    raise FileNotFoundError(f"Extracted patient directory '{patient_dir_name}' not found")


def export_ct(patient_dir, image_set_id, output_dir):
    """Copy the CT DICOM images to ``output_dir/CT``.

    Returns the list of copied file paths (empty if the image set has no DICOM
    CT, e.g. a Pinnacle-native phantom volume).
    """
    dcm_dir = os.path.join(patient_dir, f"ImageSet_{image_set_id}.DICOM")
    files = sorted(glob.glob(os.path.join(dcm_dir, "*.dcm")))
    if not files:
        return []
    ct_dir = os.path.join(output_dir, "CT")
    os.makedirs(ct_dir, exist_ok=True)
    out = []
    for f in files:
        dst = os.path.join(ct_dir, os.path.basename(f))
        if not os.path.exists(dst):
            shutil.copy2(f, dst)
        out.append(dst)
    return out


def export_plan(tar_path, patient_dir_name, plan_dir_name, output_dir, work_dir, progress=None,
                types=None):
    """Extract a patient and export one plan's selected DICOM data.

    ``plan_dir_name`` is e.g. ``"Plan_0"``.  ``types`` is an iterable of the
    data to export: any of ``"CT"``, ``"RTSTRUCT"``, ``"RTPLAN"``, ``"RTDOSE"``
    (default: all four).  Returns
    ``(rtp_path, rtd_path, ct_files, rtstruct_path, meta)``.
    """
    types = set(types) if types else {"CT", "RTSTRUCT", "RTPLAN", "RTDOSE"}

    def log(msg):
        if progress:
            progress(msg)
        else:
            print(msg)

    log(f"Extracting patient {patient_dir_name} ...")
    patient_dir = extract_patient(tar_path, patient_dir_name, work_dir)

    # determine the plan's primary CT image set
    from _make_rtdose import parse_patient_registry
    _info, _name, plans = parse_patient_registry(patient_dir)
    plan = next((p for p in plans if f"Plan_{p['plan_id']}" == plan_dir_name), None)
    image_set_id = plan["primary_ct"] if plan else 0

    rtp_path = rtd_path = rtstruct_path = None
    ct_files = []
    meta = {}

    # RTPLAN + RTDOSE (RTDOSE references the RTPLAN, so RTPLAN is produced
    # whenever RTDOSE is requested)
    if "RTPLAN" in types or "RTDOSE" in types:
        log(f"Generating RTPLAN ({plan_dir_name})...")
        from _make_rtplan import generate_rtplan
        rtp_path, rtp_uid = generate_rtplan(patient_dir, plan_dir_name, output_dir)

        if "RTDOSE" in types:
            log(f"Generating RTDOSE ({plan_dir_name})...")
            from _make_rtdose import generate_rtdose
            rtd_path, meta = generate_rtdose(patient_dir, plan_dir_name, output_dir,
                                             rtp_uid=rtp_uid)

    # CT + RTSTRUCT (RTSTRUCT references the CT, so CT is exported too when
    # RTSTRUCT is requested)
    if "CT" in types or "RTSTRUCT" in types:
        ct_files = export_ct(patient_dir, image_set_id, output_dir)
        if ct_files:
            if "CT" in types:
                log(f"Exporting CT ({len(ct_files)} slice(s))...")
            if "RTSTRUCT" in types:
                log(f"Generating RTSTRUCT ({plan_dir_name})...")
                from _make_rtstruct import generate_rtstruct
                rtstruct_path = generate_rtstruct(patient_dir, plan_dir_name, output_dir,
                                                  image_set_id=image_set_id)
        else:
            log("No DICOM CT for this plan, skipping CT/RTSTRUCT export")

    if "RTDOSE" in types:
        log(f"Done: {plan_dir_name}, max dose {meta['max_dose_gy']:.3f} Gy")
    else:
        log(f"Done: {plan_dir_name}")
    return rtp_path, rtd_path, ct_files, rtstruct_path, meta
