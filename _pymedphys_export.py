import os, sys, io, warnings

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from pymedphys.pinnacle import PinnacleExport

PATIENT_DIR = r"C:\Users\wt\OneDrive\桌面\2026082002\_patient_5407\Institution_34\Mount_0\Patient_5407"
OUT_ROOT = r"C:\Users\wt\OneDrive\桌面\2026082002\DICOM_1100122014_v2"

pinnacle = PinnacleExport(PATIENT_DIR)

print("Patient:", pinnacle.patient_info["FullName"], " MRN:", pinnacle.patient_info["MedicalRecordNumber"])
print("Images available:")
for im in pinnacle.images:
    hdr = im.image_header
    print(f"  ImageSet_{im.image['ImageSetID']}: {im.image['Modality']} {im.image['NumberOfImages']} slices, series {im.image['SeriesUID']}")

print("\nPlans:")
for plan in pinnacle.plans:
    name = plan.plan_info["PlanName"]
    print(f"  Plan_{plan.plan_info['PlanID']}: {name}  (PrimaryCT ImageSet_{plan.plan_info['PrimaryCTImageSetID']})")
    for t in plan.trials:
        print(f"     Trial: {t['Name']}")

print("\n=== Exporting RTDOSE + RTPLAN per plan ===")
for plan in pinnacle.plans:
    name = plan.plan_info["PlanName"]
    plan_id = plan.plan_info["PlanID"]
    outdir = os.path.join(OUT_ROOT, f"Plan_{plan_id}_{name}")
    os.makedirs(outdir, exist_ok=True)
    # use the first (active) trial
    if plan.trials:
        plan.active_trial = plan.trials[0]["Name"]
        print(f"\n[{name}] active trial: {plan.active_trial['Name']}")
    try:
        PinnacleExport.export_plan(plan, export_path=outdir)
        print(f"  RTPLAN exported -> {outdir}")
    except Exception as e:
        print(f"  RTPLAN export FAILED: {type(e).__name__}: {e}")
    try:
        PinnacleExport.export_dose(plan, export_path=outdir)
        print(f"  RTDOSE exported -> {outdir}")
    except Exception as e:
        print(f"  RTDOSE export FAILED: {type(e).__name__}: {e}")

print("\n=== Done. Output files ===")
for root, dirs, files in os.walk(OUT_ROOT):
    for f in files:
        p = os.path.join(root, f)
        print(f"  {os.path.relpath(p, OUT_ROOT)}  ({os.path.getsize(p):,} bytes)")
