# Pinnacle Backup DICOM Exporter
### Pinnacle 备份 DICOM 导出工具

A desktop tool (Tkinter GUI) that unpacks Pinnacle Treatment Planning System (v16.2) `.tar` backups and exports them as DICOM objects: **CT, RTSTRUCT, RTPLAN, RTDOSE**.

将 Pinnacle 治疗计划系统（v16.2）的 `.tar` 备份文件解包，并导出为 DICOM 标准文件的桌面工具（Tkinter 图形界面），导出对象包括 **CT、RTSTRUCT、RTPLAN、RTDOSE**。

---

## Features · 功能

- Read multiple `.tar` backups, scan and list all patients and their plans.
  读取多个 `.tar` 备份，扫描并列出所有患者及其计划。
- Tree view: expand a patient → check the plans to export; the patient checkbox toggles all of that patient's plans.
  树形界面：展开患者 → 勾选要导出的计划；患者前的勾选框 = 全选 / 取消全选该患者的所有计划。
- Choose which data types to export: **CT / RTSTRUCT / RTPLAN / RTDOSE**.
  可选导出数据类型：CT / RTSTRUCT / RTPLAN / RTDOSE。
- Generates from Pinnacle native formats:
  从 Pinnacle 原生格式生成：
  - **RTDOSE** — per-beam big-endian float32 dose summed and normalised at the prescription point; DICOM-compliant (matches the TPS export voxel-by-voxel, correlation 0.9999).
    **RTDOSE** — 逐射束大端 float32 剂量求和、按处方点归一化，DICOM 合规（与 TPS 导出逐体素比对，相关系数 0.9999）。
  - **RTPLAN** — beams / control points / jaw & MLC leaf positions.
    **RTPLAN** — 射束 / 控制点 / 铅门与 MLC 叶片位置。
  - **RTSTRUCT** — ROI contours + POI markers (point-by-point identical to the TPS export).
    **RTSTRUCT** — ROI 轮廓 + POI 标记点（与 TPS 导出逐点一致）。
  - **CT** — copies the ready-made DICOM CT series from the backup.
    **CT** — 直接复制备份中现成的 DICOM CT 序列。

## Dependencies · 环境依赖

- Python 3.8+
- [pydicom](https://github.com/pydicom/pydicom)
- [numpy](https://numpy.org/)

```bash
pip install pydicom numpy
```

## Usage · 使用

```bash
python pinnacle_app.py
```

1. Click **"Select .tar files"** and choose the Pinnacle backups (multiple allowed).
   点击「Select .tar files」选择 Pinnacle 备份（可多选）。
2. Click **"Scan"** to list all patients.
   点击「Scan」列出所有患者。
3. Expand a patient and check the plans to export.
   展开患者，勾选要导出的计划。
4. Check the data types to export (CT / RTSTRUCT / RTPLAN / RTDOSE).
   勾选要导出的数据类型（CT / RTSTRUCT / RTPLAN / RTDOSE）。
5. Click **"Export Selected Plans"** and pick an output directory.
   点击「Export Selected Plans」，选择输出目录。

## Modules · 模块说明

| File · 文件 | Purpose · 作用 |
|------|------|
| `pinnacle_app.py` | Tkinter front end · Tkinter 前端入口 |
| `pinnacle_core.py` | Back end: scan tar, extract patient, orchestrate export · 核心后端：扫描 tar、提取患者、导出编排 |
| `_make_rtdose.py` | RTDOSE generator (DICOM-compliant, with auto-validation) · RTDOSE 生成（含 DICOM 合规与自动校验） |
| `_make_rtplan.py` | RTPLAN generator · RTPLAN 生成 |
| `_make_rtstruct.py` | RTSTRUCT generator · RTSTRUCT 生成 |
| `_pymedphys_export.py` | Optional PyMedPhys-based reference implementation · 基于 PyMedPhys 的参考实现（可选） |

## Technical notes · 技术要点

- Pinnacle dose binary: **big-endian float32**, x-fastest, y-descending, z-ascending.
  Pinnacle 剂量二进制：**大端 float32**，x 最快、y 降序、z 升序。
- Dose normalisation: per beam, normalised by `MonitorUnitInfo.PrescriptionDose` at the prescription point, then summed.
  剂量归一化：逐射束按 `MonitorUnitInfo.PrescriptionDose` 在处方点归一化后求和。
- Coordinate conversion (HFS): DICOM = `(X*10, -Y*10, -Z*10)` mm.
  坐标转换（HFS）：DICOM = `(X*10, -Y*10, -Z*10)` mm。
- Transfer syntax: Explicit VR Little Endian (`1.2.840.10008.1.2.1`).
  DICOM 传输语法：Explicit VR Little Endian（`1.2.840.10008.1.2.1`）。

## Disclaimer · 免责声明

This tool is for **research use only**. The generated DICOM (especially the RTPLAN jaw/MLC sign conventions) **must be verified against the treatment planning system's own export before clinical use**. Users are responsible for protecting patient privacy data and must never commit patient data to a public repository.

本工具仅供**研究用途**。生成的 DICOM 文件（尤其 RTPLAN 的铅门 / MLC 符号约定）在**临床使用前需与治疗计划系统的导出结果比对验证**。使用者需自行确保患者隐私数据的安全，切勿将患者数据提交到公开仓库。

## License · 许可

[MIT License](LICENSE)
