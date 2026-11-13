# -*- coding: utf-8 -*-
"""
Pinnacle backup -> DICOM export tool (Tkinter desktop GUI).

Features:
  1. Select multiple Pinnacle backup .tar files;
  2. Scan and list all patients and their plans;
  3. Check the plans to export (patient checkbox = select/deselect all);
  4. Export the selected plans as RTPLAN / RTDOSE / RTSTRUCT / CT.
"""

import os
import queue
import tempfile
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import pinnacle_core as core


class PinnacleApp:
    def __init__(self, root):
        self.root = root
        root.title("Pinnacle Backup → DICOM Export Tool")
        root.geometry("960x700")

        self.tar_files = []
        self.patients = []
        self.msg_queue = queue.Queue()

        self._build_ui()
        self.root.after(100, self._drain_queue)

    # ---- UI -----------------------------------------------------------------
    def _build_ui(self):
        bar = ttk.Frame(self.root, padding=8)
        bar.pack(fill="x")
        ttk.Button(bar, text="Select .tar files", command=self.choose_tars).pack(side="left")
        self.tar_combo = ttk.Combobox(bar, state="readonly", width=28)
        self.tar_combo.pack(side="left", padx=8)
        self.tar_combo.bind("<<ComboboxSelected>>", self._on_tar_select)
        self._update_tar_combo()
        ttk.Button(bar, text="Scan", command=self.scan).pack(side="left", padx=4)
        ttk.Button(bar, text="Select All", command=lambda: self._set_all(True)).pack(side="left", padx=2)
        ttk.Button(bar, text="Deselect All", command=lambda: self._set_all(False)).pack(side="left", padx=2)

        tree_frame = ttk.Frame(self.root)
        tree_frame.pack(fill="both", expand=True, padx=8, pady=4)
        self.tree = ttk.Treeview(tree_frame, columns=("check", "label"), show="tree headings")
        self.tree.heading("#0", text="Expand")
        self.tree.heading("check", text="Export")
        self.tree.heading("label", text="Patient / Plan")
        self.tree.column("#0", width=55, stretch=False, anchor="center")
        self.tree.column("check", width=55, anchor="center", stretch=False)
        self.tree.column("label", width=815, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        scroll = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        scroll.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=scroll.set)
        self.tree.bind("<Button-1>", self._on_click)

        bottom = ttk.Frame(self.root, padding=8)
        bottom.pack(fill="x")
        self.export_btn = ttk.Button(bottom, text="Export Selected Plans", command=self.export_selected)
        self.export_btn.pack(side="left")
        ttk.Label(bottom, text="Export data:").pack(side="left", padx=(12, 2))
        self.var_ct = tk.BooleanVar(value=True)
        self.var_rs = tk.BooleanVar(value=True)
        self.var_rp = tk.BooleanVar(value=True)
        self.var_rd = tk.BooleanVar(value=True)
        ttk.Checkbutton(bottom, text="CT", variable=self.var_ct).pack(side="left")
        ttk.Checkbutton(bottom, text="RTSTRUCT", variable=self.var_rs).pack(side="left")
        ttk.Checkbutton(bottom, text="RTPLAN", variable=self.var_rp).pack(side="left")
        ttk.Checkbutton(bottom, text="RTDOSE", variable=self.var_rd).pack(side="left")
        self.progress = ttk.Progressbar(bottom, mode="indeterminate")
        self.progress.pack(side="left", fill="x", expand=True, padx=8)

        self.log_text = tk.Text(self.root, height=11, state="disabled", wrap="word")
        self.log_text.pack(fill="both", padx=8, pady=(0, 8))

    def log(self, msg):
        self.msg_queue.put(msg)

    def _drain_queue(self):
        try:
            while True:
                msg = self.msg_queue.get_nowait()
                self.log_text.configure(state="normal")
                self.log_text.insert("end", msg + "\n")
                self.log_text.see("end")
                self.log_text.configure(state="disabled")
        except queue.Empty:
            pass
        self.root.after(100, self._drain_queue)

    # ---- actions ------------------------------------------------------------
    def choose_tars(self):
        files = filedialog.askopenfilenames(
            title="Select Pinnacle backup .tar files",
            filetypes=[("tar files", "*.tar"), ("All files", "*.*")],
        )
        if files:
            self.tar_files = list(files)
            self._update_tar_combo()

    def _update_tar_combo(self):
        n = len(self.tar_files)
        if n:
            self.tar_combo["values"] = self.tar_files
            self.tar_combo.set(f"{n} file(s) selected (click to expand)")
        else:
            self.tar_combo["values"] = []
            self.tar_combo.set("No files selected")

    def _on_tar_select(self, _event):
        # keep the summary text; the dropdown itself already listed the files
        self._update_tar_combo()

    def scan(self):
        if not self.tar_files:
            messagebox.showwarning("Warning", "Please select .tar file(s) first")
            return
        self.log("Scanning .tar files ...")
        self.patients = []
        for t in self.tar_files:
            self.log(f"  Scanning {os.path.basename(t)} ...")
            try:
                self.patients.extend(core.scan_tar(t))
            except Exception as e:
                self.log(f"  [Error] scan failed: {e}")
        self.log(f"Found {len(self.patients)} patient(s).")
        self._populate()

    def _populate(self):
        self.tree.delete(*self.tree.get_children())
        for i, p in enumerate(self.patients):
            pid = f"p{i}"
            label = f"{p['dir_name']}    {p['name']}    (MRN {p['mrn']})    {len(p['plans'])} plan(s)"
            self.tree.insert("", "end", iid=pid, text="", values=("☐", label), open=False)
            for j, pl in enumerate(p["plans"]):
                pj = f"p{i}_pl{j}"
                self.tree.insert(pid, "end", iid=pj, text="",
                                 values=("☐", f"Plan_{pl['plan_id']}    {pl['plan_name']}"))

    def _on_click(self, event):
        iid = self.tree.identify_row(event.y)
        col = self.tree.identify_column(event.x)
        if not iid or col != "#1":
            return
        self._toggle(iid)

    def _toggle(self, iid):
        if "_pl" not in iid:
            # patient node -> toggle all children
            children = self.tree.get_children(iid)
            if not children:
                return
            all_checked = all(self.tree.set(c, "check") == "☑" for c in children)
            new = "☐" if all_checked else "☑"
            self.tree.set(iid, "check", new)
            for c in children:
                self.tree.set(c, "check", new)
        else:
            cur = self.tree.set(iid, "check")
            new = "☐" if cur == "☑" else "☑"
            self.tree.set(iid, "check", new)
            parent = self.tree.parent(iid)
            if parent:
                children = self.tree.get_children(parent)
                states = [self.tree.set(c, "check") for c in children]
                if all(s == "☑" for s in states):
                    self.tree.set(parent, "check", "☑")
                elif any(s == "☑" for s in states):
                    self.tree.set(parent, "check", "◐")
                else:
                    self.tree.set(parent, "check", "☐")

    def _set_all(self, checked):
        ch = "☑" if checked else "☐"
        for p_iid in self.tree.get_children(""):
            self.tree.set(p_iid, "check", ch)
            for c in self.tree.get_children(p_iid):
                self.tree.set(c, "check", ch)

    def _collect_checked(self):
        """Return list of (patient_idx, plan_idx) for checked plans."""
        result = []
        for p_iid in self.tree.get_children(""):
            pi = int(p_iid[1:])
            for c in self.tree.get_children(p_iid):
                if self.tree.set(c, "check") == "☑":
                    pj = int(c.split("_pl")[1])
                    result.append((pi, pj))
        return result

    def export_selected(self):
        checked = self._collect_checked()
        if not checked:
            messagebox.showinfo("Info", "Please check the plan(s) to export")
            return
        types = set()
        if self.var_ct.get():
            types.add("CT")
        if self.var_rs.get():
            types.add("RTSTRUCT")
        if self.var_rp.get():
            types.add("RTPLAN")
        if self.var_rd.get():
            types.add("RTDOSE")
        if not types:
            messagebox.showinfo("Info", "Please check at least one data type (CT/RTSTRUCT/RTPLAN/RTDOSE)")
            return
        out_dir = filedialog.askdirectory(title="Select output directory")
        if not out_dir:
            return
        self.export_btn.configure(state="disabled")
        self.progress.start(12)
        threading.Thread(target=self._export_worker, args=(checked, types, out_dir), daemon=True).start()

    def _export_worker(self, checked, types, out_dir):
        work_dir = tempfile.mkdtemp(prefix="pinnacle_export_")
        try:
            for pi, pj in checked:
                patient = self.patients[pi]
                plan = patient["plans"][pj]
                plan_dir = f"Plan_{plan['plan_id']}"
                tar = self._find_tar_for(patient["dir_name"])
                if tar is None:
                    self.log(f"[Error] tar file not found for {patient['dir_name']}")
                    continue
                self.log(f"=== Exporting {patient['dir_name']} / {plan_dir} ({plan['plan_name']}) ===")
                try:
                    core.export_plan(tar, patient["dir_name"], plan_dir, out_dir, work_dir,
                                     progress=self.log, types=types)
                except Exception as e:
                    self.log(f"[Error] export failed: {e}")
            self.log("All done.")
        finally:
            self.root.after(0, self._export_done)

    def _find_tar_for(self, patient_dir_name):
        for p in self.patients:
            if p["dir_name"] == patient_dir_name:
                return p.get("tar_path") or (self.tar_files[0] if self.tar_files else None)
        return self.tar_files[0] if self.tar_files else None

    def _export_done(self):
        self.progress.stop()
        self.export_btn.configure(state="normal")


def main():
    root = tk.Tk()
    PinnacleApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
