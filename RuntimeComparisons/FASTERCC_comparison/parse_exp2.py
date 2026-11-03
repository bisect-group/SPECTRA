#!/usr/bin/env python3
"""Tidy table of every expanded-Recon3D (experiment 2) run, parsed from logs.

Reads the SLURM logs rather than the stored .mat instances: at k=50 a single
instance is ~450,000 reactions and loading 36 of them takes longer than the
runs did. The logs carry the same numbers.

It also recovers the spectracc results from tasks SLURM marked
OUT_OF_MEMORY. bench_variants checkpoints per ARM, so when fastercc_r0
exhausted the node the spectracc result for that instance was already
written and printed.
"""
import csv
import glob
import os
import re

rows = []
for pre in ("exp2a_", "exp2k30_", "exp2k50_"):
    for f in sorted(glob.glob(f"logs/{pre}*.out")):
        txt = open(f, errors="replace").read()
        m = re.search(r"task \d+ -> k=(\d+) level=(\d+) draw=(\d+)(?: timing=(\d+))?", txt)
        if not m:
            continue
        k, lv, dr = int(m.group(1)), int(m.group(2)), int(m.group(3))
        tm = int(m.group(4) or 1)
        nr = re.search(r":\s*(\d+) rxns,", txt)
        nrx = int(nr.group(1)) if nr else ""
        errf = f.replace(".out", ".err")
        err_txt = open(errf, errors="replace").read() if os.path.exists(errf) else ""
        for arm in ("spectracc", "fastercc_r0"):
            a = re.search(rf"^  {arm}\s+([\d.]+) s.*?\|cons\|=\s*(\d+)\s+LPs=(\d+)",
                          txt, re.M)
            if a:
                rows.append(dict(k=k, level=lv, draw=dr, timing=tm, arm=arm, nRxns=nrx,
                                 runtimeSec=float(a.group(1)),
                                 nConsistent=int(a.group(2)),
                                 nLPs=int(a.group(3)), err=""))
            else:
                e = re.search(rf"^  {arm}\s+FAILED: (.*)$", err_txt, re.M)
                msg = e.group(1) if e else "no result (killed: OOM or wall clock)"
                msg = " ".join(msg.split()).replace(",", ";")[:120]
                rows.append(dict(k=k, level=lv, draw=dr, timing=tm, arm=arm, nRxns=nrx,
                                 runtimeSec="", nConsistent="", nLPs="", err=msg))

out = "exp2_expanded_all.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), quoting=csv.QUOTE_MINIMAL)
    w.writeheader()
    w.writerows(rows)
print(f"wrote {len(rows)} rows -> {out}")
