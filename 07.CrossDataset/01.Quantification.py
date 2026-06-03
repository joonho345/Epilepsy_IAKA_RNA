#!/usr/bin/env python3
"""
01.Quantification.py
Build count matrices and coldata for IAKA, IHKA, and PILO models.
IHKA/PILO are extracted from the NatCom merged matrix (Kim et al., Nat Commun, 2025).
"""
import csv

# SET THESE PATHS
BASE_DIR     = "/path/to/project"
ANIMAL_SHEET = "/path/to/Animal_sheet.csv"       # sample metadata with SRR, GEO, SIALAP, Case_Control
NATCOM_MTX   = "/path/to/NatCom/merged_matrix_M.txt"  # NatCom mouse count matrix
IAKA_MTX     = BASE_DIR + "/02.Quantification/merged_matrix.txt"
IAKA_COLDATA = BASE_DIR + "/02.Quantification/coldata.txt"
OUT_DIR      = BASE_DIR + "/07.CrossDataset/01.Quantification"

import os; os.makedirs(OUT_DIR, exist_ok=True)

# Parse animal sheet for IHKA (GSE99577, GSE213393) and PILO samples
ihka_meta, pilo_meta = {}, {}
with open(ANIMAL_SHEET, newline="") as fh:
    for row in csv.DictReader(fh):
        srr    = row["SRR"].strip()
        geo    = row["GEO"].strip()
        cc     = row["Case_Control"].strip()
        sialap = row["SIALAP"].strip()
        phase  = sialap.split("-")[-1]
        if geo in ("GSE99577", "GSE213393"):
            ihka_meta[srr] = {"geo": geo, "cc": cc, "phase": phase}
        elif sialap.startswith(("M-PIL", "M-SAL")):
            pilo_meta[srr] = {"geo": geo, "cc": cc, "phase": phase}

# IAKA — write coldata and pass-through count matrix
day_phase = {"Day3": "AC", "Day7": "IM", "Day14": "CR"}
iaka_rows = []
with open(IAKA_COLDATA) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        srr, sample_id, group, day = p[0], p[1], p[2], p[4]
        is_case = not group.startswith("IASL")
        phase   = day_phase.get(day, day)
        iaka_rows.append((srr, "IAKA",
                          f"IAKA_{phase}" if is_case else f"IASL_{phase}",
                          phase, "GSE319769", "KA" if is_case else "Saline"))

with open(f"{OUT_DIR}/coldata_iaka.txt", "w") as fh:
    fh.write("SRR\tModel\tGroup\tPhase\tBatch\tTreatment\n")
    for r in iaka_rows: fh.write("\t".join(r) + "\n")

with open(IAKA_MTX) as fi, open(f"{OUT_DIR}/merged_matrix_iaka.txt", "w") as fo:
    fo.write("gene_name\t" + fi.readline().rstrip("\n") + "\n")
    fo.writelines(fi)

# IHKA / PILO — build coldata
def make_coldata(meta, model, case_grp_fn, ctrl_grp_fn, case_treat, ctrl_treat):
    rows = []
    for srr, m in meta.items():
        is_case = m["cc"] == "Case"
        rows.append((srr, model,
                     case_grp_fn(m["phase"]) if is_case else ctrl_grp_fn(m),
                     m["phase"], m["geo"],
                     case_treat if is_case else ctrl_treat))
    return rows

ihka_rows = make_coldata(ihka_meta, "IHKA",
                          lambda p: f"IHKA_{p}", lambda m: f"IHSL_{m['phase']}",
                          "KA", "Saline")
pilo_rows = make_coldata(pilo_meta, "PILO",
                          lambda p: f"PILO_{p}",
                          lambda m: "PSLC_GSE72402" if m["geo"] == "GSE72402" else f"PSLC_{m['phase']}",
                          "PILO", "Saline")

for model, rows in [("ihka", ihka_rows), ("pilo", pilo_rows)]:
    with open(f"{OUT_DIR}/coldata_{model}.txt", "w") as fh:
        fh.write("SRR\tModel\tGroup\tPhase\tBatch\tTreatment\n")
        for r in rows: fh.write("\t".join(r) + "\n")

# Extract IHKA and PILO columns from NatCom matrix
ihka_srrs = set(ihka_meta); pilo_srrs = set(pilo_meta)
with open(NATCOM_MTX) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    ihka_cols = [(i, h) for i, h in enumerate(header) if h in ihka_srrs]
    pilo_cols = [(i, h) for i, h in enumerate(header) if h in pilo_srrs]
    ihka_cols.sort(); pilo_cols.sort()
    ihka_idx = [i for i, _ in ihka_cols]; pilo_idx = [i for i, _ in pilo_cols]
    ihka_ord = [h for _, h in ihka_cols]; pilo_ord = [h for _, h in pilo_cols]

    fo_i = open(f"{OUT_DIR}/merged_matrix_ihka.txt", "w")
    fo_p = open(f"{OUT_DIR}/merged_matrix_pilo.txt", "w")
    fo_i.write("gene_name\t" + "\t".join(ihka_ord) + "\n")
    fo_p.write("gene_name\t" + "\t".join(pilo_ord) + "\n")
    for line in fh:
        p = line.rstrip("\n").split("\t"); g = p[0]
        fo_i.write(g + "\t" + "\t".join(p[i+1] for i in ihka_idx) + "\n")
        fo_p.write(g + "\t" + "\t".join(p[i+1] for i in pilo_idx) + "\n")
    fo_i.close(); fo_p.close()

# Reorder coldata to match matrix column order
def reorder(rows, srr_order, path):
    rmap = {r[0]: r for r in rows}
    with open(path, "w") as fh:
        fh.write("SRR\tModel\tGroup\tPhase\tBatch\tTreatment\n")
        for srr in srr_order:
            if srr in rmap: fh.write("\t".join(rmap[srr]) + "\n")

reorder(ihka_rows, ihka_ord, f"{OUT_DIR}/coldata_ihka.txt")
reorder(pilo_rows, pilo_ord, f"{OUT_DIR}/coldata_pilo.txt")
print(f"Done — output: {OUT_DIR}")
