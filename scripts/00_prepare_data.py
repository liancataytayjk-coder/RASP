#!/usr/bin/env python3
"""Stage 0: tip names + area table. distribution.csv was not supplied -> areas inferred from the trailing
region letter in each file name (P_species_REGION), as in the uploaded folder."""
import csv, glob, os, re
os.makedirs("results", exist_ok=True)
valid = set("ABCDEFGH")
rows = []
for f in sorted(glob.glob("data/ingroup/*.fasta")):
    hdr = open(f).readline().strip().lstrip(">")
    stem = os.path.basename(f)[:-6]
    m = re.fullmatch(r"P_([A-Za-z]+)_([A-H])", stem)
    assert m, stem
    rows.append((os.path.basename(f), hdr, stem, m.group(2)))
for f, tip in (("P_memnon.fasta", "P_memnon"), ("P_machaon.fasta", "P_machaon")):
    hdr = open("data/outgroup/" + f).readline().strip().lstrip(">")
    rows.append((f, hdr, tip, "NA"))
assert all(re.fullmatch(r"[A-Za-z0-9_]+", r[2]) for r in rows) and len({r[2] for r in rows}) == 42
with open("results/tip_names.csv", "w", newline="") as o:
    w = csv.writer(o); w.writerow(["file", "original_header", "tip_name", "area"]); w.writerows(rows)
with open("distribution.csv", "w", newline="") as o:   # inferred (file name suffix)
    w = csv.writer(o); w.writerow(["file", "area"])
    for r in rows[:40]: w.writerow([r[0], r[3]])
print(len(rows), "tips")
