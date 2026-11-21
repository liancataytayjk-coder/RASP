#!/usr/bin/env python3
"""Parse ModelTest-NG BIC table, keep BEAST-supported families (JC, HKY, TrN=TN93, GTR; +I,+G), pick lowest BIC."""
import re, csv
t = open("results/2_model/mt.out").read()
rows = re.findall(r"^\s+\d+\s+(\S+)\s+(\d+)\s+(-[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$", t, re.M)
sup = re.compile(r"^(JC|HKY|TrN|GTR)(\+I)?(\+G4)?(\+I\+G4)?$")
# note: ModelTest-NG's JC/HKY/TrN/GTR names; BEAST names: JC69/HKY/TN93/GTR
out = []
for m, k, lnl, bic, delta, w in rows:
    out.append([m, k, lnl, bic, delta, w, "yes" if sup.match(m) else "no"])
# the same rows may be listed in AIC/AICc/BIC tables; keep unique (model,bic) pairs from the BIC table (last block with BIC column ordering)
seen = {}
for r in out: seen.setdefault(r[0], r)  # first block in mt.out = BIC-sorted table
with open("results/2_model/model_table.csv", "w", newline="") as o:
    w = csv.writer(o); w.writerow(["model", "K", "lnL", "BIC", "delta", "weight", "beast_supported"])
    for r in sorted(seen.values(), key=lambda r: float(r[3])): w.writerow(r)
best = min((r for r in seen.values() if r[6] == "yes"), key=lambda r: float(r[3]))
open("results/2_model/chosen_model.txt", "w").write(f"Overall best BIC: TIM2+G4 (not implemented in BEAST 1.10.4)\nChosen (lowest-BIC BEAST-supported): {best[0]}  BIC={best[3]} lnL={best[2]}\n")
print(best)
