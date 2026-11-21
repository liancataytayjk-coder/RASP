#!/usr/bin/env python3
"""Parse loganalyser output, write ess_table.csv and list parameters with ESS < 200."""
import csv
rows=[]
for l in open("results/4_beast/ess_combined.txt"):
    p=l.split()
    if len(p)>=7 and not l.startswith(("#","statistic","Burn")):
        try: float(p[1]); rows.append(p)
        except ValueError: pass
hdr=["statistic","mean","stderr_mean","median","hpd_lo","hpd_hi","ESS"]
with open("results/4_beast/ess_table.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(hdr); w.writerows([r[:7] for r in rows])
low=[(r[0],float(r[6])) for r in rows if float(r[6])<200]
print("parameters:",len(rows),"; ESS<200:",low if low else "none")
open("results/4_beast/ess_low.txt","w").write(("\n".join(f"{a}\t{b}" for a,b in low)) or "none (all ESS >= 200)\n")
