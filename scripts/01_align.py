#!/usr/bin/env python3
"""Stage 1: rename tips, MUSCLE v5 alignment, reading-frame/stop check (NCBI table 5), trim to shared region."""
import csv, subprocess, os, sys, collections
MUSCLE="/opt/mm/root/envs/phylo/bin/muscle"
os.makedirs("results/1_alignment", exist_ok=True)
def read(f):
    return "".join(l.strip() for l in open(f) if not l.startswith(">")).upper().replace("U","T")
tips = list(csv.DictReader(open("results/tip_names.csv")))
seqs = {}
for r in tips:
    d = "data/outgroup/" if r["area"]=="NA" else "data/ingroup/"
    seqs[r["tip_name"]] = read(d+r["file"])
with open("results/1_alignment/COI_unaligned.fasta","w") as o:
    for k,v in seqs.items(): o.write(f">{k}\n{v}\n")
subprocess.run([MUSCLE,"-align","results/1_alignment/COI_unaligned.fasta","-output","results/1_alignment/COI_aligned_raw.fasta"],check=True,capture_output=True)
al = {}; name=None
for l in open("results/1_alignment/COI_aligned_raw.fasta"):
    l=l.strip()
    if l.startswith(">"): name=l[1:]; al[name]=""
    else: al[name]+=l.upper()
# keep input order
al = {k: al[k] for k in seqs}
L = len(next(iter(al.values())))
# table 5 (invertebrate mito) stop codons: TAA, TAG only (TGA=Trp, AGA/AGG=Ser, ATA=Met)
STOP={"TAA","TAG"}
def stops(seq, frame):
    return sum(seq[i:i+3] in STOP for i in range(frame,len(seq)-2,3))
# frame in alignment coordinates: choose offset minimising total stops over ungapped codons
def stops_aln(f):
    tot=0
    for s in al.values():
        for i in range(f,L-2,3):
            c=s[i:i+3]
            if "-" not in c and c in STOP: tot+=1
    return tot
res={f:stops_aln(f) for f in range(3)}
print("stop codons per alignment frame offset (raw alignment):",res)
frame=min(res,key=res.get)
# occupancy per column
nseq=len(al)
occ=[sum(s[i]!="-" for s in al.values()) for i in range(L)]
# shared region = columns where every sequence has data -> trim ragged ends only (first/last full-occupancy col)
full=[i for i in range(L) if occ[i]==nseq]
start,end=full[0],full[-1]
# snap start to codon boundary of the chosen frame, end to codon end
while (start-frame)%3: start+=1
while (end+1-frame)%3: end-=1
print("raw alignment length",L,"; shared region",start,end+1)
trim={k:v[start:end+1] for k,v in al.items()}
T=len(next(iter(trim.values())))
# drop columns all-gap (none expected) and report internal gaps
miss={k:100*sum(c in "-N?" for c in v)/T for k,v in trim.items()}
st_final={k:stops(v.replace("-","N"),0) for k,v in trim.items()}
print("final length",T,"(codons:",T/3,")")
print("stops after trim (frame 1, table 5):",{k:x for k,x in st_final.items() if x})
print("max % missing:",max(miss.values()))
with open("results/1_alignment/COI_aligned.fasta","w") as o:
    for k,v in trim.items(): o.write(f">{k}\n{v}\n")
with open("results/1_alignment/COI_aligned.nex","w") as o:
    o.write(f"#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX={len(trim)} NCHAR={T};\n  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n  MATRIX\n")
    for k,v in trim.items(): o.write(f"  {k} {v}\n")
    o.write("  ;\nEND;\n")
with open("results/1_alignment/alignment_report.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(["tip","raw_length","pct_missing_final","stop_codons_table5"])
    for k in trim: w.writerow([k,len(seqs[k]),round(miss[k],2),st_final[k]])
open("results/1_alignment/alignment_summary.txt","w").write(f"Tool: MUSCLE v5.3 (MEGA-CC unavailable: megasoftware.net blocked)\nraw length {L}; frame offset {frame} (stops per offset {res}); trimmed region {start+1}-{end+1}; final length {T}; max %missing {max(miss.values()):.2f}\n")
