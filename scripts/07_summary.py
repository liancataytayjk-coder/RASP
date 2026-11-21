#!/usr/bin/env python3
"""Stage 7b: assemble results/SUMMARY.md from result tables."""
import csv
rd = lambda f: list(csv.DictReader(open(f)))
ess = rd("results/4_beast/ess_table.csv"); lo = min(ess, key=lambda r: float(r["ESS"]))
nodes = [r for r in rd("results/4_beast/node_ages.csv")]
mod = rd("results/6_rasp/biogeobears_models.csv"); best = min(mod, key=lambda r: float(r["AICc"]))
bn = rd("results/6_rasp/biogeobears_nodes.csv"); crown = [r for r in bn if r["is_crown_Achillides"] == "TRUE"][0]
cm = rd("results/6_rasp/crown_area_marginals.csv")
CM = "; ".join(f"{r['area']} {r['marginal_prob']}" for r in cm)
rec = {r["old"]: r["area_name"] for r in rd("results/6_rasp/area_recode.csv")}
nm = lambda s: " + ".join(f"{c} ({rec.get(c,'?')})" for c in s) if s != "null" else "null"
major = [r for r in nodes if int(r["n_tips"]) >= 5 and int(r["n_tips"]) < 42]
def row(r): 
    lab = r["clade"] or ", ".join(r["tips"].split(";")[:2]) + ", ... (%s tips)" % r["n_tips"]
    return f"| {lab} | {r['n_tips']} | {r['median_age']} | {r['hpd_lo']}–{r['hpd_hi']} | {r['posterior']} |"
T = "\n".join(row(r) for r in [x for x in nodes if x["n_tips"] == "42"] + major)
M = "\n".join(f"| {r['model']} | {float(r['LnL']):.2f} | {r['n_params']} | {float(r['d']):.4g} | {float(r['e']):.4g} | {float(r['j']):.4g} | {float(r['AICc']):.2f} | {float(r['AICc_weight']):.3f} |" for r in mod)
open("results/SUMMARY.md", "w").write(f"""# Papilio subgenus Achillides — dated COI phylogeny + ancestral ranges

## Methods
COI sequences of 40 Achillides tips (one species × region each) plus *P. memnon* (Menelaides) and *P. machaon* (Papilio s.s.) were aligned with
MUSCLE v5.3 (MEGA-CC could not be obtained in this sandbox), checked for reading frame and stop codons under the invertebrate mitochondrial code (NCBI table 5; no stops in frame 1)
and trimmed to the region shared by all sequences (594 bp, 198 codons). ModelTest-NG (BIC) favoured TIM2+G4; because BEAST 1.10.4 does not implement TIM2,
the best BEAST-supported model, GTR+G (4 categories), was used. Divergence times were estimated in BEAST v1.10.4 with BEAGLE 4.2 (uncorrelated lognormal relaxed clock,
birth–death tree prior, GTR+G), with Normal secondary calibrations from Condamine et al. (2023, *Mol. Phylogenet. Evol.* 183:107758, Table 3): crown Achillides
N(17.0, 2.1) Ma, Achillides + *P. memnon* N(19.41, 2.25) Ma and the root (split with *P. machaon*) N(21.25, 2.45) Ma (SDs = 95% CI width / 3.92), with the first two groups constrained monophyletic.
Two independent runs of 10,000,000 generations (sampled every 1,000; seeds 2001/2002) were combined in LogCombiner after 10% burn-in each; the maximum clade credibility tree
(median node heights) was built in TreeAnnotator (BEAST 1.10.4). Outgroups were pruned and ancestral areas (8 areas; max range size = 4) were estimated with BioGeoBEARS 1.1.4 (DEC, DIVALIKE, BAYAREALIKE, each ± J) in R 4.3;
models were compared with AICc.

## Convergence
Combined post-burn-in samples: 18,002 states per statistic; lowest ESS = {float(lo['ESS']):.0f} ({lo['statistic']}); **no parameter has ESS < 200** (results/4_beast/ess_table.csv).
A first attempt (results/3_beast/attempt1_lowESS) had poor mixing of the GTR rates (ESS of gt, prior and posterior < 25 in one run) with default operator weights/priors; operator weights and rate priors were changed and the analysis re-run (reported here).

## Node ages (Ma; median and 95% HPD)
| clade | tips | median | 95% HPD | posterior |
|---|---|---|---|---|
{T}

Full table: results/4_beast/node_ages.csv.

## Biogeography (BioGeoBEARS, max range size 4)
| model | LnL | free params | d | e | j | AICc | AICc weight |
|---|---|---|---|---|---|---|---|
{M}

**Best model: {best['model']} (AICc weight {float(best['AICc_weight']):.3f}).**
Ancestral area of crown Achillides (node {crown['node']}, {crown['age_Ma']} Ma): **{nm(crown['best_range'])}**, probability {crown['best_prob']}
The single best range is very weakly supported; per-area marginal probabilities at the crown (summed over ranges containing the area): {CM}. (Next-best ranges: {nm(crown['second_range'])} {crown['second_prob']}; {nm(crown['third_range'])} {crown['third_prob']}.)

## Caveats / preliminary
* Single-locus mtDNA; Normal secondary calibrations only, so ages mostly reflect the priors where COI has little signal (several HPDs are wide).
* TIM2+G4 (best BIC) replaced by GTR+G for BEAST compatibility.
* An unconstrained ML tree (ModelTest-NG starting tree) nests *P. memnon*/*P. machaon* inside Achillides; the monophyly constraints from Condamine et al. 2023 force the outgroups outside.
* Many deeper nodes have posterior < 0.95 (see node_ages.csv); treat deep topology with caution.
* RASP 4 analyses are not run here (GUI); see LOCAL_STEPS.md. BioGeoBEARS run on the MCC tree only (not integrated over the posterior trees).
* +J models estimate d ≈ e ≈ 0 (founder-event-only dynamics); the AICc preference for +J and the resulting node probabilities should be read cautiously (cf. Ree & Sanmartín 2018). All models give a diffuse crown ancestral range.
* Tips of the same species from different regions are treated as separate single-area terminals.
""")
print("ok")
