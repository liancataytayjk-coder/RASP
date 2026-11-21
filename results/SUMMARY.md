# Papilio subgenus Achillides — dated COI phylogeny + ancestral ranges

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
Combined post-burn-in samples: 18,002 states per statistic; lowest ESS = 1275 (ucld.stdev); **no parameter has ESS < 200** (results/4_beast/ess_table.csv).
A first attempt (results/3_beast/attempt1_lowESS) had poor mixing of the GTR rates (ESS of gt, prior and posterior < 25 in one run) with default operator weights/priors; operator weights and rate priors were changed and the analysis re-run (reported here).

## Node ages (Ma; median and 95% HPD)
| clade | tips | median | 95% HPD | posterior |
|---|---|---|---|---|
| root (Papilio s.s. split) | 42 | 20.59 | 17.17–24.21 | 1 |
| Achillides+memnon | 41 | 18.52 | 15.69–21.46 | 1 |
| crown Achillides | 40 | 17.3 | 14.62–20.26 | 1 |
| P_angustatus_F, P_blumei_G, ... (9 tips) | 9 | 15.25 | 11.91–18.75 | 0.706 |
| P_arcturus_A, P_arcturus_C, ... (31 tips) | 31 | 15.02 | 11.11–18.81 | 0.668 |
| P_crino_C, P_lorquinianus_G, ... (7 tips) | 7 | 13.54 | 9.35–17.69 | 0.587 |
| P_arcturus_A, P_arcturus_C, ... (24 tips) | 24 | 11.63 | 8.11–15.3 | 0.982 |
| P_lorquinianus_G, P_lorquinianus_H, ... (6 tips) | 6 | 11.33 | 7.27–15.57 | 0.902 |
| P_bianor_A, P_bianor_B, ... (17 tips) | 17 | 10.12 | 6.91–13.6 | 0.361 |
| P_lorquinianus_G, P_lorquinianus_H, ... (5 tips) | 5 | 9.98 | 5.92–13.83 | 0.545 |
| P_montrouzieri_H, P_orsippus_H, ... (5 tips) | 5 | 9.97 | 6.4–13.97 | 1 |
| P_bianor_A, P_bianor_B, ... (11 tips) | 11 | 8.47 | 5.51–11.71 | 0.995 |
| P_arcturus_A, P_arcturus_C, ... (7 tips) | 7 | 8.18 | 4.89–12.22 | 0.609 |
| P_bianor_A, P_bianor_B, ... (9 tips) | 9 | 7.47 | 4.85–10.43 | 0.431 |
| P_bianor_A, P_bianor_B, ... (7 tips) | 7 | 6.58 | 3.89–9.43 | 0.379 |
| P_hermosanus_B, P_karna_D, ... (6 tips) | 6 | 5.99 | 2.94–9.74 | 0.992 |
| P_bianor_A, P_bianor_B, ... (6 tips) | 6 | 5.45 | 2.99–8.1 | 0.545 |
| P_bianor_A, P_bianor_B, ... (5 tips) | 5 | 4.09 | 2.19–6.57 | 0.667 |
| P_hermosanus_B, P_karna_D, ... (5 tips) | 5 | 3.81 | 1.94–6.25 | 0.81 |

Full table: results/4_beast/node_ages.csv.

## Biogeography (BioGeoBEARS, max range size 4)
| model | LnL | free params | d | e | j | AICc | AICc weight |
|---|---|---|---|---|---|---|---|
| DEC | -121.80 | 2 | 0.009857 | 0.007279 | 0 | 247.93 | 0.000 |
| DEC+J | -77.86 | 3 | 1e-12 | 1e-12 | 0.1887 | 162.38 | 0.286 |
| DIVALIKE | -107.66 | 2 | 0.01215 | 1e-12 | 0 | 219.64 | 0.000 |
| DIVALIKE+J | -76.99 | 3 | 1e-12 | 1e-12 | 0.1579 | 160.65 | 0.680 |
| BAYAREALIKE | -158.77 | 2 | 0.01408 | 0.1769 | 0 | 321.86 | 0.000 |
| BAYAREALIKE+J | -79.99 | 3 | 1e-12 | 1e-12 | 0.1809 | 166.65 | 0.034 |

**Best model: DIVALIKE+J (AICc weight 0.680).**
Ancestral area of crown Achillides (node 41, 17.3 Ma): **C (Indo-Himalaya & Indochina) + F (Oceanic Philippines) + G (Wallacea) + H (Sahul & Melanesia)**, probability 0.044
The single best range is very weakly supported; per-area marginal probabilities at the crown (summed over ranges containing the area): A 0.343; B 0.252; C 0.55; D 0.337; E 0.032; F 0.497; G 0.606; H 0.516. (Next-best ranges: C (Indo-Himalaya & Indochina) + G (Wallacea) + H (Sahul & Melanesia) 0.034; C (Indo-Himalaya & Indochina) + F (Oceanic Philippines) + G (Wallacea) 0.028.)

## Caveats / preliminary
* Single-locus mtDNA; Normal secondary calibrations only, so ages mostly reflect the priors where COI has little signal (several HPDs are wide).
* TIM2+G4 (best BIC) replaced by GTR+G for BEAST compatibility.
* An unconstrained ML tree (ModelTest-NG starting tree) nests *P. memnon*/*P. machaon* inside Achillides; the monophyly constraints from Condamine et al. 2023 force the outgroups outside.
* Many deeper nodes have posterior < 0.95 (see node_ages.csv); treat deep topology with caution.
* RASP 4 analyses are not run here (GUI); see LOCAL_STEPS.md. BioGeoBEARS run on the MCC tree only (not integrated over the posterior trees).
* +J models estimate d ≈ e ≈ 0 (founder-event-only dynamics); the AICc preference for +J and the resulting node probabilities should be read cautiously (cf. Ree & Sanmartín 2018). All models give a diffuse crown ancestral range.
* Tips of the same species from different regions are treated as separate single-area terminals.
