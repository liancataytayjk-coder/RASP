# Local steps (only what needs a desktop GUI)

Everything else (alignment, model test, BEAST ×2 runs, LogCombiner, TreeAnnotator, MCC tree, figures, BioGeoBEARS) was run in the sandbox;
re-run any stage with `scripts/` (see header comments). BEAST/TreeAnnotator therefore do **not** need to be repeated locally.

## RASP 4 — ancestral range reconstruction on the BEAST posterior

Inputs (all in `results/6_rasp/`):
| file | use |
|---|---|
| `rasp_trees.trees` | 1,000 post-burn-in posterior trees, outgroups pruned (NEXUS) |
| `rasp_MCC.tree` | MCC (median-height) tree, outgroups pruned (NEXUS); `rasp_MCC.tree.nwk` is the same in Newick |
| `rasp_distribution.csv` | `ID,species name,distribution` (no header row; 40 tips; letters A–H, none dropped) |
| `area_recode.csv` | old→new area letters (identity here: all 8 areas have tips) |

1. Start RASP 4 (Java; double-click `RASP.jar` or the app launcher).
2. **Load trees**: in the *Trees* tab (top-left panel) choose *Load Trees* → `results/6_rasp/rasp_trees.trees`. When asked, set **burn-in = 0** (burn-in was already removed) and confirm the taxon names load.
3. **Condensed (final) tree**: *Load Condensed Tree* / *Final tree* → `results/6_rasp/rasp_MCC.tree`.
4. **Distribution**: *Load Distribution* → `results/6_rasp/rasp_distribution.csv`. Check that the 40 names match the tree tips and the areas A–H are recognised (if RASP expects a tab-delimited file with a header, save a copy from Excel/LibreOffice as tab-delimited text with header `ID  Species  Distribution`).
5. *Reconstruct Ancestral State* → choose **S-DIVA**, **Statistical DEC (S-DEC)**, and/or **BBM**:
   - Max areas at each node: **4** (matches BioGeoBEARS).
   - Use the 1,000 trees as the sampling set and the MCC tree as the display tree.
   - BBM: 2 chains, 500,000+ generations, sampling every 100, root distribution *Null*, model *Jukes-Cantor (JC)* + *Equal*; discard the first 25% as burn-in.
6. Run; when finished open *Results* → pies/labels on the MCC tree. Export the tree figure (PDF/SVG) and the node table (*Save results*) into `results/6_rasp/` as `rasp_<method>_tree.pdf` / `rasp_<method>_nodes.txt`.
7. Report the node for **crown Achillides** (the root of the pruned MCC tree, node of all 40 tips) and compare it with `results/6_rasp/biogeobears_nodes.csv` (`is_crown_Achillides = TRUE`).

Area legend: A Mainland East Asia · B Taiwan & Ryukyus · C Indo-Himalaya & Indochina · D Sundaland · E Palawan · F Oceanic Philippines · G Wallacea · H Sahul & Melanesia.

## Only if you want to repeat BEAST locally (not required)
* Open `results/3_beast/run1/Achillides_run1.xml` in BEAST 1.10.4 (GUI: *Choose File* → tick BEAGLE → *Run*), or
  `beast -beagle_SSE -seed 2001 Achillides_run1.xml` (run2: seed 2002).
* LogCombiner: *Log files*, burn-in 1,000,000 states, add both `.log` / both `.trees` (for trees, resample 18,000 → 1,000 trees).
* TreeAnnotator: burn-in 0 (already removed), target *Maximum clade credibility tree*, node heights *Median heights*, input `combined_full.trees`.
