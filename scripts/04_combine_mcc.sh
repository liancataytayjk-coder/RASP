#!/usr/bin/env bash
# Stage 4: LogCombiner (10% burn-in = 1,000,000 states each), ESS via loganalyser, TreeAnnotator MCC (median heights),
# thinned 1,000-tree posterior sample.
set -euo pipefail
cd "$(dirname "$0")/.."
P=/opt/mm/root/envs/phylo/bin; O=results/4_beast; R=results/3_beast
$P/logcombiner -burnin 1000000 -log $R/run1/Achillides_run1.log -log $R/run2/Achillides_run2.log -o $O/combined.log 2>&1 | tail -2
$P/logcombiner -burnin 1000000 -trees $R/run1/Achillides_run1.trees -trees $R/run2/Achillides_run2.trees -o $O/combined_full.trees 2>&1 | tail -2
$P/logcombiner -burnin 1000000 -resample 18000 -trees $R/run1/Achillides_run1.trees -trees $R/run2/Achillides_run2.trees -o $O/Achillides_posterior_1000.trees 2>&1 | tail -2
# ESS (combined, burn-in already removed) and per-run
$P/loganalyser -burnin 0 $O/combined.log > $O/ess_combined.txt 2>/dev/null
for r in 1 2; do $P/loganalyser -burnin 1000000 $R/run$r/Achillides_run$r.log > $O/ess_run$r.txt 2>/dev/null; done
python3 scripts/04_ess.py
$P/treeannotator -burninTrees 0 -heights median -hpd 0.95 $O/combined_full.trees $O/Achillides_MCC.tree 2>&1 | tail -15
