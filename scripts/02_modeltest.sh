#!/usr/bin/env bash
# Stage 2: ModelTest-NG (MEGA-CC unavailable), BIC; choose best model among those BEAST supports (JC/HKY/TN93/GTR +I/+G)
set -euo pipefail
P=/opt/mm/root/envs/phylo/bin
cd "$(dirname "$0")/.."
$P/modeltest-ng -i results/1_alignment/COI_aligned.fasta -d nt -t ml -p 4 -o results/2_model/mt -r 1 -c 2>&1 >/dev/null | tail -3 || true
