#!/usr/bin/env bash
# Stage 3: test run (100k) then 2 independent full runs (10M, log every 1000) in parallel, BEAGLE CPU/SSE.
set -euo pipefail
cd "$(dirname "$0")/.."
P=/opt/mm/root/envs/phylo/bin
python3 scripts/03_make_xml.py results/3_beast/test/test.xml 100000 1000 test
(cd results/3_beast/test && $P/beast -overwrite -seed 11 -beagle_CPU -beagle_SSE -threads 2 test.xml > test.out 2>&1)
for r in 1 2; do
  mkdir -p results/3_beast/run$r
  python3 scripts/03_make_xml.py results/3_beast/run$r/Achillides_run$r.xml 10000000 1000 Achillides_run$r
  (cd results/3_beast/run$r && nohup $P/beast -overwrite -seed $((1000+r)) -beagle_CPU -beagle_SSE -threads 2 Achillides_run$r.xml > run$r.out 2>&1 &)
done
