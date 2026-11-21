#!/usr/bin/env bash
# Stage 0: install software with micromamba (conda-forge/bioconda). github.com / CRAN are blocked in this sandbox, conda is not.
set -euo pipefail
export MAMBA_ROOT_PREFIX=/opt/mm/root
mkdir -p /opt/mm && cd /opt/mm
if [ ! -x bin/micromamba ]; then
  curl -s -O https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.10.0-0.tar.bz2
  tar xjf micromamba-2.10.0-0.tar.bz2 bin/micromamba
fi
./bin/micromamba create -y -n phylo -c conda-forge -c bioconda beast=1.10.4 beagle-lib mafft muscle iqtree modeltest-ng openjdk=11
./bin/micromamba create -y -n fig -c conda-forge -c bioconda figtree
./bin/micromamba create -y -n renv -c conda-forge -c bioconda r-base=4.3 r-ape bioconductor-treeio bioconductor-ggtree r-optimx \
  r-cladorcpp r-rexpokit r-phylobase r-fd r-snow r-gensa r-multinomialci r-gtools r-stringr r-data.table r-ggplot2 r-biogeobears
