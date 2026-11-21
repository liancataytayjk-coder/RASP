#!/usr/bin/env bash
# Stage 5: FigTree CLI PDF + ggtree PDF
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=/opt/mm/root/envs/renv/bin:$PATH
FT=/opt/mm/root/envs/fig/bin/figtree
MCC=${1:-results/4_beast/Achillides_MCC.tree}; OUT=${2:-results/5_figure}
mkdir -p $OUT
# FigTree needs a nexus with FigTree block for labels; add display settings via a block appended to a copy
cp $MCC $OUT/_figtree_in.tree
cat >> $OUT/_figtree_in.tree <<'EOT'
begin figtree;
	set appearance.backgroundColorAttribute="Default";
	set layout.expansion=0;
	set layout.spread=0;
	set nodeLabels.displayAttribute="height_median";
	set nodeLabels.fontSize=8;
	set nodeLabels.isShown=true;
	set nodeLabels.significantDigits=1;
	set nodeBars.barWidth=4.0;
	set nodeBars.displayAttribute="height_95%_HPD";
	set nodeBars.isShown=true;
	set nodeBars.significantDigits=2;
	set scaleAxis.isShown=false;
	set scaleBar.isShown=false;
	set timeScale.isShown=true;
	set rootAge.isShown=false;
	set tipLabels.fontSize=8;
	set trees.order=true;
	set trees.orderType="increasing";
	set trees.rooting=false;
	set trees.scaleFactor=1.0;
end;
EOT
if xvfb-run -a true 2>/dev/null; then RUN="xvfb-run -a"; else RUN=""; fi
$RUN $FT -graphic PDF -width 900 -height 1100 $OUT/_figtree_in.tree $OUT/chronogram_figtree.pdf || echo "FigTree CLI failed"
rm -f $OUT/_figtree_in.tree
Rscript scripts/05_ggtree.R $MCC $OUT/chronogram_ggtree.pdf
