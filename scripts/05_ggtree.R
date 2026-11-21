# Stage 5: ggtree chronogram: Ma axis, median node ages, 95% HPD bars, posterior >= 0.95 marked
args <- commandArgs(TRUE); mcc <- ifelse(length(args)>0, args[1], "results/4_beast/Achillides_MCC.tree")
out <- ifelse(length(args)>1, args[2], "results/5_figure/chronogram_ggtree.pdf")
suppressMessages({library(treeio); library(ggtree); library(ggplot2); library(ape)})
tb <- read.beast(mcc); tr <- tb@phylo
d <- ggtree(tb)$data
d$h  <- as.numeric(d$height_median); if (is.null(d$h) || all(is.na(d$h))) d$h <- as.numeric(d$height)
hp <- d[["height_0.95_HPD"]]
d$lo <- sapply(hp, function(v) if (is.null(v)||all(is.na(v))) NA else min(as.numeric(v)))
d$hi <- sapply(hp, function(v) if (is.null(v)||all(is.na(v))) NA else max(as.numeric(v)))
# ggtree x = distance from root; age of a node = tree depth - x
depth <- max(d$x); d$age <- depth - d$x
d$post <- as.numeric(d$posterior)
age_scale <- depth                                  # x position of age 0 (tips)
bars <- d[!d$isTip & !is.na(d$lo), ]
# HPD in ages (Ma) -> x positions
bars$xl <- depth - bars$hi; bars$xr <- depth - bars$lo
ints <- d[!d$isTip, ]
p <- ggtree(tb, size=0.45) +
  geom_segment(data=bars, aes(x=xl, xend=xr, y=y, yend=y), colour="#2c7fb8", linewidth=2.2, alpha=0.45, lineend="round") +
  geom_tiplab(size=2.6, offset=0.15, fontface="italic") +
  geom_text(data=ints, aes(x=x, y=y, label=sprintf("%.1f", age)), size=2.2, nudge_y=0.45, hjust=1.05, colour="grey20") +
  geom_point(data=ints[!is.na(ints$post) & ints$post >= 0.95, ], aes(x=x, y=y), shape=21, fill="black", colour="white", size=2) +
  scale_x_continuous(breaks=depth - seq(0, 30, 5), labels=seq(0, 30, 5), limits=c(0, depth*1.35)) +
  labs(x="Time (Ma)", caption="Node labels: median age (Ma); bars: 95% HPD; filled circle: posterior probability >= 0.95") +
  theme_tree2() + theme(axis.line.x=element_line(), plot.caption=element_text(size=7))
n <- Ntip(tr); ggsave(out, p, width=8.5, height=max(6, 0.22*n+1.5))
cat("wrote", out, "\n")
