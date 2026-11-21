# Stage 7a: node ages (median, 95% HPD, posterior) from the MCC tree for all internal nodes
suppressMessages({library(treeio); library(ape)})
tb <- read.beast("results/4_beast/Achillides_MCC.tree"); tr <- tb@phylo; d <- as_tibble(tb)
n <- Ntip(tr); rows <- NULL
for (nd in (n+1):(n+tr$Nnode)) {
  r <- d[d$node == nd, ]; tips <- extract.clade(tr, nd)$tip.label
  h <- unlist(r[["height_0.95_HPD"]])
  rows <- rbind(rows, data.frame(node=nd, n_tips=length(tips), median_age=round(as.numeric(r$height_median),2),
     hpd_lo=round(min(h),2), hpd_hi=round(max(h),2), posterior=round(as.numeric(r$posterior),3),
     tips=paste(tips, collapse=";")))
}
rows$clade <- ifelse(rows$n_tips==42,"root (Papilio s.s. split)", ifelse(rows$n_tips==41,"Achillides+memnon",
              ifelse(rows$n_tips==40,"crown Achillides","")))
write.csv(rows[order(-rows$median_age),], "results/4_beast/node_ages.csv", row.names=FALSE)
print(rows[order(-rows$median_age), c("node","n_tips","median_age","hpd_lo","hpd_hi","posterior","clade")][1:14,])
