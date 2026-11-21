# Stage 6a: prune outgroups; area recode; RASP inputs; BioGeoBEARS geography file
suppressMessages({library(ape)})
args <- commandArgs(TRUE)
mcc <- ifelse(length(args)>0, args[1], "results/4_beast/Achillides_MCC.tree")
post <- ifelse(length(args)>1, args[2], "results/4_beast/Achillides_posterior_1000.trees")
od <- ifelse(length(args)>2, args[3], "results/6_rasp"); dir.create(od, showWarnings=FALSE, recursive=TRUE)
out <- c("P_memnon","P_machaon")
td <- read.csv("results/tip_names.csv"); td <- td[td$area != "NA" & !td$tip_name %in% out, ]
# --- area recode (A-H -> continuous letters in use)
used <- sort(unique(td$area)); new <- LETTERS[seq_along(used)]
rec <- data.frame(old=used, new=new)
areanames <- c(A="Mainland East Asia",B="Taiwan & Ryukyus",C="Indo-Himalaya & Indochina",D="Sundaland",E="Palawan",F="Oceanic Philippines",G="Wallacea",H="Sahul & Melanesia")
rec$area_name <- areanames[rec$old]; rec$n_tips <- as.integer(table(td$area)[rec$old])
write.csv(rec, file.path(od,"area_recode.csv"), row.names=FALSE)
td$new_area <- rec$new[match(td$area, rec$old)]
message("areas without tips dropped: ", paste(setdiff(names(areanames), used), collapse=","), " (none if empty)")
# --- trees. MCC read as plain tree (annotations stripped); posterior sample
mt <- read.nexus(mcc); if (inherits(mt,"multiPhylo")) mt <- mt[[1]]
mt <- drop.tip(mt, out); stopifnot(all(mt$tip.label %in% td$tip_name), Ntip(mt)==40)
write.tree(mt, file.path(od,"rasp_MCC.tree.nwk"))
write.nexus(mt, file=file.path(od,"rasp_MCC.tree"))
pt <- read.nexus(post)
pt <- lapply(pt, drop.tip, tip=out); class(pt) <- "multiPhylo"
write.nexus(pt, file=file.path(od,"rasp_trees.trees"))
message(length(pt), " pruned posterior trees written")
# --- RASP distribution (ID, species name, distribution) ; no header
dd <- data.frame(ID=seq_len(nrow(td)), name=td$tip_name, dist=td$new_area)
dd <- dd[match(mt$tip.label, dd$name), ]; dd$ID <- seq_len(nrow(dd))
write.table(dd, file.path(od,"rasp_distribution.csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
# --- BioGeoBEARS geography file (PHYLIP-like)
k <- length(new)
g <- file.path(od,"geog_biogeobears.data")
cat(nrow(td), k, paste0("(", paste(new, collapse=" "), ")"), "\n", file=g)
for (i in seq_len(nrow(td))) { v <- rep(0, k); v[match(td$new_area[i], new)] <- 1; cat(td$tip_name[i], "\t", paste(v, collapse=""), "\n", file=g, append=TRUE, sep="") }
