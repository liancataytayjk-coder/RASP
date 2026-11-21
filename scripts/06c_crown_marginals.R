# Per-area marginal probability (sum of state probabilities containing that area) at crown Achillides, best model
suppressMessages({library(BioGeoBEARS); library(ape)})
res <- readRDS("results/6_rasp/biogeobears_results.rds"); tab <- read.csv("results/6_rasp/biogeobears_models.csv")
best <- tab$model[which.min(tab$AICc)]; pr <- res[[best]]$ML_marginal_prob_each_state_at_branch_top_AT_node
areas <- LETTERS[1:8]; sl <- cladoRcpp::rcpp_areas_list_to_states_list(areas=areas, maxareas=4, include_null_range=TRUE)
root <- pr[41,]
m <- sapply(seq_along(areas), function(a) sum(root[sapply(sl, function(s) (a-1) %in% s)]))
out <- data.frame(model=best, area=areas, area_name=read.csv("results/6_rasp/area_recode.csv")$area_name, marginal_prob=round(m,3))
write.csv(out, "results/6_rasp/crown_area_marginals.csv", row.names=FALSE); print(out)
