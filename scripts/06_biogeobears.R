# Stage 6b: BioGeoBEARS DEC/DIVALIKE/BAYAREALIKE (+/-J), max_range_size = 4, on pruned MCC tree
suppressMessages({library(ape); library(BioGeoBEARS)})
od <- Sys.getenv("OD","results/6_rasp")
tr <- read.tree(file.path(od,"rasp_MCC.tree.nwk"))
# MCC (median heights) can have tiny negative / non-ultrametric branches: repair
tr$edge.length[tr$edge.length < 1e-4] <- 1e-4
dd <- node.depth.edgelength(tr)[1:Ntip(tr)]
if (diff(range(dd)) > 1e-6) { message("non-ultrametric (range ", signif(diff(range(dd)),3), "): extending tips"); 
  tipe <- match(seq_len(Ntip(tr)), tr$edge[,2]); tr$edge.length[tipe] <- tr$edge.length[tipe] + (max(dd)-dd) }
write.tree(tr, file.path(od,"biogeobears_tree.nwk"))
geogfn <- file.path(od,"geog_biogeobears.data")
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max_range_size <- 4
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- file.path(od,"biogeobears_tree.nwk")
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$on_NaN_error <- -1e50
BioGeoBEARS_run_object$speedup <- TRUE
BioGeoBEARS_run_object$use_optimx <- TRUE   # optimx L-BFGS-B (GenSA was far too slow)
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
setp <- function(obj, nm, type=NULL, init=NULL, est=NULL) {
  pt <- obj$BioGeoBEARS_model_object@params_table
  if (!is.null(type)) pt[nm,"type"] <- type
  if (!is.null(init)) pt[nm,"init"] <- init
  if (!is.null(est))  pt[nm,"est"]  <- est
  obj$BioGeoBEARS_model_object@params_table <- pt; obj }
mk <- function(base, J) {
  o <- BioGeoBEARS_run_object
  if (base == "DIVALIKE") {
    o <- setp(o,"s","fixed",0,0); o <- setp(o,"ysv","2-j"); o <- setp(o,"ys","ysv*1/2"); o <- setp(o,"y","1-j"); o <- setp(o,"v","ysv*1/2")
    o <- setp(o,"mx01v","fixed",0.5,0.5)
  }
  if (base == "BAYAREALIKE") {
    o <- setp(o,"s","fixed",0,0); o <- setp(o,"v","fixed",0,0); o <- setp(o,"ysv","1-j"); o <- setp(o,"ys","ysv*1/1"); o <- setp(o,"y","1-j")
    o <- setp(o,"mx01y","fixed",0.9999,0.9999)
  }
  if (J) o <- setp(o,"j","free",0.01,0.01)
  o <- fix_BioGeoBEARS_params_minmax(o); o }
models <- list(DEC=c("DEC",FALSE), "DEC+J"=c("DEC",TRUE), DIVALIKE=c("DIVALIKE",FALSE), "DIVALIKE+J"=c("DIVALIKE",TRUE),
               BAYAREALIKE=c("BAYAREALIKE",FALSE), "BAYAREALIKE+J"=c("BAYAREALIKE",TRUE))
library(parallel)
fit1 <- function(m) { o <- mk(models[[m]][1], models[[m]][2]=="TRUE"); bears_optim_run(o) }
ncores <- as.integer(Sys.getenv("NCORES","6"))
res <- mclapply(names(models), fit1, mc.cores=ncores, mc.preschedule=FALSE); names(res) <- names(models)
stopifnot(!any(sapply(res, function(x) inherits(x,"try-error"))))
tab <- NULL
for (m in names(models)) {
  r <- res[[m]]; p <- r$outputs@params_table
  tab <- rbind(tab, data.frame(model=m, LnL=get_LnL_from_BioGeoBEARS_results(r)$ML_value_of_LnL,
        n_params=sum(p$type=="free"), d=p["d","est"], e=p["e","est"], j=p["j","est"]))
}
tab$AICc <- with(tab, { n <- Ntip(tr); 2*n_params - 2*LnL + (2*n_params*(n_params+1))/(n-n_params-1) })
tab$delta_AICc <- tab$AICc - min(tab$AICc)
tab$AICc_weight <- exp(-0.5*tab$delta_AICc)/sum(exp(-0.5*tab$delta_AICc))
write.csv(tab, file.path(od,"biogeobears_models.csv"), row.names=FALSE)
print(tab)
best <- tab$model[which.min(tab$AICc)]; message("best model: ", best)
saveRDS(res, file.path(od,"biogeobears_results.rds"))
rb <- res[[best]]
# --- pie chart tree
areas <- getareas_from_tipranges_object(tipranges)
states <- rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
states_list <- states  # list of 0-based area indices
statenames <- sapply(states_list, function(s) if (length(s)==0 || all(is.na(s))) "null" else paste(areas[s+1], collapse=""))
pdf(file.path(od,"biogeobears_best.pdf"), width=9, height=11)
plot_BioGeoBEARS_results(rb, analysis_titletxt=paste0(best," (max range size 4)"), addl_params=list("j"), plotwhat="pie",
   label.offset=0.5, tipcex=0.7, statecex=0.5, splitcex=0.5, titlecex=0.9, plotsplits=FALSE,
   include_null_range=TRUE, tr=tr, tipranges=tipranges)
dev.off()
# --- node table: most probable state + probability (marginal ML at node = "ML_marginal_prob_each_state_at_branch_top_AT_node")
pr <- rb$ML_marginal_prob_each_state_at_branch_top_AT_node
n <- Ntip(tr); ages <- max(node.depth.edgelength(tr)) - node.depth.edgelength(tr)
rows <- NULL
for (nd in (n+1):(n+tr$Nnode)) {
  p <- pr[nd,]; o <- order(p, decreasing=TRUE)[1:3]
  tips <- extract.clade(tr, nd)$tip.label
  rows <- rbind(rows, data.frame(node=nd, age_Ma=round(ages[nd],2), n_tips=length(tips),
     clade=paste(head(tips,3), collapse="|"), best_range=statenames[o[1]], best_prob=round(p[o[1]],3),
     second_range=statenames[o[2]], second_prob=round(p[o[2]],3), third_range=statenames[o[3]], third_prob=round(p[o[3]],3)))
}
rows$is_crown_Achillides <- rows$node == n+1
write.csv(rows, file.path(od,"biogeobears_nodes.csv"), row.names=FALSE)
cat("crown Achillides (root):\n"); print(rows[rows$node==n+1,])
