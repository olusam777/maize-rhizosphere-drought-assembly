cat("=== Loading packages ===\n")
suppressPackageStartupMessages({
  library(phyloseq)
  library(iCAMP)
  library(ape)
  library(vegan)
})

## Reload the same filtered data
ps <- readRDS("ps_final_filtered.rds")
meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
tree <- phy_tree(ps)

planted_idx <- meta$condition == "Planted"
ps_planted <- prune_taxa(taxa_sums(prune_samples(planted_idx, ps)) > 0,
                         prune_samples(planted_idx, ps))
comm_full <- as.matrix(otu_table(ps_planted))
if (taxa_are_rows(ps_planted)) comm_full <- t(comm_full)
meta_planted <- meta[planted_idx, , drop = FALSE]

## Same filtering as before
prevalence <- colSums(comm_full > 0)
total_counts <- colSums(comm_full)
keep <- prevalence >= 4 & total_counts >= 50
comm <- comm_full[, keep, drop = FALSE]

shared_tips <- intersect(colnames(comm), tree$tip.label)
comm <- comm[, shared_tips, drop = FALSE]

cat("Filtered:", ncol(comm), "ASVs,", nrow(comm), "samples\n")

## RC-bray WITHOUT memo.size.GB
cat("=== Computing RC-bray ===\n")
cat("Started:", format(Sys.time()), "\n")

set.seed(42)
rcbray <- iCAMP::RC.pc(comm = comm, rand = 999, nworker = 20)

cat("RC-bray done:", format(Sys.time()), "\n")

## Load βNTI from the successful run (it's in memory format, need to recompute)
## Actually easier: reload the pd and recompute βNTI too since it's fast
pd <- cophenetic(drop.tip(tree, setdiff(tree$tip.label, colnames(comm))))

cat("Recomputing βNTI for safety...\n")
set.seed(42)
bnti <- iCAMP::bNTIn.p(
  comm = comm, dis = pd,
  nworker = 20, weighted = TRUE, rand = 999,
  output.bMNTD = TRUE, sig.index = "SES"
)
cat("βNTI done:", format(Sys.time()), "\n")

## Assemble results
snames <- rownames(comm)
n <- length(snames)
results_list <- vector("list", n*(n-1)/2)
idx <- 0

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    idx <- idx + 1
    s1 <- snames[i]; s2 <- snames[j]
    bv <- bnti$index[s1, s2]
    rv <- rcbray$index[s1, s2]
    
    process <- if (!is.na(bv) && abs(bv) > 1.96) {
      if (bv > 1.96) "Variable.Selection" else "Homogeneous.Selection"
    } else if (!is.na(rv)) {
      if (rv > 0.95) "Dispersal.Limitation"
      else if (rv < -0.95) "Homogenizing.Dispersal"
      else "Drift"
    } else NA_character_
    
    results_list[[idx]] <- data.frame(
      sample1=s1, sample2=s2, bNTI=bv, RC_bray=rv, Process=process,
      stringsAsFactors=FALSE
    )
  }
}

proc_df <- do.call(rbind, results_list)
proc_df$treatment1 <- as.character(meta_planted[proc_df$sample1, "treatment"])
proc_df$treatment2 <- as.character(meta_planted[proc_df$sample2, "treatment"])
proc_df$trait1     <- as.character(meta_planted[proc_df$sample1, "trait"])
proc_df$trait2     <- as.character(meta_planted[proc_df$sample2, "trait"])
proc_df$harvest1   <- as.character(meta_planted[proc_df$sample1, "harvest"])
proc_df$harvest2   <- as.character(meta_planted[proc_df$sample2, "harvest"])

proc_df$comparison_type <- ifelse(
  proc_df$treatment1 == proc_df$treatment2 & proc_df$treatment1 == "Drought",
  "Within_Drought",
  ifelse(proc_df$treatment1 == proc_df$treatment2 & proc_df$treatment1 == "Watered",
         "Within_Watered", "Between_treatments"))

## Save
write.csv(proc_df, "iCAMP_process_results_actual.csv", row.names=FALSE)
write.csv(as.data.frame(as.matrix(bnti$index)), "bNTI_matrix.csv")
write.csv(as.data.frame(as.matrix(rcbray$index)), "RC_bray_matrix.csv")
write.csv(data.frame(
  Metric=c("ASVs used","Reads retained (%)","Pairwise comparisons"),
  Value=c(ncol(comm), round(100*sum(comm)/sum(comm_full),1), nrow(proc_df))
), "iCAMP_filtering_report.csv", row.names=FALSE)

cat("\n=== RESULTS ===\n")
cat("Pairs:", nrow(proc_df), "\n")
cat("\nOverall:\n")
print(round(prop.table(table(proc_df$Process)), 3))

det <- sum(proc_df$Process %in% c("Variable.Selection","Homogeneous.Selection"), na.rm=TRUE)
sto <- sum(proc_df$Process %in% c("Dispersal.Limitation","Homogenizing.Dispersal","Drift"), na.rm=TRUE)
cat(sprintf("\nDeterministic: %.1f%% | Stochastic: %.1f%%\n", 100*det/(det+sto), 100*sto/(det+sto)))

for (ct in c("Within_Drought","Within_Watered","Between_treatments")) {
  sub <- proc_df[proc_df$comparison_type == ct, ]
  if (nrow(sub) > 0) {
    cat(sprintf("\n%s (n=%d):\n", ct, nrow(sub)))
    print(round(prop.table(table(sub$Process)), 3))
  }
}

bvals <- proc_df$bNTI[!is.na(proc_df$bNTI)]
cat(sprintf("\nbetaNTI: mean=%.3f, median=%.3f\n", mean(bvals), median(bvals)))
cat(sprintf("|betaNTI|>1.96: %.1f%%\n", 100*mean(abs(bvals)>1.96)))

cat("\n=== Complete:", format(Sys.time()), "===\n")
