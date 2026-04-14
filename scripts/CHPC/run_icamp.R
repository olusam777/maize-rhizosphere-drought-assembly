## =============================================================================
##  run_icamp.R 
##  iCAMP βNTI + RC-bray on actual drought dataset
## =============================================================================

cat("=== Loading packages ===\n")
suppressPackageStartupMessages({
  library(phyloseq)
  library(iCAMP)
  library(ape)
  library(vegan)
})

cat("=== Loading data ===\n")
ps <- readRDS("ps_final_filtered.rds")
cat("Phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")

## ── Extract and subset to planted samples ───────────────────────────────────

meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
tree <- phy_tree(ps)

planted_idx <- meta$condition == "Planted"
ps_planted <- prune_taxa(taxa_sums(prune_samples(planted_idx, ps)) > 0,
                         prune_samples(planted_idx, ps))

comm_full <- as.matrix(otu_table(ps_planted))
if (taxa_are_rows(ps_planted)) comm_full <- t(comm_full)

meta_planted <- meta[planted_idx, , drop = FALSE]

cat("Planted samples:", nrow(comm_full), "\n")
cat("Total ASVs in planted samples:", ncol(comm_full), "\n")

## ── CRITICAL: Filter ASVs for iCAMP ────────────────────────────────────────
## 140K ASVs creates a 140K x 140K distance matrix (~148 GB) which overflows
## R integer limit. We filter to prevalent ASVs — standard practice for
## betaNTI. Rare singletons are phylogenetically uninformative for null
## models.

cat("\n=== Filtering ASVs for iCAMP ===\n")

## Prevalence: keep ASVs in >=5% of planted samples
n_planted <- nrow(comm_full)
min_prev <- ceiling(0.05 * n_planted)
prevalence <- colSums(comm_full > 0)
prev_pass <- prevalence >= min_prev
cat(sprintf("Prevalence filter (>=%d samples): %d / %d pass\n",
            min_prev, sum(prev_pass), ncol(comm_full)))

## Abundance: keep ASVs with >=50 total reads
total_counts <- colSums(comm_full)
abund_pass <- total_counts >= 50
cat(sprintf("Abundance filter (>=50 reads): %d / %d pass\n",
            sum(abund_pass), ncol(comm_full)))

## Combined
keep_asvs <- prev_pass & abund_pass
comm <- comm_full[, keep_asvs, drop = FALSE]
cat(sprintf("Combined: %d ASVs retained (%.1f%%)\n",
            ncol(comm), 100 * ncol(comm) / ncol(comm_full)))
cat(sprintf("Reads retained: %.1f%%\n",
            100 * sum(comm) / sum(comm_full)))

## Safety check: cap at 15,000 ASVs if still too large
n_taxa <- ncol(comm)
matrix_gb <- as.double(n_taxa)^2 * 8 / 1e9
cat(sprintf("Distance matrix size: %d x %d = %.1f GB\n", n_taxa, n_taxa, matrix_gb))

if (n_taxa > 15000) {
  cat("Still too many ASVs — taking top 15,000 by prevalence x abundance\n")
  score <- prevalence[keep_asvs] * log1p(total_counts[keep_asvs])
  top_n <- names(sort(score, decreasing = TRUE))[1:15000]
  comm <- comm[, top_n, drop = FALSE]
  n_taxa <- ncol(comm)
  matrix_gb <- as.double(n_taxa)^2 * 8 / 1e9
  cat(sprintf("Reduced to %d ASVs (%.1f GB matrix)\n", n_taxa, matrix_gb))
}

## ── Align tree ─────────────────────────────────────────────────────────────

shared_tips <- intersect(colnames(comm), tree$tip.label)
comm <- comm[, shared_tips, drop = FALSE]
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared_tips))

cat(sprintf("\nFinal: %d samples x %d ASVs\n", nrow(comm), ncol(comm)))
cat("Tree tips:", length(tree_pruned$tip.label), "\n")

## ── Phylogenetic distance matrix ───────────────────────────────────────────

cat("\n=== Computing phylogenetic distance matrix ===\n")
cat("Started:", format(Sys.time()), "\n")

pd <- cophenetic(tree_pruned)
cat("Done:", nrow(pd), "x", ncol(pd), "\n")
cat("Finished:", format(Sys.time()), "\n")

## ── Groups ─────────────────────────────────────────────────────────────────

groups <- data.frame(
  sample = rownames(comm),
  treatment = as.character(meta_planted[rownames(comm), "treatment"]),
  stringsAsFactors = FALSE
)
cat("\nSamples per treatment:\n")
print(table(groups$treatment))

## ── betaNTI ────────────────────────────────────────────────────────────────

cat("\n=== Computing betaNTI (expect several hours) ===\n")
cat("Started:", format(Sys.time()), "\n")

set.seed(42)
bnti <- tryCatch(
  iCAMP::bNTIn.p(
    comm = comm, dis = pd,
    nworker = 20, memo.size.GB = 90,
    weighted = TRUE, rand = 999,
    output.bMNTD = TRUE, sig.index = "SES"
  ),
  error = function(e) {
    cat("bNTIn.p error:", e$message, "\n")
    cat("Retrying with 499 randomizations...\n")
    tryCatch(
      iCAMP::bNTIn.p(
        comm = comm, dis = pd,
        nworker = 20, memo.size.GB = 90,
        weighted = TRUE, rand = 499,
        output.bMNTD = TRUE, sig.index = "SES"
      ),
      error = function(e2) { cat("Failed again:", e2$message, "\n"); NULL }
    )
  }
)

if (is.null(bnti)) {
  cat("FATAL: betaNTI failed\n")
  quit(save = "no", status = 1)
}
cat("betaNTI done:", format(Sys.time()), "\n")

## ── RC-bray ────────────────────────────────────────────────────────────────

cat("\n=== Computing RC-bray ===\n")
cat("Started:", format(Sys.time()), "\n")

set.seed(42)
rcbray <- tryCatch(
  iCAMP::RC.pc(
    comm = comm, rand = 999,
    nworker = 20, memo.size.GB = 90
  ),
  error = function(e) {
    cat("RC.pc error:", e$message, "\n")
    tryCatch(
      iCAMP::RC.pc(comm = comm, rand = 499, nworker = 20, memo.size.GB = 90),
      error = function(e2) { cat("Failed again:", e2$message, "\n"); NULL }
    )
  }
)

if (is.null(rcbray)) {
  cat("FATAL: RC-bray failed\n")
  quit(save = "no", status = 1)
}
cat("RC-bray done:", format(Sys.time()), "\n")

## ── Assemble process fractions ─────────────────────────────────────────────

cat("\n=== Classifying assembly processes ===\n")

snames <- rownames(comm)
n <- length(snames)

results_list <- vector("list", n * (n - 1) / 2)
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
      sample1 = s1, sample2 = s2,
      bNTI = bv, RC_bray = rv, Process = process,
      stringsAsFactors = FALSE
    )
  }
}

proc_df <- do.call(rbind, results_list)

## Add metadata
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
         "Within_Watered", "Between_treatments")
)

## ── Save everything ────────────────────────────────────────────────────────

write.csv(proc_df, "iCAMP_process_results_actual.csv", row.names = FALSE)
write.csv(as.data.frame(as.matrix(bnti$index)), "bNTI_matrix.csv")
write.csv(as.data.frame(as.matrix(rcbray$index)), "RC_bray_matrix.csv")

## Filtering report
write.csv(data.frame(
  Metric = c("Total ASVs (planted)", "ASVs after filtering",
             "Reads retained (%)", "Pairwise comparisons"),
  Value = c(ncol(comm_full), ncol(comm),
            round(100 * sum(comm) / sum(comm_full), 1), nrow(proc_df))
), "iCAMP_filtering_report.csv", row.names = FALSE)

## ── Summary ────────────────────────────────────────────────────────────────

cat("\n=== RESULTS ===\n")
cat("Pairs:", nrow(proc_df), "\n\n")

cat("Overall assembly fractions:\n")
print(round(prop.table(table(proc_df$Process)), 3))

det <- sum(proc_df$Process %in% c("Variable.Selection","Homogeneous.Selection"), na.rm=TRUE)
sto <- sum(proc_df$Process %in% c("Dispersal.Limitation","Homogenizing.Dispersal","Drift"), na.rm=TRUE)
cat(sprintf("\nDeterministic: %.1f%% | Stochastic: %.1f%%\n",
            100*det/(det+sto), 100*sto/(det+sto)))

for (ct in c("Within_Drought","Within_Watered","Between_treatments")) {
  sub <- proc_df[proc_df$comparison_type == ct, ]
  if (nrow(sub) > 0) {
    cat(sprintf("\n%s (n=%d):\n", ct, nrow(sub)))
    print(round(prop.table(table(sub$Process)), 3))
  }
}

bvals <- proc_df$bNTI[!is.na(proc_df$bNTI)]
cat(sprintf("\nbetaNTI: mean=%.3f, median=%.3f, SD=%.3f\n",
            mean(bvals), median(bvals), sd(bvals)))
cat(sprintf("|betaNTI|>1.96: %.1f%%  |betaNTI|<=1.96: %.1f%%\n",
            100*mean(abs(bvals)>1.96), 100*mean(abs(bvals)<=1.96)))

cat("\n=== Complete:", format(Sys.time()), "===\n")
writeLines(capture.output(sessionInfo()), "icamp_sessionInfo.txt")
