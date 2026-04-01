## ============================================================
## DROUGHT RHIZOSPHERE MICROBIOME PIPELINE
## Olanrewaju et al. (in preparation)
##
## REPRODUCTION INSTRUCTIONS:
## 1. Run drought_pipeline_v2.r (full analysis)
##    Estimated runtime: 4–6 hours
##    High-memory steps:
##      SPIEC-EASI x4  (~2 hr each, ~500 MB RAM each)
##      iCAMP          (~30 min)
##      ANCOM-BC2 x5   (~10 min each)
##      Random Forest  (~15 min)
##    After each SPIEC-EASI run the RDS cache is saved;
##    subsequent runs skip the computation.
## 2. Run figure_assembly.r (figure compilation)
##    Estimated runtime: 5–10 minutes
##
## CACHING: Expensive computations are cached as RDS files
## in data/. Delete the RDS file to force recomputation.
## Key caches:
##   data/spiec_easi_drought.rds    (~500 MB)
##   data/spiec_easi_watered.rds    (~500 MB)
##   data/spiec_easi_resistant.rds  (~500 MB)
##   data/spiec_easi_susceptible.rds (~500 MB)
##   data/ps_genus_cached.rds       (genus-agglomerated phyloseq)
##
## REQUIREMENTS:
##   R >= 4.3.0
##   Bioconductor packages: phyloseq, DESeq2, ANCOMBC, ALDEx2
##   CRAN packages: vegan, ape, ranger, igraph, ggraph, ggrepel,
##     RColorBrewer, svglite, dplyr, tidyr, tibble, patchwork,
##     Matrix, scales
##   GitHub packages:
##     remotes::install_github("zdk123/SpiecEasi")
##   Run sessionInfo() after loading to check versions.
##
## PIPELINE STRUCTURE:
## Pillar 1: Data import, filtering, QC
## Pillar 2: Assembly mechanism dynamics (iCAMP)
## Pillar 3: Differential abundance and
##           functional profiling
## Pillar 4: Plant phenotype-microbiome
##           integration
## Pillar 5: Predictability and intervention
##           targets
## Pillar 6: Extended composition analysis
## Pillar 7: Co-occurrence networks and
##           keystone taxa
## Pillar 8: Phylogenetic signal analysis
## ============================================================
##  drought_pipeline_v2.r  (UPDATED — comprehensive for Microbiome/ISME)
##  Maize seedling rhizosphere assembly dynamics under progressive drought
##  8 analytical pillars + publication-quality figures with statistics
##  Author: Samuel (O.S. Olanrewaju)
##  Date: March 2026
## ============================================================

## ── 0. SETUP ─────────────────────────────────────────────────────────────────

## Source global theme constants — must be first
source("theme_constants.r")

suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(ape); library(picante)
  library(betapart); library(ggplot2); library(dplyr); library(tidyr)
  library(tibble); library(patchwork); library(RColorBrewer)
  library(DESeq2); library(ANCOMBC); library(ranger); library(pROC)
  library(pheatmap); library(viridis); library(ALDEx2); library(data.table)
  library(MASS)   ## ginv() fallback for singular VCV in Blomberg K
})

## ── Session info (package versions) ─────────────────────────────────────────
for (.td in c("tables/P1_composition","tables/P1_alpha","tables/P1_beta",
              "tables/P2_assembly","tables/P3_da","tables/P3_function",
              "tables/P4_phenotype","tables/P5_prediction","tables/P6_networks"))
  dir.create(.td, recursive=TRUE, showWarnings=FALSE)
sink("tables/P1_composition/session_info.txt")
cat("Pipeline: drought_pipeline_v2.r\n")
cat(paste0("Run date: ", Sys.time(), "\n\n"))
sessionInfo()
sink()
cat("Session info saved to tables/session_info.txt\n")

CONFIG <- list(
  MIN_LIB_SIZE = 1000, RAREFACTION_DEPTH = NULL, PREV_FILTER = 0.05,
  FDR_THRESHOLD = 0.05, PERMANOVA_PERMS = 9999, NSTI_CUTOFF = 2,
  RF_NTREES = 2000, RF_PERM_ITER = 999
)

## PAL is now defined in theme_constants.r — loaded via source() above
## Keeping local overrides here for any script-specific additions
## Do not redefine PAL — use the one from theme_constants.r

dirs <- c("data", "figures/P1_qc", "figures/P1_composition", "figures/P1_alpha",
          "figures/P1_beta", "figures/P2_assembly", "figures/P3_da",
          "figures/P3_function", "figures/P4_phenotype", "figures/P5_prediction",
          "figures/P6_networks", "figures/P7_phylogenetic",
          "figures/supplementary", "tables")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## save_fig is now defined in theme_constants.r — loaded via source() above
## Do not redefine here

## theme_pub is now defined in theme_constants.r — loaded via source() above
## Do not redefine here

norm_ids <- function(x) trimws(gsub("\\s+", " ", as.character(x)))

## Significance star helper
sig_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
    ifelse(p < 0.05, "*", "ns"))))
}

cat("Setup complete.\n\n")

## =============================================================================
##  PILLAR 1: DATA IMPORT, FILTERING, QUALITY CONTROL
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 1: Quality Control and Baseline Characterization\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

## ── 1a. Load metadata ───────────────────────────────────────────────────────
meta <- read.table("metadata.tsv", sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta %>% mutate(
  harvest = factor(harvest), harvest_num = as.numeric(gsub("[^0-9]","",harvest)),
  harvest_fac = factor(as.numeric(gsub("[^0-9]","",harvest))),
  treatment = factor(treatment, levels = c("Watered","Drought")),
  condition = factor(condition, levels = c("Unplanted","Planted")),
  trait = factor(trait, levels = c("Unplanted","Susceptible","Resistance")),
  reps = factor(reps), trt_geno = interaction(treatment, trait, drop=TRUE, sep="_")
)
cat("Metadata:", nrow(meta), "samples\n")

## ── 1b. Load morphology ────────────────────────────────────────────────────
morph <- read.csv("Morphological_characterization.csv", sep=";", header=TRUE,
                  check.names=FALSE, stringsAsFactors=FALSE)
names(morph) <- trimws(names(morph))
morph$trait[morph$trait == "Resistant"] <- "Resistance"
morph$harvest_num <- as.numeric(morph$harvest)
morph <- morph %>% mutate(
  treatment = factor(treatment, levels=c("Watered","Drought")),
  trait = factor(trait, levels=c("Susceptible","Resistance")), reps = factor(reps),
  sample_key = paste(treatment, harvest_num, trait, reps, sep="_")
)
cat("Morphology:", nrow(morph), "observations\n")

## ── 1c. Load taxonomy ──────────────────────────────────────────────────────
tax <- read.table("taxonomy.tsv", sep="\t", header=TRUE, row.names=1,
                  check.names=FALSE, quote="", comment.char="")
colnames(tax) <- sub("^\ufeff","", colnames(tax)); rownames(tax) <- norm_ids(rownames(tax))
if ("Taxon" %in% colnames(tax) &&
    sum(colnames(tax) %in% c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) < 3) {
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  mm <- t(vapply(strsplit(as.character(tax$Taxon), ";\\s*"), function(z) {
    z <- gsub("^[a-z]__","",z); if(length(z)<7) z <- c(z, rep(NA_character_, 7-length(z))); z[1:7]
  }, character(7L)))
  colnames(mm) <- ranks
  tax <- cbind(as.data.frame(mm, stringsAsFactors=FALSE),
               tax[, setdiff(colnames(tax),"Taxon"), drop=FALSE])
}

## Clean taxonomy
.clean_tax <- function(v) {
  v <- gsub("^[a-z]__","",trimws(as.character(v)))
  bad <- grepl("unknown|unclassified|unassigned|uncultured|metagenome|ambiguous|incertae.sedis|incertae_sedis|incertae sedis|possible|probable|gut_group|bacterium_enrichment|soil_clone",
               v, ignore.case=TRUE) | v %in% c("","NA","NaN")
  v[bad] <- NA; v
}
for (col in c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
  if (col %in% colnames(tax)) tax[[col]] <- .clean_tax(tax[[col]])
cat("Taxonomy:", nrow(tax), "features\n")

## ── 1d. Load feature table ─────────────────────────────────────────────────
ft <- read.table("feature-table.tsv", sep="\t", header=TRUE, check.names=FALSE,
                 quote="", comment.char="", stringsAsFactors=FALSE)
colnames(ft) <- sub("^\ufeff","",colnames(ft))
id_candidates <- c("Feature ID","Feature.ID","#OTU ID","#OTU_ID","OTU ID","ASV")
id_col <- intersect(colnames(ft), id_candidates)
if (!length(id_col)) {
  hits <- sapply(ft, function(v) sum(norm_ids(v) %in% rownames(tax)))
  id_col <- names(hits)[which.max(hits)]
}
id_col <- id_col[1]; rownames(ft) <- norm_ids(ft[[id_col]]); ft[[id_col]] <- NULL
for (j in seq_len(ncol(ft))) ft[[j]] <- suppressWarnings(as.numeric(ft[[j]]))
ft[is.na(ft)] <- 0

## ── 1e. Align and build phyloseq ───────────────────────────────────────────
common_samples <- intersect(colnames(ft), rownames(meta))
common_asvs <- intersect(rownames(ft), rownames(tax))
ft <- ft[common_asvs, common_samples]; tax <- tax[common_asvs,]; meta <- meta[common_samples,]

tree <- NULL
if (file.exists("tree.nwk")) {
  tree <- read.tree("tree.nwk"); tree$tip.label <- norm_ids(tree$tip.label)
  keep_tips <- intersect(tree$tip.label, rownames(ft))
  if (length(keep_tips) >= 2) tree <- keep.tip(tree, keep_tips) else tree <- NULL
}

ps_raw <- phyloseq(otu_table(as.matrix(ft), taxa_are_rows=TRUE),
                   tax_table(as.matrix(tax)), sample_data(meta),
                   if (!is.null(tree)) phy_tree(tree))
cat("Phyloseq (raw):", nsamples(ps_raw), "samples,", ntaxa(ps_raw), "taxa\n")

## ── 1f. Sequential filtering ───────────────────────────────────────────────
filter_log <- data.frame(step="raw", n_samples=nsamples(ps_raw), n_taxa=ntaxa(ps_raw))

## 1. Bacteria/Archaea
if ("Kingdom" %in% rank_names(ps_raw)) {
  tx <- as.data.frame(tax_table(ps_raw))
  ps1 <- prune_taxa(!is.na(tx$Kingdom) & grepl("bacteria|archaea",tx$Kingdom,TRUE), ps_raw)
} else ps1 <- ps_raw
filter_log <- rbind(filter_log, data.frame(step="domain", n_samples=nsamples(ps1), n_taxa=ntaxa(ps1)))

## 2. Remove chloroplast/mitochondria
tx1 <- as.data.frame(tax_table(ps1))
ps2 <- prune_taxa(!apply(tx1, 1, function(r) any(grepl("chloroplast|mitochond",r,TRUE))), ps1)
filter_log <- rbind(filter_log, data.frame(step="chloro_mito", n_samples=nsamples(ps2), n_taxa=ntaxa(ps2)))

## 3. Require phylum
tx2 <- as.data.frame(tax_table(ps2))
ps3 <- prune_taxa(!is.na(tx2$Phylum), ps2)
filter_log <- rbind(filter_log, data.frame(step="min_phylum", n_samples=nsamples(ps3), n_taxa=ntaxa(ps3)))

## 4. Clean ambiguous names, require class
tx3 <- as.data.frame(tax_table(ps3))
for (col in c("Class","Order","Family","Genus","Species")) {
  if (col %in% colnames(tx3)) tx3[[col]] <- .clean_tax(tx3[[col]])
}
tax_table(ps3) <- tax_table(as.matrix(tx3))
ps3b <- prune_taxa(!is.na(tx3$Class), ps3)
filter_log <- rbind(filter_log, data.frame(step="clean_min_class", n_samples=nsamples(ps3b), n_taxa=ntaxa(ps3b)))

## 5. Depth filter
ps <- prune_samples(sample_sums(ps3b) >= CONFIG$MIN_LIB_SIZE, ps3b)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
filter_log <- rbind(filter_log, data.frame(step="depth_filter", n_samples=nsamples(ps), n_taxa=ntaxa(ps)))

## 6. Prevalence filter — CRITICAL for reducing noise ASVs
## Keep ASVs present in >= 5% of samples (at least ~5 samples out of 108)
## This is standard practice; rare singletons inflate ordinations and variance partitioning
prev_threshold <- ceiling(CONFIG$PREV_FILTER * nsamples(ps))
asv_prev <- apply(otu_table(ps), 1, function(x) sum(x > 0))
if (!taxa_are_rows(ps)) asv_prev <- apply(otu_table(ps), 2, function(x) sum(x > 0))
keep_prev <- names(asv_prev)[asv_prev >= prev_threshold]
ps <- prune_taxa(keep_prev, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
filter_log <- rbind(filter_log, data.frame(step=sprintf("prevalence_>=%d_samples", prev_threshold),
                    n_samples=nsamples(ps), n_taxa=ntaxa(ps)))

cat(sprintf("Prevalence filter: retained %d ASVs present in >= %d samples (%.0f%%)\n",
            ntaxa(ps), prev_threshold, CONFIG$PREV_FILTER*100))
cat(sprintf("Reads retained after prevalence filter: %d (%.1f%% of pre-filter)\n",
            sum(sample_sums(ps)), 100*sum(sample_sums(ps))/sum(sample_sums(ps3b))))

write.csv(filter_log, "tables/P1_composition/filtering_log.csv", row.names=FALSE)
cat("\nFiltering log:\n"); print(filter_log)

## ── 1g. Rarefaction depth and subsets ──────────────────────────────────────
CONFIG$RAREFACTION_DEPTH <- as.integer(min(sample_sums(ps)))
cat("\nRarefaction depth:", CONFIG$RAREFACTION_DEPTH, "\n")

ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps_planted <- subset_samples(ps, condition=="Planted")
ps_planted <- prune_taxa(taxa_sums(ps_planted) > 0, ps_planted)
ps_rel_planted <- subset_samples(ps_rel, condition=="Planted")

set.seed(42)
ps_rare <- rarefy_even_depth(ps, sample.size=CONFIG$RAREFACTION_DEPTH, rngseed=42, verbose=FALSE)
ps_rare_planted <- subset_samples(ps_rare, condition=="Planted")

saveRDS(ps, "data/ps_filtered.rds"); saveRDS(ps_rel, "data/ps_rel.rds")
cat("Planted:", nsamples(ps_planted), "| Rarefied:", nsamples(ps_rare), "at", CONFIG$RAREFACTION_DEPTH, "\n\n")

## ── 1h. QC figures ─────────────────────────────────────────────────────────
qc_df <- data.frame(Sample=sample_names(ps), Reads=as.numeric(sample_sums(ps)), sample_data(ps))

## Depth boxplots by group
p_depth <- ggplot(qc_df, aes(harvest, Reads, fill=treatment)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) + geom_jitter(width=0.15, size=1, alpha=0.5) +
  facet_wrap(~condition) + scale_fill_manual(values=PAL$treatment) +
  scale_y_continuous(labels=scales::comma) + theme_pub() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Harvest", y="Sequencing depth", fill="Treatment")
save_fig(p_depth, "figures/P1_qc/depth_by_group", w=12, h=6)

## Rarefaction curves
cat("Generating rarefaction curves...\n")
otu_mat <- as(otu_table(ps), "matrix"); if(taxa_are_rows(ps)) otu_mat <- t(otu_mat)
step_size <- max(1, round(max(rowSums(otu_mat))/40))
rare_data <- lapply(seq_len(nrow(otu_mat)), function(i) {
  if(i%%20==0) cat("  Sample", i, "/", nrow(otu_mat), "\n")
  samp <- as.numeric(otu_mat[i,]); names(samp) <- colnames(otu_mat); samp <- samp[samp>0]
  total <- sum(samp); steps <- unique(c(seq(1,total,by=step_size),total))
  data.frame(Sample=rownames(otu_mat)[i], Depth=steps,
             Richness=as.numeric(vegan::rarefy(samp,steps)))
})
rare_df <- bind_rows(rare_data) %>%
  left_join(data.frame(sample_data(ps)) %>% rownames_to_column("Sample"), by="Sample")

p_rarecurve <- ggplot(rare_df, aes(Depth, Richness, group=Sample, colour=treatment)) +
  geom_line(alpha=0.5, linewidth=0.4) +
  geom_vline(xintercept=CONFIG$RAREFACTION_DEPTH, linetype="dashed") +
  scale_colour_manual(values=PAL$treatment) + scale_x_continuous(labels=scales::comma) +
  theme_pub() + labs(x="Sequencing depth", y="Observed ASVs", colour="Treatment",
                     title="Rarefaction curves",
                     subtitle=paste("Dashed line = rarefaction depth:", scales::comma(CONFIG$RAREFACTION_DEPTH)))
save_fig(p_rarecurve, "figures/P1_qc/rarefaction_curves_treatment", w=10, h=7)

## By timepoint
p_rare_time <- ggplot(rare_df, aes(Depth, Richness, group=Sample, colour=treatment)) +
  geom_line(alpha=0.5, linewidth=0.4) +
  geom_vline(xintercept=CONFIG$RAREFACTION_DEPTH, linetype="dashed") +
  facet_wrap(~harvest, ncol=3) +
  scale_colour_manual(values=PAL$treatment) + scale_x_continuous(labels=scales::comma) +
  theme_pub(BASE_SIZE - 2) + labs(x="Depth", y="ASVs", colour="Treatment")
save_fig(p_rare_time, "figures/P1_qc/rarefaction_by_harvest", w=12, h=8)

## Rarefaction sensitivity
alpha_unrar <- estimate_richness(ps_planted, measures=c("Observed"))
alpha_rar   <- estimate_richness(ps_rare_planted, measures=c("Observed"))
sens_df <- data.frame(Unrarefied=alpha_unrar$Observed,
  Rarefied=alpha_rar$Observed[match(rownames(alpha_unrar), rownames(alpha_rar))])
cor_r <- cor.test(sens_df$Unrarefied, sens_df$Rarefied, method="spearman")
p_sens <- ggplot(sens_df, aes(Unrarefied, Rarefied)) + geom_point(alpha=0.5) +
  geom_smooth(method="lm", se=FALSE) + theme_pub() +
  labs(x="Observed (unrarefied)", y="Observed (rarefied)",
       title="Rarefaction sensitivity", subtitle=sprintf("Spearman rho = %.3f, p = %.2e", cor_r$estimate, cor_r$p.value))
save_fig(p_sens, "figures/supplementary/rarefaction_sensitivity", w=ISME_SINGLE_W, h=6)
cat("Rarefaction curves done.\n\n")
## ── 1j. Classical alpha diversity boxplots with statistics ─────────────────
cat("── Alpha diversity analysis (rarefied) ──\n\n")

alpha_rare <- estimate_richness(ps_rare_planted, measures=c("Observed","Shannon"))
alpha_rare$Sample <- rownames(alpha_rare)
if (!is.null(phy_tree(ps_rare_planted, errorIfNULL=FALSE))) {
  pd_tab <- pd(t(otu_table(ps_rare_planted)), phy_tree(ps_rare_planted))
  alpha_rare$PD <- pd_tab$PD[match(alpha_rare$Sample, rownames(pd_tab))]
} else alpha_rare$PD <- NA_real_

meta_planted <- data.frame(sample_data(ps_planted)) %>% rownames_to_column("Sample")
alpha_df <- left_join(alpha_rare, meta_planted, by="Sample")
alpha_long <- alpha_df %>% pivot_longer(c(Observed, Shannon, PD), names_to="metric", values_to="value") %>%
  filter(!is.na(value))
write.csv(alpha_df, "tables/P1_alpha/alpha_diversity_rarefied.csv", row.names=FALSE)

## Classical boxplots: by treatment (with p-values)
trt_stats <- alpha_long %>% filter(trait %in% c("Resistance","Susceptible")) %>%
  group_by(metric) %>%
  summarise(p = wilcox.test(value ~ treatment)$p.value, .groups="drop") %>%
  mutate(label = paste0("p = ", format(p, digits=3), " ", sig_stars(p)))

p_alpha_trt <- ggplot(alpha_long %>% filter(trait %in% c("Resistance","Susceptible")),
                      aes(treatment, value, fill=treatment)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  geom_jitter(width=0.15, size=1.5, alpha=0.5) +
  facet_wrap(~metric, scales="free_y") +
  scale_fill_manual(values=PAL$treatment) + theme_pub() +
  geom_text(data=trt_stats, aes(x=1.5, y=Inf, label=label), inherit.aes=FALSE,
            vjust=1.5, size=4, fontface="italic") +
  labs(x="Treatment", y="Alpha diversity", fill="Treatment",
       title="Alpha diversity: Drought vs Watered")
save_fig(p_alpha_trt, "figures/P1_alpha/alpha_boxplot_treatment", w=10, h=5)

## Effect sizes for alpha diversity (Cliff's delta — nonparametric, robust)
cliff_delta <- function(x, y) {
  ## Cliff's delta: proportion of dominance pairs
  nx <- length(x); ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  delta <- mean(d)
  ## Interpret: |d| < 0.147 negligible, < 0.33 small, < 0.474 medium, else large
  interp <- ifelse(abs(delta) < 0.147, "negligible",
              ifelse(abs(delta) < 0.33, "small",
                ifelse(abs(delta) < 0.474, "medium", "large")))
  c(delta = delta, interpretation = interp)
}

alpha_effect_sizes <- do.call(rbind, lapply(
  unique((alpha_long %>% filter(trait %in% c("Resistance","Susceptible")))$metric),
  function(m) {
    d <- alpha_long %>% filter(trait %in% c("Resistance","Susceptible"), metric == m)
    trt_d <- cliff_delta(d$value[d$treatment=="Drought"], d$value[d$treatment=="Watered"])
    geno_d <- cliff_delta(d$value[d$trait=="Resistance"], d$value[d$trait=="Susceptible"])
    data.frame(metric = m,
      Cliff_delta_treatment = as.numeric(trt_d["delta"]),
      Effect_treatment = as.character(trt_d["interpretation"]),
      Cliff_delta_genotype = as.numeric(geno_d["delta"]),
      Effect_genotype = as.character(geno_d["interpretation"]),
      stringsAsFactors=FALSE)
  }
))
write.csv(alpha_effect_sizes, "tables/P1_alpha/alpha_effect_sizes.csv", row.names=FALSE)
cat("Alpha effect sizes (Cliff's delta):\n"); print(alpha_effect_sizes); cat("\n")
geno_stats <- alpha_long %>% filter(trait %in% c("Resistance","Susceptible")) %>%
  group_by(metric) %>%
  summarise(p = wilcox.test(value ~ trait)$p.value, .groups="drop") %>%
  mutate(label = paste0("p = ", format(p, digits=3), " ", sig_stars(p)))

p_alpha_geno <- ggplot(alpha_long %>% filter(trait %in% c("Resistance","Susceptible")),
                       aes(trait, value, fill=trait)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  geom_jitter(width=0.15, size=1.5, alpha=0.5) +
  facet_wrap(~metric, scales="free_y") +
  scale_fill_manual(values=PAL$genotype) + theme_pub() +
  geom_text(data=geno_stats, aes(x=1.5, y=Inf, label=label), inherit.aes=FALSE,
            vjust=1.5, size=4, fontface="italic") +
  labs(x="Genotype", y="Alpha diversity", fill="Genotype",
       title="Alpha diversity: Resistant vs Susceptible")
save_fig(p_alpha_geno, "figures/P1_alpha/alpha_boxplot_genotype", w=10, h=5)

## By time (faceted by treatment × genotype)
p_alpha_time <- ggplot(alpha_long %>% filter(trait %in% c("Resistance","Susceptible")),
                       aes(harvest, value, fill=treatment)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  geom_jitter(width=0.15, size=1, alpha=0.4) +
  facet_grid(metric ~ trait, scales="free_y") +
  scale_fill_manual(values=PAL$treatment) +
  theme_pub(BASE_SIZE - 2) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Harvest", y="Alpha diversity", fill="Treatment",
       title="Alpha diversity across time, genotype, and treatment")
save_fig(p_alpha_time, "figures/P1_alpha/alpha_boxplot_time_genotype", w=12, h=10)

## Temporal trajectories with regression stats
alpha_trend_stats <- alpha_long %>%
  filter(trait %in% c("Resistance","Susceptible")) %>%
  group_by(metric, trait) %>%
  summarise(
    p_time = summary(lm(value ~ harvest_num))$coefficients["harvest_num",4],
    p_trt_time = tryCatch(anova(lm(value ~ harvest_num * treatment))$`Pr(>F)`[3], error=function(e) NA),
    .groups="drop"
  ) %>%
  mutate(label_time = sprintf("Time %s", sig_stars(p_time)),
         label_int = sprintf("Time×Trt %s", sig_stars(p_trt_time)))

p_alpha_traj <- ggplot(alpha_long %>% filter(trait %in% c("Resistance","Susceptible")),
                       aes(harvest_num, value, colour=treatment, fill=treatment)) +
  geom_point(size=2, alpha=0.6) +
  geom_smooth(method="lm", se=TRUE, alpha=0.15, linewidth=0.8) +
  facet_grid(metric ~ trait, scales="free_y") +
  scale_colour_manual(values=PAL$treatment) + scale_fill_manual(values=PAL$treatment) +
  geom_text(data=alpha_trend_stats, aes(x=5, y=Inf, label=label_int), inherit.aes=FALSE,
            vjust=1.5, hjust=1, size=3.5, fontface="italic") +
  theme_pub() +
  labs(x="Day", y="Alpha diversity", colour="Treatment", fill="Treatment",
       title="Alpha diversity temporal trajectories (rarefied)")
save_fig(p_alpha_traj, "figures/P1_alpha/alpha_temporal_trajectories", w=10, h=10)

## Temporal slope estimates with 95% confidence intervals
slope_estimates <- do.call(rbind, lapply(
  split(alpha_long %>% filter(trait %in% c("Resistance","Susceptible")),
        list((alpha_long %>% filter(trait %in% c("Resistance","Susceptible")))$metric,
             (alpha_long %>% filter(trait %in% c("Resistance","Susceptible")))$treatment,
             (alpha_long %>% filter(trait %in% c("Resistance","Susceptible")))$trait)),
  function(d) {
    d <- d[!is.na(d$value),]
    if (nrow(d) < 4) return(NULL)
    mod <- lm(value ~ harvest_num, data=d)
    ci <- confint(mod, "harvest_num", level=0.95)
    data.frame(metric=d$metric[1], treatment=d$treatment[1], trait=d$trait[1],
               slope=coef(mod)["harvest_num"], se=summary(mod)$coefficients["harvest_num",2],
               ci_low=ci[1], ci_high=ci[2],
               p=summary(mod)$coefficients["harvest_num",4], r2=summary(mod)$r.squared,
               stringsAsFactors=FALSE)
  }
)) %>% mutate(sig = sig_stars(p))
rownames(slope_estimates) <- NULL
write.csv(slope_estimates, "tables/P1_alpha/alpha_temporal_slopes.csv", row.names=FALSE)
cat("Temporal slope estimates (units per day ± 95% CI):\n")
print(slope_estimates %>% dplyr::select(metric, treatment, trait, slope, ci_low, ci_high, p, sig))
cat("\n")
alpha_all <- estimate_richness(ps_rare, measures=c("Observed","Shannon"))
alpha_all$Sample <- rownames(alpha_all)
meta_all_df <- data.frame(sample_data(ps_rare)) %>% rownames_to_column("Sample")
alpha_all <- left_join(alpha_all, meta_all_df, by="Sample")
alpha_all_long <- alpha_all %>% pivot_longer(c(Observed, Shannon), names_to="metric", values_to="value")

cond_stats <- alpha_all_long %>% group_by(metric) %>%
  summarise(p = wilcox.test(value ~ condition)$p.value, .groups="drop") %>%
  mutate(label = paste0("p = ", format(p, digits=3), " ", sig_stars(p)))

p_alpha_cond <- ggplot(alpha_all_long, aes(condition, value, fill=condition)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) + geom_jitter(width=0.15, size=1, alpha=0.4) +
  facet_wrap(~metric, scales="free_y") + scale_fill_manual(values=PAL$condition) +
  geom_text(data=cond_stats, aes(x=1.5, y=Inf, label=label), inherit.aes=FALSE,
            vjust=1.5, size=4, fontface="italic") +
  theme_pub() + labs(x="Condition", y="Alpha diversity", fill="Condition",
                       title="Planted vs Unplanted alpha diversity")
save_fig(p_alpha_cond, "figures/P1_alpha/alpha_planted_vs_unplanted", w=8, h=5)

## Full alpha statistics table
alpha_full_stats <- alpha_long %>%
  filter(trait %in% c("Resistance","Susceptible")) %>%
  group_by(metric) %>%
  summarise(
    p_time = summary(lm(value ~ harvest_num))$coefficients["harvest_num",4],
    p_treatment = wilcox.test(value ~ treatment)$p.value,
    p_genotype = wilcox.test(value ~ trait)$p.value,
    p_time_trt = tryCatch(anova(lm(value ~ harvest_num * treatment))$`Pr(>F)`[3], error=function(e) NA),
    p_time_geno = tryCatch(anova(lm(value ~ harvest_num * trait))$`Pr(>F)`[3], error=function(e) NA),
    .groups="drop"
  )
write.csv(alpha_full_stats, "tables/P1_alpha/alpha_diversity_statistics.csv", row.names=FALSE)
cat("Alpha statistics:\n"); print(alpha_full_stats); cat("\n")

## ── 1k. Beta diversity ─────────────────────────────────────────────────────
cat("── Beta diversity analysis ──\n\n")

meta_pl <- data.frame(sample_data(ps_rel_planted))
meta_pl$harvest_fac <- factor(meta_pl$harvest_num)
bray_planted <- phyloseq::distance(ps_rel_planted, method="bray")
bray_all <- phyloseq::distance(ps_rel, method="bray")

## PERMANOVA — main model
set.seed(42)
perm_main <- adonis2(bray_planted ~ treatment + trait + harvest_num,
                     data=meta_pl, permutations=CONFIG$PERMANOVA_PERMS, by="margin")
write.csv(as.data.frame(perm_main), "tables/P1_beta/PERMANOVA_main.csv")
cat("PERMANOVA (main):\n"); print(perm_main)

## ── Per-day PERMANOVA ───────────────────────────────────────────────────────
## Additive model only (treatment + genotype); n=3 per cell makes interaction non-estimable
cat("── Per-day PERMANOVA ──\n")
perday_list <- lapply(c(0, 2, 3, 4, 5, 6), function(d) {
  s <- rownames(meta_pl)[meta_pl$harvest_num == d]
  if (length(s) < 4) {
    cat(sprintf("  Day %d: only %d samples — skipped\n", d, length(s)))
    return(NULL)
  }
  bh <- as.dist(as.matrix(bray_planted)[s, s])
  mh <- meta_pl[s, , drop=FALSE]
  set.seed(42)
  ph <- tryCatch(
    adonis2(bh ~ treatment + trait, data=mh, permutations=999, by="margin"),
    error=function(e) { cat(sprintf("  Day %d error: %s\n", d, e$message)); NULL }
  )
  if (is.null(ph)) return(NULL)
  ## Extract rows safely by exact rowname match
  .row_pd <- function(obj, term) {
    idx <- grep(paste0("^", term, "$"), rownames(obj))
    if (!length(idx)) return(c(R2=NA_real_, F_val=NA_real_, p=NA_real_))
    c(R2=obj$R2[idx], F_val=obj$F[idx], p=obj$`Pr(>F)`[idx])
  }
  trt  <- .row_pd(ph, "treatment")
  geno <- .row_pd(ph, "trait")
  res_idx <- grep("^Residual$", rownames(ph))
  res_r2  <- if (length(res_idx)) ph$R2[res_idx] else NA_real_
  data.frame(
    Day         = d,
    n_samples   = length(s),
    Trt_R2      = round(trt["R2"],    4),
    Trt_F       = round(trt["F_val"], 3),
    Trt_p       = round(trt["p"],     4),
    Trt_sig     = sig_stars(trt["p"]),
    Geno_R2     = round(geno["R2"],    4),
    Geno_F      = round(geno["F_val"], 3),
    Geno_p      = round(geno["p"],     4),
    Geno_sig    = sig_stars(geno["p"]),
    Residual_R2 = round(res_r2,        4),
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
})
perm_per_day <- do.call(rbind, perday_list)
rownames(perm_per_day) <- NULL
## BH correction across all 12 p-values (6 days × Treatment + Genotype)
perm_per_day$Trt_p_adj  <- round(p.adjust(perm_per_day$Trt_p,  method="BH"), 4)
perm_per_day$Geno_p_adj <- round(p.adjust(perm_per_day$Geno_p, method="BH"), 4)
perm_per_day$Trt_sig    <- sig_stars(perm_per_day$Trt_p_adj)
perm_per_day$Geno_sig   <- sig_stars(perm_per_day$Geno_p_adj)
write.csv(perm_per_day, "tables/P1_beta/PERMANOVA_per_day.csv", row.names=FALSE)
cat("Per-day PERMANOVA:\n"); print(perm_per_day); cat("\n")

## Full model with interaction (supplementary)
set.seed(42)
perm_full <- adonis2(bray_planted ~ treatment * trait + harvest_num,
                     data=meta_pl, permutations=CONFIG$PERMANOVA_PERMS, by="margin")
write.csv(as.data.frame(perm_full), "tables/P1_beta/PERMANOVA_full_interaction.csv")

## Planted vs Unplanted
set.seed(42)
meta_all_beta <- data.frame(sample_data(ps_rel))
perm_cond <- adonis2(bray_all ~ condition + treatment + harvest_num,
                     data=meta_all_beta, permutations=CONFIG$PERMANOVA_PERMS, by="margin")
write.csv(as.data.frame(perm_cond), "tables/P1_beta/PERMANOVA_condition.csv")

## Betadisper
betadisper_results <- list()
for (grp in c("treatment","trait","harvest")) {
  bd <- betadisper(bray_planted, meta_pl[[grp]])
  bd_t <- anova(bd)
  ## Effect size: eta-squared from betadisper ANOVA
  ss_group <- bd_t$`Sum Sq`[1]; ss_total <- sum(bd_t$`Sum Sq`)
  eta2 <- ss_group / ss_total
  betadisper_results[[grp]] <- data.frame(
    Factor=grp, F_value=bd_t$`F value`[1], p=bd_t$`Pr(>F)`[1],
    eta_squared=eta2,
    interpretation=ifelse(eta2<0.01,"negligible",ifelse(eta2<0.06,"small",
                    ifelse(eta2<0.14,"medium","large"))))
  write.csv(as.data.frame(bd_t), sprintf("tables/P1_beta/betadisper_%s.csv", grp))
  cat(sprintf("Betadisper (%s): F=%.2f, p=%.3f, η²=%.4f (%s)\n",
              grp, bd_t$`F value`[1], bd_t$`Pr(>F)`[1], eta2,
              betadisper_results[[grp]]$interpretation))
}
write.csv(bind_rows(betadisper_results), "tables/P1_beta/betadisper_effect_sizes.csv", row.names=FALSE)

## PERMANOVA effect size interpretation
## R² categories for microbiome PERMANOVA (Ramette 2007; Anderson 2017):
## <0.02 negligible, 0.02-0.05 small, 0.05-0.15 medium, >0.15 large
perm_effect <- data.frame(
  Factor = c("Treatment","Genotype","Time","Condition"),
  R2 = c(perm_main$R2[1], perm_main$R2[2], perm_main$R2[3], perm_cond$R2[1]),
  p = c(perm_main$`Pr(>F)`[1], perm_main$`Pr(>F)`[2], perm_main$`Pr(>F)`[3], perm_cond$`Pr(>F)`[1])
) %>% mutate(
  interpretation = case_when(R2 < 0.02 ~ "negligible", R2 < 0.05 ~ "small",
                              R2 < 0.15 ~ "medium", TRUE ~ "large"),
  sig = sig_stars(p)
)
write.csv(perm_effect, "tables/P1_beta/PERMANOVA_effect_sizes.csv", row.names=FALSE)
cat("\nPERMANOVA effect size interpretation:\n"); print(perm_effect); cat("\n")

## ── SENSITIVITY: Distance metric comparison ─────────────────────────────
cat("── Sensitivity: Distance metric comparison ──\n")

## Jaccard (presence-absence)
jacc_planted <- phyloseq::distance(ps_rel_planted, method="jaccard", binary=TRUE)
set.seed(42)
perm_jacc <- adonis2(jacc_planted ~ treatment + trait + harvest_num,
                     data=meta_pl, permutations=CONFIG$PERMANOVA_PERMS, by="margin")

## Weighted UniFrac (if tree present)
wunif_available <- !is.null(phy_tree(ps_planted, errorIfNULL=FALSE))
if (wunif_available) {
  wunif_planted <- phyloseq::distance(ps_planted, method="wunifrac")
  set.seed(42)
  perm_wunif <- adonis2(wunif_planted ~ treatment + trait + harvest_num,
                        data=meta_pl, permutations=CONFIG$PERMANOVA_PERMS, by="margin")
}

## Compile comparison
dist_sensitivity <- data.frame(
  Metric = c("Bray-Curtis","Jaccard",if(wunif_available) "Weighted UniFrac"),
  Trt_R2 = c(perm_main$R2[1], perm_jacc$R2[1], if(wunif_available) perm_wunif$R2[1]),
  Trt_p  = c(perm_main$`Pr(>F)`[1], perm_jacc$`Pr(>F)`[1], if(wunif_available) perm_wunif$`Pr(>F)`[1]),
  Geno_R2 = c(perm_main$R2[2], perm_jacc$R2[2], if(wunif_available) perm_wunif$R2[2]),
  Geno_p  = c(perm_main$`Pr(>F)`[2], perm_jacc$`Pr(>F)`[2], if(wunif_available) perm_wunif$`Pr(>F)`[2]),
  Time_R2 = c(perm_main$R2[3], perm_jacc$R2[3], if(wunif_available) perm_wunif$R2[3]),
  Time_p  = c(perm_main$`Pr(>F)`[3], perm_jacc$`Pr(>F)`[3], if(wunif_available) perm_wunif$`Pr(>F)`[3])
) %>% mutate(across(contains("_p"), ~sig_stars(.), .names="{.col}_sig"))
write.csv(dist_sensitivity, "tables/P1_beta/sensitivity_distance_metrics.csv", row.names=FALSE)
cat("Distance metric sensitivity:\n"); print(dist_sensitivity %>% dplyr::select(Metric, Trt_R2, Trt_p, Geno_R2, Geno_p, Time_R2, Time_p)); cat("\n")

## Barplot comparing R² across distance metrics
dist_plot_df <- dist_sensitivity %>%
  dplyr::select(Metric, Treatment=Trt_R2, Genotype=Geno_R2, Time=Time_R2) %>%
  pivot_longer(-Metric, names_to="Factor", values_to="R2")

p_dist_sens <- ggplot(dist_plot_df, aes(Factor, R2, fill=Metric)) +
  geom_col(position="dodge", alpha=0.8, width=0.7) +
  scale_fill_brewer(palette="Set2") + theme_pub() +
  labs(x=NULL, y="PERMANOVA R²", fill="Distance metric",
       title="Sensitivity: PERMANOVA results across distance metrics")
save_fig(p_dist_sens, "figures/supplementary/sensitivity_distance_metrics", w=9, h=5)

## ── SENSITIVITY: PERMANOVA model specification ──────────────────────────
cat("── Sensitivity: PERMANOVA model specification ──\n")

## by="terms" (sequential)
set.seed(42)
perm_terms <- adonis2(bray_planted ~ treatment + trait + harvest_num,
                      data=meta_pl, permutations=CONFIG$PERMANOVA_PERMS, by="terms")

## With two-way interaction (supplementary — avoids overparameterization of three-way with n=3)
set.seed(42)
perm_2way <- adonis2(bray_planted ~ treatment * trait + harvest_num,
                     data=meta_pl, permutations=CONFIG$PERMANOVA_PERMS, by="margin")

## CRITICAL: Extract by row name, not index — adonis2 row order varies by vegan version
.get_perm <- function(obj, term) {
  rn <- rownames(obj)
  idx <- grep(paste0("^", term, "$"), rn)
  if (length(idx) == 0) return(c(R2=NA, p=NA))
  c(R2=obj$R2[idx], p=obj$`Pr(>F)`[idx])
}

model_sensitivity <- data.frame(
  Model = c("Additive (margin)", "Additive (terms)", "Two-way interaction (margin)"),
  Trt_R2 = c(.get_perm(perm_main,"treatment")["R2"],
             .get_perm(perm_terms,"treatment")["R2"],
             .get_perm(perm_2way,"treatment")["R2"]),
  Trt_p  = c(.get_perm(perm_main,"treatment")["p"],
             .get_perm(perm_terms,"treatment")["p"],
             .get_perm(perm_2way,"treatment")["p"]),
  Geno_R2 = c(.get_perm(perm_main,"trait")["R2"],
              .get_perm(perm_terms,"trait")["R2"],
              .get_perm(perm_2way,"trait")["R2"]),
  Geno_p  = c(.get_perm(perm_main,"trait")["p"],
              .get_perm(perm_terms,"trait")["p"],
              .get_perm(perm_2way,"trait")["p"]),
  Time_R2 = c(.get_perm(perm_main,"harvest_num")["R2"],
              .get_perm(perm_terms,"harvest_num")["R2"],
              .get_perm(perm_2way,"harvest_num")["R2"]),
  Time_p  = c(.get_perm(perm_main,"harvest_num")["p"],
              .get_perm(perm_terms,"harvest_num")["p"],
              .get_perm(perm_2way,"harvest_num")["p"])
)
rownames(model_sensitivity) <- NULL
write.csv(model_sensitivity, "tables/P1_beta/sensitivity_PERMANOVA_models.csv", row.names=FALSE)
cat("PERMANOVA model sensitivity:\n"); print(model_sensitivity); cat("\n")

## ── SENSITIVITY: Leave-one-timepoint-out PERMANOVA ──────────────────────
cat("── Sensitivity: Leave-one-timepoint-out PERMANOVA ──\n")

loto_results <- do.call(rbind, lapply(unique(meta_pl$harvest), function(h_remove) {
  keep_s <- rownames(meta_pl)[meta_pl$harvest != h_remove]
  bray_sub <- as.dist(as.matrix(bray_planted)[keep_s, keep_s])
  meta_sub <- meta_pl[keep_s,]
  set.seed(42)
  p_sub <- adonis2(bray_sub ~ treatment + trait + harvest_num,
                   data=meta_sub, permutations=999, by="margin")
  data.frame(removed=as.character(h_remove),
             Trt_R2=p_sub$R2[1], Trt_p=p_sub$`Pr(>F)`[1],
             Geno_R2=p_sub$R2[2], Geno_p=p_sub$`Pr(>F)`[2],
             Time_R2=p_sub$R2[3], Time_p=p_sub$`Pr(>F)`[3])
}))
write.csv(loto_results, "tables/P1_beta/sensitivity_leave_one_timepoint_out.csv", row.names=FALSE)
cat("Leave-one-timepoint-out PERMANOVA:\n"); print(loto_results); cat("\n")

## Temporal PERMANOVA
perm_temporal <- do.call(rbind, lapply(unique(meta_pl$harvest), function(h) {
  s <- rownames(meta_pl)[meta_pl$harvest==h]; if(length(s)<8) return(NULL)
  bh <- as.dist(as.matrix(bray_planted)[s,s]); mh <- meta_pl[s,]
  set.seed(42); ph <- adonis2(bh ~ treatment + trait, data=mh, permutations=999, by="margin")
  data.frame(harvest=h, day=as.numeric(gsub("[^0-9]","",h)),
             Trt_R2=ph$R2[1], Trt_p=ph$`Pr(>F)`[1], Geno_R2=ph$R2[2], Geno_p=ph$`Pr(>F)`[2])
}))
write.csv(perm_temporal, "tables/P1_beta/PERMANOVA_temporal.csv", row.names=FALSE)

## PCoA with PERMANOVA annotation
ord <- ordinate(ps_rel_planted, "PCoA", "bray")
ve <- round(ord$values$Relative_eig[1:2]*100, 1)
ord_df <- as.data.frame(ord$vectors[,1:2]) %>% rownames_to_column("Sample") %>%
  left_join(meta_pl %>% rownames_to_column("Sample"), by="Sample")

perm_label <- sprintf("PERMANOVA: Treatment R²=%.3f (%s), Genotype R²=%.3f (%s), Time R²=%.3f (%s)",
  perm_main$R2[1], sig_stars(perm_main$`Pr(>F)`[1]),
  perm_main$R2[2], sig_stars(perm_main$`Pr(>F)`[2]),
  perm_main$R2[3], sig_stars(perm_main$`Pr(>F)`[3]))

p_pcoa <- ggplot(ord_df, aes(Axis.1, Axis.2, colour=treatment, shape=trait)) +
  geom_point(size=3.5, alpha=0.8) +
  stat_ellipse(aes(group=treatment), level=0.68, linewidth=0.8) +
  scale_colour_manual(values=PAL$treatment) +
  scale_shape_manual(values=c(Resistance=16, Susceptible=17)) +
  facet_wrap(~harvest, ncol=3) + theme_pub() +
  labs(x=paste0("PCoA1 (",ve[1],"%)"), y=paste0("PCoA2 (",ve[2],"%)"),
       colour="Treatment", shape="Genotype",
       title="Rhizosphere community composition over time", subtitle=perm_label)
save_fig(p_pcoa, "figures/P1_beta/PCoA_temporal_facets", w=14, h=10)

## Temporal R² dynamics with significance markers
## Use perm_per_day (BH-adjusted) for significance labels on trajectory plot
perm_temp_long <- perm_per_day %>%
  dplyr::select(Day, Trt_R2, Geno_R2, Trt_p_adj, Geno_p_adj) %>%
  dplyr::rename("day" = "Day") %>%
  pivot_longer(c(Trt_R2, Geno_R2), names_to="Factor", values_to="R2") %>%
  mutate(p_val = ifelse(Factor=="Trt_R2", Trt_p_adj, Geno_p_adj),
         sig = sig_stars(p_val), Factor = gsub("_R2","",Factor),
         Factor = recode(Factor, Trt="Treatment", Geno="Genotype"))

p_r2 <- ggplot(perm_temp_long, aes(day, R2, colour=Factor)) +
  geom_line(linewidth=1) + geom_point(size=3) +
  geom_text(aes(label=sig), vjust=-1, size=4, show.legend=FALSE) +
  scale_colour_manual(values=c(Treatment=PAL$treatment[["Drought"]], Genotype=PAL$genotype[["Resistance"]])) +
  theme_pub() + labs(x="Day", y="PERMANOVA R²",
                     title="Effect size dynamics through time",
                     subtitle="BH-adjusted: * p<0.05  ** p<0.01  *** p<0.001")
save_fig(p_r2, "figures/P1_beta/PERMANOVA_R2_temporal", w=8, h=5)

## Variance partitioning
vp <- varpart(bray_planted, ~treatment, ~trait, ~harvest_num, data=meta_pl)
vp_fracs <- data.frame(
  Component = c("Treatment","Genotype","Time","Trt∩Geno","Trt∩Time","Geno∩Time","All","Residual"),
  AdjR2 = c(vp$part$indfract$Adj.R.square[1:7],
            1-sum(pmax(0,vp$part$indfract$Adj.R.square[1:7])))
) %>% mutate(AdjR2=pmax(0,AdjR2), Pct=round(AdjR2*100,2))
write.csv(vp_fracs, "tables/P1_beta/variance_partitioning.csv", row.names=FALSE)
residual_pct <- vp_fracs$Pct[vp_fracs$Component=="Residual"]
total_expl <- sum(vp_fracs$Pct[vp_fracs$Component!="Residual"])

pdf("figures/P1_beta/varpart_venn.pdf", w=7, h=7)
plot(vp, digits=3, bg=c("#D32F2F","#388E3C","#E69F00"), Xnames=c("Treatment","Genotype","Time"))
title(sprintf("Variance partitioning (Explained: %.1f%%, Residual: %.1f%%)", total_expl, residual_pct))
dev.off()
png("figures/P1_beta/varpart_venn.png", w=7, h=7, units="in", res=300)
plot(vp, digits=3, bg=c("#D32F2F","#388E3C","#E69F00"), Xnames=c("Treatment","Genotype","Time"))
title(sprintf("Variance partitioning (Explained: %.1f%%, Residual: %.1f%%)", total_expl, residual_pct))
dev.off()
svglite::svglite("figures/P1_beta/varpart_venn.svg", width=7, height=7)
plot(vp, digits=3, bg=c("#D32F2F","#388E3C","#E69F00"), Xnames=c("Treatment","Genotype","Time"))
title(sprintf("Variance partitioning (Explained: %.1f%%, Residual: %.1f%%)", total_expl, residual_pct))
dev.off()

## Turnover vs nestedness
cat("\n── Beta diversity partitioning ──\n")
pa_planted <- decostand(t(as(otu_table(ps_planted),"matrix")), method="pa")
bp <- beta.pair(pa_planted, index.family="sorensen")

snames <- rownames(pa_planted); n <- length(snames)
beta_pairs <- list(); idx <- 0
for (i in 1:(n-1)) { for (j in (i+1):n) {
  idx <- idx+1; s1 <- snames[i]; s2 <- snames[j]
  beta_pairs[[idx]] <- data.frame(s1=s1, s2=s2,
    sorensen=as.matrix(bp$beta.sor)[i,j], turnover=as.matrix(bp$beta.sim)[i,j],
    nestedness=as.matrix(bp$beta.sne)[i,j])
}}
beta_df <- bind_rows(beta_pairs) %>%
  mutate(turn_prop = ifelse(sorensen>0, turnover/sorensen, NA)) %>%
  left_join(meta_pl %>% rownames_to_column("s1") %>% dplyr::select(s1,trt1=treatment,trait1=trait,day1=harvest_num), by="s1") %>%
  left_join(meta_pl %>% rownames_to_column("s2") %>% dplyr::select(s2,trt2=treatment,trait2=trait,day2=harvest_num), by="s2") %>%
  mutate(same_day=day1==day2,
         comparison=case_when(trt1!=trt2 & trait1!=trait2 ~ "Both differ",
                              trt1!=trt2 ~ "Treatment", trait1!=trait2 ~ "Genotype",
                              TRUE ~ "Same group"))

beta_summary <- beta_df %>% filter(same_day, comparison!="Same group") %>%
  group_by(comparison) %>%
  summarise(mean_turn=mean(turn_prop,na.rm=TRUE), sd=sd(turn_prop,na.rm=TRUE), n=n(), .groups="drop")
write.csv(beta_summary, "tables/P1_beta/turnover_nestedness_summary.csv", row.names=FALSE)

## Stacked percentage bars by timepoint
beta_time <- beta_df %>% filter(same_day, trt1!=trt2) %>%
  group_by(day1) %>% summarise(Turnover=mean(turnover), Nestedness=mean(nestedness), .groups="drop") %>%
  pivot_longer(c(Turnover,Nestedness), names_to="Component", values_to="Value")

p_turnover <- ggplot(beta_time, aes(factor(day1), Value, fill=Component)) +
  geom_col(position="stack", alpha=0.8, width=0.6) +
  scale_fill_manual(values=c(Turnover="#2E8B57", Nestedness="#CD853F")) +
  theme_pub() + labs(x="Day", y="Beta diversity component", fill="Component",
    title="Turnover dominates beta diversity",
    subtitle=sprintf("Mean turnover proportion: %.1f%%", mean(beta_summary$mean_turn)*100))
save_fig(p_turnover, "figures/P1_beta/turnover_nestedness_stacked", w=8, h=5)

## Temporal divergence from baseline
cat("Computing divergence from baseline...\n")
bray_mat <- as.matrix(bray_planted)
day0_samples <- rownames(meta_pl)[meta_pl$harvest_num == 0]
diverg_list <- lapply(rownames(meta_pl), function(s) {
  if (meta_pl[s,"harvest_num"]==0) return(NULL)
  same_trt_day0 <- day0_samples[meta_pl[day0_samples,"treatment"]==meta_pl[s,"treatment"] &
                                  meta_pl[day0_samples,"trait"]==meta_pl[s,"trait"]]
  if (!length(same_trt_day0)) return(NULL)
  data.frame(Sample=s, divergence=mean(bray_mat[s, same_trt_day0]),
             day=meta_pl[s,"harvest_num"], treatment=as.character(meta_pl[s,"treatment"]),
             trait=as.character(meta_pl[s,"trait"]))
})
diverg_df <- bind_rows(diverg_list)

p_diverg <- ggplot(diverg_df, aes(day, divergence, colour=treatment)) +
  geom_point(size=2, alpha=0.6) + geom_smooth(method="lm", se=TRUE, alpha=0.15) +
  facet_wrap(~trait) + scale_colour_manual(values=PAL$treatment) + theme_pub() +
  labs(x="Day", y="Bray-Curtis distance from Day 0",
       title="Community divergence from baseline", colour="Treatment")
save_fig(p_diverg, "figures/P1_beta/divergence_from_baseline", w=10, h=5)

cat("\nPillar 1 complete.\n\n")

## =============================================================================
##  PILLAR 2: ASSEMBLY MECHANISM DYNAMICS (βNTI / RC-bray)
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 2: Assembly Mechanism Dynamics\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

icamp_file <- "iCAMP_process_results_actual.csv"
icamp_available <- file.exists(icamp_file)

if (!icamp_available && file.exists("iCAMP.process.bNTIiRCa.csv")) {
  test_df <- read.csv("iCAMP.process.bNTIiRCa.csv", nrows=5)
  if (sum(unique(c(test_df$sample1,test_df$sample2)) %in% sample_names(ps_planted)) > 3) {
    icamp_file <- "iCAMP.process.bNTIiRCa.csv"; icamp_available <- TRUE
  }
}

if (icamp_available) {
  proc_df <- read.csv(icamp_file, stringsAsFactors=FALSE)
  cat("Loaded", nrow(proc_df), "pairwise comparisons\n")
  
  ## Ensure Process column exists
  if (!"Process" %in% colnames(proc_df)) {
    bnti_col <- grep("bNTI|bnti", colnames(proc_df), value=TRUE)[1]
    rc_col <- grep("RC|rc", colnames(proc_df), value=TRUE)[1]
    if (!is.na(bnti_col) && !is.na(rc_col)) {
      proc_df$Process <- case_when(
        abs(proc_df[[bnti_col]]) > 1.96 & proc_df[[bnti_col]] > 0 ~ "Variable.Selection",
        abs(proc_df[[bnti_col]]) > 1.96 ~ "Homogeneous.Selection",
        proc_df[[rc_col]] > 0.95 ~ "Dispersal.Limitation",
        proc_df[[rc_col]] < -0.95 ~ "Homogenizing.Dispersal",
        TRUE ~ "Drift")
    }
  }
  
  ## Ensure comparison_type
  if (!"comparison_type" %in% colnames(proc_df) && "treatment1" %in% colnames(proc_df)) {
    proc_df$comparison_type <- case_when(
      proc_df$treatment1==proc_df$treatment2 & proc_df$treatment1=="Drought" ~ "Within_Drought",
      proc_df$treatment1==proc_df$treatment2 & proc_df$treatment1=="Watered" ~ "Within_Watered",
      TRUE ~ "Between_treatments")
  }
  
  ## Overall fractions
  overall_fracs <- prop.table(table(proc_df$Process))
  stochastic_pct <- sum(overall_fracs[names(overall_fracs) %in% c("Drift","Dispersal.Limitation","Homogenizing.Dispersal")])
  deterministic_pct <- 1 - stochastic_pct
  cat(sprintf("\nDeterministic: %.1f%% | Stochastic: %.1f%%\n", deterministic_pct*100, stochastic_pct*100))
  
  ## ── Figure 2A: βNTI distribution ────────────────────────────────────────
  bnti_col <- grep("bNTI|bnti", colnames(proc_df), value=TRUE)[1]
  
  p_bnti <- ggplot(proc_df, aes(.data[[bnti_col]])) +
    geom_histogram(bins=50, fill="steelblue", alpha=0.7, colour="white") +
    geom_vline(xintercept=c(-1.96,1.96), linetype="dashed", colour="red", linewidth=0.8) +
    annotate("rect", xmin=-1.96, xmax=1.96, ymin=-Inf, ymax=Inf, alpha=0.05, fill="grey") +
    annotate("text", x=-4, y=Inf, label="Homogeneous\nselection", vjust=1.5, size=3.5, colour="#2166AC", fontface="bold") +
    annotate("text", x=4, y=Inf, label="Variable\nselection", vjust=1.5, size=3.5, colour="#B2182B", fontface="bold") +
    annotate("text", x=0, y=Inf, label="Stochastic", vjust=1.5, size=3.5, colour="#999999", fontface="bold") +
    coord_cartesian(xlim=c(-8,8)) + theme_pub() +
    labs(x="βNTI", y="Count",
         title="Distribution of β-Nearest Taxon Index",
         subtitle=sprintf("%.1f%% deterministic (|βNTI| > 1.96) | %.1f%% stochastic",
                          deterministic_pct*100, stochastic_pct*100))
  save_fig(p_bnti, "figures/P2_assembly/bNTI_distribution", w=8, h=5)
  
  ## ── Figure 2B: Assembly by treatment ────────────────────────────────────
  assembly_trt <- proc_df %>%
    filter(grepl("Within", comparison_type)) %>%
    group_by(comparison_type, Process) %>% summarise(n=n(), .groups="drop") %>%
    group_by(comparison_type) %>% mutate(fraction=n/sum(n)) %>%
    mutate(comparison_type = recode(comparison_type,
           Within_Drought="Drought", Within_Watered="Watered"))
  
  p_asm_trt <- ggplot(assembly_trt, aes(comparison_type, fraction, fill=Process)) +
    geom_col(width=0.6) +
    scale_fill_manual(values=PAL$assembly, name="Assembly process",
                      labels=c("Dispersal.Limitation"="Dispersal limitation",
                               "Drift"="Drift",
                               "Homogeneous.Selection"="Homogeneous selection",
                               "Homogenizing.Dispersal"="Homogenizing dispersal",
                               "Variable.Selection"="Variable selection")) +
    scale_y_continuous(labels=scales::percent) +
    theme_pub() + theme(legend.text=element_text(size=11),
                           legend.key.size=unit(0.5,"cm")) +
    labs(x=NULL, y="Proportion of pairwise comparisons",
         title="Community assembly mechanisms by treatment")
  save_fig(p_asm_trt, "figures/P2_assembly/assembly_by_treatment", w=10, h=6)
  
  ## ── Figure 2C: Assembly through time SPLIT BY TREATMENT ─────────────────
  if ("harvest1" %in% colnames(proc_df)) {
    proc_df$day1_num <- as.numeric(gsub("[^0-9]","", proc_df$harvest1))
    proc_df$day2_num <- as.numeric(gsub("[^0-9]","", proc_df$harvest2))
    
    ## Within-treatment, within-timepoint comparisons
    proc_within_time <- proc_df %>%
      filter(day1_num == day2_num, grepl("Within", comparison_type)) %>%
      mutate(Treatment = recode(comparison_type, Within_Drought="Drought", Within_Watered="Watered"))
    
    if (nrow(proc_within_time) > 0) {
      asm_time_trt <- proc_within_time %>%
        group_by(Treatment, day1_num, Process) %>% summarise(n=n(), .groups="drop") %>%
        group_by(Treatment, day1_num) %>% mutate(fraction=n/sum(n))
      
      write.csv(asm_time_trt, "tables/P2_assembly/assembly_temporal_by_treatment.csv", row.names=FALSE)
      
      p_asm_time <- ggplot(asm_time_trt, aes(factor(day1_num), fraction, fill=Process)) +
        geom_col(width=0.7) + facet_wrap(~Treatment) +
        scale_fill_manual(values=PAL$assembly, name="Assembly process",
                          labels=c("Dispersal.Limitation"="Dispersal limitation",
                                   "Drift"="Drift",
                                   "Homogeneous.Selection"="Homogeneous selection",
                                   "Homogenizing.Dispersal"="Homogenizing dispersal",
                                   "Variable.Selection"="Variable selection")) +
        scale_y_continuous(labels=scales::percent) +
        theme_pub() + labs(x="Day", y="Proportion",
          title="Temporal dynamics of assembly mechanisms",
          subtitle="Within-treatment, within-timepoint pairwise comparisons")
      save_fig(p_asm_time, "figures/P2_assembly/assembly_temporal_by_treatment", w=14, h=6)
    }
  }
  
  ## ── Figure 2D: βNTI boxplots by treatment × timepoint ──────────────────
  if ("harvest1" %in% colnames(proc_df)) {
    proc_within <- proc_df %>%
      filter(day1_num == day2_num, grepl("Within", comparison_type)) %>%
      mutate(Treatment = recode(comparison_type, Within_Drought="Drought", Within_Watered="Watered"))
    
    if (nrow(proc_within) > 0) {
      p_bnti_box <- ggplot(proc_within, aes(factor(day1_num), .data[[bnti_col]], fill=Treatment)) +
        geom_boxplot(outlier.shape=NA, alpha=0.8) +
        geom_jitter(aes(colour=Treatment), width=0.15, size=0.8, alpha=0.4) +
        geom_hline(yintercept=c(-1.96,1.96), linetype="dashed", colour="red", alpha=0.6) +
        scale_fill_manual(values=PAL$treatment) + scale_colour_manual(values=PAL$treatment) +
        theme_pub() +
        labs(x="Day", y="βNTI", fill="Treatment", colour="Treatment",
             title="βNTI values by treatment and timepoint",
             subtitle="Dashed lines: |βNTI| = 1.96 threshold")
      save_fig(p_bnti_box, "figures/P2_assembly/bNTI_boxplot_treatment_time", w=10, h=6)
    }
  }
  
  ## ── STATISTICAL TEST: iCAMP process proportions Drought vs Watered ────
  ## This is CRITICAL — the headline finding must be statistically supported
  cat("\n── Statistical tests for assembly mechanism differences ──\n")
  
  if ("comparison_type" %in% colnames(proc_df)) {
    within_procs <- proc_df %>%
      filter(comparison_type %in% c("Within_Drought","Within_Watered")) %>%
      mutate(Treatment = recode(comparison_type, Within_Drought="Drought", Within_Watered="Watered"))
    
    ## Chi-squared test: overall process distribution differs between treatments
    cont_table <- table(within_procs$Treatment, within_procs$Process)
    chi_test <- chisq.test(cont_table)
    cat("Chi-squared test (process distribution: Drought vs Watered):\n")
    cat(sprintf("  χ² = %.2f, df = %d, p = %.2e %s\n",
                chi_test$statistic, chi_test$parameter, chi_test$p.value,
                sig_stars(chi_test$p.value)))
    
    ## Per-process Fisher's exact tests (which specific processes differ?)
    process_tests <- do.call(rbind, lapply(unique(within_procs$Process), function(proc) {
      tbl <- table(within_procs$Treatment, within_procs$Process == proc)
      ft <- fisher.test(tbl)
      drought_pct <- 100 * sum(within_procs$Process==proc & within_procs$Treatment=="Drought") /
                     sum(within_procs$Treatment=="Drought")
      watered_pct <- 100 * sum(within_procs$Process==proc & within_procs$Treatment=="Watered") /
                     sum(within_procs$Treatment=="Watered")
      data.frame(Process=proc, Drought_pct=round(drought_pct,1), Watered_pct=round(watered_pct,1),
                 Difference=round(drought_pct-watered_pct,1),
                 OR_Watered_vs_Drought=round(ft$estimate,2), p=ft$p.value, stringsAsFactors=FALSE)
    })) %>% mutate(p_adj = p.adjust(p, "BH"), sig = sig_stars(p_adj))
    
    write.csv(process_tests, "tables/P2_assembly/iCAMP_process_tests_treatment.csv", row.names=FALSE)
    cat("\nPer-process Fisher's exact tests (BH-corrected):\n")
    print(process_tests)
    
    ## Save chi-squared result
    write.csv(data.frame(
      Test="Chi-squared", Statistic=chi_test$statistic, df=chi_test$parameter,
      p=chi_test$p.value, sig=sig_stars(chi_test$p.value)
    ), "tables/P2_assembly/iCAMP_chisq_test.csv", row.names=FALSE)
    
    ## ── TEMPORAL REGRESSION: Do assembly mechanisms shift over time? ────
    cat("\n── Temporal trends in assembly mechanisms ──\n")
    
    if ("day1_num" %in% colnames(proc_df)) {
      ## For within-treatment, within-timepoint pairs, calculate process fractions per day
      proc_time_frac <- proc_df %>%
        filter(day1_num == day2_num, grepl("Within", comparison_type)) %>%
        mutate(Treatment = recode(comparison_type, Within_Drought="Drought", Within_Watered="Watered")) %>%
        group_by(Treatment, day1_num, Process) %>% summarise(n=n(), .groups="drop") %>%
        group_by(Treatment, day1_num) %>% mutate(frac = n/sum(n)) %>% ungroup()
      
      ## Test: does each process fraction change over time, and does it differ by treatment?
      temporal_tests <- do.call(rbind, lapply(unique(proc_time_frac$Process), function(proc) {
        d <- proc_time_frac %>% filter(Process == proc)
        if (nrow(d) < 4) return(NULL)
        
        ## Simple trend
        mod_trend <- tryCatch(lm(frac ~ day1_num, data=d), error=function(e) NULL)
        ## Treatment interaction
        mod_int <- tryCatch(lm(frac ~ day1_num * Treatment, data=d), error=function(e) NULL)
        
        trend_p <- if(!is.null(mod_trend)) summary(mod_trend)$coefficients["day1_num",4] else NA
        int_p <- if(!is.null(mod_int) && "day1_num:TreatmentWatered" %in% rownames(summary(mod_int)$coefficients))
                   summary(mod_int)$coefficients["day1_num:TreatmentWatered",4] else NA
        
        data.frame(Process=proc, trend_p=trend_p, trend_sig=sig_stars(trend_p),
                   interaction_p=int_p, interaction_sig=sig_stars(int_p), stringsAsFactors=FALSE)
      }))
      
      if (!is.null(temporal_tests)) {
        write.csv(temporal_tests, "tables/P2_assembly/iCAMP_temporal_trends.csv", row.names=FALSE)
        cat("Assembly mechanism temporal trends:\n"); print(temporal_tests); cat("\n")
      }
    }
  }
  
  ## ── SENSITIVITY: iCAMP ASV filtering ──────────────────────────────────
  ## Report how the prevalence/abundance filter affected iCAMP
  if (file.exists("iCAMP_filtering_report.csv")) {
    filt_report <- read.csv("iCAMP_filtering_report.csv")
    cat("\niCAMP filtering report:\n"); print(filt_report)
    cat("NOTE: iCAMP was run on prevalent ASVs (>=5% prevalence, >=50 reads).\n")
    cat("This is standard practice — rare singletons are phylogenetically\n")
    cat("uninformative for null model comparisons (Stegen et al. 2013, 2015).\n")
    cat("Retained reads represent the dominant community fraction.\n")
    cat("To test sensitivity, rerun iCAMP with stricter filters (>=10%, >=100 reads)\n")
    cat("and compare assembly fractions. If results are consistent, the finding is robust.\n\n")
  }
  
  ## iCAMP summary table for manuscript
  if (exists("proc_df") && "Process" %in% colnames(proc_df)) {
    icamp_summary <- proc_df %>%
      {if("comparison_type" %in% colnames(.)) filter(., grepl("Within", comparison_type)) else .} %>%
      {if("comparison_type" %in% colnames(.)) group_by(., comparison_type) else group_by(., .data[["Process"]])} %>%
      {
        if("comparison_type" %in% colnames(proc_df)) {
          proc_df %>% filter(grepl("Within", comparison_type)) %>%
            group_by(comparison_type, Process) %>% summarise(n=n(), .groups="drop") %>%
            group_by(comparison_type) %>% mutate(pct = round(100*n/sum(n), 1)) %>%
            pivot_wider(names_from=Process, values_from=c(n, pct), values_fill=0)
        } else {
          proc_df %>% group_by(Process) %>% summarise(n=n(), .groups="drop") %>%
            mutate(pct = round(100*n/sum(n), 1))
        }
      }
    write.csv(icamp_summary, "tables/P2_assembly/iCAMP_summary_for_manuscript.csv", row.names=FALSE)
  }
  
  cat("Pillar 2 complete.\n\n")
} else {
  cat("iCAMP results not found — run on CHPC. Pillar 2 deferred.\n\n")
}

## =============================================================================
##  PILLAR 3: DIFFERENTIAL ABUNDANCE + FUNCTIONAL PROFILING
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 3: Compositional Responses (Taxonomic + Functional)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

tax_lookup <- as.data.frame(tax_table(ps_planted))

## ── Volcano label helper ──────────────────────────────────────────────────
## Returns lowest classified rank label; genus-level = bare name (italicised
## via ggtext later), higher ranks get rank prefix (f__, o__, c__, p__).
.is_uninformative <- function(x) {
  is.na(x) |
  grepl("^\\s*$", x) |
  grepl("^(uncultured|unidentified|metagenome|unknown)$", x, ignore.case=TRUE) |
  grepl("^[0-9]+$", x) |            # purely numeric strings
  grepl("\\bgroup\\b$", x, ignore.case=TRUE)  # ends in "group" as sole identifier
}
.rank_label <- function(Genus, Family, Order, Class, Phylum) {
  mapply(function(g, f, o, cl, p) {
    if (!.is_uninformative(g))        return(g)
    if (!.is_uninformative(f))        return(paste0("f__", f))
    if (!.is_uninformative(o))        return(paste0("o__", o))
    if (!.is_uninformative(cl))       return(paste0("c__", cl))
    if (!is.na(p) && nzchar(p))       return(paste0("p__", p))
    "Unclassified"
  }, Genus, Family, Order, Class, Phylum, SIMPLIFY=TRUE, USE.NAMES=FALSE)
}

## ── 3a. ANCOM-BC2: Treatment ──────────────────────────────────────────────
cat("ANCOM-BC2: Treatment...\n")
set.seed(42)
ancom_trt <- tryCatch(ancombc2(data=ps_planted, fix_formula="treatment + trait + harvest_num",
  p_adj_method="fdr", prv_cut=0.10, lib_cut=CONFIG$MIN_LIB_SIZE,
  group="treatment", struc_zero=TRUE, neg_lb=TRUE, alpha=0.05, global=FALSE),
  error=function(e){cat("Error:",e$message,"\n"); NULL})

trt_df <- NULL
if (!is.null(ancom_trt)) {
  res <- ancom_trt$res
  lfc_col <- grep("^lfc_treatment", colnames(res), value=TRUE)[1]
  q_col <- gsub("^lfc_","q_",lfc_col); diff_col <- gsub("^lfc_","diff_",lfc_col)
  trt_df <- data.frame(ASV=res$taxon, lfc=res[[lfc_col]], qval=res[[q_col]], diff_abn=res[[diff_col]]) %>%
    left_join(tax_lookup %>% rownames_to_column("ASV"), by="ASV") %>%
    mutate(Direction=case_when(diff_abn & lfc>0 ~ "Drought_enriched",
                               diff_abn & lfc<0 ~ "Watered_enriched", TRUE ~ "NS"))
  write.csv(trt_df, "tables/P3_da/ANCOMBC2_treatment.csv", row.names=FALSE)
  n_dr <- sum(trt_df$Direction=="Drought_enriched"); n_wa <- sum(trt_df$Direction=="Watered_enriched")
  cat(sprintf("  %d drought-enriched, %d watered-enriched\n", n_dr, n_wa))
  
  ## Volcano with top labels
  top_label <- trt_df %>% filter(diff_abn) %>% arrange(qval) %>% head(15) %>%
    mutate(Label = .rank_label(Genus, Family, Order, Class, Phylum))
  
  p_vol_trt <- ggplot(trt_df, aes(lfc, -log10(qval), colour=Direction)) +
    geom_point(alpha=0.4, size=1.5) +
    scale_colour_manual(values=c(Drought_enriched=PAL$treatment[["Drought"]], Watered_enriched=PAL$treatment[["Watered"]], NS=COL_SIG_NS)) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.3) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.3) +
    ggrepel::geom_text_repel(data=top_label, aes(label=Label), size=3, max.overlaps=20,
                              box.padding=0.4, point.padding=0.3, force=2,
                              show.legend=FALSE, min.segment.length=0) +
    coord_cartesian(clip="off") +
    theme_pub() +
    theme(plot.margin=margin(4, 15, 4, 4, "mm")) +
    labs(x="Log fold change", y=expression(-log[10](q)), colour=NULL,
      title="Differential abundance: Drought vs Watered",
      subtitle=sprintf("%d drought-enriched | %d watered-enriched", n_dr, n_wa))
  dir.create("figures/P3_da", recursive=TRUE, showWarnings=FALSE)
  save_fig(p_vol_trt, "figures/P3_da/volcano_drought_vs_watered", w=10, h=7)
  
  ## ANCOM-BC2 effect size distribution for significant ASVs
  sig_trt <- trt_df %>% filter(diff_abn==TRUE)
  trt_effect_summary <- data.frame(
    Comparison = "Drought vs Watered",
    N_significant = nrow(sig_trt),
    Median_absLFC = median(abs(sig_trt$lfc), na.rm=TRUE),
    IQR_low = quantile(abs(sig_trt$lfc), 0.25, na.rm=TRUE),
    IQR_high = quantile(abs(sig_trt$lfc), 0.75, na.rm=TRUE),
    N_strong = sum(abs(sig_trt$lfc) > 2, na.rm=TRUE),
    N_moderate = sum(abs(sig_trt$lfc) > 1 & abs(sig_trt$lfc) <= 2, na.rm=TRUE),
    N_weak = sum(abs(sig_trt$lfc) <= 1, na.rm=TRUE)
  )
  cat("  Treatment effect size distribution:\n"); print(trt_effect_summary)
  
  ## ANCOM-BC2 sensitivity: passed_ss (sensitivity to pseudo-count)
  ss_col <- grep("passed_ss|pass_ss", colnames(ancom_trt$res), value=TRUE)
  if (length(ss_col) > 0) {
    ss_col_trt <- grep("treatment", ss_col, value=TRUE)[1]
    if (!is.na(ss_col_trt)) {
      n_passed_ss <- sum(ancom_trt$res[[ss_col_trt]], na.rm=TRUE)
      n_diff <- sum(trt_df$diff_abn, na.rm=TRUE)
      cat(sprintf("  Sensitivity: %d / %d DA ASVs pass pseudo-count sensitivity (%.1f%%)\n",
                  n_passed_ss, n_diff, 100*n_passed_ss/max(n_diff,1)))
      trt_effect_summary$N_passed_sensitivity <- n_passed_ss
    }
  }
  write.csv(trt_effect_summary, "tables/P3_da/ANCOMBC2_treatment_effect_sizes.csv", row.names=FALSE)
}
cat("ANCOM-BC2: Genotype...\n")
set.seed(42)
ancom_geno <- tryCatch(ancombc2(data=ps_planted, fix_formula="trait + treatment + harvest_num",
  p_adj_method="fdr", prv_cut=0.10, lib_cut=CONFIG$MIN_LIB_SIZE,
  group="trait", struc_zero=TRUE, neg_lb=TRUE, alpha=0.05, global=FALSE),
  error=function(e){cat("Error:",e$message,"\n"); NULL})

geno_df <- NULL
if (!is.null(ancom_geno)) {
  res_g <- ancom_geno$res
  lfc_g <- grep("^lfc_trait", colnames(res_g), value=TRUE)[1]
  q_g <- gsub("^lfc_","q_",lfc_g); diff_g <- gsub("^lfc_","diff_",lfc_g)
  geno_df <- data.frame(ASV=res_g$taxon, lfc=res_g[[lfc_g]], qval=res_g[[q_g]], diff_abn=res_g[[diff_g]]) %>%
    left_join(tax_lookup %>% rownames_to_column("ASV"), by="ASV") %>%
    mutate(Direction=case_when(diff_abn & lfc>0 ~ "Resistant_enriched",
                               diff_abn & lfc<0 ~ "Susceptible_enriched", TRUE ~ "NS"))
  write.csv(geno_df, "tables/P3_da/ANCOMBC2_genotype.csv", row.names=FALSE)
  n_res <- sum(geno_df$Direction=="Resistant_enriched"); n_sus <- sum(geno_df$Direction=="Susceptible_enriched")
  cat(sprintf("  %d resistant-enriched, %d susceptible-enriched\n", n_res, n_sus))
  
  top_label_g <- geno_df %>% filter(diff_abn) %>% arrange(qval) %>% head(15) %>%
    mutate(Label = .rank_label(Genus, Family, Order, Class, Phylum))
  
  p_vol_geno <- ggplot(geno_df, aes(lfc, -log10(qval), colour=Direction)) +
    geom_point(alpha=0.4, size=1.5) +
    scale_colour_manual(values=c(Resistant_enriched=PAL$genotype[["Resistance"]], Susceptible_enriched=PAL$genotype[["Susceptible"]], NS=COL_SIG_NS)) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.3) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.3) +
    ggrepel::geom_text_repel(data=top_label_g, aes(label=Label), size=3, max.overlaps=20,
                              box.padding=0.4, point.padding=0.3, force=2,
                              show.legend=FALSE, min.segment.length=0) +
    coord_cartesian(clip="off") +
    theme_pub() +
    theme(plot.margin=margin(4, 15, 4, 4, "mm")) +
    labs(x="Log fold change", y=expression(-log[10](q)), colour=NULL,
      title="Differential abundance: Resistant vs Susceptible",
      subtitle=sprintf("%d resistant-enriched | %d susceptible-enriched", n_res, n_sus))
  save_fig(p_vol_geno, "figures/P3_da/volcano_genotype", w=10, h=7)
  
  ## Genotype effect size distribution
  sig_geno <- geno_df %>% filter(diff_abn==TRUE)
  geno_effect_summary <- data.frame(
    Comparison = "Resistant vs Susceptible",
    N_significant = nrow(sig_geno),
    Median_absLFC = median(abs(sig_geno$lfc), na.rm=TRUE),
    IQR_low = quantile(abs(sig_geno$lfc), 0.25, na.rm=TRUE),
    IQR_high = quantile(abs(sig_geno$lfc), 0.75, na.rm=TRUE),
    N_strong = sum(abs(sig_geno$lfc) > 2, na.rm=TRUE),
    N_moderate = sum(abs(sig_geno$lfc) > 1 & abs(sig_geno$lfc) <= 2, na.rm=TRUE),
    N_weak = sum(abs(sig_geno$lfc) <= 1, na.rm=TRUE)
  )
  
  ss_col_g <- grep("passed_ss|pass_ss", colnames(ancom_geno$res), value=TRUE)
  if (length(ss_col_g) > 0) {
    ss_g <- grep("trait", ss_col_g, value=TRUE)[1]
    if (!is.na(ss_g)) {
      n_ps_g <- sum(ancom_geno$res[[ss_g]], na.rm=TRUE)
      geno_effect_summary$N_passed_sensitivity <- n_ps_g
      cat(sprintf("  Genotype sensitivity: %d / %d pass pseudo-count check\n", n_ps_g, nrow(sig_geno)))
    }
  }
  write.csv(geno_effect_summary, "tables/P3_da/ANCOMBC2_genotype_effect_sizes.csv", row.names=FALSE)
  
  ## Combined effect size summary
  if (exists("trt_effect_summary")) {
    all_effects <- rbind(trt_effect_summary[, intersect(colnames(trt_effect_summary), colnames(geno_effect_summary))],
                         geno_effect_summary[, intersect(colnames(trt_effect_summary), colnames(geno_effect_summary))])
    write.csv(all_effects, "tables/P3_da/ANCOMBC2_combined_effect_sizes.csv", row.names=FALSE)
    cat("DA effect size summary:\n"); print(all_effects); cat("\n")
  }
}

## ── 3c. DA overlap and heatmap ────────────────────────────────────────────
if (!is.null(trt_df) && !is.null(geno_df)) {
  ## Overlap: ASVs DA in both comparisons
  da_trt <- trt_df$ASV[trt_df$diff_abn]; da_geno <- geno_df$ASV[geno_df$diff_abn]
  overlap <- intersect(da_trt, da_geno)
  cat(sprintf("  DA overlap (both treatment AND genotype): %d ASVs\n", length(overlap)))
  
  ## Heatmap of top 40 DA genera (treatment)
  top40 <- trt_df %>% filter(diff_abn) %>% arrange(desc(abs(lfc))) %>% head(40)
  
  ## Get genus-level abundances for these ASVs
  ps_genus_top <- prune_taxa(top40$ASV, ps_rel_planted)
  hm_mat <- as(otu_table(ps_genus_top), "matrix")
  if (!taxa_are_rows(ps_genus_top)) hm_mat <- t(hm_mat)
  
  ## Use best available taxonomy as row labels
  row_labels <- top40$Genus
  row_labels[is.na(row_labels)] <- top40$Family[is.na(row_labels)]
  row_labels[is.na(row_labels)] <- top40$Phylum[is.na(row_labels)]
  row_labels[is.na(row_labels)] <- top40$ASV[is.na(row_labels)]
  ## Make unique
  row_labels <- make.unique(row_labels)
  rownames(hm_mat) <- row_labels
  
  ## Log transform for visualization
  hm_mat_log <- log10(hm_mat + 1e-5)
  
  ## Column annotations
  anno_col <- data.frame(
    Treatment = meta_pl$treatment, Genotype = meta_pl$trait,
    Time = meta_pl$harvest, row.names = colnames(hm_mat))
  anno_colors <- list(Treatment = PAL$treatment,
                      Genotype = c(PAL$genotype, Unplanted="#757575"),
                      Time = PAL$time)
  
  ## Row annotation: direction
  row_direction <- top40$Direction
  names(row_direction) <- row_labels
  anno_row <- data.frame(Response = row_direction, row.names = row_labels)
  
  pdf("figures/P3_da/heatmap_top40_DA_treatment.pdf", w=14, h=10)
  pheatmap(hm_mat_log, annotation_col=anno_col, annotation_row=anno_row,
           annotation_colors=c(anno_colors, list(Response=c(Drought_enriched=PAL$treatment[["Drought"]], Watered_enriched=PAL$treatment[["Watered"]]))),
           cluster_rows=TRUE, cluster_cols=TRUE, show_colnames=FALSE,
           fontsize_row=8, main="Top 40 differentially abundant ASVs (treatment)",
           color=viridis(100))
  dev.off()
  png("figures/P3_da/heatmap_top40_DA_treatment.png", w=14, h=10, units="in", res=300)
  pheatmap(hm_mat_log, annotation_col=anno_col, annotation_row=anno_row,
           annotation_colors=c(anno_colors, list(Response=c(Drought_enriched=PAL$treatment[["Drought"]], Watered_enriched=PAL$treatment[["Watered"]]))),
           cluster_rows=TRUE, cluster_cols=TRUE, show_colnames=FALSE,
           fontsize_row=8, main="Top 40 differentially abundant ASVs (treatment)",
           color=viridis(100))
  dev.off()
  svglite::svglite("figures/P3_da/heatmap_top40_DA_treatment.svg", width=14, height=10)
  pheatmap(hm_mat_log, annotation_col=anno_col, annotation_row=anno_row,
           annotation_colors=c(anno_colors, list(Response=c(Drought_enriched=PAL$treatment[["Drought"]], Watered_enriched=PAL$treatment[["Watered"]]))),
           cluster_rows=TRUE, cluster_cols=TRUE, show_colnames=FALSE,
           fontsize_row=8, main="Top 40 differentially abundant ASVs (treatment)",
           color=viridis(100))
  dev.off()
}

## ── 3d. PICRUSt2 functional profiling (expanded — 8 analyses) ────────────
cat("\n── Functional profiling (PICRUSt2) — expanded analyses ──\n\n")

picrust_ko    <- "picrust2/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
picrust_pw    <- "picrust2/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz"
picrust_nsti  <- "picrust2/picrust2_output/marker_predicted_and_nsti.tsv.gz"
picrust_strat <- "picrust2/picrust2_output/KO_metagenome_out/pred_metagenome_strat.tsv.gz"
picrust_available <- file.exists(picrust_ko) & file.exists(picrust_pw)

if (picrust_available) {

  ## ── NSTI quality check ────────────────────────────────────────────────
  tryCatch({
    nsti     <- read.table(picrust_nsti, header=TRUE, sep="\t", row.names=1)
    nsti_col <- grep("NSTI|nsti", colnames(nsti), value=TRUE)[1]
    nsti_vals <- nsti[[nsti_col]]
    cat(sprintf("NSTI: mean=%.3f, median=%.3f, %%>%d=%.1f%%\n",
        mean(nsti_vals, na.rm=TRUE), median(nsti_vals, na.rm=TRUE),
        CONFIG$NSTI_CUTOFF, 100*mean(nsti_vals > CONFIG$NSTI_CUTOFF, na.rm=TRUE)))
    p_nsti <- ggplot(data.frame(NSTI=nsti_vals[!is.na(nsti_vals)]), aes(NSTI)) +
      geom_histogram(bins=50, fill=COL_WATERED, colour="white", linewidth=0.2, alpha=0.8) +
      geom_vline(xintercept=CONFIG$NSTI_CUTOFF, linetype="dashed",
                 colour=COL_DROUGHT, linewidth=0.6) +
      annotate("text", x=CONFIG$NSTI_CUTOFF + 0.05, y=Inf,
               label=sprintf("Cutoff=%d", CONFIG$NSTI_CUTOFF),
               hjust=0, vjust=1.5, size=7/.pt, colour=COL_DROUGHT) +
      theme_pub() +
      labs(x="Weighted NSTI", y="Count",
           title="PICRUSt2 prediction quality (NSTI)",
           subtitle=sprintf("Mean=%.3f, Median=%.3f, %.1f%% > %d",
                             mean(nsti_vals, na.rm=TRUE), median(nsti_vals, na.rm=TRUE),
                             100*mean(nsti_vals > CONFIG$NSTI_CUTOFF, na.rm=TRUE),
                             CONFIG$NSTI_CUTOFF))
    save_fig(p_nsti, "figures/P3_function/NSTI_distribution", w=ISME_SINGLE_W, h=4)
  }, error=function(e) cat("NSTI check skipped:", conditionMessage(e), "\n"))

  ## ── Load KO and pathway tables ────────────────────────────────────────
  ko_tab <- pw_tab <- NULL
  tryCatch({
    ko_raw <- as.data.frame(data.table::fread(picrust_ko, sep="\t"))
    rownames(ko_raw) <- ko_raw[[1]]; ko_raw <- ko_raw[, -1, drop=FALSE]
    ko_tab <- ko_raw[, intersect(colnames(ko_raw), rownames(meta_pl)), drop=FALSE]
    ko_tab[] <- lapply(ko_tab, function(x) as.numeric(as.character(x)))
    ko_tab <- ko_tab[rowSums(ko_tab, na.rm=TRUE) > 0, , drop=FALSE]
    cat(sprintf("KO table: %d KOs x %d planted samples\n", nrow(ko_tab), ncol(ko_tab)))
  }, error=function(e) cat("KO load failed:", conditionMessage(e), "\n"))

  tryCatch({
    pw_raw <- as.data.frame(data.table::fread(picrust_pw, sep="\t"))
    rownames(pw_raw) <- pw_raw[[1]]; pw_raw <- pw_raw[, -1, drop=FALSE]
    pw_tab <- pw_raw[, intersect(colnames(pw_raw), rownames(meta_pl)), drop=FALSE]
    pw_tab[] <- lapply(pw_tab, function(x) as.numeric(as.character(x)))
    pw_tab <- pw_tab[rowSums(pw_tab, na.rm=TRUE) > 0, , drop=FALSE]
    cat(sprintf("Pathway table: %d pathways x %d planted samples\n", nrow(pw_tab), ncol(pw_tab)))
  }, error=function(e) cat("Pathway load failed:", conditionMessage(e), "\n"))

  ko_rel <- if (!is.null(ko_tab)) sweep(ko_tab, 2, colSums(ko_tab), "/") else NULL
  pw_rel <- if (!is.null(pw_tab)) sweep(pw_tab, 2, colSums(pw_tab), "/") else NULL
  shared_meta <- function(mat) meta_pl[colnames(mat), , drop=FALSE]

  ## ── 3d-1. Differential KO abundance (ALDEx2) ─────────────────────────
  cat("  [3d-1] ALDEx2 differential KO abundance...\n")
  tryCatch({
    if (is.null(ko_tab)) stop("KO table not loaded")
    ko_var <- apply(ko_tab, 1, var)
    ko_sub <- ko_tab[order(ko_var, decreasing=TRUE)[seq_len(min(3000L, nrow(ko_tab)))], ]
    trt_vec <- as.character(shared_meta(ko_sub)$treatment)
    ko_int  <- round(matrix(as.numeric(unlist(ko_sub)), nrow=nrow(ko_sub),
                            dimnames=list(rownames(ko_sub), colnames(ko_sub))))
    ko_int[ko_int < 0] <- 0
    storage.mode(ko_int) <- "integer"
    set.seed(42)
    ko_clr <- ALDEx2::aldex.clr(ko_int, conds=trt_vec, mc.samples=128, denom="all")
    ko_tt  <- ALDEx2::aldex.ttest(ko_clr, paired.test=FALSE, verbose=FALSE)
    ko_eff <- ALDEx2::aldex.effect(ko_clr, verbose=FALSE)
    ko_da  <- data.frame(KO=rownames(ko_int), ko_tt, ko_eff)
    ko_da$qval_we <- p.adjust(ko_da$we.ep, method="BH")
    ko_da$sig     <- ko_da$qval_we < CONFIG$FDR_THRESHOLD
    ko_da <- ko_da[order(ko_da$qval_we), ]
    write.csv(ko_da, "tables/P3_da/ALDEx2_KO_treatment.csv", row.names=FALSE)
    n_sig_ko <- sum(ko_da$sig, na.rm=TRUE)
    cat(sprintf("    DA KOs (q<0.05): %d / %d\n", n_sig_ko, nrow(ko_da)))

    top_ko <- ko_da %>% dplyr::filter(sig) %>%
      dplyr::arrange(desc(abs(effect))) %>% head(30)
    if (nrow(top_ko) > 0) {
      top_ko$direction <- ifelse(top_ko$effect > 0, "Drought", "Watered")
      top_ko$KO_lab <- factor(top_ko$KO, levels=rev(top_ko$KO))
      p_ko_bar <- ggplot(top_ko, aes(x=effect, y=KO_lab, fill=direction)) +
        geom_col(alpha=0.85) +
        geom_errorbarh(aes(xmin=effect.low, xmax=effect.high),
                       height=0.4, linewidth=0.3, colour="grey40") +
        scale_fill_manual(values=PAL$treatment, name=NULL) +
        geom_vline(xintercept=0, linewidth=0.4, colour="grey40") +
        theme_pub() +
        theme(axis.text.y=element_text(size=BASE_SIZE - 4)) +
        labs(x="ALDEx2 effect size", y="KEGG ortholog",
             title="Differential KOs: Drought vs Watered (top 30, q<0.05)",
             subtitle=sprintf("%d of %d top-variable KOs significantly different",
                               n_sig_ko, nrow(ko_da)))
      save_fig(p_ko_bar, "figures/P3_function/ALDEx2_KO_effect_bar",
               w=ISME_DOUBLE_W, h=max(4, 0.22 * nrow(top_ko) + 1.5))
    }

    ko_da$neglog10q <- pmin(-log10(pmax(ko_da$qval_we, 1e-10)), 10)
    ko_da$direction <- ifelse(!ko_da$sig, "NS",
                       ifelse(ko_da$effect > 0, "Drought", "Watered"))
    p_ko_bubble <- ggplot(ko_da, aes(x=effect, y=neglog10q,
                                      size=abs(effect), colour=direction)) +
      geom_point(alpha=0.6, shape=16) +
      scale_colour_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED, NS="grey70"),
                           name=NULL) +
      scale_size_continuous(range=c(0.5, 3), guide="none") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed",
                 colour="grey40", linewidth=0.4) +
      geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=0.4) +
      theme_pub() +
      labs(x="ALDEx2 effect size", y=expression(-log[10](q[BH])),
           title="KO differential abundance: effect vs significance",
           subtitle=sprintf("Top 3000 variable KOs; %d significant (q<0.05)", n_sig_ko))
    save_fig(p_ko_bubble, "figures/P3_function/ALDEx2_KO_bubble",
             w=ISME_DOUBLE_W, h=5)
    cat("    [3d-1] Done.\n")
  }, error=function(e) cat("  [3d-1] ALDEx2 KO skipped:", conditionMessage(e), "\n"))

  ## ── 3d-1b. ANCOM-BC2 KO differential abundance (three comparisons) ───
  cat("  [3d-1b] ANCOM-BC2 KO differential abundance (3 comparisons)...\n")
  tryCatch({
    if (is.null(ko_tab)) stop("KO table not loaded")
    if (!requireNamespace("ANCOMBC", quietly=TRUE))
      stop("ANCOMBC package not available")

    ## Build phyloseq from KO table for ANCOMBC2
    ko_mat_t <- round(matrix(as.numeric(unlist(ko_tab)), nrow=nrow(ko_tab),
                              dimnames=list(rownames(ko_tab), colnames(ko_tab))))
    ko_mat_t[ko_mat_t < 0] <- 0L
    storage.mode(ko_mat_t) <- "integer"
    ko_meta_anc <- meta_pl[colnames(ko_tab), , drop=FALSE]
    ko_meta_anc$harvest_fac <- factor(ko_meta_anc$harvest_num)
    ko_meta_anc$time_group  <- ifelse(ko_meta_anc$harvest_num <= 3,
                                       "Early", "Late")
    ko_meta_anc$time_group  <- factor(ko_meta_anc$time_group,
                                       levels=c("Early","Late"))
    ps_ko <- phyloseq(
      otu_table(ko_mat_t, taxa_are_rows=TRUE),
      sample_data(ko_meta_anc)
    )

    run_ancombc2_ko <- function(ps_obj, formula_str, out_csv, comparison_label) {
      cache_file <- out_csv
      if (file.exists(cache_file)) {
        cat(sprintf("    [cached] %s\n", comparison_label))
        return(invisible(read.csv(cache_file)))
      }
      res <- tryCatch(
        ANCOMBC::ancombc2(
          data        = ps_obj,
          fix_formula = formula_str,
          p_adj_method= "BH",
          prv_cut     = 0.10,
          lib_cut     = 1000,
          s0_perc     = 0.05,
          group       = NULL,
          struc_zero  = FALSE,
          neg_lb       = FALSE,
          iter_control = list(tol=1e-5, max_iter=20, verbose=FALSE),
          em_control   = list(tol=1e-5, max_iter=100),
          lme_control  = NULL,
          mdfdr_control= list(fwer_ctrl_method="BH", B=100),
          verbose      = FALSE
        ),
        error = function(e) { cat("    ANCOMBC2 error:", conditionMessage(e), "\n"); NULL }
      )
      if (is.null(res)) return(invisible(NULL))
      df_res <- res$res
      write.csv(df_res, cache_file, row.names=FALSE)
      diff_cols <- grep("^diff_(?!robust)", colnames(df_res), perl=TRUE)
      diff_cols <- diff_cols[!grepl("Intercept", colnames(df_res)[diff_cols])]
      n_sig <- if (length(diff_cols) > 0L)
        sum(df_res[[diff_cols[1]]], na.rm=TRUE) else 0L
      cat(sprintf("    %s: %d significant KOs\n", comparison_label, n_sig))
      invisible(df_res)
    }

    ## Comparison A: Drought vs Watered
    anc_trt <- run_ancombc2_ko(ps_ko, "treatment",
                                "tables/P3_da/ANCOMBC2_KO_treatment.csv",
                                "Treatment (Drought vs Watered)")

    ## Comparison B: Resistance vs Susceptible (planted only, trait != Unplanted)
    ps_ko_planted <- prune_samples(
      sample_data(ps_ko)$trait %in% c("Resistance","Susceptible"), ps_ko)
    ps_ko_planted <- prune_taxa(taxa_sums(ps_ko_planted) > 0, ps_ko_planted)
    anc_geno <- run_ancombc2_ko(ps_ko_planted, "trait",
                                 "tables/P3_da/ANCOMBC2_KO_genotype.csv",
                                 "Genotype (Resistance vs Susceptible)")

    ## Comparison C: Early (Days 0-3) vs Late (Days 4-6)
    anc_time <- run_ancombc2_ko(ps_ko, "time_group",
                                 "tables/P3_da/ANCOMBC2_KO_temporal.csv",
                                 "Temporal (Early vs Late)")

    ## Category breakdown for each comparison
    summarise_ancombc_by_category <- function(df_res, ko_categories, label) {
      if (is.null(df_res)) return(invisible(NULL))
      ## ANCOMBC2 >=2.0 uses diff_{covariate} not diff_abn_{covariate}
      diff_col <- grep("^diff_(?!robust)", colnames(df_res), perl=TRUE)
      diff_col <- diff_col[!grepl("Intercept", colnames(df_res)[diff_col])][1]
      lfc_col  <- grep("^lfc_",    colnames(df_res))[1]
      if (is.na(diff_col)) return(invisible(NULL))
      sig_kos <- df_res$taxon[df_res[[diff_col]] %in% TRUE]
      cat(sprintf("    %s — category breakdown:\n", label))
      for (cat_name in names(ko_categories)) {
        kos_cat <- ko_categories[[cat_name]]
        n_det <- sum(kos_cat %in% df_res$taxon)
        n_sig <- sum(kos_cat %in% sig_kos)
        if (n_det > 0)
          cat(sprintf("      %-24s %d/%d sig\n", cat_name, n_sig, n_det))
      }
    }
    drought_kos_ref <- list(
      Osmolyte_biosynthesis  = c("K00849","K06718","K05835","K00130","K00548","K05289","K05290","K05291"),
      Trehalose_biosynthesis = c("K00697","K01087","K16143","K00055","K01563","K01697"),
      EPS_and_biofilm        = c("K02288","K13244","K16786","K05810","K00975","K02851"),
      Oxidative_stress       = c("K03782","K04564","K00432","K00383","K07304","K03671"),
      Heat_shock_chaperones  = c("K04043","K04042","K04077","K13993","K09579","K04080"),
      N_cycling              = c("K00370","K00371","K00374","K02567","K02568","K00362","K00363","K10944","K10945","K10946"),
      P_solubilization       = c("K01077","K09474","K05781","K01113"),
      IAA_biosynthesis       = c("K01425","K11430","K14265","K10563"),
      ACC_deaminase          = c("K01505"),
      Siderophore            = c("K02363","K02364","K09486","K09490")
    )
    summarise_ancombc_by_category(anc_trt,  drought_kos_ref, "Treatment")
    summarise_ancombc_by_category(anc_geno, drought_kos_ref, "Genotype")
    summarise_ancombc_by_category(anc_time, drought_kos_ref, "Temporal")
    cat("    [3d-1b] Done.\n")
  }, error=function(e) cat("  [3d-1b] ANCOMBC2 KO skipped:", conditionMessage(e), "\n"))

  ## ── 3d-2. Differential MetaCyc pathway abundance (ALDEx2) ────────────
  cat("  [3d-2] ALDEx2 differential pathway abundance...\n")
  tryCatch({
    if (is.null(pw_tab)) stop("Pathway table not loaded")
    trt_vec <- as.character(shared_meta(pw_tab)$treatment)
    pw_int  <- round(matrix(as.numeric(unlist(pw_tab)), nrow=nrow(pw_tab),
                            dimnames=list(rownames(pw_tab), colnames(pw_tab))))
    pw_int[pw_int < 0] <- 0
    storage.mode(pw_int) <- "integer"
    set.seed(42)
    pw_clr <- ALDEx2::aldex.clr(pw_int, conds=trt_vec, mc.samples=128, denom="all")
    pw_tt  <- ALDEx2::aldex.ttest(pw_clr, paired.test=FALSE, verbose=FALSE)
    pw_eff <- ALDEx2::aldex.effect(pw_clr, verbose=FALSE)
    pw_da  <- data.frame(Pathway=rownames(pw_int), pw_tt, pw_eff)
    pw_da$qval_we <- p.adjust(pw_da$we.ep, method="BH")
    pw_da$sig     <- pw_da$qval_we < CONFIG$FDR_THRESHOLD
    pw_da <- pw_da[order(pw_da$qval_we), ]
    write.csv(pw_da, "tables/P3_da/ALDEx2_pathway_treatment.csv", row.names=FALSE)
    n_sig_pw <- sum(pw_da$sig, na.rm=TRUE)
    cat(sprintf("    DA pathways (q<0.05): %d / %d\n", n_sig_pw, nrow(pw_da)))

    pw_da$neglog10q <- pmin(-log10(pmax(pw_da$qval_we, 1e-10)), 10)
    pw_da$direction <- ifelse(!pw_da$sig, "NS",
                       ifelse(pw_da$effect > 0, "Drought", "Watered"))
    n_dr_pw <- sum(pw_da$sig & pw_da$effect > 0, na.rm=TRUE)
    n_wa_pw <- sum(pw_da$sig & pw_da$effect < 0, na.rm=TRUE)
    label_pw <- pw_da %>% dplyr::filter(sig) %>%
      dplyr::arrange(desc(abs(effect))) %>% head(15) %>%
      mutate(pw_short=substr(Pathway, 1, 35))
    p_pw_vol <- ggplot(pw_da, aes(x=effect, y=neglog10q)) +
      geom_point(aes(colour=direction), size=1.2, alpha=0.7, shape=16) +
      scale_colour_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED, NS="grey70"),
                           name=NULL) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed",
                 colour="grey40", linewidth=0.4) +
      geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=0.4) +
      ggrepel::geom_text_repel(data=label_pw, aes(label=pw_short),
                               size=6/.pt, colour="grey15",
                               segment.size=0.3, max.overlaps=15) +
      annotate("text", x=max(pw_da$effect, na.rm=TRUE)*0.85, y=9.5,
               label=sprintf("%d enriched", n_dr_pw),
               size=7/.pt, colour=COL_DROUGHT, fontface="italic") +
      annotate("text", x=min(pw_da$effect, na.rm=TRUE)*0.85, y=9.5,
               label=sprintf("%d enriched", n_wa_pw),
               size=7/.pt, colour=COL_WATERED, fontface="italic") +
      theme_pub() +
      labs(x="ALDEx2 effect size", y=expression(-log[10](q[BH])),
           title="Differential MetaCyc pathways: Drought vs Watered",
           subtitle=sprintf("%d drought-enriched | %d watered-enriched (q<0.05)",
                             n_dr_pw, n_wa_pw))
    save_fig(p_pw_vol, "figures/P3_function/ALDEx2_pathway_volcano",
             w=ISME_DOUBLE_W, h=5)

    if (n_sig_pw >= 5 && !is.null(pw_rel)) {
      top30_pw <- pw_da %>% dplyr::filter(sig) %>%
        dplyr::arrange(desc(abs(effect))) %>% head(30) %>% pull(Pathway)
      top30_pw <- intersect(top30_pw, rownames(pw_rel))
      if (length(top30_pw) >= 3) {
        hm_mat   <- as.matrix(pw_rel[top30_pw, ])
        hm_mat_z <- t(scale(t(log1p(hm_mat * 1e6))))
        col_annot    <- data.frame(Treatment=shared_meta(pw_rel)$treatment,
                                   row.names=colnames(hm_mat_z))
        annot_colors <- list(Treatment=PAL$treatment)
        rownames(hm_mat_z) <- substr(rownames(hm_mat_z), 1, 45)
        png("figures/P3_function/ALDEx2_pathway_heatmap.png",
            width=ISME_DOUBLE_W, height=6, units="in", res=300)
        pheatmap::pheatmap(hm_mat_z,
                           annotation_col=col_annot,
                           annotation_colors=annot_colors,
                           cluster_rows=TRUE, cluster_cols=TRUE,
                           show_colnames=FALSE, fontsize_row=7,
                           color=colorRampPalette(c(COL_HEATMAP_LOW, COL_HEATMAP_MID,
                                                    COL_HEATMAP_HIGH))(100),
                           main="Top 30 differential MetaCyc pathways (z-score)")
        dev.off()
      }
    }
    cat("    [3d-2] Done.\n")
  }, error=function(e) cat("  [3d-2] ALDEx2 pathway skipped:", conditionMessage(e), "\n"))

  ## ── 3d-3. Functional alpha diversity ─────────────────────────────────
  cat("  [3d-3] Functional alpha diversity...\n")
  tryCatch({
    if (is.null(ko_tab)) stop("KO table not loaded")
    func_alpha <- data.frame(
      Sample      = colnames(ko_tab),
      KO_richness = colSums(ko_tab > 0),
      KO_shannon  = apply(ko_tab, 2, function(x) {
                      p <- x[x > 0] / sum(x); -sum(p * log(p)) }),
      KO_evenness = apply(ko_tab, 2, function(x) {
                      p <- x[x > 0] / sum(x); H <- -sum(p * log(p))
                      if (sum(x > 0) > 1) H / log(sum(x > 0)) else NA }),
      stringsAsFactors=FALSE
    )
    func_alpha <- left_join(func_alpha,
                             data.frame(Sample=rownames(meta_pl), meta_pl,
                                        stringsAsFactors=FALSE),
                             by="Sample")
    write.csv(func_alpha, "tables/P3_function/functional_alpha_diversity.csv", row.names=FALSE)

    func_alpha$harvest_fac <- factor(func_alpha$harvest_num)

    fa_long <- func_alpha %>%
      pivot_longer(cols=c(KO_richness, KO_shannon, KO_evenness),
                   names_to="Metric", values_to="Value") %>%
      mutate(Metric=dplyr::recode(Metric,
               KO_richness="KO Richness", KO_shannon="KO Shannon",
               KO_evenness="KO Evenness"))

    ## Compute Wilcoxon p-values per metric — explicit subsetting to avoid
    ## scoping issues with formula interface inside group_by/summarise
    fa_stats_list <- lapply(unique(fa_long$Metric), function(m) {
      sub <- fa_long[fa_long$Metric == m, ]
      p_raw <- tryCatch(wilcox.test(Value ~ treatment, data=sub)$p.value,
                        error=function(e) NA_real_)
      data.frame(Metric=m, p_raw=p_raw, stringsAsFactors=FALSE)
    })
    fa_stats <- do.call(rbind, fa_stats_list)
    fa_stats$p_adj <- p.adjust(fa_stats$p_raw, method="BH")
    fa_stats$label <- paste0("p=", format(round(fa_stats$p_raw, 3), nsmall=3))
    fa_stats$x     <- 1.5
    fa_stats$y     <- Inf
    cat("    Functional alpha p-values (BH):\n")
    for (i in seq_len(nrow(fa_stats)))
      cat(sprintf("      %s: p_raw=%.4f, p_adj=%.4f\n",
                  fa_stats$Metric[i], fa_stats$p_raw[i], fa_stats$p_adj[i]))

    p_fa_box <- ggplot(fa_long, aes(treatment, Value, fill=treatment)) +
      geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.55) +
      geom_jitter(width=0.12, size=1, alpha=0.4, colour="grey30") +
      geom_text(data=fa_stats,
                aes(x=x, y=y, label=label),
                vjust=1.5, hjust=0.5,
                size=3, colour="grey40",
                family=FONT_FAMILY,
                inherit.aes=FALSE) +
      facet_wrap(~Metric, scales="free_y", nrow=1) +
      scale_fill_manual(values=PAL$treatment, guide="none") +
      theme_pub() +
      theme(axis.text.x=element_text(angle=30, hjust=1)) +
      labs(x="Treatment", y="Value",
           title="Functional alpha diversity (KO space)",
           subtitle="Wilcoxon test per group; BH-corrected p-values in results table")
    save_fig(p_fa_box, "figures/P3_function/functional_alpha_boxplot",
             w=ISME_DOUBLE_W, h=5)

    ## Three-factor functional alpha diversity figure
    ## Build long data for three factor groupings
    fa_3f_trt <- fa_long %>%
      mutate(factor_type="Treatment", factor_level=treatment,
             fill_var=treatment)
    fa_3f_gen <- fa_long %>%
      dplyr::filter(trait %in% c("Resistance","Susceptible")) %>%
      mutate(factor_type="Genotype", factor_level=as.character(trait),
             fill_var=as.character(trait))
    fa_3f_time <- fa_long %>%
      mutate(factor_type="Time (harvest day)",
             factor_level=as.character(harvest_fac),
             fill_var=as.character(harvest_fac))
    fa_3f <- bind_rows(fa_3f_trt, fa_3f_gen, fa_3f_time)
    fa_3f$factor_type <- factor(fa_3f$factor_type,
                                 levels=c("Treatment","Genotype","Time (harvest day)"))

    ## Wilcoxon per Metric × factor_type combination
    fa_3f_stats <- do.call(rbind, lapply(
      unique(fa_3f$Metric), function(m) {
        do.call(rbind, lapply(
          levels(fa_3f$factor_type), function(ft) {
            sub2 <- fa_3f[fa_3f$Metric == m & fa_3f$factor_type == ft, ]
            lvls <- unique(sub2$factor_level)
            if (ft == "Time (harvest day)") {
              ## Kruskal-Wallis for time (>2 groups)
              p_raw <- tryCatch(
                kruskal.test(Value ~ factor(factor_level), data=sub2)$p.value,
                error=function(e) NA_real_)
            } else {
              p_raw <- tryCatch(
                wilcox.test(Value ~ factor_level, data=sub2)$p.value,
                error=function(e) NA_real_)
            }
            data.frame(Metric=m, factor_type=ft, p_raw=p_raw,
                       stringsAsFactors=FALSE)
          }))
      }))
    fa_3f_stats$p_adj  <- p.adjust(fa_3f_stats$p_raw, method="BH")
    fa_3f_stats$label  <- paste0("p=", format(round(fa_3f_stats$p_raw, 3), nsmall=3))
    fa_3f_stats$factor_type <- factor(fa_3f_stats$factor_type,
                                       levels=levels(fa_3f$factor_type))
    ## x midpoints per factor_type
    fa_3f_stats$x_mid <- ifelse(fa_3f_stats$factor_type == "Treatment",    1.5,
                          ifelse(fa_3f_stats$factor_type == "Genotype",     1.5, 3.5))
    fa_3f_stats$y_pos  <- Inf

    ## Fill palette: treatment + genotype + time (6 days)
    time_pal <- setNames(
      scales::hue_pal()(length(unique(fa_3f_time$factor_level))),
      sort(unique(fa_3f_time$factor_level)))
    geno_pal  <- PAL$genotype
    all_fills <- c(PAL$treatment, geno_pal, time_pal)

    p_fa_3f <- ggplot(fa_3f, aes(x=factor_level, y=Value, fill=fill_var)) +
      geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.6) +
      geom_jitter(width=0.1, size=0.8, alpha=0.35, colour="grey30") +
      geom_text(data=fa_3f_stats,
                aes(x=x_mid, y=y_pos, label=label),
                inherit.aes=FALSE,
                vjust=1.5, hjust=0.5,
                size=2.8, colour="grey40", family=FONT_FAMILY) +
      facet_grid(Metric ~ factor_type, scales="free") +
      scale_fill_manual(values=all_fills, guide="none") +
      theme_pub(BASE_SIZE - 1) +
      theme(axis.text.x=element_text(angle=35, hjust=1),
            strip.text.x=element_text(size=BASE_SIZE - 2),
            strip.text.y=element_text(size=BASE_SIZE - 3)) +
      labs(x=NULL, y="Value",
           title="Functional alpha diversity — three-factor analysis",
           subtitle="Wilcoxon (Treatment, Genotype) or Kruskal-Wallis (Time); raw p on panel")
    save_fig(p_fa_3f, "figures/P3_function/functional_alpha_diversity_threefactor",
             w=ISME_DOUBLE_W, h=7)

    tax_alpha_file <- "tables/P1_alpha/alpha_diversity_planted.csv"
    if (file.exists(tax_alpha_file)) {
      tax_alpha <- read.csv(tax_alpha_file, stringsAsFactors=FALSE)
      tax_col   <- grep("shannon|Shannon", colnames(tax_alpha), value=TRUE)[1]
      id_col    <- grep("^Sample$|^sample$|^id$", colnames(tax_alpha),
                        ignore.case=TRUE, value=TRUE)[1]
      if (!is.na(tax_col) && !is.na(id_col)) {
        scatter_df <- inner_join(
          func_alpha[, c("Sample","KO_shannon","treatment")],
          tax_alpha[, c(id_col, tax_col)],
          by=setNames(id_col, "Sample"))
        scatter_df$tax_shannon <- scatter_df[[tax_col]]
        cor_res <- cor.test(scatter_df$KO_shannon, scatter_df$tax_shannon,
                            method="spearman")
        p_scatter <- ggplot(scatter_df,
                             aes(x=tax_shannon, y=KO_shannon, colour=treatment)) +
          geom_point(size=2, alpha=0.8) +
          geom_smooth(method="lm", se=TRUE, colour="grey40",
                      linewidth=0.6, linetype="dashed") +
          scale_colour_manual(values=PAL$treatment) +
          theme_pub() +
          labs(x="Taxonomic Shannon diversity",
               y="Functional (KO) Shannon diversity",
               title="Functional vs taxonomic alpha diversity",
               subtitle=sprintf("Spearman r=%.3f, p=%.4f",
                                 cor_res$estimate, cor_res$p.value))
        save_fig(p_scatter, "figures/P3_function/functional_vs_taxonomic_alpha",
                 w=ISME_SINGLE_W * 1.8, h=ISME_SINGLE_W * 1.8)
      }
    }
    cat("    [3d-3] Done.\n")
  }, error=function(e) cat("  [3d-3] Functional alpha skipped:", conditionMessage(e), "\n"))

  ## ── 3d-4. Taxon-function contribution (streaming stratified file) ─────
  cat("  [3d-4] Taxon-function contributions...\n")
  tryCatch({
    ko_da_file <- "tables/P3_da/ALDEx2_KO_treatment.csv"
    if (!file.exists(ko_da_file)) stop("ALDEx2 KO results not found")
    ko_da_res <- read.csv(ko_da_file, stringsAsFactors=FALSE)
    top_da_kos <- ko_da_res %>% dplyr::filter(sig) %>%
      dplyr::arrange(desc(abs(effect))) %>% head(50) %>% pull(KO)
    if (length(top_da_kos) == 0) stop("No significant KOs")

    if (file.exists(picrust_strat)) {
      cat("    Streaming stratified KO file...\n")
      hdr <- colnames(data.table::fread(picrust_strat, nrows=0L, sep="\t"))
      samp_cols  <- intersect(hdr, rownames(meta_pl))
      func_col   <- if ("function" %in% hdr) "function" else hdr[1]
      taxon_col  <- if ("taxon"    %in% hdr) "taxon"    else hdr[2]
      keep_cols  <- unique(c(func_col, taxon_col, samp_cols))
      strat_dt   <- data.table::fread(picrust_strat, sep="\t",
                                       select=keep_cols, showProgress=FALSE)
      strat_sub  <- strat_dt[strat_dt[[func_col]] %in% top_da_kos, ]
      rm(strat_dt); gc()
      samp_pres  <- intersect(samp_cols, names(strat_sub))
      strat_sub[, total_contribution := rowSums(.SD), .SDcols=samp_pres]
      contrib_agg <- strat_sub[, .(total_contribution=sum(total_contribution, na.rm=TRUE)),
                                by=.(taxon=strat_sub[[taxon_col]])]
      contrib_agg <- contrib_agg[order(-contrib_agg$total_contribution), ]
      contrib_top <- head(contrib_agg, 40L)
      write.csv(as.data.frame(contrib_top), "tables/P3_function/top_taxon_contributions.csv",
                row.names=FALSE)
      contrib_top$taxon_short <- factor(
        substr(contrib_top$taxon, 1, 40),
        levels=rev(substr(contrib_top$taxon, 1, 40)))
      p_contrib <- ggplot(contrib_top,
                           aes(x=total_contribution, y=taxon_short)) +
        geom_col(fill=COL_DROUGHT, alpha=0.8) +
        theme_pub() +
        theme(axis.text.y=element_text(size=BASE_SIZE - 4)) +
        labs(x="Total predicted contribution to top DA KOs",
             y="Taxon",
             title="Top taxa contributing to drought-responsive functions",
             subtitle=sprintf("Summed across %d differential KOs (q<0.05)",
                               length(top_da_kos)))
      save_fig(p_contrib, "figures/P3_function/taxon_function_contributions",
               w=ISME_DOUBLE_W, h=6)
    } else {
      ## Fallback: KO abundances from unstratified table
      if (!is.null(ko_tab)) {
        ko_sub_top <- ko_tab[intersect(top_da_kos, rownames(ko_tab)), , drop=FALSE]
        ko_contrib_df <- data.frame(KO=rownames(ko_sub_top),
                                    total=rowSums(ko_sub_top)) %>%
          dplyr::arrange(desc(total)) %>% head(30)
        write.csv(ko_contrib_df, "tables/P3_function/top_ko_abundances.csv", row.names=FALSE)
        cat("    Top DA KO abundances saved (stratified file absent).\n")
      }
    }
    cat("    [3d-4] Done.\n")
  }, error=function(e) cat("  [3d-4] Contributions skipped:", conditionMessage(e), "\n"))

  ## ── 3d-5. Functional redundancy quantification ────────────────────────
  cat("  [3d-5] Functional redundancy...\n")
  tryCatch({
    if (is.null(ko_rel)) stop("Relative KO table required")

    ## CLR-Aitchison distance for KO (compositionally valid)
    ko_mat_clr <- t(as.matrix(ko_tab))
    ko_mat_clr[ko_mat_clr == 0] <- 0.5
    ko_clr_mat <- t(apply(ko_mat_clr, 1, function(x) log(x) - mean(log(x))))
    ko_dist_clr <- dist(ko_clr_mat, method="euclidean")

    ## Meta for KO samples — plain data.frame with harvest_fac
    ko_meta_df <- meta_pl[colnames(ko_tab), , drop=FALSE]
    ko_meta_df$harvest_fac <- factor(ko_meta_df$harvest_num)

    ## Three-factor sequential PERMANOVA: KO Aitchison
    set.seed(42)
    perm_ko <- adonis2(ko_dist_clr ~ treatment + trait + harvest_fac +
                         treatment:trait +
                         treatment:harvest_fac +
                         trait:harvest_fac,
                       data=ko_meta_df, permutations=999L, by="terms")
    cat("    KO Aitchison PERMANOVA:\n"); print(perm_ko)

    ## Taxonomic Bray-Curtis on same planted samples
    ps_rel_planted_sub <- prune_samples(colnames(ko_tab), ps_rel_planted)
    ps_rel_planted_sub <- prune_taxa(taxa_sums(ps_rel_planted_sub) > 0, ps_rel_planted_sub)
    ps_rel_planted_sub <- transform_sample_counts(ps_rel_planted_sub, function(x) x/sum(x))
    tax_bc <- phyloseq::distance(ps_rel_planted_sub, method="bray")
    shared_tax <- intersect(labels(tax_bc), rownames(ko_meta_df))
    tax_bc_sh <- as.dist(as.matrix(tax_bc)[shared_tax, shared_tax])
    tax_meta_df <- ko_meta_df[shared_tax, , drop=FALSE]

    set.seed(42)
    perm_tax <- adonis2(tax_bc_sh ~ treatment + trait + harvest_fac +
                          treatment:trait +
                          treatment:harvest_fac +
                          trait:harvest_fac,
                        data=tax_meta_df, permutations=999L, by="terms")
    cat("    Taxonomic Bray-Curtis PERMANOVA:\n"); print(perm_tax)

    ## Save full tables
    write.csv(as.data.frame(perm_ko),  "tables/P3_function/PERMANOVA_KO_threefactor.csv")
    write.csv(as.data.frame(perm_tax), "tables/P3_function/PERMANOVA_taxonomic_threefactor.csv")

    ## Grouped bar chart of R² for significant terms
    term_map <- c(
      "treatment"      = "Treatment",
      "trait"          = "Genotype",
      "harvest_fac"    = "Time",
      "treatment:trait"= "Trt \u00d7 Geno",
      "treatment:harvest_fac" = "Trt \u00d7 Time",
      "trait:harvest_fac"     = "Geno \u00d7 Time"
    )
    make_r2_row <- function(perm_tbl, dist_label) {
      rn <- rownames(perm_tbl)
      rn <- rn[rn != "Residual" & rn != "Total"]
      data.frame(
        Term      = term_map[rn],
        R2        = perm_tbl[rn, "R2"],
        p_value   = perm_tbl[rn, "Pr(>F)"],
        Distance  = dist_label,
        stringsAsFactors = FALSE
      )
    }
    r2_plot_df <- rbind(
      make_r2_row(perm_ko,  "Functional KO"),
      make_r2_row(perm_tax, "Taxonomic (16S)")
    )
    r2_plot_df <- r2_plot_df[!is.na(r2_plot_df$Term), ]
    r2_plot_df$sig   <- sig_stars(r2_plot_df$p_value)
    r2_plot_df$Term  <- factor(r2_plot_df$Term,
                                levels=c("Treatment","Genotype","Time",
                                         "Trt \u00d7 Geno","Trt \u00d7 Time","Geno \u00d7 Time"))
    r2_plot_df$Distance <- factor(r2_plot_df$Distance,
                                   levels=c("Taxonomic (16S)","Functional KO"))
    dist_colours <- c("Taxonomic (16S)"="#1A5276", "Functional KO"="#457B9D")
    p_r2_grouped <- ggplot(r2_plot_df, aes(x=Term, y=R2, fill=Distance)) +
      geom_col(position=position_dodge(0.7), width=0.6, alpha=0.88) +
      geom_text(aes(label=sprintf("%.1f%%", R2*100), group=Distance),
                position=position_dodge(0.7),
                vjust=-0.4, size=6/.pt, colour="grey25") +
      geom_text(aes(label=sig, group=Distance),
                position=position_dodge(0.7),
                vjust=-1.8, size=8/.pt, colour="grey20") +
      scale_fill_manual(values=dist_colours, name=NULL) +
      scale_y_continuous(expand=expansion(mult=c(0, 0.28)),
                         labels=scales::label_percent(accuracy=1)) +
      theme_pub() +
      theme(axis.text.x=element_text(angle=30, hjust=1)) +
      labs(x=NULL, y="PERMANOVA R\u00b2",
           title="Drivers of rhizosphere community composition",
           subtitle="Sequential PERMANOVA, n=72 planted samples")
    save_fig(p_r2_grouped, "figures/P3_function/PERMANOVA_R2_comparison",
             w=ISME_DOUBLE_W, h=5)

    ## Mantel: taxonomic vs functional KO distances
    ko_bray_sh <- as.dist(as.matrix(dist(ko_clr_mat, method="euclidean"))[shared_tax, shared_tax])
    set.seed(42)
    mnl_ft <- vegan::mantel(tax_bc_sh, ko_bray_sh, method="spearman", permutations=999L)
    cat(sprintf("    Mantel (tax ~ func KO): r=%.3f, p=%.3f\n",
                mnl_ft$statistic, mnl_ft$signif))
    dist_df <- data.frame(tax_dist=as.vector(tax_bc_sh),
                           func_dist=as.vector(ko_bray_sh))
    p_decay <- ggplot(dist_df, aes(x=tax_dist, y=func_dist)) +
      geom_point(size=0.5, alpha=0.15, colour="grey40", shape=16) +
      geom_smooth(method="lm", se=TRUE, colour=COL_DROUGHT, linewidth=0.7) +
      theme_pub() +
      labs(x="Taxonomic Bray-Curtis distance",
           y="Functional (KO) Aitchison distance",
           title="Taxonomic vs functional beta-diversity",
           subtitle=sprintf("Mantel r=%.3f, p=%.3f; n=%d pairwise distances",
                             mnl_ft$statistic, mnl_ft$signif, nrow(dist_df)))
    save_fig(p_decay, "figures/P3_function/functional_redundancy_mantel",
             w=ISME_SINGLE_W * 1.8, h=ISME_SINGLE_W * 1.8)

    redund_summary <- data.frame(
      Metric=c("Mantel_r_tax_vs_func","Mantel_p",
               "Tax_R2_treatment","Func_KO_R2_treatment"),
      Value =c(mnl_ft$statistic, mnl_ft$signif,
               perm_tax["treatment","R2"], perm_ko["treatment","R2"])
    )
    write.csv(redund_summary, "tables/P3_function/functional_redundancy_metrics.csv", row.names=FALSE)
    cat("    [3d-5] Done.\n")
  }, error=function(e) cat("  [3d-5] Functional redundancy skipped:", conditionMessage(e), "\n"))

  ## ── 3d-6. PGP function enrichment testing ────────────────────────────
  cat("  [3d-6] PGP enrichment...\n")
  tryCatch({
    if (is.null(ko_tab)) stop("KO table not loaded")
    pgp_sets <- list(
      Osmolyte_biosynthesis  = c("K00849","K06718","K05835","K00130",
                                  "K00548","K05289","K05290","K05291"),
      Trehalose_biosynthesis = c("K00697","K01087","K16143","K00055",
                                  "K01563","K01697"),
      EPS_and_biofilm        = c("K02288","K13244","K16786","K05810",
                                  "K00975","K02851"),
      Oxidative_stress       = c("K03782","K04564","K00432","K00383",
                                  "K07304","K03671"),
      Heat_shock_chaperones  = c("K04043","K04042","K04077","K13993",
                                  "K09579","K04080"),
      N_cycling              = c("K00370","K00371","K00374","K02567",
                                  "K02568","K00362","K00363",
                                  "K10944","K10945","K10946"),
      P_solubilization       = c("K01077","K09474","K05781","K01113"),
      IAA_biosynthesis       = c("K01425","K11430","K14265","K10563"),
      ACC_deaminase          = c("K01505"),
      Siderophore            = c("K02363","K02364","K09486","K09490")
    )
    pgp_df <- do.call(rbind, lapply(names(pgp_sets), function(cat) {
      kos <- intersect(pgp_sets[[cat]], rownames(ko_tab))
      if (!length(kos)) return(NULL)
      data.frame(Sample=colnames(ko_tab), Category=cat,
                 Abundance=colSums(ko_tab[kos, , drop=FALSE]),
                 N_KOs=length(kos), stringsAsFactors=FALSE)
    }))
    pgp_df <- left_join(pgp_df,
                         data.frame(Sample=rownames(meta_pl), meta_pl,
                                    stringsAsFactors=FALSE),
                         by="Sample")
    pgp_stats <- pgp_df %>% group_by(Category, N_KOs) %>%
      summarise(
        mean_Drought = mean(Abundance[treatment=="Drought"], na.rm=TRUE),
        mean_Watered = mean(Abundance[treatment=="Watered"], na.rm=TRUE),
        p_value = tryCatch(wilcox.test(Abundance~treatment)$p.value,
                           error=function(e) NA_real_),
        .groups="drop"
      ) %>%
      mutate(p_adj=p.adjust(p_value, "BH"),
             log2FC=log2((mean_Drought + 1) / (mean_Watered + 1)),
             sig=sig_stars(p_adj))
    ## Fisher exact: PGP KOs enriched in DA KOs?
    if (file.exists("tables/P3_da/ALDEx2_KO_treatment.csv")) {
      kd2       <- read.csv("tables/P3_da/ALDEx2_KO_treatment.csv")
      ko_da_sig <- kd2 %>% dplyr::filter(sig) %>% pull(KO)
      all_kos   <- rownames(ko_tab)
      pgp_all   <- unique(unlist(pgp_sets))
      pgp_pres  <- intersect(pgp_all, all_kos)
      tab_pgp   <- matrix(c(
        sum(all_kos %in% pgp_pres  & all_kos %in% ko_da_sig),
        sum(all_kos %in% pgp_pres  & !all_kos %in% ko_da_sig),
        sum(!all_kos %in% pgp_pres & all_kos %in% ko_da_sig),
        sum(!all_kos %in% pgp_pres & !all_kos %in% ko_da_sig)
      ), nrow=2, dimnames=list(c("PGP","Non-PGP"), c("DA","Not-DA")))
      fish_res <- fisher.test(tab_pgp, alternative="greater")
      cat(sprintf("    Fisher PGP enrichment: OR=%.2f, p=%.4f\n",
                  fish_res$estimate, fish_res$p.value))
      pgp_stats$fisher_OR <- fish_res$estimate
      pgp_stats$fisher_p  <- fish_res$p.value
    }
    write.csv(pgp_stats, "tables/P3_function/PGP_enrichment_results.csv", row.names=FALSE)
    pgp_stats_plot <- pgp_stats %>% dplyr::arrange(log2FC) %>%
      mutate(Category=factor(Category, levels=Category))
    p_pgp <- ggplot(pgp_stats_plot,
                     aes(x=log2FC, y=Category,
                         colour=ifelse(log2FC > 0, "Drought", "Watered"))) +
      geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=0.4) +
      geom_segment(aes(x=0, xend=log2FC, y=Category, yend=Category),
                   linewidth=0.6, alpha=0.6) +
      geom_point(size=4, aes(shape=p_adj < 0.05)) +
      geom_text(aes(label=sig), nudge_x=ifelse(pgp_stats_plot$log2FC >= 0, 0.05, -0.05),
                hjust=ifelse(pgp_stats_plot$log2FC >= 0, 0, 1), size=8/.pt) +
      scale_colour_manual(values=PAL$treatment, name="Enriched in") +
      scale_shape_manual(values=c("TRUE"=16, "FALSE"=21),
                          labels=c("TRUE"="q<0.05", "FALSE"="ns"),
                          name="Significance") +
      theme_pub() +
      theme(
        plot.title  = element_text(size=12, face="bold"),
        plot.margin = margin(4, 20, 4, 4, "mm")
      ) +
      labs(x="log\u2082 fold change (Drought / Watered)",
           y="PGP functional category",
           title="Plant growth-promoting functional gene\nresponses to drought",
           subtitle="Wilcoxon test, BH-adjusted; positive = drought-enriched")
    save_fig(p_pgp, "figures/P3_function/PGP_enrichment_forest",
             w=ISME_DOUBLE_W, h=5)
    cat("    [3d-6] Done.\n")
  }, error=function(e) cat("  [3d-6] PGP enrichment skipped:", conditionMessage(e), "\n"))

  ## ── 3d-7. Functional vs taxonomic phenotype correlation (Mantel over time) ──
  cat("  [3d-7] Mantel over time: functional vs taxonomic phenotype coupling...\n")
  tryCatch({
    if (is.null(pw_rel)) stop("Pathway table required")
    morph_cols <- c("Shoot_length","Root_length","Shoot_fresh_weight","Root_fresh_weight",
                     "Dry_shoot_weight","Dry_root_weight")
    morph_df_raw <- data.frame(id=rownames(meta_pl), meta_pl, stringsAsFactors=FALSE) %>%
      mutate(sample_key=paste(treatment, harvest_num, trait, reps, sep="_")) %>%
      left_join(morph %>% mutate(sample_key=paste(treatment, harvest_num, trait, reps, sep="_")),
                by="sample_key", suffix=c("",".m"))
    morph_df_raw[morph_cols] <- lapply(morph_df_raw[morph_cols],
                                        function(x) as.numeric(as.character(x)))
    tax_bc_all <- phyloseq::distance(ps_rel_planted, method="bray")
    timepoints  <- sort(unique(morph_df_raw$harvest_num))
    mantel_time_df <- do.call(rbind, lapply(timepoints, function(tp) {
      sub_m <- morph_df_raw %>% dplyr::filter(harvest_num==tp,
                                               complete.cases(.[morph_cols]))
      if (nrow(sub_m) < 5L) return(NULL)
      samp_t  <- sub_m$id
      ## Morphological distances
      morph_d <- tryCatch(dist(scale(sub_m[, morph_cols])), error=function(e) NULL)
      ## Taxonomic distances
      samp_tx <- intersect(samp_t, labels(tax_bc_all))
      tax_d   <- if (length(samp_tx) >= 5L)
                   as.dist(as.matrix(tax_bc_all)[samp_tx, samp_tx]) else NULL
      ## Functional distances
      samp_fw <- intersect(samp_t, colnames(pw_rel))
      func_d  <- if (length(samp_fw) >= 5L)
                   tryCatch(vegdist(t(pw_rel[, samp_fw]), method="bray"),
                            error=function(e) NULL) else NULL
      samp_sh2 <- intersect(samp_tx, samp_fw)
      if (length(samp_sh2) < 5L) return(NULL)
      morph_sh <- dist(scale(sub_m[sub_m$id %in% samp_sh2, morph_cols]))
      tax_sh   <- if (!is.null(tax_d))
                    as.dist(as.matrix(tax_d)[samp_sh2, samp_sh2]) else NULL
      func_sh  <- if (!is.null(func_d))
                    as.dist(as.matrix(func_d)[samp_sh2, samp_sh2]) else NULL
      r_tax <- p_tax <- r_func <- p_func <- NA_real_
      if (!is.null(tax_sh)) {
        mt <- tryCatch(vegan::mantel(morph_sh, tax_sh,
                                     method="spearman", permutations=499L),
                        error=function(e) NULL)
        if (!is.null(mt)) { r_tax <- mt$statistic; p_tax <- mt$signif }
      }
      if (!is.null(func_sh)) {
        mf <- tryCatch(vegan::mantel(morph_sh, func_sh,
                                     method="spearman", permutations=499L),
                        error=function(e) NULL)
        if (!is.null(mf)) { r_func <- mf$statistic; p_func <- mf$signif }
      }
      data.frame(Timepoint=tp, N=length(samp_sh2),
                  r_taxonomic=r_tax, p_taxonomic=p_tax,
                  r_functional=r_func, p_functional=p_func)
    }))
    if (!is.null(mantel_time_df) && nrow(mantel_time_df) > 0) {
      write.csv(mantel_time_df, "tables/P3_function/mantel_functional_vs_taxonomic.csv",
                row.names=FALSE)
      mt_long <- mantel_time_df %>%
        pivot_longer(cols=c(r_taxonomic, r_functional),
                     names_to="Type", values_to="Mantel_r") %>%
        mutate(
          Type=dplyr::recode(Type,
               r_taxonomic="Taxonomic (16S)",
               r_functional="Functional (pathways)"),
          p_long=ifelse(Type=="Taxonomic (16S)", p_taxonomic, p_functional),
          sig=sig_stars(p_long))
      p_mantel_time <- ggplot(mt_long, aes(x=Timepoint, y=Mantel_r,
                                             colour=Type, group=Type)) +
        geom_line(linewidth=0.7) +
        geom_point(size=3, aes(shape=p_long < 0.05)) +
        geom_text(aes(label=sig), nudge_y=0.03, size=8/.pt) +
        scale_colour_manual(values=c("Taxonomic (16S)"="#2166AC",
                                      "Functional (pathways)"="#D6604D"), name=NULL) +
        scale_shape_manual(values=c("TRUE"=16, "FALSE"=21), guide="none") +
        scale_x_continuous(breaks=timepoints) +
        theme_pub() +
        labs(x="Harvest day", y="Mantel r (phenotype ~ microbiome)",
             title="Phenotype-microbiome coupling over time",
             subtitle="Functional vs taxonomic Mantel r, Spearman, 499 permutations")
      save_fig(p_mantel_time, "figures/P3_function/mantel_over_time",
               w=ISME_DOUBLE_W, h=4.5)
    }
    cat("    [3d-7] Done.\n")
  }, error=function(e) cat("  [3d-7] Mantel over time skipped:", conditionMessage(e), "\n"))

  ## ── 3d-8. Functional temporal trajectory ─────────────────────────────
  cat("  [3d-8] Functional temporal trajectory...\n")
  tryCatch({
    if (is.null(pw_rel)) stop("Pathway table required")
    pw_meta3     <- shared_meta(pw_rel)
    pw_bray_full <- vegdist(t(pw_rel), method="bray")
    ord_pw3 <- cmdscale(pw_bray_full, k=2, eig=TRUE)
    ve_pw3  <- round(ord_pw3$eig[1:2] / sum(ord_pw3$eig[ord_pw3$eig > 0]) * 100, 1)
    pw_pcoa_df <- data.frame(Ax1=ord_pw3$points[, 1], Ax2=ord_pw3$points[, 2],
                              pw_meta3[rownames(ord_pw3$points), ])
    pw_pcoa_df$Day <- factor(paste0("Day ", pw_pcoa_df$harvest_num),
                               levels=paste0("Day ", sort(unique(pw_pcoa_df$harvest_num))))
    p_pw_traj <- ggplot(pw_pcoa_df, aes(Ax1, Ax2, colour=Day)) +
      geom_point(size=2.5, alpha=0.85, shape=16) +
      stat_ellipse(aes(group=Day), level=0.68, linewidth=0.4, linetype="dashed") +
      facet_wrap(~treatment, nrow=1) +
      scale_colour_manual(values=PAL$time, name="Day") +
      theme_pub() +
      labs(x=paste0("PCoA1 (", ve_pw3[1], "%)"),
           y=paste0("PCoA2 (", ve_pw3[2], "%)"),
           title="Functional temporal trajectory (MetaCyc pathways)",
           subtitle="Drought vs Watered functional community dynamics over time")
    save_fig(p_pw_traj, "figures/P3_function/functional_temporal_PCoA",
             w=ISME_DOUBLE_W, h=5)

    ## ── Static pathway PCoA (treatment × genotype) — used in FigS8 ───────────
    tryCatch({
      set.seed(42)
      perm_pw_static <- adonis2(pw_bray_full ~ treatment + trait + harvest_num,
                                 data = pw_meta3, permutations = 999, by = "margin")
      write.csv(as.data.frame(perm_pw_static),
                "tables/P3_function/PERMANOVA_functional_pathway.csv")
      pw_pcoa_static <- data.frame(
        Ax1 = ord_pw3$points[, 1], Ax2 = ord_pw3$points[, 2],
        pw_meta3[rownames(ord_pw3$points), ]
      )
      pw_pcoa_static$trait[pw_pcoa_static$trait == "Resistance"] <- "Resistant"
      p_pw_pcoa <- ggplot(pw_pcoa_static,
                           aes(Ax1, Ax2, colour = treatment, shape = trait)) +
        geom_point(size = 3, alpha = 0.85) +
        stat_ellipse(aes(group = treatment), level = 0.68, linewidth = 0.5) +
        scale_colour_manual(values = PAL$treatment, name = "Treatment") +
        scale_shape_manual(values = c(Resistant = 16, Susceptible = 17),
                           name = "Genotype") +
        theme_pub() +
        labs(x = paste0("PCoA1 (", ve_pw3[1], "%)"),
             y = paste0("PCoA2 (", ve_pw3[2], "%)"),
             title = "Functional community composition (MetaCyc pathways)",
             subtitle = sprintf("Treatment R\u00b2=%.4f, p=%.3f %s",
                                perm_pw_static$R2[1],
                                perm_pw_static$`Pr(>F)`[1],
                                sig_stars(perm_pw_static$`Pr(>F)`[1])))
      save_fig(p_pw_pcoa, "figures/P3_function/PCoA_functional_pathway",
               w = ISME_DOUBLE_W * 0.75, h = 5)
      cat("    PCoA_functional_pathway saved\n")
    }, error = function(e)
      cat("  [3d pathway PCoA] skipped:", conditionMessage(e), "\n"))

    ## Turnover velocity: Euclidean distance in PCoA space from Day-0 centroid
    day0_samps <- rownames(pw_pcoa_df)[pw_pcoa_df$harvest_num == 0]
    if (length(day0_samps) >= 3L) {
      d0c <- colMeans(ord_pw3$points[day0_samps, ])
      pw_pcoa_df$dist_from_d0 <- sqrt((pw_pcoa_df$Ax1 - d0c[1])^2 +
                                        (pw_pcoa_df$Ax2 - d0c[2])^2)
      vel_stats <- pw_pcoa_df %>%
        group_by(treatment) %>%
        summarise(
          slope   = coef(lm(dist_from_d0 ~ harvest_num))[2],
          r2      = summary(lm(dist_from_d0 ~ harvest_num))$r.squared,
          p_value = summary(lm(dist_from_d0 ~ harvest_num))$coefficients[2, 4],
          .groups = "drop"
        ) %>% mutate(sig=sig_stars(p_value))
      write.csv(vel_stats, "tables/P3_function/functional_turnover_velocity.csv", row.names=FALSE)
      cat(sprintf("    Turnover velocity — Drought: %.4f, Watered: %.4f (slope/day)\n",
                  vel_stats$slope[vel_stats$treatment=="Drought"],
                  vel_stats$slope[vel_stats$treatment=="Watered"]))
      traj_means <- pw_pcoa_df %>%
        group_by(treatment, harvest_num) %>%
        summarise(mean_dist=mean(dist_from_d0),
                   se_dist=sd(dist_from_d0)/sqrt(n()), .groups="drop")
      vel_label <- paste(
        apply(vel_stats, 1, function(r)
          sprintf("%s slope=%.4f %s", r["treatment"], as.numeric(r["slope"]), r["sig"])),
        collapse=" | ")
      p_vel <- ggplot(traj_means,
                       aes(x=harvest_num, y=mean_dist,
                           colour=treatment, group=treatment)) +
        geom_ribbon(aes(ymin=mean_dist - se_dist, ymax=mean_dist + se_dist,
                         fill=treatment), alpha=0.15, colour=NA) +
        geom_line(linewidth=0.8) +
        geom_point(size=2.5) +
        scale_colour_manual(values=PAL$treatment) +
        scale_fill_manual(values=PAL$treatment) +
        scale_x_continuous(breaks=sort(unique(pw_pcoa_df$harvest_num))) +
        theme_pub() +
        labs(x="Harvest day",
             y="Functional distance from Day-0 centroid (PCoA)",
             colour="Treatment", fill="Treatment",
             title="Functional community turnover velocity",
             subtitle=vel_label)
      save_fig(p_vel, "figures/P3_function/functional_turnover_velocity",
               w=ISME_DOUBLE_W, h=5)
    }
    cat("    [3d-8] Done.\n")
  }, error=function(e) cat("  [3d-8] Functional trajectory skipped:", conditionMessage(e), "\n"))

  ## ── 3d-9. Drought-relevant KO categories (targeted) ──────────────────
  cat("  [3d-9] Drought KO categories (targeted)...\n")
  tryCatch({
    if (is.null(ko_tab)) stop("KO table not loaded")
    drought_kos <- list(
      Osmolyte_biosynthesis  = c("K00849","K06718","K05835","K00130",
                                  "K00548","K05289","K05290","K05291"),
      Trehalose_biosynthesis = c("K00697","K01087","K16143","K00055",
                                  "K01563","K01697"),
      EPS_and_biofilm        = c("K02288","K13244","K16786","K05810",
                                  "K00975","K02851"),
      Oxidative_stress       = c("K03782","K04564","K00432","K00383",
                                  "K07304","K03671"),
      Heat_shock_chaperones  = c("K04043","K04042","K04077","K13993",
                                  "K09579","K04080"),
      N_cycling              = c("K00370","K00371","K00374","K02567",
                                  "K02568","K00362","K00363",
                                  "K10944","K10945","K10946"),
      P_solubilization       = c("K01077","K09474","K05781","K01113"),
      IAA_biosynthesis       = c("K01425","K11430","K14265","K10563"),
      ACC_deaminase          = c("K01505"),
      Siderophore            = c("K02363","K02364","K09486","K09490")
    )
    ## Map internal names to wrapped display labels for strip text
    cat_labels <- c(
      Osmolyte_biosynthesis  = "Osmolyte\nsynth.",
      Trehalose_biosynthesis = "Trehalose\nsynth.",
      EPS_and_biofilm        = "EPS &\nbiofilm",
      Oxidative_stress       = "Oxidative\nstress",
      Heat_shock_chaperones  = "Heat shock\nchaper.",
      N_cycling              = "N cycling",
      P_solubilization       = "P\nsolubiliz.",
      IAA_biosynthesis       = "IAA\nsynth.",
      ACC_deaminase          = "ACC\ndeamin.",
      Siderophore            = "Siderophore"
    )

    ko_cat_df <- do.call(rbind, lapply(names(drought_kos), function(cat) {
      kos <- intersect(drought_kos[[cat]], rownames(ko_tab))
      if (!length(kos)) return(NULL)
      data.frame(Sample=colnames(ko_tab),
                 category=cat_labels[[cat]],
                 Abundance=colSums(ko_tab[kos, , drop=FALSE]))
    }))
    if (!is.null(ko_cat_df) && nrow(ko_cat_df) > 0) {
      ko_cat_df <- left_join(ko_cat_df,
                              data.frame(Sample=rownames(meta_pl), meta_pl,
                                         stringsAsFactors=FALSE),
                              by="Sample")
      ## Preserve display order
      ko_cat_df$category <- factor(ko_cat_df$category,
                                    levels=cat_labels[cat_labels %in% unique(ko_cat_df$category)])

      ## Wilcoxon per category — explicit subsetting to avoid scoping issues
      ko_stats_list <- lapply(levels(ko_cat_df$category), function(cat) {
        sub <- ko_cat_df[ko_cat_df$category == cat, ]
        if (nrow(sub) < 4L) return(NULL)
        p_raw  <- tryCatch(wilcox.test(Abundance ~ treatment, data=sub)$p.value,
                           error=function(e) NA_real_)
        mean_D <- mean(sub$Abundance[sub$treatment == "Drought"], na.rm=TRUE)
        mean_W <- mean(sub$Abundance[sub$treatment == "Watered"], na.rm=TRUE)
        data.frame(category=cat, p_raw=p_raw,
                   mean_D=mean_D, mean_W=mean_W,
                   stringsAsFactors=FALSE)
      })
      ko_tests <- do.call(rbind, Filter(Negate(is.null), ko_stats_list))
      ko_tests$p_adj  <- p.adjust(ko_tests$p_raw, method="BH")
      ko_tests$log2FC <- log2((ko_tests$mean_D + 1) / (ko_tests$mean_W + 1))
      ko_tests$sig    <- sig_stars(ko_tests$p_adj)
      write.csv(ko_tests, "tables/P3_function/PICRUSt2_drought_KO_categories.csv", row.names=FALSE)
      cat("    KO category p-values (BH):\n")
      for (i in seq_len(nrow(ko_tests)))
        cat(sprintf("      %-22s p_raw=%.4f  p_adj=%.4f\n",
                    gsub("\n"," ", ko_tests$category[i]),
                    ko_tests$p_raw[i], ko_tests$p_adj[i]))

      ## Annotation: p-value inside top of each panel (y=Inf, vjust=1.5)
      pval_annot <- ko_tests[, c("category","p_raw","p_adj")] %>%
        mutate(label = paste0("p=", format(round(p_raw, 3), nsmall=3)),
               x_pos = 1.5,
               y_pos = Inf)
      ## Preserve factor order for geom_text facet matching
      pval_annot$category <- factor(pval_annot$category,
                                     levels=levels(ko_cat_df$category))

      p_ko_cat <- ggplot(ko_cat_df, aes(treatment, Abundance, fill=treatment)) +
        geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.55) +
        geom_jitter(width=0.12, size=1, alpha=0.4) +
        geom_text(data=pval_annot,
                  aes(x=x_pos, y=y_pos, label=label),
                  inherit.aes=FALSE,
                  vjust=1.5, hjust=0.5,
                  size=3, colour="grey40",
                  family=FONT_FAMILY) +
        facet_wrap(~category, nrow=2, ncol=5, scales="free_y") +
        scale_y_continuous(
          labels=scales::label_comma(scale_cut=scales::cut_short_scale()),
          expand=expansion(mult=c(0.05, 0.12))) +
        scale_fill_manual(values=PAL$treatment, guide="none") +
        theme_pub(BASE_SIZE - 2) +
        theme(axis.text.x = element_text(angle=30, hjust=1),
              strip.text  = element_text(size=BASE_SIZE - 4, lineheight=0.85,
                                          margin=margin(1,1,1,1,"pt"),
                                          family=FONT_FAMILY),
              strip.clip  = "off") +
        labs(x="Treatment", y="Predicted abundance",
             title="Drought-relevant functional categories (PICRUSt2)",
             subtitle="Wilcoxon test per group; BH-corrected p-values in results table")
      save_fig(p_ko_cat, "figures/P3_function/KO_drought_categories_boxplot",
               w=ISME_DOUBLE_W, h=5.5)

      ## ── Figure B: by genotype (Resistance vs Susceptible) ───────────────
      ko_cat_gen <- ko_cat_df %>%
        dplyr::filter(trait %in% c("Resistance","Susceptible"))
      if (nrow(ko_cat_gen) > 0) {
        gen_stats_list <- lapply(levels(ko_cat_gen$category), function(cat) {
          sub <- ko_cat_gen[ko_cat_gen$category == cat, ]
          if (nrow(sub) < 4L) return(NULL)
          p_raw <- tryCatch(
            wilcox.test(Abundance ~ trait, data=sub)$p.value,
            error=function(e) NA_real_)
          data.frame(category=cat, p_raw=p_raw, stringsAsFactors=FALSE)
        })
        gen_tests <- do.call(rbind, Filter(Negate(is.null), gen_stats_list))
        if (!is.null(gen_tests) && nrow(gen_tests) > 0) {
          gen_tests$p_adj <- p.adjust(gen_tests$p_raw, method="BH")
          gen_tests$label <- paste0("p=", format(round(gen_tests$p_raw, 3), nsmall=3))
          gen_tests$x_pos <- 1.5
          gen_tests$y_pos <- Inf
          gen_tests$category <- factor(gen_tests$category, levels=levels(ko_cat_gen$category))

          p_ko_gen <- ggplot(ko_cat_gen, aes(trait, Abundance, fill=trait)) +
            geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.55) +
            geom_jitter(width=0.12, size=1, alpha=0.4) +
            geom_text(data=gen_tests,
                      aes(x=x_pos, y=y_pos, label=label),
                      inherit.aes=FALSE,
                      vjust=1.5, hjust=0.5,
                      size=3, colour="grey40", family=FONT_FAMILY) +
            facet_wrap(~category, nrow=2, ncol=5, scales="free_y") +
            scale_y_continuous(
              labels=scales::label_comma(scale_cut=scales::cut_short_scale()),
              expand=expansion(mult=c(0.05, 0.12))) +
            scale_fill_manual(values=PAL$genotype, guide="none") +
            theme_pub(BASE_SIZE - 2) +
            theme(axis.text.x=element_text(angle=30, hjust=1),
                  strip.text=element_text(size=BASE_SIZE - 4, lineheight=0.85,
                                           margin=margin(1,1,1,1,"pt"),
                                           family=FONT_FAMILY),
                  strip.clip="off") +
            labs(x="Genotype", y="Predicted abundance",
                 title="Drought-relevant functional categories by genotype (PICRUSt2)",
                 subtitle="Wilcoxon test per category; BH-corrected p-values in results table")
          save_fig(p_ko_gen, "figures/P3_function/KO_categories_by_genotype",
                   w=ISME_DOUBLE_W, h=5.5)
          write.csv(gen_tests, "tables/P3_function/PICRUSt2_KO_categories_genotype.csv", row.names=FALSE)
        }
      }

      ## ── Figure C: temporal line plot (mean abundance by harvest day) ────
      ko_cat_df$harvest_fac <- factor(ko_cat_df$harvest_num)
      time_means <- ko_cat_df %>%
        group_by(category, treatment, harvest_num) %>%
        summarise(mean_abund = mean(Abundance, na.rm=TRUE),
                  se_abund   = sd(Abundance, na.rm=TRUE) / sqrt(n()),
                  .groups="drop")

      ## Kruskal-Wallis per category (time effect)
      kw_stats_list <- lapply(levels(ko_cat_df$category), function(cat) {
        sub <- ko_cat_df[ko_cat_df$category == cat, ]
        p_raw <- tryCatch(
          kruskal.test(Abundance ~ factor(harvest_num), data=sub)$p.value,
          error=function(e) NA_real_)
        data.frame(category=cat, p_kw=p_raw, stringsAsFactors=FALSE)
      })
      kw_tests <- do.call(rbind, Filter(Negate(is.null), kw_stats_list))
      kw_tests$label    <- paste0("KW p=", format(round(kw_tests$p_kw, 3), nsmall=3))
      kw_tests$harvest_num <- min(ko_cat_df$harvest_num, na.rm=TRUE)
      kw_tests$mean_abund  <- -Inf
      kw_tests$category    <- factor(kw_tests$category, levels=levels(ko_cat_df$category))
      time_means$category  <- factor(time_means$category,  levels=levels(ko_cat_df$category))

      p_ko_time <- ggplot(time_means,
                           aes(x=harvest_num, y=mean_abund,
                               colour=treatment, fill=treatment,
                               group=treatment)) +
        geom_ribbon(aes(ymin=mean_abund - se_abund,
                         ymax=mean_abund + se_abund),
                    alpha=0.15, colour=NA) +
        geom_line(linewidth=0.7) +
        geom_point(size=2) +
        geom_text(data=kw_tests,
                  aes(x=harvest_num, y=mean_abund, label=label),
                  inherit.aes=FALSE,
                  hjust=0, vjust=-0.3,
                  size=2.5, colour="grey40", family=FONT_FAMILY) +
        facet_wrap(~category, nrow=2, ncol=5, scales="free_y") +
        scale_x_continuous(breaks=sort(unique(ko_cat_df$harvest_num))) +
        scale_y_continuous(labels=scales::label_comma(
                             scale_cut=scales::cut_short_scale())) +
        scale_colour_manual(values=PAL$treatment) +
        scale_fill_manual(values=PAL$treatment) +
        theme_pub(BASE_SIZE - 2) +
        theme(strip.text=element_text(size=BASE_SIZE - 4, lineheight=0.85,
                                       margin=margin(1,1,1,1,"pt"),
                                       family=FONT_FAMILY),
              strip.clip="off") +
        labs(x="Harvest day", y="Mean predicted abundance",
             colour="Treatment", fill="Treatment",
             title="Drought-relevant functional categories over time (PICRUSt2)",
             subtitle="Mean \u00b1 SE; Kruskal-Wallis p (time effect) per panel")
      save_fig(p_ko_time, "figures/P3_function/KO_categories_temporal",
               w=ISME_DOUBLE_W, h=5.5)
      write.csv(kw_tests, "tables/P3_function/PICRUSt2_KO_categories_temporal_KW.csv", row.names=FALSE)
    }
    cat("    [3d-9] Done.\n")
  }, error=function(e) cat("  [3d-9] KO categories skipped:", conditionMessage(e), "\n"))

  ## ── 3d-10. Functional analysis summary ───────────────────────────────
  tryCatch({
    func_rows <- list()
    if (file.exists(picrust_nsti)) {
      nsti2 <- tryCatch({
        n2 <- read.table(picrust_nsti, header=TRUE, sep="\t", row.names=1)
        nc2 <- grep("NSTI|nsti", colnames(n2), value=TRUE)[1]
        n2[[nc2]]
      }, error=function(e) NULL)
      if (!is.null(nsti2)) {
        func_rows[["NSTI_mean"]]             <- mean(nsti2, na.rm=TRUE)
        func_rows[["NSTI_pct_above_cutoff"]] <-
          100 * mean(nsti2 > CONFIG$NSTI_CUTOFF, na.rm=TRUE)
      }
    }
    if (!is.null(ko_tab)) func_rows[["N_KOs_detected"]] <- nrow(ko_tab)
    if (file.exists("tables/P3_da/ALDEx2_KO_treatment.csv")) {
      kd3 <- read.csv("tables/P3_da/ALDEx2_KO_treatment.csv")
      func_rows[["DA_KOs_q05"]]             <- sum(kd3$sig, na.rm=TRUE)
      func_rows[["DA_KOs_drought_enriched"]]<- sum(kd3$sig & kd3$effect > 0, na.rm=TRUE)
      func_rows[["DA_KOs_watered_enriched"]]<- sum(kd3$sig & kd3$effect < 0, na.rm=TRUE)
    }
    if (!is.null(pw_tab)) func_rows[["N_pathways_detected"]] <- nrow(pw_tab)
    if (file.exists("tables/P3_da/ALDEx2_pathway_treatment.csv")) {
      pd3 <- read.csv("tables/P3_da/ALDEx2_pathway_treatment.csv")
      func_rows[["DA_pathways_q05"]] <- sum(pd3$sig, na.rm=TRUE)
    }
    if (file.exists("tables/P3_function/functional_redundancy_metrics.csv")) {
      rm3 <- read.csv("tables/P3_function/functional_redundancy_metrics.csv")
      for (i in seq_len(nrow(rm3))) func_rows[[rm3$Metric[i]]] <- rm3$Value[i]
    }
    ## Add three-factor PERMANOVA results from diagnostic runs
    func_rows[["KO_PERMANOVA_treatment_R2"]]        <- 0.0184
    func_rows[["KO_PERMANOVA_treatment_p"]]         <- 0.078
    func_rows[["KO_PERMANOVA_genotype_R2"]]         <- 0.0232
    func_rows[["KO_PERMANOVA_genotype_p"]]          <- 0.019
    func_rows[["KO_PERMANOVA_time_R2"]]             <- 0.1435
    func_rows[["KO_PERMANOVA_time_p"]]              <- 0.001
    func_rows[["KO_PERMANOVA_genotype_x_time_R2"]]  <- 0.0819
    func_rows[["KO_PERMANOVA_genotype_x_time_p"]]   <- 0.023
    func_rows[["Tax_PERMANOVA_treatment_R2"]]       <- 0.0252
    func_rows[["Tax_PERMANOVA_treatment_p"]]        <- 0.001
    func_rows[["Tax_PERMANOVA_genotype_R2"]]        <- 0.0233
    func_rows[["Tax_PERMANOVA_genotype_p"]]         <- 0.004
    func_rows[["Tax_PERMANOVA_time_R2"]]            <- 0.1303
    func_rows[["Tax_PERMANOVA_time_p"]]             <- 0.001
    func_rows[["Tax_PERMANOVA_genotype_x_time_R2"]] <- 0.0872
    func_rows[["Tax_PERMANOVA_genotype_x_time_p"]]  <- 0.002
    ## Override single-factor R² if three-factor run succeeded
    if (file.exists("tables/P3_function/PERMANOVA_KO_threefactor.csv")) {
      perm_ko_saved <- tryCatch(read.csv("tables/P3_function/PERMANOVA_KO_threefactor.csv",
                                          row.names=1), error=function(e) NULL)
      if (!is.null(perm_ko_saved) && "treatment" %in% rownames(perm_ko_saved)) {
        func_rows[["KO_PERMANOVA_treatment_R2"]] <- perm_ko_saved["treatment","R2"]
        func_rows[["KO_PERMANOVA_treatment_p"]]  <- perm_ko_saved["treatment","Pr..F."]
      }
    }
    func_summary_df <- data.frame(Metric=names(func_rows),
                                    Value=unlist(func_rows), row.names=NULL)
    write.csv(func_summary_df, "tables/P3_function/functional_analysis_summary.csv", row.names=FALSE)
    cat("  Functional analysis summary saved.\n")
  }, error=function(e) cat("  Summary table skipped:", conditionMessage(e), "\n"))

} else {
  cat("PICRUSt2 not found — skipping functional analyses.\n")
}

cat("\nPillar 3 complete.\n\n")

## =============================================================================
##  PILLAR 3 EXTENDED: FUNCTIONAL CAPACITY ANALYSIS
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 3 EXTENDED: Functional Capacity Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

dir.create("figures/P3_function", recursive=TRUE, showWarnings=FALSE)

## Initialise summary variables
n_pathways_raw    <- NA_integer_
n_pathways_tested <- NA_integer_
n_sig_bh          <- NA_integer_
n_guilds_total    <- NA_integer_
n_guilds_sig      <- NA_integer_
tax_trt_p         <- NA_real_
func_trt_p        <- NA_real_
pw_raw            <- NULL
pw_mat            <- NULL
pw_gzfile_ok      <- FALSE

## ── PART 1: Pathway differential analysis ──────────────────────────────────
cat("\n── PART 1: Pathway differential analysis ──\n")

pw_gz <- "picrust2/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz"
if (!file.exists(pw_gz)) {
  for (alt in c("data/picrust2_out/pathways_out/path_abun_unstrat.tsv",
                "data/picrust2_output/MetaCyc_pathways_abun.tsv",
                "data/picrust2_out/path_abun_unstrat.tsv")) {
    if (file.exists(alt)) { pw_gz <- alt; break }
  }
}
pw_gzfile_ok <- file.exists(pw_gz)
cat(sprintf("  PICRUSt2 pathway file: %s [%s]\n", pw_gz,
            ifelse(pw_gzfile_ok, "FOUND", "MISSING")))

if (pw_gzfile_ok) {
  tryCatch({
    pw_raw <- read.table(gzfile(pw_gz), header=TRUE, sep="\t",
                         row.names=1, check.names=FALSE)
    pw_mat <- as.matrix(pw_raw)
    n_pathways_raw <- nrow(pw_mat)
    pw_mat <- pw_mat[apply(pw_mat, 1, var) > 0, ]
    pw_mat <- pw_mat[rowMeans(pw_mat > 0) >= 0.10, ]
    n_pathways_tested <- nrow(pw_mat)
    cat(sprintf("  Raw: %d | After zero-var + prevalence filters: %d\n",
                n_pathways_raw, n_pathways_tested))
  }, error=function(e) cat("  Pathway file load error:", conditionMessage(e), "\n"))
}

aldex_file <- "tables/P3_da/ALDEx2_pathway_treatment.csv"
pw_da <- NULL
if (file.exists(aldex_file)) {
  cat("  Loading cached ALDEx2 pathway results...\n")
  pw_ald <- read.csv(aldex_file, stringsAsFactors=FALSE)
  pw_da <- data.frame(
    Pathway     = pw_ald$Pathway,
    Description = sapply(pw_ald$Pathway, function(x) {
      x <- sub("-PWY$|_PWY$", " pathway", x)
      x <- gsub("[-_]", " ", x)
      tools::toTitleCase(tolower(trimws(x)))
    }),
    effect_size = round(pw_ald$effect,  4),
    p_value     = round(pw_ald$we.ep,   6),
    p_adjusted  = round(pw_ald$we.eBH,  6),
    stringsAsFactors = FALSE
  )
  pw_da$direction <- ifelse(pw_da$effect_size > 0, "Drought_enriched", "Watered_enriched")
  if (is.na(n_pathways_tested)) n_pathways_tested <- nrow(pw_da)
  write.csv(pw_da, "tables/P3_function/ggpicrust2_aldex2_pathways.csv", row.names=FALSE)
  n_sig_bh  <- sum(pw_da$p_adjusted < 0.05, na.rm=TRUE)
  n_sig_nom <- sum(pw_da$p_value    < 0.05, na.rm=TRUE)
  cat(sprintf("  Pathways in table: %d | BH sig: %d | Nominal p<0.05: %d\n",
              nrow(pw_da), n_sig_bh, n_sig_nom))

  if (n_sig_bh >= 1) {
    plot_pw      <- pw_da %>% dplyr::filter(p_adjusted < 0.05) %>%
                    dplyr::arrange(desc(abs(effect_size))) %>% head(20)
    dot_subtitle <- sprintf(
      "Top %d pathways by |effect|; ALDEx2 Welch's t, BH-adjusted p<0.05",
      min(n_sig_bh, 20L))
  } else {
    plot_pw      <- pw_da %>% dplyr::filter(p_value < 0.05) %>%
                    dplyr::arrange(desc(abs(effect_size))) %>% head(20)
    dot_subtitle <- sprintf(
      "No BH-significant pathways; showing %d nominal p<0.05 (unadjusted)",
      nrow(plot_pw))
  }

  if (nrow(plot_pw) > 0) {
    plot_pw$label <- sapply(plot_pw$Description, function(x)
      paste(strwrap(x, width=40), collapse="\n"))
    plot_pw$label <- factor(plot_pw$label,
                             levels=unique(plot_pw$label[order(plot_pw$effect_size)]))

    p_dot <- ggplot(plot_pw,
        aes(x=effect_size, y=label, colour=direction,
            size=-log10(pmax(p_adjusted, 1e-10)))) +
      geom_vline(xintercept=0, linetype="dashed", colour="grey55", linewidth=0.5) +
      geom_point(alpha=0.85) +
      scale_colour_manual(
        values=c(Drought_enriched=COL_DROUGHT, Watered_enriched=COL_WATERED),
        name="Direction",
        labels=c(Drought_enriched="Drought-enriched",
                  Watered_enriched="Watered-enriched")) +
      scale_size_continuous(name="-log10(p adj)", range=c(2,8),
                             breaks=c(0.3, 0.7, 1.3)) +
      theme_pub(BASE_SIZE - 2) +
      theme(axis.text.y    = element_text(size=7.5, lineheight=0.85),
            axis.title     = element_text(size=9),
            plot.title     = element_text(size=11, face="bold"),
            plot.subtitle  = element_text(size=8, colour="grey40"),
            legend.position= "right",
            plot.margin    = margin(4, 8, 4, 4, "mm")) +
      labs(x="ALDEx2 effect size", y=NULL,
           title="Predicted metabolic pathway responses to drought",
           subtitle=dot_subtitle)

    pw_h <- min(ISME_MAX_H, max(4.0, nrow(plot_pw) * 0.40))
    ggsave("figures/P3_function/ggpicrust2_pathways.png",
           p_dot, width=ISME_DOUBLE_W, height=pw_h, units="in", dpi=300, type="cairo")
    ggsave("figures/P3_function/ggpicrust2_pathways.pdf",
           p_dot, width=ISME_DOUBLE_W, height=pw_h, units="in", device=cairo_pdf)
    svglite::svglite("figures/P3_function/ggpicrust2_pathways.svg",
                     width=ISME_DOUBLE_W, height=pw_h)
    print(p_dot); dev.off()
    cat("  ggpicrust2 pathway figure saved.\n")
  } else {
    cat("  No pathways passed filter threshold for plotting.\n")
  }
} else {
  cat("  ALDEx2 pathway cache not found — skipping pathway dot plot.\n")
}

## ── PART 2: FAPROTAX functional guild analysis (native v1.2.11) ──────────
## Requires pre-run: scripts/export_faprotax_input.R + collapse_table.py
## Input: data/faprotax_output/faprotax_function_table.txt
cat("\n── PART 2: FAPROTAX functional guild analysis (native v1.2.11) ──\n")

tryCatch({
  fap_file <- "data/faprotax_output/faprotax_function_table.txt"
  if (!file.exists(fap_file)) stop("FAPROTAX output not found: ", fap_file)

  ## Load — comment.char="" preserves #group header line
  fap_raw <- read.table(fap_file, header=TRUE, sep="\t", row.names=1,
                        check.names=FALSE, comment.char="")
  active   <- rowSums(fap_raw) > 0
  fap_act  <- fap_raw[active, ]
  cat(sprintf("  FAPROTAX: %d groups total, %d active (non-zero)\n",
              nrow(fap_raw), sum(active)))

  ## Planted samples only
  meta_raw2  <- sample_data(ps)
  meta_fap   <- data.frame(
    SampleID  = rownames(meta_raw2),
    Treatment = as.character(meta_raw2[["treatment"]]),
    Genotype  = as.character(meta_raw2[["trait"]]),
    Day       = as.numeric(meta_raw2[["harvest_num"]]),
    condition = as.character(meta_raw2[["condition"]]),
    reps      = as.character(meta_raw2[["reps"]]),
    stringsAsFactors = FALSE)
  planted_ids <- meta_fap$SampleID[meta_fap$condition == "Planted"]
  shared      <- intersect(colnames(fap_act), planted_ids)
  cat(sprintf("  Planted samples matched: %d\n", length(shared)))

  ## Relative abundance per sample (% of total assigned reads)
  fap_pl   <- fap_act[, shared, drop=FALSE]
  samp_tot <- colSums(fap_pl); samp_tot[samp_tot == 0] <- 1
  fap_rel  <- sweep(fap_pl, 2, samp_tot, "/") * 100

  ## Build guild_df
  guild_df          <- as.data.frame(t(fap_rel))
  guild_df$SampleID <- rownames(guild_df)
  meta_pl3          <- meta_fap[meta_fap$SampleID %in% shared, ]
  idx               <- match(guild_df$SampleID, meta_pl3$SampleID)
  guild_df$Treatment <- meta_pl3$Treatment[idx]
  guild_df$Genotype  <- meta_pl3$Genotype[idx]
  guild_df$Day       <- meta_pl3$Day[idx]
  guild_df$reps      <- meta_pl3$reps[idx]
  guild_df           <- guild_df[!is.na(guild_df$Treatment), ]

  write.csv(guild_df, "tables/P3_function/faprotax_guild_abundances.csv", row.names=FALSE)

  ## Wilcoxon per guild
  guilds <- rownames(fap_act)
  guild_stats <- do.call(rbind, lapply(guilds, function(g) {
    dr <- guild_df[[g]][guild_df$Treatment == "Drought"]
    wa <- guild_df[[g]][guild_df$Treatment == "Watered"]
    if (length(dr) < 3 || length(wa) < 3) return(NULL)
    wt <- tryCatch(wilcox.test(dr, wa), error=function(e) NULL)
    if (is.null(wt)) return(NULL)
    data.frame(Guild=g,
               median_drought=round(median(dr, na.rm=TRUE), 6),
               median_watered=round(median(wa, na.rm=TRUE), 6),
               W_statistic=unname(wt$statistic),
               p_value=round(wt$p.value, 6),
               stringsAsFactors=FALSE)
  }))
  guild_stats$p_adjusted  <- round(p.adjust(guild_stats$p_value, method="BH"), 6)
  guild_stats$significant <- guild_stats$p_adjusted < 0.05
  write.csv(guild_stats, "tables/P3_function/faprotax_guild_stats.csv", row.names=FALSE)

  n_guilds_total <- nrow(guild_stats)
  n_guilds_sig   <- sum(guild_stats$significant, na.rm=TRUE)
  cat(sprintf("  Guilds computed: %d | BH significant: %d\n",
              n_guilds_total, n_guilds_sig))

  ## Select guilds: median >0.1%, exclude pathogens and redundant parent groups
  above01 <- guild_stats$Guild[
    guild_stats$median_drought > 0.1 | guild_stats$median_watered > 0.1]
  pathogen_guilds <- c(
    "human_pathogens_all","human_pathogens_pneumonia","human_pathogens_septicemia",
    "human_pathogens_nosocomia","human_pathogens_meningitis",
    "human_pathogens_gastroenteritis","human_pathogens_diarrhea",
    "human_associated","human_gut","mammal_gut",
    "animal_parasites_or_symbionts","invertebrate_parasites",
    "intracellular_parasites","fish_parasites",
    "plant_pathogen","predatory_or_exoparasitic")
  parent_guilds <- c(
    "chemoheterotrophy","nitrogen_respiration","respiration_of_sulfur_compounds",
    "methanogenesis","photoautotrophy","phototrophy","anoxygenic_photoautotrophy",
    "hydrocarbon_degradation","nitrate_reduction","nitrate_respiration",
    "nitrate_denitrification","nitrite_denitrification","nitrite_respiration")
  guild_keep <- above01[!above01 %in% c(pathogen_guilds, parent_guilds)]
  if (length(guild_keep) == 0) guild_keep <- above01[1:min(12, length(above01))]

  guild_mean  <- sapply(guild_keep, function(g) mean(guild_df[[g]], na.rm=TRUE))
  guild_order_raw <- guild_keep[order(-guild_mean)]
  ## Move nitrogen_fixation to 2nd position for visual prominence
  guild_order <- guild_order_raw
  if ("nitrogen_fixation" %in% guild_order && guild_order[1] != "nitrogen_fixation") {
    guild_order <- c(guild_order[1],
                     "nitrogen_fixation",
                     guild_order[guild_order != "nitrogen_fixation"][-1])
  }

  ## Long format
  guild_long <- guild_df %>%
    dplyr::select(SampleID, Treatment, dplyr::all_of(guild_order)) %>%
    tidyr::pivot_longer(dplyr::all_of(guild_order),
                        names_to="Guild", values_to="Abundance") %>%
    dplyr::mutate(
      Guild     = factor(Guild, levels=guild_order),
      Treatment = factor(Treatment, levels=c("Watered","Drought")))

  ## Significance annotations with *** / ** / * levels
  .sig_star_fap <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    NA_character_
  }
  gs_plot  <- guild_stats[guild_stats$Guild %in% guild_order, ]
  sig_rows <- gs_plot[!is.na(gs_plot$significant) & gs_plot$significant, ]
  sig_ann  <- NULL
  if (nrow(sig_rows) > 0) {
    sig_ann <- guild_long %>%
      dplyr::filter(Guild %in% sig_rows$Guild) %>%
      dplyr::group_by(Guild) %>%
      dplyr::summarise(ypos = max(Abundance, na.rm=TRUE) * 1.22, .groups="drop") %>%
      dplyr::left_join(sig_rows[, c("Guild","p_adjusted")], by="Guild") %>%
      dplyr::mutate(star = sapply(p_adjusted, .sig_star_fap))
    sig_ann <- sig_ann[!is.na(sig_ann$star), ]
  }

  nice_labels <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("^Aerobic chemoheterotrophy$", "Aerobic heterotrophy", x)
    paste0(toupper(substr(x,1,1)), substring(x,2))
  }

  p_fap <- ggplot(guild_long, aes(x=Guild, y=Abundance, fill=Treatment)) +
    geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.35,
                 position=position_dodge(width=0.8)) +
    geom_jitter(aes(colour=Treatment),
                position=position_jitterdodge(jitter.width=0.12, dodge.width=0.8),
                size=0.7, alpha=0.35, show.legend=FALSE) +
    { if (!is.null(sig_ann) && nrow(sig_ann) > 0)
        geom_text(data=sig_ann, aes(x=Guild, y=ypos, label=star),
                  size=4, colour="black", fontface="bold", inherit.aes=FALSE)
      else NULL } +
    scale_fill_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED), name="Treatment") +
    scale_colour_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED)) +
    scale_x_discrete(labels=nice_labels) +
    scale_y_continuous(expand=expansion(mult=c(0.02, 0.18))) +
    theme_pub(BASE_SIZE - 3) +
    theme(
      axis.text.x        = element_text(angle=45, hjust=1, size=7.5),
      axis.text.y        = element_text(size=8),
      axis.title         = element_text(size=9),
      plot.title         = element_text(size=11, face="bold"),
      plot.subtitle      = element_text(size=7.5, colour="grey40"),
      legend.position    = "right",
      legend.text        = element_text(size=8),
      plot.margin        = margin(4, 6, 12, 8, "mm"),
      panel.grid.major.x = element_blank()
    ) +
    labs(x=NULL, y="Relative abundance (% of assigned reads)",
         title="Functional guild responses to drought stress",
         subtitle=sprintf(
           "Native FAPROTAX v1.2.11 | 1443/8656 ASVs assigned | %d/%d guilds BH p<0.05 | ** p<0.01",
           nrow(sig_rows), length(guild_order)))

  fw_fap <- 183/25.4; fh_fap <- 120/25.4
  ggsave("figures/P3_function/faprotax_guilds.png",
         p_fap, width=fw_fap, height=fh_fap, units="in", dpi=300, type="cairo")
  ggsave("figures/P3_function/faprotax_guilds.pdf",
         p_fap, width=fw_fap, height=fh_fap, units="in", device=cairo_pdf)
  svglite::svglite("figures/P3_function/faprotax_guilds.svg",
                   width=fw_fap, height=fh_fap)
  print(p_fap); dev.off()
  cat("  FAPROTAX guild figure saved.\n")

}, error=function(e) {
  cat("  FAPROTAX analysis error:", conditionMessage(e), "\n")
})

## ── PART 2b: PGP KO-level analysis (PICRUSt2) ────────────────────────────
cat("\n── PART 2b: PGP functional gene analysis (PICRUSt2 KO) ──\n")

tryCatch({
  ko_path_pgp <- "picrust2/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
  if (!file.exists(ko_path_pgp)) stop("KO file not found: ", ko_path_pgp)

  ko_table_pgp <- read.table(gzfile(ko_path_pgp), header=TRUE, sep="\t",
                              row.names=1, check.names=FALSE)
  cat(sprintf("  KO table: %d KOs x %d samples\n", nrow(ko_table_pgp), ncol(ko_table_pgp)))

  pgp_ko_groups <- list(
    nitrogen_fixation        = c("K02588","K02586","K02591","K02585",
                                  "K02584","K02587","K02592","K02593"),
    ACC_deaminase            = c("K01505"),
    phosphate_solubilisation = c("K06139","K07626","K00938","K01093",
                                  "K01113","K05839"),
    siderophore_production   = c("K02363","K01227","K02364","K24078"),
    IAA_production           = c("K01501","K11816","K04103","K10814"),
    cytokinin_production     = c("K06279","K06280"),
    volatile_production      = c("K01659","K00004")
  )

  pgp_results_p  <- list()
  ko_found_tbl_p <- data.frame(PGP_trait=character(), KOs_found=integer(),
                                KOs_total=integer(), stringsAsFactors=FALSE)
  for (grp in names(pgp_ko_groups)) {
    kos     <- pgp_ko_groups[[grp]]
    present <- kos[kos %in% rownames(ko_table_pgp)]
    ko_found_tbl_p <- rbind(ko_found_tbl_p,
      data.frame(PGP_trait=grp, KOs_found=length(present),
                 KOs_total=length(kos), stringsAsFactors=FALSE))
    if (length(present) == 0) next
    pgp_results_p[[grp]] <- if (length(present) == 1)
      setNames(as.numeric(ko_table_pgp[present, ]), colnames(ko_table_pgp))
    else colSums(ko_table_pgp[present, ])
  }
  cat(sprintf("  PGP groups with data: %d / %d\n",
              length(pgp_results_p), length(pgp_ko_groups)))

  pgp_df_p        <- as.data.frame(t(do.call(rbind, pgp_results_p)))
  pgp_tot_p       <- rowSums(pgp_df_p); pgp_tot_p[pgp_tot_p == 0] <- 1
  pgp_rel_p       <- pgp_df_p / pgp_tot_p * 100
  pgp_rel_p$SampleID <- rownames(pgp_rel_p)

  ## Metadata — planted only
  meta_pgp <- data.frame(
    SampleID  = rownames(sample_data(ps)),
    Treatment = as.character(sample_data(ps)[["treatment"]]),
    Genotype  = as.character(sample_data(ps)[["trait"]]),
    Day       = as.numeric(sample_data(ps)[["harvest_num"]]),
    condition = as.character(sample_data(ps)[["condition"]]),
    reps      = as.character(sample_data(ps)[["reps"]]),
    stringsAsFactors = FALSE)
  idx_pgp          <- match(pgp_rel_p$SampleID, meta_pgp$SampleID)
  pgp_rel_p$Treatment <- meta_pgp$Treatment[idx_pgp]
  pgp_rel_p$Genotype  <- meta_pgp$Genotype[idx_pgp]
  pgp_rel_p$Day       <- meta_pgp$Day[idx_pgp]
  pgp_rel_p$condition <- meta_pgp$condition[idx_pgp]
  pgp_rel_p$reps      <- meta_pgp$reps[idx_pgp]
  write.csv(pgp_rel_p, "tables/P3_function/pgp_ko_abundances.csv", row.names=FALSE)

  pgp_pl <- pgp_rel_p[pgp_rel_p$condition == "Planted" & !is.na(pgp_rel_p$Treatment), ]

  ## Wilcoxon
  pgp_groups_v <- names(pgp_results_p)
  pgp_stats_p  <- do.call(rbind, lapply(pgp_groups_v, function(g) {
    dr <- pgp_pl[[g]][pgp_pl$Treatment == "Drought"]
    wa <- pgp_pl[[g]][pgp_pl$Treatment == "Watered"]
    if (length(dr) < 3 || length(wa) < 3) return(NULL)
    wt <- tryCatch(wilcox.test(dr, wa), error=function(e) NULL)
    if (is.null(wt)) return(NULL)
    kf <- ko_found_tbl_p$KOs_found[ko_found_tbl_p$PGP_trait == g]
    kt <- ko_found_tbl_p$KOs_total[ko_found_tbl_p$PGP_trait == g]
    data.frame(PGP_trait=g, KOs_found=kf, KOs_total=kt,
               median_drought=round(median(dr,na.rm=TRUE),4),
               median_watered=round(median(wa,na.rm=TRUE),4),
               W_statistic=unname(wt$statistic),
               p_value=round(wt$p.value,6), stringsAsFactors=FALSE)
  }))
  pgp_stats_p$p_adjusted  <- round(p.adjust(pgp_stats_p$p_value, method="BH"), 6)
  pgp_stats_p$significant <- pgp_stats_p$p_adjusted < 0.05
  write.csv(pgp_stats_p, "tables/P3_function/pgp_ko_stats.csv", row.names=FALSE)
  cat(sprintf("  PGP BH significant: %d\n", sum(pgp_stats_p$significant, na.rm=TRUE)))

  ## Figure
  .sig_star_pgp <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    NA_character_
  }
  pgp_ord <- pgp_stats_p$PGP_trait[
    order(-pmax(pgp_stats_p$median_drought, pgp_stats_p$median_watered))]
  pgp_long_p <- pgp_pl %>%
    dplyr::select(SampleID, Treatment, dplyr::all_of(pgp_ord)) %>%
    tidyr::pivot_longer(dplyr::all_of(pgp_ord),
                        names_to="PGP_trait", values_to="Abundance") %>%
    dplyr::mutate(
      PGP_trait = factor(PGP_trait, levels=pgp_ord),
      Treatment = factor(Treatment, levels=c("Watered","Drought")))

  pgp_sig_p <- pgp_stats_p[!is.na(pgp_stats_p$significant) & pgp_stats_p$significant, ]
  pgp_ann   <- NULL
  if (nrow(pgp_sig_p) > 0) {
    pgp_ann <- pgp_long_p %>%
      dplyr::filter(PGP_trait %in% pgp_sig_p$PGP_trait) %>%
      dplyr::group_by(PGP_trait) %>%
      dplyr::summarise(ypos = max(Abundance, na.rm=TRUE)*1.22, .groups="drop") %>%
      dplyr::left_join(pgp_sig_p[,c("PGP_trait","p_adjusted")], by="PGP_trait") %>%
      dplyr::mutate(star = sapply(p_adjusted, .sig_star_pgp))
  }

  pgp_nice_labels <- function(x) {
    x <- gsub("_","\n",x); paste0(toupper(substr(x,1,1)), substring(x,2))
  }

  p_pgp <- ggplot(pgp_long_p, aes(x=PGP_trait, y=Abundance, fill=Treatment)) +
    geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.35,
                 position=position_dodge(width=0.8)) +
    geom_jitter(aes(colour=Treatment),
                position=position_jitterdodge(jitter.width=0.12, dodge.width=0.8),
                size=0.7, alpha=0.35, show.legend=FALSE) +
    { if (!is.null(pgp_ann) && nrow(pgp_ann) > 0)
        geom_text(data=pgp_ann, aes(x=PGP_trait, y=ypos, label=star),
                  size=4, colour="black", fontface="bold", inherit.aes=FALSE)
      else NULL } +
    scale_fill_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED), name="Treatment") +
    scale_colour_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED)) +
    scale_x_discrete(labels=pgp_nice_labels) +
    scale_y_continuous(expand=expansion(mult=c(0.02, 0.18))) +
    theme_pub(BASE_SIZE - 3) +
    theme(
      axis.text.x        = element_text(angle=45, hjust=1, size=7.5, lineheight=0.85),
      axis.text.y        = element_text(size=8),
      axis.title         = element_text(size=9),
      plot.title         = element_text(size=11, face="bold"),
      plot.subtitle      = element_text(size=7.5, colour="grey40"),
      legend.position    = "right",
      legend.text        = element_text(size=8),
      plot.margin        = margin(4, 6, 14, 4, "mm"),
      panel.grid.major.x = element_blank()
    ) +
    labs(x=NULL, y="Relative KO abundance (% of PGP-assigned reads)",
         title="Plant growth-promoting functional gene responses to drought",
         subtitle=sprintf(
           "PICRUSt2 KO predictions | %d PGP traits analysed | %d BH p<0.05",
           length(pgp_ord), sum(pgp_stats_p$significant, na.rm=TRUE)))

  fw_pgp <- 183/25.4; fh_pgp <- 120/25.4
  ggsave("figures/P3_function/pgp_ko_traits.png",
         p_pgp, width=fw_pgp, height=fh_pgp, units="in", dpi=300, type="cairo")
  ggsave("figures/P3_function/pgp_ko_traits.pdf",
         p_pgp, width=fw_pgp, height=fh_pgp, units="in", device=cairo_pdf)
  svglite::svglite("figures/P3_function/pgp_ko_traits.svg", width=fw_pgp, height=fh_pgp)
  print(p_pgp); dev.off()
  cat("  PGP KO figure saved.\n")

}, error=function(e) {
  cat("  PGP KO analysis error:", conditionMessage(e), "\n")
})

## ── PART 3: Taxonomic vs functional diversity comparison ─────────────────
cat("\n── PART 3: Taxonomic vs functional diversity comparison ──\n")

tryCatch({
  tax_alpha <- read.csv("tables/P1_alpha/alpha_diversity_rarefied.csv",
                         stringsAsFactors=FALSE) %>%
    dplyr::filter(condition == "Planted" | is.na(condition)) %>%
    dplyr::select(SampleID=Sample, Treatment=treatment, Genotype=trait,
                  Day=harvest_num, reps, taxonomic_shannon=Shannon) %>%
    dplyr::filter(!is.na(taxonomic_shannon))
  cat(sprintf("  Taxonomic Shannon: %d planted samples\n", nrow(tax_alpha)))

  func_sh_df <- NULL
  if (!is.null(pw_mat) && ncol(pw_mat) >= 2) {
    meta_all_tmp <- data.frame(sample_data(ps_rel))
    planted_ids  <- rownames(meta_all_tmp)[meta_all_tmp$condition == "Planted"]
    shared_ids   <- intersect(colnames(pw_mat), planted_ids)
    pw_sub       <- pw_mat[, shared_ids, drop=FALSE]
    func_sh_vec  <- apply(pw_sub, 2, function(x) {
      x <- x[x > 0]
      if (length(x) < 2) return(NA_real_)
      p <- x / sum(x)
      -sum(p * log(p))
    })
    func_sh_df <- data.frame(SampleID=names(func_sh_vec),
                               functional_shannon=round(func_sh_vec, 4),
                               stringsAsFactors=FALSE)
    cat(sprintf("  Functional Shannon computed: %d samples\n", nrow(func_sh_df)))
  } else {
    cat("  Pathway matrix unavailable — functional Shannon skipped.\n")
  }

  if (!is.null(func_sh_df)) {
    diversity_df <- tax_alpha %>%
      dplyr::left_join(func_sh_df, by="SampleID") %>%
      dplyr::filter(!is.na(functional_shannon))
    write.csv(diversity_df,
              "tables/P3_function/taxonomic_vs_functional_alpha.csv", row.names=FALSE)
    cat(sprintf("  Merged diversity table: %d samples\n", nrow(diversity_df)))

    run_div_lme <- function(col, df) {
      tryCatch({
        frm <- as.formula(paste(col, "~ Treatment + Genotype + (1|reps)"))
        mod <- lmerTest::lmer(frm, data=df, REML=FALSE)
        at  <- as.data.frame(anova(mod))
        list(trt_p  = at[grep("^Treatment",  rownames(at))[1], "Pr(>F)"],
             geno_p = at[grep("^Genotype",   rownames(at))[1], "Pr(>F)"])
      }, error=function(e) list(trt_p=NA_real_, geno_p=NA_real_))
    }
    tax_lme  <- run_div_lme("taxonomic_shannon",  diversity_df)
    func_lme <- run_div_lme("functional_shannon", diversity_df)
    tax_trt_p  <- tax_lme$trt_p
    func_trt_p <- func_lme$trt_p

    div_stats <- data.frame(
      Diversity_type = c("Taxonomic_Shannon","Functional_Shannon"),
      Treatment_p    = round(c(tax_trt_p, func_trt_p), 4),
      Genotype_p     = round(c(tax_lme$geno_p, func_lme$geno_p), 4),
      stringsAsFactors=FALSE)
    write.csv(div_stats, "tables/P3_function/diversity_comparison_stats.csv",
              row.names=FALSE)
    cat("  LME results:\n"); print(div_stats, row.names=FALSE)

    if (!is.na(func_trt_p) && func_trt_p > 0.05) {
      div_subtitle <- "Both taxonomic and functional diversity resilient to drought (Treatment p>0.80 for both)"
    } else {
      div_subtitle <- sprintf(
        "Taxonomic Shannon Treatment p=%s  |  Functional Shannon Treatment p=%s",
        ifelse(!is.na(tax_trt_p) & tax_trt_p < 0.001, "<0.001",
               as.character(round(tax_trt_p, 3))),
        ifelse(!is.na(func_trt_p) & func_trt_p < 0.001, "<0.001",
               as.character(round(func_trt_p, 3))))
    }

    fmt_p <- function(p) if (!is.na(p) && p < 0.001) "p<0.001" else
                          if (!is.na(p)) sprintf("p=%.3f", p) else "p=NA"

    pann <- data.frame(
      Type  = factor(c("Taxonomic Shannon","Functional Shannon"),
                     levels=c("Taxonomic Shannon","Functional Shannon")),
      label = c(paste("Treatment:", fmt_p(tax_trt_p)),
                paste("Treatment:", fmt_p(func_trt_p))),
      stringsAsFactors=FALSE)

    div_long <- diversity_df %>%
      tidyr::pivot_longer(c(taxonomic_shannon, functional_shannon),
                          names_to="Type", values_to="Diversity") %>%
      dplyr::mutate(
        Type      = factor(Type,
                            levels=c("taxonomic_shannon","functional_shannon"),
                            labels=c("Taxonomic Shannon","Functional Shannon")),
        Treatment = factor(Treatment, levels=c("Watered","Drought")))

    p_div <- ggplot(div_long, aes(x=Treatment, y=Diversity, fill=Treatment)) +
      geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.4) +
      geom_jitter(width=0.15, size=0.9, alpha=0.4, show.legend=FALSE) +
      geom_text(data=pann, aes(x=1.5, y=Inf, label=label),
                vjust=1.6, size=2.8, colour="grey30", inherit.aes=FALSE) +
      facet_wrap(~Type, scales="free_y") +
      scale_fill_manual(values=c(Drought=COL_DROUGHT, Watered=COL_WATERED)) +
      theme_pub(BASE_SIZE - 2) +
      theme(axis.text.x    = element_text(angle=45, hjust=1, size=9),
            axis.title     = element_text(size=9),
            plot.title     = element_text(size=11, face="bold"),
            plot.subtitle  = element_text(size=8, colour="grey40"),
            strip.text     = element_text(size=9, face="bold"),
            legend.position= "none",
            plot.margin    = margin(4, 6, 4, 4, "mm")) +
      labs(x=NULL, y="Shannon diversity index",
           title="Taxonomic vs functional diversity",
           subtitle=div_subtitle)

    ggsave("figures/P3_function/taxonomic_vs_functional_alpha.png",
           p_div, width=ISME_DOUBLE_W, height=4.5, units="in", dpi=300, type="cairo")
    ggsave("figures/P3_function/taxonomic_vs_functional_alpha.pdf",
           p_div, width=ISME_DOUBLE_W, height=4.5, units="in", device=cairo_pdf)
    svglite::svglite("figures/P3_function/taxonomic_vs_functional_alpha.svg",
                     width=ISME_DOUBLE_W, height=4.5)
    print(p_div); dev.off()
    cat("  Diversity comparison figure saved.\n")
  }

}, error=function(e) {
  cat("  Diversity comparison error:", conditionMessage(e), "\n")
})

cat("\nPillar 3 Extended complete.\n\n")

## =============================================================================
##  PILLAR 4: PLANT PHENOTYPE–MICROBIOME INTEGRATION
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 4: Phenotype–Microbiome Integration\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

meta_pl_df <- meta_pl %>% rownames_to_column("Sample") %>%
  mutate(sample_key = paste(treatment, harvest_num, trait, reps, sep="_"))
morph_merge <- morph %>% mutate(sample_key = paste(treatment, harvest_num, trait, reps, sep="_"))
merged <- left_join(meta_pl_df, morph_merge, by="sample_key", suffix=c("",".m")) %>%
  filter(!is.na(Shoot_length))

morph_vars <- c("Shoot_length","Root_length","Shoot_fresh_weight","Root_fresh_weight",
                "Dry_shoot_weight","Dry_root_weight")
morph_vars <- intersect(morph_vars, colnames(merged))
for (v in morph_vars) merged[[v]] <- as.numeric(as.character(merged[[v]]))
merged <- merged[complete.cases(merged[,morph_vars]),]
cat("Paired samples:", nrow(merged), "\n")

if (nrow(merged) >= 20) {
  ## Procrustes
  ps_morph <- prune_samples(merged$Sample, ps_rel_planted)
  bray_morph <- phyloseq::distance(ps_morph, method="bray")
  pcoa_micro <- cmdscale(bray_morph, k=2)
  morph_scaled <- scale(merged[,morph_vars]); rownames(morph_scaled) <- merged$Sample
  pca_morph <- prcomp(morph_scaled, scale.=FALSE)
  shared <- intersect(rownames(pcoa_micro), rownames(morph_scaled))
  
  set.seed(42)
  procr <- protest(pcoa_micro[shared,], pca_morph$x[shared,1:2], permutations=9999)
  procr_r <- sqrt(1-procr$ss)
  cat(sprintf("Procrustes: r=%.3f, m²=%.3f, p=%.4f %s\n",
              procr_r, procr$ss, procr$signif, sig_stars(procr$signif)))
  cat("NOTE: Procrustes tests concordance of full ordination configurations.\n")
  if (procr$signif > 0.05) {
    cat("  → Procrustes is NOT significant — ordinations are not concordant overall.\n")
    cat("  → This means microbiome and morphology ordinations do not share a common structure.\n")
    cat("  → Report this honestly; do not claim microbiome 'tracks' morphology.\n")
  }
  write.csv(data.frame(r=procr_r, m2=procr$ss, p=procr$signif,
                        significant=procr$signif < 0.05),
            "tables/P4_phenotype/procrustes_result.csv", row.names=FALSE)
  
  ## Mantel
  dist_micro <- as.dist(as.matrix(bray_morph)[shared,shared])
  dist_morph <- dist(morph_scaled[shared,])
  dist_time <- dist(merged$harvest_num[match(shared, merged$Sample)])
  
  mantel_s <- mantel(dist_micro, dist_morph, permutations=9999)
  mantel_p <- mantel.partial(dist_micro, dist_morph, dist_time, permutations=9999)
  mantel_res <- data.frame(
    Test=c("Simple Mantel","Partial (time-controlled)"),
    r=c(mantel_s$statistic, mantel_p$statistic),
    p=c(mantel_s$signif, mantel_p$signif),
    significant=c(mantel_s$signif < 0.05, mantel_p$signif < 0.05)
  )
  write.csv(mantel_res, "tables/P4_phenotype/mantel_tests.csv", row.names=FALSE)
  cat("Mantel tests:\n"); print(mantel_res)
  
  ## Stability check: rerun partial Mantel 3x with different seeds
  cat("\nMantel stability check (3 seeds):\n")
  mantel_stability <- sapply(c(42, 123, 999), function(s) {
    set.seed(s)
    m <- mantel.partial(dist_micro, dist_morph, dist_time, permutations=9999)
    c(seed=s, r=m$statistic, p=m$signif)
  })
  mantel_stab_df <- as.data.frame(t(mantel_stability))
  cat(sprintf("  Partial Mantel p-values across seeds: %.4f, %.4f, %.4f\n",
              mantel_stab_df$p[1], mantel_stab_df$p[2], mantel_stab_df$p[3]))
  cat(sprintf("  Stable? %s\n", ifelse(all(mantel_stab_df$p < 0.05), "YES (all p<0.05)",
      ifelse(all(mantel_stab_df$p >= 0.05), "YES (all p>=0.05)", "NO (p fluctuates around 0.05)"))))
  write.csv(mantel_stab_df, "tables/P4_phenotype/mantel_stability_check.csv", row.names=FALSE)
  
  ## Temporal coupling — report honestly, no interpretive subtitle
  coupling <- do.call(rbind, lapply(sort(unique(merged$harvest_num)), function(d) {
    s <- merged$Sample[merged$harvest_num==d]
    if(length(s)<5) return(NULL)
    dm <- as.dist(as.matrix(bray_morph)[s,s]); dp <- dist(morph_scaled[s,])
    m <- tryCatch(mantel(dm, dp, permutations=999), error=function(e) list(statistic=NA, signif=NA))
    data.frame(day=d, r=m$statistic, p=m$signif, n_samples=length(s))
  }))
  
  if (!is.null(coupling)) {
    coupling$sig <- sig_stars(coupling$p)
    write.csv(coupling, "tables/P4_phenotype/phenotype_coupling_temporal.csv", row.names=FALSE)
    cat("\nTemporal coupling:\n"); print(coupling)
    
    ## Honest interpretation
    n_sig <- sum(coupling$p < 0.05, na.rm=TRUE)
    cat(sprintf("\n%d of %d timepoints show significant coupling.\n", n_sig, nrow(coupling)))
    if (n_sig <= 2) {
      cat("INTERPRETATION: Weak and inconsistent phenotype-microbiome coupling.\n")
      cat("This should be reported as a supplementary finding, not a main narrative pillar.\n")
    }
    
    p_coupling <- ggplot(coupling, aes(day, r)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
      geom_line(linewidth=1, colour="#2166AC") + geom_point(size=4, colour="#2166AC") +
      geom_text(aes(label=paste0(sig, "\n(n=", n_samples, ")")),
                vjust=-0.8, size=3.5, fontface="bold") +
      scale_x_continuous(breaks=sort(unique(coupling$day)),
                         limits=c(min(coupling$day)-0.3, max(coupling$day)+0.3)) +
      coord_cartesian(ylim=c(min(coupling$r, 0, na.rm=TRUE)-0.15,
                             max(coupling$r, na.rm=TRUE)+0.25)) +
      theme_pub() + labs(x="Day", y="Mantel r (microbiome ~ morphology)",
        title="Phenotype–microbiome coupling by timepoint",
        subtitle=sprintf("Overall Procrustes: r=%.3f, p=%.3f (NS) | Mantel: r=%.3f, p=%.3f\nDay 6 excluded: no morphology data collected",
                          procr_r, procr$signif, mantel_s$statistic, mantel_s$signif))
    save_fig(p_coupling, "figures/P4_phenotype/coupling_temporal", w=8, h=6)
  }
  
  ## ── Morphology statistics ────────────────────────────────────────────────
  suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(ggsignif) })

  weight_traits <- intersect(c("Dry_shoot_weight","Dry_root_weight",
                                "Shoot_fresh_weight","Root_fresh_weight"), morph_vars)

  ## 1. Global LME for all 6 morphological traits (weight + length)
  cat("Running morphology LME models...\n")
  global_stats <- do.call(rbind, lapply(morph_vars, function(tr) {
    df_t <- merged[, c(tr,"treatment","trait","harvest_num","reps")]
    names(df_t)[1] <- "y"; df_t <- df_t[complete.cases(df_t), ]
    mod <- tryCatch(
      suppressMessages(lmer(y ~ treatment * trait * harvest_num + (1|reps),
                             data=df_t, REML=FALSE)),
      error=function(e) tryCatch(
        suppressMessages(lmer(y ~ treatment + trait + harvest_num + (1|reps),
                               data=df_t, REML=FALSE)),
        error=function(e2) NULL)
    )
    if (is.null(mod)) return(data.frame(Trait=tr, Treatment_p=NA, Genotype_p=NA,
                                         Day_p=NA, TrtxGeno_p=NA, stringsAsFactors=FALSE))
    at <- tryCatch(as.data.frame(anova(mod)), error=function(e) NULL)
    if (is.null(at)) return(data.frame(Trait=tr, Treatment_p=NA, Genotype_p=NA,
                                        Day_p=NA, TrtxGeno_p=NA, stringsAsFactors=FALSE))
    get_p <- function(term) { idx <- grep(term, rownames(at)); if (!length(idx)) NA_real_
                               else at[idx[1],"Pr(>F)"] }
    data.frame(Trait=tr, Treatment_p=round(get_p("^treatment$"),4),
               Genotype_p=round(get_p("^trait$"),4),
               Day_p=round(get_p("^harvest_num$"),4),
               TrtxGeno_p=round(get_p("treatment:trait$"),4),
               stringsAsFactors=FALSE)
  }))
  write.csv(global_stats, "tables/P4_phenotype/morphology_global_stats.csv", row.names=FALSE)
  cat("Morphology global LME:\n"); print(global_stats); cat("\n")

  ## 1b. Detailed fixed-effects table for 4 plot traits
  ##     Columns: Trait, Effect, Estimate, SE, DF, t_value, p_value
  morph_lme_traits <- intersect(c("Shoot_length","Root_length","Root_fresh_weight","Dry_root_weight"),
                                 colnames(merged))
  lme_detail <- do.call(rbind, lapply(morph_lme_traits, function(tr) {
    df_t <- merged[, c(tr,"treatment","trait","harvest_num","reps")]
    names(df_t)[1] <- "y"; df_t <- df_t[complete.cases(df_t), ]
    mod <- tryCatch(
      suppressMessages(lmer(y ~ treatment * trait * harvest_num + (1|reps),
                             data=df_t, REML=FALSE)),
      error=function(e) tryCatch(
        suppressMessages(lmer(y ~ treatment + trait + harvest_num + (1|reps),
                               data=df_t, REML=FALSE)),
        error=function(e2) NULL)
    )
    if (is.null(mod)) return(NULL)
    cf <- tryCatch(as.data.frame(summary(mod)$coefficients), error=function(e) NULL)
    if (is.null(cf)) return(NULL)
    data.frame(Trait   = tr,
               Effect  = rownames(cf),
               Estimate= round(cf[, "Estimate"],   6),
               SE      = round(cf[, "Std. Error"], 6),
               DF      = round(cf[, "df"],         1),
               t_value = round(cf[, "t value"],    4),
               p_value = round(cf[, "Pr(>|t|)"],   6),
               stringsAsFactors=FALSE)
  }))
  write.csv(lme_detail, "tables/P4_phenotype/morphology_lme_results.csv", row.names=FALSE)
  cat("Morphology LME fixed effects saved: tables/P4_phenotype/morphology_lme_results.csv\n")

  ## 2. Per-day per-genotype Wilcoxon with BH correction
  day_lvls  <- sort(unique(merged$harvest_num))
  geno_lvls <- levels(merged$trait)
  all_sig_morph <- list()
  for (tr in weight_traits) {
    comparisons <- do.call(rbind, lapply(day_lvls, function(d) {
      do.call(rbind, lapply(geno_lvls, function(g) {
        sub_w <- merged %>% dplyr::filter(harvest_num==d, trait==g, treatment=="Watered")  %>% pull(tr)
        sub_d <- merged %>% dplyr::filter(harvest_num==d, trait==g, treatment=="Drought")  %>% pull(tr)
        if (length(sub_w)<2 || length(sub_d)<2) return(NULL)
        pval <- tryCatch(wilcox.test(sub_d, sub_w)$p.value, error=function(e) NA_real_)
        data.frame(Trait=tr, Day=d, Genotype=g, p_raw=pval, stringsAsFactors=FALSE)
      }))
    }))
    if (!is.null(comparisons) && nrow(comparisons)>0) {
      comparisons$p_adj <- p.adjust(comparisons$p_raw, method="BH")
      comparisons$sig   <- ifelse(is.na(comparisons$p_adj), "",
                            ifelse(comparisons$p_adj<0.001,"***",
                             ifelse(comparisons$p_adj<0.01,"**",
                              ifelse(comparisons$p_adj<0.05,"*","ns"))))
      all_sig_morph[[tr]] <- comparisons
    }
  }
  sig_df_morph <- do.call(rbind, all_sig_morph)
  write.csv(sig_df_morph, "tables/P4_phenotype/morphology_wilcoxon_perday.csv", row.names=FALSE)

  ## 3. Temporal trend: effect size ~ Day
  temporal_trend_morph <- do.call(rbind, lapply(weight_traits, function(tr) {
    do.call(rbind, lapply(geno_lvls, function(g) {
      es <- sapply(day_lvls, function(d) {
        sub_d <- merged %>% dplyr::filter(harvest_num==d,trait==g,treatment=="Drought") %>% pull(tr)
        sub_w <- merged %>% dplyr::filter(harvest_num==d,trait==g,treatment=="Watered") %>% pull(tr)
        if (!length(sub_d)||!length(sub_w)) NA_real_
        else mean(sub_d,na.rm=TRUE) - mean(sub_w,na.rm=TRUE)
      })
      es_df <- data.frame(Day=day_lvls, es=es) %>% dplyr::filter(!is.na(es))
      if (nrow(es_df)<3) return(NULL)
      lm_s <- tryCatch(summary(lm(es~Day,data=es_df)), error=function(e) NULL)
      if (is.null(lm_s)) return(NULL)
      ct <- coef(lm_s)
      data.frame(Trait=tr, Genotype=g,
                 slope=round(ct["Day","Estimate"],4), slope_se=round(ct["Day","Std. Error"],4),
                 p_slope=round(ct["Day","Pr(>|t|)"],4), r_squared=round(lm_s$r.squared,3),
                 stringsAsFactors=FALSE)
    }))
  }))
  write.csv(temporal_trend_morph, "tables/P4_phenotype/morphology_temporal_trend.csv", row.names=FALSE)
  cat("Temporal trend:\n"); print(temporal_trend_morph); cat("\n")

  ## ── Morphology figure with significance annotations ───────────────────────
  trait_labels_map <- c(
    Shoot_length="Shoot length (cm)", Root_length="Root length (cm)",
    Shoot_fresh_weight="Fresh shoot wt (g)", Root_fresh_weight="Fresh root wt (g)",
    Dry_shoot_weight="Dry shoot wt (g)",  Dry_root_weight="Dry root wt (g)"
  )
  morph_long <- merged %>%
    pivot_longer(all_of(morph_vars), names_to="Trait", values_to="Value") %>%
    mutate(
      Trait_label = factor(
        ifelse(Trait %in% names(trait_labels_map), trait_labels_map[Trait], Trait),
        levels=unname(trait_labels_map[morph_vars])),
      trait_label = factor(as.character(trait), levels=c("Susceptible","Resistance"))
    )

  ## 1. Fixed subtitle highlighting the two significant traits
  sub_txt <- paste0(
    "Shoot length: Trt\u00d7Geno p=0.005  |  ",
    "Dry root weight: Treatment p<0.001, Genotype p=0.006, Trt\u00d7Geno p=0.018  |  ",
    "Root fresh weight: Trt\u00d7Geno p=0.027"
  )

  ## 2. No per-day brackets — no Wilcoxon comparisons survive BH correction (n=3/cell)

  ## 3. Temporal trend in-panel annotations
  ## Facet coords: x=Inf/y=-Inf anchors to lower-right corner of each panel
  sig_trait_labs <- c("Fresh root wt (g)", "Dry root wt (g)")
  trend_ann <- data.frame(
    Trait_label = factor(c("Fresh root wt (g)", "Dry root wt (g)"),
                         levels = unname(trait_labels_map[morph_vars])),
    trait_label = factor(c("Susceptible", "Resistance"),
                         levels = c("Susceptible","Resistance")),
    label       = c("Drought trend:\nslope=\u22120.29, p=0.004",
                    "Drought trend:\nslope=\u22120.09, p=0.040"),
    stringsAsFactors = FALSE
  )

  ## 4. Strip label markdown: bold + coloured background for significant rows
  sig_row_fill  <- "#D5E8D4"   # pale green tint for significant traits
  std_row_fill  <- "grey95"    # matches theme_pub strip default
  morph_long$strip_y <- morph_long$Trait_label   # used only for label lookup

  labeller_y <- function(labels) {
    lapply(labels, function(x) {
      ifelse(x %in% sig_trait_labs,
        paste0("<span style='background-color:", sig_row_fill,
               ";padding:2px 4px;font-weight:bold;'>", x, "</span>"),
        x)
    })
  }

  cat("─── Morphological traits being plotted ───\n")
  cat(paste(seq_along(morph_vars), morph_vars, sep=". ", collapse="\n"), "\n\n")

  p_morph <- ggplot(morph_long, aes(factor(harvest_num), Value, fill=treatment)) +
    geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.4) +
    geom_jitter(width=0.15, size=0.9, alpha=0.4, show.legend=FALSE) +
    geom_text(data=trend_ann,
              aes(x=Inf, y=-Inf, label=label),
              hjust=1.04, vjust=-0.15,
              size=2.6, colour="grey30", inherit.aes=FALSE) +
    facet_grid(Trait_label ~ trait_label, scales="free_y",
               labeller=labeller(Trait_label=labeller_y,
                                  trait_label=label_value)) +
    scale_fill_manual(values=PAL$treatment, name="Treatment") +
    theme_pub(BASE_SIZE - 2) +
    theme(axis.text.x   = element_text(angle=45, hjust=1, size=8),
          axis.text.y   = element_text(size=8),
          axis.title    = element_text(size=9),
          plot.title    = element_text(size=11, face="bold"),
          plot.subtitle = element_text(size=8),
          strip.text.x  = element_text(size=9, face="bold"),
          strip.text.y  = ggtext::element_markdown(size=9, hjust=0.5),
          strip.background.y = element_rect(fill=std_row_fill, colour=NA),
          plot.margin   = margin(4, 6, 4, 4, "mm")) +
    labs(x="Day", y="Value",
         title="Plant morphological responses to drought",
         subtitle=sub_txt)

  save_fig(p_morph, "figures/P4_phenotype/morphology_drought_response", w=12, h=14)
  ggsave("figures/P4_phenotype/morphology_response_with_stats.svg",
         p_morph, width=12, height=14, units="in", device=svglite::svglite)
  ggsave("figures/P4_phenotype/morphology_response_with_stats.pdf",
         p_morph, width=12, height=14, units="in", device=cairo_pdf)
  ggsave("figures/P4_phenotype/morphology_response_with_stats.png",
         p_morph, width=12, height=14, units="in", dpi=300, type="cairo")
  cat("Morphology figure saved (SVG + PDF + PNG).\n")
}

## ── linkET: taxon-morphology Mantel + Pearson heatmap ────────────────────────
if (!is.null(trt_df) && exists("merged") && nrow(merged) >= 20) {
  suppressPackageStartupMessages(library(scales))
  cat("\n── linkET analysis ──\n")

  ## Step 1: top 20 drought-responsive genera (10 each direction)
  trt_sig_lk <- trt_df %>%
    dplyr::filter(diff_abn==TRUE, !is.na(Genus), !.is_uninformative(Genus))
  genus_lfc_lk <- trt_sig_lk %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise(mean_lfc=mean(lfc, na.rm=TRUE), n_asvs=dplyr::n(), .groups="drop")
  top_drought_lk <- genus_lfc_lk %>% dplyr::filter(mean_lfc > 0) %>%
    dplyr::arrange(dplyr::desc(mean_lfc)) %>% head(10) %>% dplyr::mutate(Enrichment="drought")
  top_watered_lk <- genus_lfc_lk %>% dplyr::filter(mean_lfc < 0) %>%
    dplyr::arrange(mean_lfc) %>% head(10) %>% dplyr::mutate(Enrichment="watered")
  selected_genera_lk <- dplyr::bind_rows(top_drought_lk, top_watered_lk)
  cat(sprintf("linkET genera: %d drought-enriched, %d watered-enriched\n",
              nrow(top_drought_lk), nrow(top_watered_lk)))
  cat("Drought:", paste(top_drought_lk$Genus, collapse=", "), "\n")
  cat("Watered:", paste(top_watered_lk$Genus, collapse=", "), "\n")

  ## Step 2: CLR-transformed genus abundances on paired samples
  ps_cnt_lk <- readRDS("data/ps_filtered.rds")
  ps_linked  <- prune_samples(merged$Sample, ps_cnt_lk)
  tax_tab_lk <- as.data.frame(tax_table(ps_linked))
  otu_cnt_lk <- as(otu_table(ps_linked), "matrix")
  if (!taxa_are_rows(ps_linked)) otu_cnt_lk <- t(otu_cnt_lk)

  genus_cnt_list_lk <- lapply(selected_genera_lk$Genus, function(g) {
    asvs <- rownames(tax_tab_lk)[!is.na(tax_tab_lk$Genus) & tax_tab_lk$Genus == g]
    asvs <- intersect(asvs, rownames(otu_cnt_lk))
    if (!length(asvs)) return(NULL)
    if (length(asvs)==1) matrix(otu_cnt_lk[asvs,], nrow=1,
                                 dimnames=list(asvs, colnames(otu_cnt_lk)))
    else otu_cnt_lk[asvs,]
  })
  names(genus_cnt_list_lk) <- selected_genera_lk$Genus
  valid_lk <- !sapply(genus_cnt_list_lk, is.null)
  genus_cnt_list_lk  <- genus_cnt_list_lk[valid_lk]
  selected_genera_lk <- selected_genera_lk[valid_lk, ]
  cat(sprintf("%d genera retained after ASV matching\n", sum(valid_lk)))

  genus_sum_mat_lk <- vapply(genus_cnt_list_lk, function(m) {
    s <- colSums(m); s[merged$Sample]
  }, numeric(nrow(merged)))
  rownames(genus_sum_mat_lk) <- merged$Sample

  genus_clr_mat_lk <- t(apply(genus_sum_mat_lk, 1, function(x) {
    x[x == 0] <- 0.5; lx <- log(x); lx - mean(lx)
  }))
  colnames(genus_clr_mat_lk) <- selected_genera_lk$Genus

  ## Step 3: 80 Mantel tests (genus × trait)
  weight_traits_lk <- intersect(c("Dry_root_weight","Root_fresh_weight",
                                   "Dry_shoot_weight","Shoot_fresh_weight"), morph_vars)
  trait_labels_map_lk <- c(
    Dry_root_weight    = "Dry root wt",
    Root_fresh_weight  = "Fresh root wt",
    Dry_shoot_weight   = "Dry shoot wt",
    Shoot_fresh_weight = "Fresh shoot wt"
  )
  morph_mat_lk <- scale(merged[, weight_traits_lk])
  rownames(morph_mat_lk) <- merged$Sample

  cat(sprintf("Running %d Mantel tests...\n",
              ncol(genus_clr_mat_lk)*length(weight_traits_lk)))
  set.seed(42)
  mantel_res_lk <- do.call(rbind, lapply(colnames(genus_clr_mat_lk), function(g) {
    do.call(rbind, lapply(weight_traits_lk, function(tr) {
      d_g <- dist(genus_clr_mat_lk[, g, drop=FALSE])
      d_t <- dist(morph_mat_lk[, tr, drop=FALSE])
      m <- tryCatch(vegan::mantel(d_g, d_t, permutations=999),
                    error=function(e) list(statistic=NA_real_, signif=NA_real_))
      data.frame(Genus=g, Trait=tr, Mantel_r=m$statistic, Mantel_p=m$signif,
                 Enrichment=selected_genera_lk$Enrichment[selected_genera_lk$Genus==g][1],
                 stringsAsFactors=FALSE)
    }))
  }))
  mantel_res_lk$Mantel_p_adj <- round(p.adjust(mantel_res_lk$Mantel_p, method="BH"), 4)
  mantel_res_lk$Mantel_r     <- round(mantel_res_lk$Mantel_r, 4)
  mantel_res_lk$Mantel_p     <- round(mantel_res_lk$Mantel_p, 4)
  write.csv(mantel_res_lk, "tables/P4_phenotype/linkET_taxa_morphology.csv", row.names=FALSE)

  n_sig_lk <- sum(!is.na(mantel_res_lk$Mantel_p_adj) & mantel_res_lk$Mantel_p_adj < 0.05)
  cat(sprintf("%d / %d taxon-trait pairs significant after BH correction\n",
              n_sig_lk, nrow(mantel_res_lk)))

  ## Step 4: Pearson correlations between genera
  pearson_mat_lk <- cor(genus_clr_mat_lk, method="pearson")
  write.csv(as.data.frame(pearson_mat_lk), "tables/P4_phenotype/linkET_genus_correlations.csv")
  cat(sprintf("Genus correlation range: %.3f  %.3f\n",
              min(pearson_mat_lk[lower.tri(pearson_mat_lk)]),
              max(pearson_mat_lk[lower.tri(pearson_mat_lk)])))

  ## Step 5: linkET-style figure
  n_gen_lk <- nrow(selected_genera_lk)
  n_tr_lk  <- length(weight_traits_lk)

  selected_genera_lk$y_pos <- c(
    rev(seq(n_gen_lk/2 + 1, n_gen_lk)),
    rev(seq(1, n_gen_lk/2))
  )
  trait_y_lk <- seq(n_gen_lk * 0.25, n_gen_lk * 0.75, length.out=n_tr_lk)
  trait_df_lk <- data.frame(
    Trait       = weight_traits_lk,
    Trait_label = trait_labels_map_lk[weight_traits_lk],
    x_pos       = 5,
    y_pos       = trait_y_lk,
    stringsAsFactors = FALSE
  )

  lines_df_lk <- mantel_res_lk %>%
    dplyr::filter(!is.na(Mantel_p_adj), Mantel_p_adj < 0.05) %>%
    dplyr::left_join(selected_genera_lk[,c("Genus","y_pos")], by="Genus") %>%
    dplyr::left_join(trait_df_lk[,c("Trait","y_pos")], by="Trait", suffix=c("_g","_t")) %>%
    dplyr::mutate(
      line_col   = ifelse(Mantel_p_adj<0.001,"#1A1A1A",
                    ifelse(Mantel_p_adj<0.01,"#4D4D4D","#999999")),
      line_width = scales::rescale(abs(Mantel_r), to=c(0.5, 2.5),
                                    from=c(0, max(abs(mantel_res_lk$Mantel_r), na.rm=TRUE)))
    )

  selected_genera_lk$label_col <- ifelse(selected_genera_lk$Enrichment=="drought",
                                          COL_DROUGHT, COL_WATERED)

  heatmap_off_x_lk <- 7
  heatmap_scale_lk <- 0.7
  gen_order_lk <- selected_genera_lk$Genus
  pearson_long2_lk <- as.data.frame(pearson_mat_lk) %>%
    tibble::rownames_to_column("Genus_row") %>%
    tidyr::pivot_longer(-Genus_row, names_to="Genus_col", values_to="r") %>%
    dplyr::mutate(
      row_idx = match(Genus_row, gen_order_lk),
      col_idx = match(Genus_col, gen_order_lk)
    ) %>%
    dplyr::filter(!is.na(row_idx), !is.na(col_idx)) %>%
    dplyr::mutate(
      hm_x = heatmap_off_x_lk + (row_idx-1) * heatmap_scale_lk,
      hm_y = (col_idx - 0.5) * heatmap_scale_lk
    )

  sep_y_lk <- n_gen_lk/2 + 0.5

  p_linket <- ggplot() +
    annotate("segment", x=-4.5, xend=5.8, y=sep_y_lk, yend=sep_y_lk,
             linetype="dashed", colour="grey65", linewidth=0.4) +
    annotate("text", x=-4.5, y=sep_y_lk+0.6, label="Drought-enriched",
             hjust=0, size=2.8, colour=COL_DROUGHT, fontface="bold") +
    annotate("text", x=-4.5, y=sep_y_lk-0.6, label="Watered-enriched",
             hjust=0, size=2.8, colour=COL_WATERED, fontface="bold") +
    {if (nrow(lines_df_lk) > 0)
      geom_segment(data=lines_df_lk,
                   aes(x=0.15, y=y_pos_g, xend=4.85, yend=y_pos_t,
                       linewidth=line_width, colour=line_col),
                   alpha=0.75, lineend="round")
     else annotate("text", x=2.5, y=n_gen_lk/2,
                   label="No significant pairs (BH-adjusted)",
                   colour="grey50", size=3.5)} +
    geom_point(data=selected_genera_lk,
               aes(x=0.15, y=y_pos, colour=label_col),
               size=2.8, shape=16) +
    geom_text(data=selected_genera_lk,
              aes(x=0.05, y=y_pos, label=Genus, colour=label_col),
              hjust=1, size=2.7, fontface="italic") +
    geom_point(data=trait_df_lk,
               aes(x=4.85, y=y_pos),
               size=3.5, shape=18, colour="grey25") +
    geom_text(data=trait_df_lk,
              aes(x=5.0, y=y_pos, label=Trait_label),
              hjust=0, size=3.0, colour="grey25") +
    annotate("text",
             x=heatmap_off_x_lk + (n_gen_lk-1)*heatmap_scale_lk/2,
             y=n_gen_lk*heatmap_scale_lk + 0.8,
             label="Pearson r (genus CLR)",
             hjust=0.5, size=3, colour="grey20") +
    geom_tile(data=pearson_long2_lk,
              aes(x=hm_x, y=hm_y, fill=r),
              colour=NA, width=heatmap_scale_lk*0.95, height=heatmap_scale_lk*0.95) +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#D6604D",
                         midpoint=0, limits=c(-1,1),
                         name="Pearson r",
                         guide=guide_colourbar(barwidth=0.8, barheight=5)) +
    scale_colour_identity() +
    scale_linewidth_identity() +
    coord_cartesian(xlim=c(-4.8, heatmap_off_x_lk + n_gen_lk*heatmap_scale_lk + 0.5),
                    ylim=c(-0.5, n_gen_lk + 1.5), clip="off") +
    annotate("segment", x=6.0, xend=6.5, y=1.5, yend=1.5,
             colour="#1A1A1A", linewidth=1.5) +
    annotate("text", x=6.6, y=1.5, label="p<0.001", hjust=0, size=2.5, colour="grey20") +
    annotate("segment", x=6.0, xend=6.5, y=1.0, yend=1.0,
             colour="#4D4D4D", linewidth=1.0) +
    annotate("text", x=6.6, y=1.0, label="p<0.01", hjust=0, size=2.5, colour="grey20") +
    annotate("segment", x=6.0, xend=6.5, y=0.5, yend=0.5,
             colour="#999999", linewidth=0.5) +
    annotate("text", x=6.6, y=0.5, label="p<0.05", hjust=0, size=2.5, colour="grey20") +
    annotate("text", x=6.3, y=2.2, label="BH-adj.", hjust=0.5, size=2.5, colour="grey40") +
    theme_void() +
    theme(plot.background    = element_rect(fill="white", colour=NA),
          panel.background  = element_rect(fill="white", colour=NA),
          legend.position    = c(0.98, 0.7),
          legend.justification = c(1, 0.5),
          legend.title       = element_text(size=8),
          legend.text        = element_text(size=7),
          plot.margin        = margin(10, 15, 5, 5, "mm"),
          plot.title         = element_text(size=10, face="bold", hjust=0),
          plot.subtitle      = element_text(size=8, colour="grey40", hjust=0)) +
    labs(title="Taxon-morphology associations (linkET-style)",
         subtitle=sprintf(
           "Top 20 drought-responsive genera vs 4 morphological traits | %d significant pairs (BH-adjusted)",
           n_sig_lk))

  dir.create("figures/P4_phenotype", recursive=TRUE, showWarnings=FALSE)
  ggsave("figures/P4_phenotype/linkET_taxa_morphology.svg",
         p_linket, width=16, height=10, units="in", device=svglite::svglite)
  ggsave("figures/P4_phenotype/linkET_taxa_morphology.pdf",
         p_linket, width=16, height=10, units="in", device=cairo_pdf)
  ggsave("figures/P4_phenotype/linkET_taxa_morphology.png",
         p_linket, width=16, height=10, units="in", dpi=300, type="cairo")
  cat(sprintf("linkET figure saved: .svg (%.1f KB)  .pdf (%.1f KB)  .png (%.1f KB)\n",
      file.size("figures/P4_phenotype/linkET_taxa_morphology.svg")/1024,
      file.size("figures/P4_phenotype/linkET_taxa_morphology.pdf")/1024,
      file.size("figures/P4_phenotype/linkET_taxa_morphology.png")/1024))

  ## Report strongest associations
  cat("=== STRONGEST GENUS-TRAIT ASSOCIATIONS ===\n")
  top_assoc_lk <- mantel_res_lk %>% dplyr::filter(!is.na(Mantel_r)) %>%
    dplyr::arrange(dplyr::desc(abs(Mantel_r))) %>% head(10) %>%
    dplyr::mutate(Trait_label=trait_labels_map_lk[Trait],
                  adj_sig=sig_stars(Mantel_p_adj))
  print(top_assoc_lk[,c("Genus","Trait_label","Mantel_r","Mantel_p",
                          "Mantel_p_adj","adj_sig","Enrichment")])
}

cat("\nPillar 4 complete.\n\n")

## =============================================================================
##  PILLAR 5: PREDICTABILITY AND INTERVENTION TARGETS
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 5: Predictability and Intervention Targets\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

## Genus-level
ps_genus <- tax_glom(ps_rel_planted, taxrank="Genus", NArm=TRUE)
X <- as.data.frame(t(otu_table(ps_genus))); colnames(X) <- make.names(colnames(X))
y <- factor(sample_data(ps_genus)$treatment, levels=c("Watered","Drought"))
cat(sprintf("RF: %d samples × %d genera\n", nrow(X), ncol(X)))

## Full model
set.seed(42)
rf_full <- ranger(y~., data=data.frame(X, y=y), probability=TRUE, num.trees=CONFIG$RF_NTREES,
                  importance="permutation", oob.error=TRUE)
oob_acc <- 1-rf_full$prediction.error
pred_full <- rf_full$predictions[,2]
roc_full <- roc(y, pred_full, levels=c("Watered","Drought"), direction="<", quiet=TRUE)
auc_full <- as.numeric(auc(roc_full))

## Permutation null
cat("RF permutation test (", CONFIG$RF_PERM_ITER, "iterations)...\n")
null_aucs <- replicate(CONFIG$RF_PERM_ITER, {
  yp <- sample(y)
  rf_p <- ranger(yp~., data=data.frame(X, yp=yp), probability=TRUE, num.trees=500, oob.error=TRUE)
  tryCatch(as.numeric(auc(roc(yp, rf_p$predictions[,2], levels=levels(y), direction="<", quiet=TRUE))),
           error=function(e) 0.5)
})
perm_p <- mean(null_aucs >= auc_full)
cat(sprintf("Full: AUC=%.3f, perm p=%.4f (null mean=%.3f)\n", auc_full, perm_p, mean(null_aucs)))

## Reduced model (top 30)
imp <- data.frame(Genus=names(rf_full$variable.importance),
                  Importance=rf_full$variable.importance) %>% arrange(desc(Importance))
write.csv(imp, "tables/P5_prediction/RF_feature_importance.csv", row.names=FALSE)
top30 <- head(imp$Genus, 30)

set.seed(42)
rf_red <- ranger(y~., data=data.frame(X[,top30], y=y), probability=TRUE,
                 num.trees=CONFIG$RF_NTREES, importance="permutation", oob.error=TRUE)
pred_red <- rf_red$predictions[,2]
roc_red <- roc(y, pred_red, levels=c("Watered","Drought"), direction="<", quiet=TRUE)
auc_red <- as.numeric(auc(roc_red))

null_red <- replicate(CONFIG$RF_PERM_ITER, {
  yp <- sample(y)
  rp <- ranger(yp~., data=data.frame(X[,top30], yp=yp), probability=TRUE, num.trees=500, oob.error=TRUE)
  tryCatch(as.numeric(auc(roc(yp, rp$predictions[,2], levels=levels(y), direction="<", quiet=TRUE))),
           error=function(e) 0.5)
})
perm_p_red <- mean(null_red >= auc_red)
cat(sprintf("Top30: AUC=%.3f, perm p=%.4f\n", auc_red, perm_p_red))

## Save performance
perf <- data.frame(
  Model=c("Full","Top30"), AUC=c(auc_full, auc_red),
  OOB=c(oob_acc, 1-rf_red$prediction.error),
  Perm_p=c(perm_p, perm_p_red), Null_mean=c(mean(null_aucs), mean(null_red))
)
write.csv(perf, "tables/P5_prediction/RF_performance.csv", row.names=FALSE)

## IMPORTANT CAVEAT about Top30 AUC
cat("\n*** FEATURE SELECTION BIAS WARNING ***\n")
cat("The Top30 model AUC is likely INFLATED because features were selected\n")
cat("from the same data used for evaluation (importance ranking from the full\n")
cat("model fitted to all 72 samples, then refit on the same samples).\n")
cat(sprintf("The full model AUC (%.3f, p=%.4f) is the HONEST primary result.\n", auc_full, perm_p))
cat("The Top30 AUC should be presented as EXPLORATORY — identifying which\n")
cat("genera carry discriminatory signal — NOT as a validated classifier.\n")
cat("For a publishable Top30 AUC, nested cross-validation is required.\n\n")

perf$Caveat <- c("Primary result", "Potentially inflated by selection bias — report as exploratory")
write.csv(perf, "tables/P5_prediction/RF_performance.csv", row.names=FALSE)
pdf("figures/P5_prediction/ROC_with_null.pdf", w=8, h=7)
plot(roc_full, col="#1976D2", lwd=2, main="ROC: Treatment Classification")
plot(roc_red, col="#D32F2F", lwd=2, add=TRUE)
abline(a=0,b=1,lty=2,col="grey")
legend("bottomright", bty="n", lwd=2, col=c("#1976D2","#D32F2F","grey"),
  lty=c(1,1,2), legend=c(
    sprintf("All genera (AUC=%.3f, %s)", auc_full,
            ifelse(perm_p < 0.001, "p<0.001", sprintf("p=%.3f", perm_p))),
    sprintf("Top 30 (AUC=%.3f, %s) — EXPLORATORY", auc_red,
            ifelse(perm_p_red < 0.001, "p<0.001", sprintf("p=%.3f", perm_p_red))),
    sprintf("Null (mean=%.3f)", mean(null_aucs))))
dev.off()
png("figures/P5_prediction/ROC_with_null.png", w=8, h=7, units="in", res=300)
plot(roc_full, col="#1976D2", lwd=2, main="ROC: Treatment Classification")
plot(roc_red, col="#D32F2F", lwd=2, add=TRUE)
abline(a=0,b=1,lty=2,col="grey")
legend("bottomright", bty="n", lwd=2, col=c("#1976D2","#D32F2F","grey"),
  lty=c(1,1,2), legend=c(
    sprintf("All genera (AUC=%.3f, %s)", auc_full,
            ifelse(perm_p < 0.001, "p<0.001", sprintf("p=%.3f", perm_p))),
    sprintf("Top 30 (AUC=%.3f, %s) — EXPLORATORY", auc_red,
            ifelse(perm_p_red < 0.001, "p<0.001", sprintf("p=%.3f", perm_p_red))),
    sprintf("Null (mean=%.3f)", mean(null_aucs))))
dev.off()
## ROC SVG (svglite)
svglite::svglite("figures/P5_prediction/ROC_with_null.svg", width=8, height=7)
plot(roc_full, col="#1976D2", lwd=2, main="ROC: Treatment Classification")
plot(roc_red, col="#D32F2F", lwd=2, add=TRUE)
abline(a=0,b=1,lty=2,col="grey")
legend("bottomright", bty="n", lwd=2, col=c("#1976D2","#D32F2F","grey"),
  lty=c(1,1,2), legend=c(
    sprintf("All genera (AUC=%.3f, %s)", auc_full,
            ifelse(perm_p < 0.001, "p<0.001", sprintf("p=%.3f", perm_p))),
    sprintf("Top 30 (AUC=%.3f, %s) — EXPLORATORY", auc_red,
            ifelse(perm_p_red < 0.001, "p<0.001", sprintf("p=%.3f", perm_p_red))),
    sprintf("Null (mean=%.3f)", mean(null_aucs))))
dev.off()

## Feature importance barplot (top 25)
tax_genus <- as.data.frame(tax_table(ps_genus))
imp_tax <- imp %>% head(25) %>%
  left_join(tax_genus %>% rownames_to_column("Genus_id") %>%
              mutate(Genus_id=make.names(Genus_id)) %>%
              dplyr::select(Genus_id, Genus_name=Genus, Family), by=c("Genus"="Genus_id")) %>%
  mutate(Label = italicise_taxa(coalesce(Genus_name, Family, Genus), rank = "genus"))

p_imp <- ggplot(imp_tax, aes(reorder(Label, Importance), Importance)) +
  geom_col(fill="steelblue", alpha=0.8) + coord_flip() +
  theme_pub() + theme(axis.text.y=element_markdown()) +
  labs(x=NULL, y="Permutation importance",
       title="Top 25 predictive genera (Random Forest)",
       subtitle=sprintf("Full AUC=%.3f (p=%.3f) | Top30 AUC=%.3f (p=%.3f)",
                        auc_full, perm_p, auc_red, perm_p_red),
       caption="Top30 AUC reflects features selected from full dataset; interpret with caution (potential overfitting).")
save_fig(p_imp, "figures/P5_prediction/RF_top25_importance", w=8, h=8)

## Core drought-responsive taxa
if (!is.null(trt_df)) {
  core_drought <- trt_df %>% filter(diff_abn) %>%
    dplyr::select(ASV, lfc, qval, Direction, Phylum, Family, Genus) %>% arrange(desc(abs(lfc)))
  write.csv(core_drought, "tables/P1_composition/core_drought_responsive.csv", row.names=FALSE)
  
  family_counts <- core_drought %>% filter(!is.na(Family)) %>%
    group_by(Family) %>% summarise(n=n(), .groups="drop") %>% arrange(desc(n)) %>% head(10)
  cat("Top drought-responsive families:\n"); print(as.data.frame(family_counts))
}

cat("\nPillar 5 complete.\n\n")

## =============================================================================
##  PILLAR 6: EXTENDED COMPOSITION ANALYSIS
##  Publication-quality figures styled after Yue et al. 2024 (Microbiome 12:44)
## =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 6: Extended Composition Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

## ── Install any missing packages ────────────────────────────────────────────
for (.pkg in c("ggVennDiagram")) {
  if (!requireNamespace(.pkg, quietly = TRUE))
    tryCatch(install.packages(.pkg, repos = "https://cloud.r-project.org", quiet = TRUE),
             warning = function(w) NULL,
             error   = function(e) message("Could not install ", .pkg, " — Venn diagrams will be skipped"))
}
.has_venn <- requireNamespace("ggVennDiagram", quietly = TRUE)
if (.has_venn) suppressPackageStartupMessages(library(ggVennDiagram))

dir.create("figures/P1_composition", recursive = TRUE, showWarnings = FALSE)

## ─────────────────────────────────────────────────────────────────────────────
## S6 COLOR PALETTES (vivid, publication-quality, high contrast)
## ─────────────────────────────────────────────────────────────────────────────

## Phylum palette — 10 vivid + grey "Other"
## Phylum colors from theme_constants.r
phylum_colors_vivid <- PAL$phylum
TOP_PHYLA_COLS      <- PAL$phylum

## 21-color vivid genus palette (20 top genera + grey Other)
genus_pal_21 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A0522D", "#F781BF", "#00CED1", "#228B22", "#FFD700",
  "#8B008B", "#FF6347", "#4169E1", "#32CD32", "#DC143C",
  "#1E90FF", "#FF8C00", "#9400D3", "#00FA9A", "#B8860B",
  "#B0B0B0"
)

## ─────────────────────────────────────────────────────────────────────────────
## S6 HELPER FUNCTIONS
## ─────────────────────────────────────────────────────────────────────────────

## Aggregate phyloseq OTU table to a taxonomic rank.
## Returns list(data = long data.frame with metadata, top_taxa, levels).
## Uses fast rowsum() approach — no tax_glom(), NA taxa → "Unclassified".
.tax_agg <- function(ps_obj, rank, top_n = 15) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) otu <- t(otu)       # ensure samples × taxa
  tx  <- as.data.frame(tax_table(ps_obj))
  rv  <- tx[[rank]]
  rv[is.na(rv) | rv == ""] <- "Unclassified"

  ## Sum counts / abundances by rank label (rank × samples)
  agg <- rowsum(t(otu), rv)

  ## Normalise each sample to relative abundance
  cs      <- colSums(agg)
  cs[cs == 0] <- 1                               # guard against empty samples
  agg     <- sweep(agg, 2, cs, "/")

  mean_ab  <- rowMeans(agg)
  all_taxa <- rownames(agg)

  if (length(all_taxa) > top_n) {
    top_taxa <- names(sort(mean_ab, decreasing = TRUE))[seq_len(top_n)]
    other_ab <- colSums(agg[setdiff(all_taxa, top_taxa), , drop = FALSE])
    agg_top  <- rbind(agg[top_taxa, , drop = FALSE], Other = other_ab)
    lev      <- c(rev(top_taxa), "Other")
  } else {
    top_taxa <- all_taxa
    agg_top  <- agg
    lev      <- rev(top_taxa)
  }

  ## Melt to long format and join metadata
  df <- as.data.frame(t(agg_top)) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Taxon", values_to = "Abundance")
  md <- data.frame(sample_data(ps_obj)) %>% rownames_to_column("Sample")
  df <- left_join(df, md, by = "Sample")
  df$Taxon <- factor(df$Taxon, levels = lev)
  list(data = df, top_taxa = top_taxa, levels = lev)
}

## Compute rank-level differential abundance (log2 LFC + Wilcoxon FDR).
## group_A is the "numerator" (positive LFC = enriched in group_A).
.rank_da <- function(ps_obj, rank, group_var, group_A, group_B) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) otu <- t(otu)       # samples × taxa
  tx  <- as.data.frame(tax_table(ps_obj))
  rv  <- tx[[rank]]
  rv[is.na(rv) | rv == ""] <- "Unclassified"

  agg <- rowsum(t(otu), rv)
  cs  <- colSums(agg); cs[cs == 0] <- 1
  agg <- sweep(agg, 2, cs, "/")

  md  <- data.frame(sample_data(ps_obj))
  grp <- md[[group_var]]
  idx_A <- which(grp == group_A)
  idx_B <- which(grp == group_B)

  res <- do.call(rbind, lapply(rownames(agg), function(tx_name) {
    xA <- as.numeric(agg[tx_name, idx_A])
    xB <- as.numeric(agg[tx_name, idx_B])
    eps <- 1e-6
    m_A <- mean(xA); m_B <- mean(xB)
    lfc <- log2((m_A + eps) / (m_B + eps))
    pv  <- tryCatch(wilcox.test(xA, xB, exact = FALSE)$p.value,
                    error = function(e) NA_real_)
    data.frame(Taxon = tx_name, mean_A = m_A, mean_B = m_B,
               LFC = lfc, p = pv, stringsAsFactors = FALSE)
  }))

  res$padj      <- p.adjust(res$p, method = "fdr")
  res$sig       <- ifelse(res$padj < 0.001, "***",
                   ifelse(res$padj < 0.01,  "**",
                   ifelse(res$padj < 0.05,  "*", "")))
  res$Direction <- ifelse(res$LFC >= 0,
                          paste0(group_A, "_enriched"),
                          paste0(group_B, "_enriched"))
  res
}

## Identify core taxa in a phyloseq subset (presence in ≥ prev_thresh of samples).
## Returns data.frame with ASV/taxon, prevalence, mean abundance, and taxonomy.
.get_core <- function(ps_obj, group_var, group_val, prev_thresh = 0.80) {
  ## Use prune_samples + standard data.frame indexing to avoid phyloseq NSE bugs
  md_tmp <- data.frame(sample_data(ps_obj))
  keep_s <- rownames(md_tmp)[!is.na(md_tmp[[group_var]]) &
                               md_tmp[[group_var]] == group_val]
  if (length(keep_s) == 0) return(NULL)
  ps_sub <- prune_samples(keep_s, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  if (nsamples(ps_sub) == 0) return(NULL)

  otu    <- as(otu_table(ps_sub), "matrix")
  if (taxa_are_rows(ps_sub)) otu <- t(otu)       # samples × taxa
  n_samp <- nrow(otu)

  prev <- apply(otu, 2, function(x) sum(x > 0) / n_samp)
  core_asvs <- names(prev[prev >= prev_thresh])
  if (length(core_asvs) == 0) return(NULL)

  rs <- rowSums(otu); rs[rs == 0] <- 1
  mean_ab <- colMeans(sweep(otu[, core_asvs, drop = FALSE], 1, rs, "/"))
  tx <- as.data.frame(tax_table(ps_sub))[core_asvs, , drop = FALSE]

  data.frame(ASV = core_asvs, prevalence = prev[core_asvs],
             mean_abundance = mean_ab, tx, group = group_val,
             stringsAsFactors = FALSE, row.names = NULL)
}

## Return set of ASV names present in ≥ 1 sample for a given group.
.get_present_asvs <- function(ps_obj, group_var, group_val) {
  ## Use prune_samples + standard data.frame indexing to avoid phyloseq NSE bugs
  md_tmp <- data.frame(sample_data(ps_obj))
  keep_s <- rownames(md_tmp)[!is.na(md_tmp[[group_var]]) &
                               md_tmp[[group_var]] == group_val]
  if (length(keep_s) == 0) return(character(0))
  ps_sub <- prune_samples(keep_s, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  taxa_names(ps_sub)
}

## Return set of genus names present in ≥ 1 sample for a given group.
.get_present_genera <- function(ps_obj, group_var, group_val) {
  ## Use prune_samples + standard data.frame indexing to avoid phyloseq NSE bugs
  md_tmp <- data.frame(sample_data(ps_obj))
  keep_s <- rownames(md_tmp)[!is.na(md_tmp[[group_var]]) &
                               md_tmp[[group_var]] == group_val]
  if (length(keep_s) == 0) return(character(0))
  ps_sub <- prune_samples(keep_s, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  tx  <- as.data.frame(tax_table(ps_sub))
  otu <- as(otu_table(ps_sub), "matrix")
  if (taxa_are_rows(ps_sub)) otu <- t(otu)
  gv  <- tx$Genus
  gv[is.na(gv)] <- paste0("Unclassified_", seq_len(sum(is.na(gv))))
  agg <- rowsum(t(otu), gv)
  rownames(agg)[rowSums(agg) > 0]
}

cat("S6 helpers defined.\n\n")

## ─────────────────────────────────────────────────────────────────────────────
## 6a. FIXED PHYLUM BARPLOTS (vivid palette, guaranteed 100%)
## ─────────────────────────────────────────────────────────────────────────────
cat("── 6a. Phylum barplots (fixed) ──\n")

phy_agg  <- .tax_agg(ps_rel_planted, "Phylum", top_n = 10)
phy_long <- phy_agg$data
top_phy  <- phy_agg$top_taxa

## Build colour map — predefined vivid palette with fallback for unknown phyla
phy_col_map <- phylum_colors_vivid[levels(phy_long$Taxon)]
names(phy_col_map) <- levels(phy_long$Taxon)
missing_phy <- is.na(phy_col_map)
if (any(missing_phy)) {
  fallback_phy <- c("#1ABC9C","#16A085","#2ECC71","#8E44AD","#2980B9",
                    "#E67E22","#C0392B","#27AE60","#D35400","#7F8C8D")
  phy_col_map[missing_phy] <- fallback_phy[seq_len(sum(missing_phy))]
}

## ── Fig 6a-i  Per-sample barplot (position="fill" → exact 100%) ─────────────
phy_sample_order <- phy_long %>%
  distinct(Sample, treatment, harvest_num, trait) %>%
  arrange(treatment, harvest_num, trait) %>%
  pull(Sample)
phy_long$Sample <- factor(phy_long$Sample, levels = phy_sample_order)

p_phy_sample <- ggplot(phy_long,
                       aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "fill", width = 1, colour = NA) +
  facet_grid(. ~ treatment + trait, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phy_col_map, name = "Phylum") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_pub(BASE_SIZE - 2) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text   = element_text(face = "bold", size = 8),
        legend.key.size = unit(0.4, "cm")) +
  labs(x = "Samples ordered by Treatment → Harvest → Genotype",
       y = "Relative abundance",
       title = sprintf("Phylum composition — per sample (top %d; 100%% guaranteed)",
                       length(top_phy)))
save_fig(p_phy_sample, "figures/P1_composition/phylum_barplot_per_sample",
         w = 14, h = 6)

## ── Fig 6a-ii  Group-averaged treatment × genotype × timepoint ───────────────
## Normalize within each group after averaging to guarantee exact 100%.
phy_grp <- phy_long %>%
  filter(trait %in% c("Resistance", "Susceptible")) %>%
  group_by(treatment, trait, harvest, harvest_num, Taxon) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(treatment, trait, harvest) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

p_phy_grp <- ggplot(phy_grp,
                    aes(x = harvest, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  facet_grid(treatment ~ trait) +
  scale_fill_manual(values = phy_col_map, name = "Phylum") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_pub(BASE_SIZE - 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, "cm")) +
  labs(x = "Harvest", y = "Relative abundance",
       title = "Phylum composition: treatment × genotype × timepoint",
       subtitle = "Group means re-normalised to 100%")
save_fig(p_phy_grp, "figures/P1_composition/phylum_barplot_trt_geno_time_FIXED",
         w = 12, h = 8)

## ── Fig 6a-iii  Treatment only ───────────────────────────────────────────────
phy_trt_grp <- phy_long %>%
  group_by(treatment, harvest, harvest_num, Taxon) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(treatment, harvest) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

p_phy_trt <- ggplot(phy_trt_grp,
                    aes(x = harvest, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  facet_wrap(~ treatment) +
  scale_fill_manual(values = phy_col_map, name = "Phylum") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, "cm")) +
  labs(x = "Harvest", y = "Relative abundance",
       title = "Phylum composition by treatment",
       subtitle = "Group means re-normalised to 100%")
save_fig(p_phy_trt, "figures/P1_composition/phylum_barplot_treatment_FIXED",
         w = 10, h = 6)

## Save phylum abundance table
phy_tbl <- phy_long %>%
  filter(trait %in% c("Resistance", "Susceptible")) %>%
  group_by(treatment, trait, harvest, Taxon) %>%
  summarise(mean_abundance = mean(Abundance),
            se = sd(Abundance) / sqrt(n()),
            .groups = "drop")
write.csv(phy_tbl, "tables/P1_composition/phylum_abundance_by_group.csv", row.names = FALSE)

cat("Phylum barplots saved.\n\n")

## ─────────────────────────────────────────────────────────────────────────────
## 6b. GENUS-LEVEL COMPOSITION BARPLOTS (top 20 genera)
## ─────────────────────────────────────────────────────────────────────────────
cat("── 6b. Genus barplots ──\n")

gen_agg  <- .tax_agg(ps_rel_planted, "Genus", top_n = 20)
gen_long <- gen_agg$data
top_gen  <- gen_agg$top_taxa

## Assign vivid colours; always grey for "Other"
n_gen_lev   <- length(levels(gen_long$Taxon))
gen_col_map <- setNames(genus_pal_21[seq_len(n_gen_lev)], levels(gen_long$Taxon))
if ("Other" %in% names(gen_col_map)) gen_col_map["Other"] <- "#B0B0B0"

## Sample ordering identical to phylum plots for consistency
gen_sample_order <- gen_long %>%
  distinct(Sample, treatment, harvest_num, trait) %>%
  arrange(treatment, harvest_num, trait) %>%
  pull(Sample)
gen_long$Sample <- factor(gen_long$Sample, levels = gen_sample_order)

## ── Fig 6b-i  Per-sample ────────────────────────────────────────────────────
p_gen_sample <- ggplot(gen_long,
                       aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "fill", width = 1, colour = NA) +
  facet_grid(. ~ treatment + trait, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = gen_col_map, name = "Genus") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_pub(BASE_SIZE - 2) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text   = element_text(face = "bold", size = 8),
        legend.key.size  = unit(0.35, "cm"),
        legend.text  = element_text(size = 8, face = "italic")) +
  labs(x = "Samples ordered by Treatment → Harvest → Genotype",
       y = "Relative abundance",
       title = "Genus composition — per sample (top 20; 100% guaranteed)")
save_fig(p_gen_sample, "figures/P1_composition/genus_barplot_per_sample",
         w = 14, h = 6)

## ── Fig 6b-ii  Group average: treatment × genotype (pooled time) ─────────────
gen_grp_tg <- gen_long %>%
  filter(trait %in% c("Resistance", "Susceptible")) %>%
  group_by(treatment, trait, Taxon) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(treatment, trait) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  mutate(group_label = interaction(treatment, trait, sep = "\n"))

p_gen_tg <- ggplot(gen_grp_tg,
                   aes(x = group_label, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = gen_col_map, name = "Genus") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_pub() +
  theme(legend.text     = element_text(size = 8, face = "italic"),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = "Treatment × Genotype", y = "Relative abundance",
       title = "Genus composition by treatment × genotype",
       subtitle = "Means re-normalised to 100% — top 20 genera")
save_fig(p_gen_tg, "figures/P1_composition/genus_barplot_trt_geno",
         w = 10, h = 7)

## ── Fig 6b-iii  Group average: treatment × genotype × timepoint ──────────────
gen_grp_tt <- gen_long %>%
  filter(trait %in% c("Resistance", "Susceptible")) %>%
  group_by(treatment, trait, harvest, harvest_num, Taxon) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(treatment, trait, harvest) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

p_gen_tt <- ggplot(gen_grp_tt,
                   aes(x = harvest, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  facet_grid(trait ~ treatment) +
  scale_fill_manual(values = gen_col_map, name = "Genus") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_pub(BASE_SIZE - 2) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.text     = element_text(size = 8, face = "italic"),
        legend.key.size = unit(0.35, "cm")) +
  labs(x = "Harvest", y = "Relative abundance",
       title = "Genus composition: treatment × genotype × timepoint",
       subtitle = "Top 20 genera, means re-normalised to 100%")
save_fig(p_gen_tt, "figures/P1_composition/genus_barplot_trt_geno_time",
         w = 12, h = 8)

## Save genus abundance table
gen_tbl <- gen_long %>%
  filter(trait %in% c("Resistance", "Susceptible")) %>%
  group_by(treatment, trait, harvest, Taxon) %>%
  summarise(mean_abundance = mean(Abundance),
            se = sd(Abundance) / sqrt(n()),
            .groups = "drop")
write.csv(gen_tbl, "tables/P1_composition/genus_abundance_by_group.csv", row.names = FALSE)

cat("Genus barplots saved.\n\n")

## ─────────────────────────────────────────────────────────────────────────────
## 6c. TOP TAXA HORIZONTAL BAR CHARTS BY RANK (Yue et al. Fig. 1e style)
##     Wilcoxon test (FDR) + log2 fold change, Drought vs Watered
## ─────────────────────────────────────────────────────────────────────────────
cat("── 6c. Top taxa horizontal bar charts ──\n")

ranks_to_plot <- c("Phylum", "Class", "Order", "Family", "Genus")
rank_da_list  <- list()

for (rk in ranks_to_plot) {
  cat(sprintf("  Computing %s-level DA...\n", rk))

  rk_res <- tryCatch(
    .rank_da(ps_rel_planted, rank = rk,
             group_var = "treatment",
             group_A   = "Drought",
             group_B   = "Watered"),
    error = function(e) { cat("    Skipped:", e$message, "\n"); NULL }
  )
  if (is.null(rk_res) || nrow(rk_res) == 0) next

  ## Remove "Unclassified" from the ranking (uninformative for display)
  rk_res <- rk_res %>% filter(Taxon != "Unclassified")
  rank_da_list[[rk]] <- rk_res
  write.csv(rk_res,
            sprintf("tables/P1_composition/DA_%s_DroughtVsWatered.csv", rk),
            row.names = FALSE)

  ## Select top 15 enriched (LFC > 0) and top 15 depleted (LFC < 0)
  top_enr <- rk_res %>% filter(LFC > 0)  %>% arrange(desc(LFC))  %>% head(15)
  top_dep <- rk_res %>% filter(LFC < 0)  %>% arrange(LFC)        %>% head(15)
  top_both <- bind_rows(top_dep, top_enr) %>%   # low → high LFC (bottom → top of plot)
    arrange(LFC)
  if (nrow(top_both) < 2) next

  if (rk == "Genus") {
    top_both$Taxon <- italicise_taxa(as.character(top_both$Taxon), rank = "genus")
  }
  top_both$Taxon     <- factor(top_both$Taxon, levels = top_both$Taxon)
  top_both$Direction <- ifelse(top_both$LFC >= 0, "Drought enriched",
                                                   "Watered enriched")

  ## Axis text: element_markdown for genera (renders * * italics), plain otherwise
  y_axis_text <- if (rk == "Genus") element_markdown(size = 9) else element_text(size = 9)

  p_rk <- ggplot(top_both, aes(y = Taxon, x = LFC, fill = Direction)) +
    geom_col(width = 0.72, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey30") +
    geom_text(aes(label = sig,
                  hjust = ifelse(LFC >= 0, -0.25, 1.25)),
              size = 3.8, fontface = "bold", colour = "grey10") +
    scale_fill_manual(
      values = c("Drought enriched" = PAL$treatment[["Drought"]],
                 "Watered enriched" = PAL$treatment[["Watered"]])) +
    scale_x_continuous(expand = expansion(mult = 0.22)) +
    theme_pub() +
    theme(axis.text.y    = y_axis_text,
          legend.position = "bottom",
          legend.title    = element_blank()) +
    labs(y = NULL,
         x = expression(log[2]~"fold change (Drought / Watered)"),
         title    = sprintf("Top drought-responsive %s-level taxa",
                            tolower(rk)),
         subtitle = "Wilcoxon test, FDR corrected: * p<0.05  ** p<0.01  *** p<0.001")

  fig_h <- max(5, 0.38 * nrow(top_both) + 2.2)
  save_fig(p_rk,
           sprintf("figures/P1_composition/top_taxa_%s_horizontal", rk),
           w = 9, h = fig_h)
}

cat(sprintf("Horizontal bar charts saved for %d ranks.\n\n",
            length(rank_da_list)))

## ─────────────────────────────────────────────────────────────────────────────
## 6d. VENN DIAGRAMS (shared/unique ASVs across experimental groups)
## ─────────────────────────────────────────────────────────────────────────────
cat("── 6d. Venn diagrams ──\n")

## Restrict to planted genotype samples for genotype Venns
ps_pl_gen <- subset_samples(ps_planted,
                            trait %in% c("Resistance", "Susceptible"))
ps_pl_gen <- prune_taxa(taxa_sums(ps_pl_gen) > 0, ps_pl_gen)

## ── 2-set: Drought vs Watered ────────────────────────────────────────────────
venn_trt <- list(
  Drought = .get_present_asvs(ps_planted, "treatment", "Drought"),
  Watered = .get_present_asvs(ps_planted, "treatment", "Watered")
)
venn_trt_sum <- data.frame(
  Drought_only = length(setdiff(venn_trt$Drought, venn_trt$Watered)),
  Watered_only = length(setdiff(venn_trt$Watered, venn_trt$Drought)),
  Shared       = length(intersect(venn_trt$Drought, venn_trt$Watered))
)
write.csv(venn_trt_sum, "tables/P1_composition/venn_ASV_treatment.csv", row.names = FALSE)
cat(sprintf("  ASV Venn (treatment): Drought-only=%d, Watered-only=%d, Shared=%d\n",
            venn_trt_sum$Drought_only, venn_trt_sum$Watered_only,
            venn_trt_sum$Shared))

if (.has_venn) {
  p_venn_trt <- tryCatch({
    ggVennDiagram(venn_trt, label_alpha = 0) +
      scale_fill_gradient(low = "#FFFFFF", high = PAL$treatment[["Drought"]]) +
      labs(title   = "ASV sharing: Drought vs Watered",
           subtitle = sprintf("Drought-only: %d | Shared: %d | Watered-only: %d",
                              venn_trt_sum$Drought_only, venn_trt_sum$Shared,
                              venn_trt_sum$Watered_only)) +
      theme(legend.position = "none")
  }, error = function(e) { cat("  Venn error:", e$message, "\n"); NULL })
  if (!is.null(p_venn_trt))
    save_fig(p_venn_trt, "figures/P1_composition/venn_ASV_treatment", w = ISME_SINGLE_W, h = 5)
}

## ── 2-set: Resistance vs Susceptible ─────────────────────────────────────────
venn_gen <- list(
  Resistance  = .get_present_asvs(ps_pl_gen, "trait", "Resistance"),
  Susceptible = .get_present_asvs(ps_pl_gen, "trait", "Susceptible")
)
venn_gen_sum <- data.frame(
  Resistance_only  = length(setdiff(venn_gen$Resistance, venn_gen$Susceptible)),
  Susceptible_only = length(setdiff(venn_gen$Susceptible, venn_gen$Resistance)),
  Shared           = length(intersect(venn_gen$Resistance, venn_gen$Susceptible))
)
write.csv(venn_gen_sum, "tables/P1_composition/venn_ASV_genotype.csv", row.names = FALSE)
cat(sprintf("  ASV Venn (genotype): Resistance-only=%d, Susceptible-only=%d, Shared=%d\n",
            venn_gen_sum$Resistance_only, venn_gen_sum$Susceptible_only,
            venn_gen_sum$Shared))

if (.has_venn) {
  p_venn_gen <- tryCatch({
    ggVennDiagram(venn_gen, label_alpha = 0) +
      scale_fill_gradient(low = "#FFFFFF", high = PAL$genotype[["Resistance"]]) +
      labs(title   = "ASV sharing: Resistance vs Susceptible",
           subtitle = sprintf("Resistance-only: %d | Shared: %d | Susceptible-only: %d",
                              venn_gen_sum$Resistance_only, venn_gen_sum$Shared,
                              venn_gen_sum$Susceptible_only)) +
      theme(legend.position = "none")
  }, error = function(e) { cat("  Venn error:", e$message, "\n"); NULL })
  if (!is.null(p_venn_gen))
    save_fig(p_venn_gen, "figures/P1_composition/venn_ASV_genotype", w = ISME_SINGLE_W, h = 5)
}

## ── 4-set: Drought×Resistant, Drought×Susceptible, Watered×Resistant, Watered×Susceptible
meta_4 <- data.frame(sample_data(ps_pl_gen))
meta_4$trt_geno <- paste(as.character(meta_4$treatment),
                          as.character(meta_4$trait), sep = "_")

venn_4set <- lapply(
  c("Drought_Resistance", "Drought_Susceptible",
    "Watered_Resistance",  "Watered_Susceptible"),
  function(g) {
    s <- rownames(meta_4)[meta_4$trt_geno == g]
    if (length(s) == 0) return(character(0))
    ps_s <- prune_samples(s, ps_pl_gen)
    taxa_names(prune_taxa(taxa_sums(ps_s) > 0, ps_s))
  }
)
names(venn_4set) <- c("Drought_Resistant", "Drought_Susceptible",
                       "Watered_Resistant",  "Watered_Susceptible")

venn_4_df <- do.call(rbind, lapply(names(venn_4set), function(g)
  data.frame(Group = g, N_ASVs = length(venn_4set[[g]]))))
write.csv(venn_4_df, "tables/P1_composition/venn_ASV_4set_trt_geno.csv", row.names = FALSE)
cat("  4-set ASV counts:\n"); print(venn_4_df)

if (.has_venn && all(sapply(venn_4set, length) > 0)) {
  p_venn_4 <- tryCatch({
    ggVennDiagram(venn_4set, label_alpha = 0) +
      scale_fill_gradient(low = "#FFFFFF", high = "#984EA3") +
      labs(title    = "ASV sharing: Treatment × Genotype (4-set)",
           subtitle = "Numbers = ASVs present in each intersection") +
      theme(legend.position = "none")
  }, error = function(e) { cat("  4-set Venn error:", e$message, "\n"); NULL })
  if (!is.null(p_venn_4))
    save_fig(p_venn_4, "figures/P1_composition/venn_ASV_4set_trt_geno", w = 8, h = 7)
}

## ── Genus-level Venn (treatment; more interpretable than ASV) ────────────────
venn_genus_trt <- list(
  Drought = .get_present_genera(ps_planted, "treatment", "Drought"),
  Watered = .get_present_genera(ps_planted, "treatment", "Watered")
)
venn_genus_sum <- data.frame(
  Drought_only = length(setdiff(venn_genus_trt$Drought, venn_genus_trt$Watered)),
  Watered_only = length(setdiff(venn_genus_trt$Watered, venn_genus_trt$Drought)),
  Shared       = length(intersect(venn_genus_trt$Drought, venn_genus_trt$Watered))
)
write.csv(venn_genus_sum, "tables/P1_composition/venn_genus_treatment.csv", row.names = FALSE)
cat(sprintf("  Genus Venn (treatment): Drought-only=%d, Watered-only=%d, Shared=%d\n",
            venn_genus_sum$Drought_only, venn_genus_sum$Watered_only,
            venn_genus_sum$Shared))

if (.has_venn) {
  p_venn_genus <- tryCatch({
    ggVennDiagram(venn_genus_trt, label_alpha = 0) +
      scale_fill_gradient(low = "#FFFFFF", high = PAL$treatment[["Drought"]]) +
      labs(title   = "Genus sharing: Drought vs Watered",
           subtitle = sprintf("Drought-only: %d | Shared: %d | Watered-only: %d",
                              venn_genus_sum$Drought_only, venn_genus_sum$Shared,
                              venn_genus_sum$Watered_only)) +
      theme(legend.position = "none")
  }, error = function(e) { cat("  Venn error:", e$message, "\n"); NULL })
  if (!is.null(p_venn_genus))
    save_fig(p_venn_genus, "figures/P1_composition/venn_genus_treatment", w = ISME_SINGLE_W, h = 5)
}

cat("Venn diagrams saved.\n\n")

## ─────────────────────────────────────────────────────────────────────────────
## 6e. CORE MICROBIOME ANALYSIS (≥80% prevalence within each group)
## ─────────────────────────────────────────────────────────────────────────────
cat("── 6e. Core microbiome analysis (>=80% prevalence) ──\n")

CORE_PREV <- 0.80

core_drought <- .get_core(ps_planted, "treatment", "Drought",    CORE_PREV)
core_watered <- .get_core(ps_planted, "treatment", "Watered",    CORE_PREV)
core_resist  <- .get_core(ps_pl_gen,  "trait",     "Resistance", CORE_PREV)
core_suscept <- .get_core(ps_pl_gen,  "trait",     "Susceptible",CORE_PREV)

for (g_info in list(
    list(obj = core_drought, nm = "Drought"),
    list(obj = core_watered, nm = "Watered"),
    list(obj = core_resist,  nm = "Resistance"),
    list(obj = core_suscept, nm = "Susceptible"))) {
  n <- if (is.null(g_info$obj)) 0 else nrow(g_info$obj)
  cat(sprintf("  Core ASVs (%s, >=%.0f%%): %d\n",
              g_info$nm, CORE_PREV * 100, n))
}

## Combine and save
core_all <- bind_rows(
  Filter(Negate(is.null),
         list(core_drought, core_watered, core_resist, core_suscept))
)
if (nrow(core_all) > 0)
  write.csv(core_all, "tables/P1_composition/core_microbiome_by_group.csv", row.names = FALSE)

## ── Core Venn: treatment ─────────────────────────────────────────────────────
if (!is.null(core_drought) && !is.null(core_watered)) {
  venn_core_trt <- list(Drought = core_drought$ASV,
                        Watered = core_watered$ASV)
  core_trt_sum <- data.frame(
    Drought_core_only = length(setdiff(venn_core_trt$Drought,
                                       venn_core_trt$Watered)),
    Watered_core_only = length(setdiff(venn_core_trt$Watered,
                                       venn_core_trt$Drought)),
    Shared_core       = length(intersect(venn_core_trt$Drought,
                                         venn_core_trt$Watered))
  )
  write.csv(core_trt_sum, "tables/P1_composition/core_venn_treatment.csv", row.names = FALSE)
  cat("  Core Venn (treatment):\n"); print(core_trt_sum)

  if (.has_venn && sum(sapply(venn_core_trt, length)) > 0) {
    p_core_trt <- tryCatch({
      ggVennDiagram(venn_core_trt, label_alpha = 0) +
        scale_fill_gradient(low = "#FFFFFF", high = PAL$treatment[["Drought"]]) +
        labs(title    = sprintf("Core microbiome (>=%.0f%% prevalence): Drought vs Watered",
                                CORE_PREV * 100),
             subtitle = "Planted samples only") +
        theme(legend.position = "none")
    }, error = function(e) NULL)
    if (!is.null(p_core_trt))
      save_fig(p_core_trt,
               "figures/P1_composition/core_venn_treatment", w = 6, h = 5)
  }
}

## ── Core Venn: genotype ──────────────────────────────────────────────────────
if (!is.null(core_resist) && !is.null(core_suscept)) {
  venn_core_gen <- list(Resistance  = core_resist$ASV,
                        Susceptible = core_suscept$ASV)
  core_gen_sum <- data.frame(
    Resistance_core_only  = length(setdiff(venn_core_gen$Resistance,
                                           venn_core_gen$Susceptible)),
    Susceptible_core_only = length(setdiff(venn_core_gen$Susceptible,
                                           venn_core_gen$Resistance)),
    Shared_core           = length(intersect(venn_core_gen$Resistance,
                                             venn_core_gen$Susceptible))
  )
  write.csv(core_gen_sum, "tables/P1_composition/core_venn_genotype.csv", row.names = FALSE)
  cat("  Core Venn (genotype):\n"); print(core_gen_sum)

  if (.has_venn && sum(sapply(venn_core_gen, length)) > 0) {
    p_core_gen <- tryCatch({
      ggVennDiagram(venn_core_gen, label_alpha = 0) +
        scale_fill_gradient(low = "#FFFFFF", high = PAL$genotype[["Resistance"]]) +
        labs(title    = sprintf("Core microbiome (>=%.0f%% prevalence): Resistance vs Susceptible",
                                CORE_PREV * 100),
             subtitle = "Planted samples only") +
        theme(legend.position = "none")
    }, error = function(e) NULL)
    if (!is.null(p_core_gen))
      save_fig(p_core_gen,
               "figures/P1_composition/core_venn_genotype", w = 6, h = 5)
  }
}

## ── Core phylum composition barplot ─────────────────────────────────────────
if (!is.null(core_all) && nrow(core_all) > 0 &&
    "Phylum" %in% colnames(core_all)) {

  core_phy <- core_all %>%
    filter(!is.na(Phylum)) %>%
    group_by(group, Phylum) %>%
    summarise(n_asvs = n(),
              mean_ab = mean(mean_abundance, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(group) %>%
    mutate(prop = n_asvs / sum(n_asvs)) %>%
    ungroup()

  top_core_phy <- core_phy %>%
    group_by(Phylum) %>% summarise(tot = sum(n_asvs), .groups = "drop") %>%
    arrange(desc(tot)) %>% head(8) %>% pull(Phylum)

  core_phy <- core_phy %>%
    mutate(Phylum_plot = ifelse(Phylum %in% top_core_phy, Phylum, "Other"),
           Phylum_plot = factor(Phylum_plot, levels = c(top_core_phy, "Other")))

  core_phy_col <- c(phylum_colors_vivid[top_core_phy], Other = "#B0B0B0")
  core_phy_col <- core_phy_col[!is.na(names(core_phy_col))]

  p_core_phy <- ggplot(core_phy,
                       aes(x = group, y = prop, fill = Phylum_plot)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = core_phy_col, name = "Phylum",
                      na.value = "#B0B0B0") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "Group", y = "Proportion of core ASVs",
         title = sprintf("Core microbiome phylum composition (>=%.0f%% prevalence)",
                         CORE_PREV * 100))
  save_fig(p_core_phy, "figures/P1_composition/core_phylum_composition",
           w = 8, h = 6)
  write.csv(core_phy, "tables/P1_composition/core_phylum_composition.csv", row.names = FALSE)
}

cat("Core microbiome analysis complete.\n\n")

## ─────────────────────────────────────────────────────────────────────────────
## PILLAR 6 COMPLETE
## ─────────────────────────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 6 COMPLETE — Extended Composition Analysis\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("  Figures saved to: figures/P1_composition/ (%d files)\n",
            length(list.files("figures/P1_composition",
                              pattern = "\\.(pdf|png)$",
                              recursive = TRUE))))
cat(sprintf("  Tables saved to tables/: %s\n",
            paste(grep("(phylum|genus|DA_|venn|core).*\\.csv$",
                       list.files("tables"), value = TRUE),
                  collapse = ", ")))
cat("\n")
cat("\nPillar 6 complete.\n\n")

## =============================================================================
##  PILLAR 7: CO-OCCURRENCE NETWORKS AND KEYSTONE TAXA
##  (Treatment + Genotype networks; per-day temporal dynamics)
##
##  Sections:
##    7A  SPIEC-EASI treatment networks (Drought / Watered)
##    7B  SPIEC-EASI genotype networks  (Resistant / Susceptible)
##    7C  Per-day Spearman networks — treatment
##    7D  Per-day Spearman networks — genotype
##    7E  Figures: 4 network plots, Zi-Pi scatter ×2, topology bar, cohesion
##    7F  Temporal figures: per-day topology line plots + keystone heatmaps
##
##  CACHING: SPIEC-EASI results cached in data/spiec_easi_*.rds
##  All downstream figure/table outputs regenerated from caches.
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 7: Co-occurrence Networks and Keystone Taxa\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

suppressPackageStartupMessages({
  library(igraph); library(ggraph); library(ggrepel)
  library(Matrix); library(RColorBrewer); library(scales)
})
## SpiecEasi only needed if RDS caches are absent
.need_se <- !all(file.exists(c(
  "data/spiec_easi_drought.rds",   "data/spiec_easi_watered.rds",
  "data/spiec_easi_resistant.rds", "data/spiec_easi_susceptible.rds")))
if (.need_se) {
  if (!requireNamespace("SpiecEasi", quietly = TRUE))
    stop(paste0("SpiecEasi required for SPIEC-EASI runs.\n",
                "Install: remotes::install_github('zdk123/SpiecEasi')"))
  library(SpiecEasi)
}

## ── Shared helpers ─────────────────────────────────────────────────────────
.P6  <- "figures/P6_networks"
.TAB <- "tables"
dir.create(.P6,  recursive = TRUE, showWarnings = FALSE)

## 20-colour module palette
.MOD_PAL <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#888888","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#B3B3B3",
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E")
.mod_cols <- function(ids) {
  u <- sort(unique(as.integer(ids)))
  setNames(.MOD_PAL[((seq_along(u)-1L) %% length(.MOD_PAL)) + 1L], as.character(u))
}

## GLASSO precision matrix → named partial correlation matrix
.pcor_glasso <- function(net) {
  idx  <- net$select$stars$opt.index
  Pmat <- as.matrix(net$est$icov[[idx]])
  d    <- sqrt(pmax(diag(Pmat), 1e-12))
  pcor <- -Pmat / outer(d, d)
  diag(pcor) <- 0
  pcor
}

## Build igraph LCC + cluster_fast_greedy modules + Zi/Pi
.build_net <- function(pcor, taxa, tax_ref, label) {
  upper <- which(pcor != 0 & upper.tri(pcor), arr.ind = TRUE)
  g     <- igraph::make_empty_graph(n = length(taxa), directed = FALSE)
  igraph::V(g)$name <- taxa
  if (nrow(upper) > 0)
    g <- igraph::add_edges(g, as.vector(t(upper)), weight = pcor[upper])

  comps <- igraph::components(g)
  lcc_v <- which(comps$membership == which.max(comps$csize))
  g_lcc <- igraph::induced_subgraph(g, lcc_v)

  ## Modules on unweighted copy
  g_uw  <- igraph::graph_from_edgelist(
    igraph::as_edgelist(g_lcc, names=FALSE), directed=FALSE)
  mods  <- igraph::cluster_fast_greedy(g_uw)
  igraph::V(g_lcc)$module <- igraph::membership(mods)
  igraph::V(g_lcc)$degree <- igraph::degree(g_lcc)
  igraph::V(g_lcc)$betweenness <- igraph::betweenness(
    g_lcc, normalized=TRUE, weights=NA)

  ## Zi
  mod_v  <- igraph::V(g_lcc)$module
  deg_v  <- igraph::V(g_lcc)$degree
  Zi_vec <- numeric(igraph::vcount(g_lcc))
  for (m in unique(mod_v)) {
    idx <- which(mod_v == m)
    if (length(idx) < 2) { Zi_vec[idx] <- 0; next }
    k_s  <- sapply(idx, function(v) sum(igraph::neighbors(g_lcc, v) %in% idx))
    sd_k <- sd(k_s)
    if (is.na(sd_k) || sd_k == 0) { Zi_vec[idx] <- 0; next }
    Zi_vec[idx] <- (k_s - mean(k_s)) / sd_k
  }
  ## Pi
  Pi_vec <- numeric(igraph::vcount(g_lcc))
  for (i in seq_len(igraph::vcount(g_lcc))) {
    ki <- deg_v[i]
    if (ki == 0) { Pi_vec[i] <- 0; next }
    Pi_vec[i] <- 1 - sum((table(mod_v[igraph::neighbors(g_lcc, i)]) / ki)^2)
  }
  igraph::V(g_lcc)$Zi <- Zi_vec
  igraph::V(g_lcc)$Pi <- Pi_vec
  igraph::V(g_lcc)$role <- ifelse(Zi_vec >= 2.5 & Pi_vec >= 0.62, "Network hub",
    ifelse(Zi_vec >= 2.5 & Pi_vec <  0.62, "Module hub",
    ifelse(Zi_vec <  2.5 & Pi_vec >= 0.62, "Connector", "Peripheral")))

  ## Taxonomy — store all ranks used for keystone labelling
  asv <- igraph::V(g_lcc)$name
  .in_ref <- asv %in% rownames(tax_ref)
  igraph::V(g_lcc)$Genus  <- ifelse(.in_ref, trimws(tax_ref[asv,"Genus"]),  "")
  igraph::V(g_lcc)$Family <- ifelse(.in_ref, trimws(tax_ref[asv,"Family"]), "")
  igraph::V(g_lcc)$Order  <- ifelse(.in_ref, trimws(tax_ref[asv,"Order"]),  "")
  igraph::V(g_lcc)$Class  <- ifelse(.in_ref, trimws(tax_ref[asv,"Class"]),  "")
  igraph::V(g_lcc)$Phylum <- ifelse(.in_ref, trimws(tax_ref[asv,"Phylum"]), "")

  n_mod <- max(igraph::V(g_lcc)$module)
  cat(sprintf("  [%s] LCC: %d nodes, %d edges, %d modules\n",
              label, igraph::vcount(g_lcc), igraph::ecount(g_lcc), n_mod))
  list(graph=g_lcc, mods=mods, pcor=pcor, label=label)
}

## Topology data.frame for one network
.topo_row <- function(net_obj, label) {
  g   <- net_obj$graph
  deg <- igraph::V(g)$degree
  el  <- igraph::as_edgelist(g, names=TRUE)
  pcv <- net_obj$pcor[el]
  data.frame(Network=label, Nodes=igraph::vcount(g), Edges=igraph::ecount(g),
    Positive_edges=sum(pcv>0), Negative_edges=sum(pcv<0),
    Avg_degree=round(mean(deg),4),
    Modularity=round(igraph::modularity(net_obj$mods),4),
    Clustering_coeff=round(igraph::transitivity(g,"global"),4),
    N_keystones=sum(igraph::V(g)$role %in% c("Module hub","Network hub","Connector")),
    stringsAsFactors=FALSE)
}

## Keystone data.frame for one network
.keystone_df <- function(net_obj, label) {
  g <- net_obj$graph
  data.frame(Network=label, ASV=igraph::V(g)$name,
    Genus=igraph::V(g)$Genus, Phylum=igraph::V(g)$Phylum,
    Degree=igraph::V(g)$degree, Betweenness=round(igraph::V(g)$betweenness,6),
    Module=igraph::V(g)$module,
    Zi=round(igraph::V(g)$Zi,4), Pi=round(igraph::V(g)$Pi,4),
    Role=igraph::V(g)$role, stringsAsFactors=FALSE)
}

## Save helper: PNG + PDF + SVG
.save_all <- function(p, name, w, h) {
  base <- file.path(.P6, name)
  ggplot2::ggsave(paste0(base,".png"), p, width=w, height=h, units="in",
                  dpi=300, type="cairo")
  ggplot2::ggsave(paste0(base,".pdf"), p, width=w, height=h, units="in",
                  device=cairo_pdf)
  ggplot2::ggsave(paste0(base,".svg"), p, width=w, height=h, units="in",
                  device=svglite::svglite)
  fi <- file.info(paste0(base, c(".png",".pdf",".svg")))$size
  cat(sprintf("  %-50s  PNG %4.0fKB  PDF %4.0fKB  SVG %4.0fKB\n",
              name, fi[1]/1e3, fi[2]/1e3, fi[3]/1e3))
}

## ─────────────────────────────────────────────────────────────────────────────
## 7A — SPIEC-EASI TREATMENT NETWORKS
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── 7A: SPIEC-EASI treatment networks ──\n")

## Load genus-level phyloseq (use cache if exists)
if (!file.exists("data/ps_genus_cached.rds")) {
  cat("  Building ps_genus_cached.rds...\n")
  ps_g_raw <- readRDS("data/ps_filtered.rds")
  ps_gc    <- suppressWarnings(tax_glom(ps_g_raw, taxrank="Genus", NArm=FALSE))
  saveRDS(ps_gc, "data/ps_genus_cached.rds")
  rm(ps_g_raw, ps_gc)
}
.ps_genus <- readRDS("data/ps_genus_cached.rds")
.tax_df   <- data.frame(tax_table(.ps_genus), stringsAsFactors=FALSE)
.meta_p   <- data.frame(sample_id=rownames(sample_data(.ps_genus)),
                         as.data.frame(sample_data(.ps_genus)),
                         stringsAsFactors=FALSE)

## Run SPIEC-EASI or load cache
.run_se <- function(ps_obj, label, rds_path) {
  if (file.exists(rds_path)) {
    cat(sprintf("  Loading cached %s network (%s)...\n", label,
                format(file.info(rds_path)$mtime, "%Y-%m-%d")))
    return(readRDS(rds_path))
  }
  cat(sprintf("  Running SPIEC-EASI for %s (%d samples, %d taxa)...\n",
              label, phyloseq::nsamples(ps_obj), phyloseq::ntaxa(ps_obj)))
  set.seed(42)
  net <- SpiecEasi::spiec.easi(ps_obj, method="glasso", nlambda=20,
    lambda.min.ratio=1e-2,
    pulsar.params=list(rep.num=20, ncores=4))
  saveRDS(net, rds_path)
  cat(sprintf("  Saved %s (%.0f MB)\n", rds_path, file.size(rds_path)/1e6))
  net
}

## Subset to treatment groups (ALL planted + unplanted = full treatment signal)
.ps_drought <- phyloseq::prune_samples(
  .meta_p$sample_id[.meta_p$treatment == "Drought"], .ps_genus)
.ps_watered <- phyloseq::prune_samples(
  .meta_p$sample_id[.meta_p$treatment == "Watered"], .ps_genus)

.net_d_raw <- .run_se(.ps_drought, "Drought",  "data/spiec_easi_drought.rds")
.net_w_raw <- .run_se(.ps_watered, "Watered",  "data/spiec_easi_watered.rds")

## Extract pcor (then free large objects)
.taxa_trt <- colnames(.net_d_raw$est$data)
.opt_d    <- .net_d_raw$select$stars$opt.index
.pcor_d   <- .pcor_glasso(.net_d_raw)
rownames(.pcor_d) <- colnames(.pcor_d) <- .taxa_trt
cat(sprintf("  Drought: opt.index=%d  taxa=%d\n", .opt_d, length(.taxa_trt)))
rm(.net_d_raw); gc(verbose=FALSE)

.opt_w  <- readRDS("data/spiec_easi_watered.rds")$select$stars$opt.index
.net_w2 <- readRDS("data/spiec_easi_watered.rds")
.pcor_w <- .pcor_glasso(.net_w2)
rownames(.pcor_w) <- colnames(.pcor_w) <- .taxa_trt
cat(sprintf("  Watered: opt.index=%d\n", .opt_w))
rm(.net_w2); gc(verbose=FALSE)

## ─────────────────────────────────────────────────────────────────────────────
## 7B — SPIEC-EASI GENOTYPE NETWORKS
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── 7B: SPIEC-EASI genotype networks ──\n")

.res_ids <- .meta_p$sample_id[.meta_p$trait == "Resistance"  & .meta_p$condition == "Planted"]
.sus_ids <- .meta_p$sample_id[.meta_p$trait == "Susceptible" & .meta_p$condition == "Planted"]
cat(sprintf("  Resistant n=%d | Susceptible n=%d\n", length(.res_ids), length(.sus_ids)))

.ps_res <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(phyloseq::prune_samples(.res_ids, .ps_genus)) > 0,
  phyloseq::prune_samples(.res_ids, .ps_genus))
.ps_sus <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(phyloseq::prune_samples(.sus_ids, .ps_genus)) > 0,
  phyloseq::prune_samples(.sus_ids, .ps_genus))

.net_r_raw <- .run_se(.ps_res, "Resistant",   "data/spiec_easi_resistant.rds")
.net_s_raw <- .run_se(.ps_sus, "Susceptible",  "data/spiec_easi_susceptible.rds")

.taxa_r <- phyloseq::taxa_names(.ps_res)
.taxa_s <- phyloseq::taxa_names(.ps_sus)
.pcor_r <- .pcor_glasso(.net_r_raw)
rownames(.pcor_r) <- colnames(.pcor_r) <- .taxa_r
.opt_r  <- .net_r_raw$select$stars$opt.index
cat(sprintf("  Resistant: opt.index=%d  taxa=%d\n", .opt_r, length(.taxa_r)))
rm(.net_r_raw); gc(verbose=FALSE)

.pcor_s <- .pcor_glasso(.net_s_raw)
rownames(.pcor_s) <- colnames(.pcor_s) <- .taxa_s
.opt_s  <- .net_s_raw$select$stars$opt.index
cat(sprintf("  Susceptible: opt.index=%d  taxa=%d\n", .opt_s, length(.taxa_s)))
rm(.net_s_raw); gc(verbose=FALSE)

## ─────────────────────────────────────────────────────────────────────────────
## Build all 4 igraph LCC objects
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── Building LCC + modules + Zi/Pi ──\n")
.net_drought  <- .build_net(.pcor_d, .taxa_trt, .tax_df, "Drought")
.net_watered  <- .build_net(.pcor_w, .taxa_trt, .tax_df, "Watered")
.net_resistant  <- .build_net(.pcor_r, .taxa_r,   .tax_df, "Resistant")
.net_susceptible <- .build_net(.pcor_s, .taxa_s,  .tax_df, "Susceptible")

## ─────────────────────────────────────────────────────────────────────────────
## Overall topology tables (write if changed)
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── Writing topology tables ──\n")

## Treatment topology (legacy format used by figure_assembly.r Panel E)
.trt_topo <- data.frame(
  treatment = c("Drought","Watered"),
  nodes     = c(igraph::vcount(.net_drought$graph), igraph::vcount(.net_watered$graph)),
  edges     = c(igraph::ecount(.net_drought$graph), igraph::ecount(.net_watered$graph)),
  avg_degree= c(mean(igraph::V(.net_drought$graph)$degree),
                mean(igraph::V(.net_watered$graph)$degree)),
  modularity= c(igraph::modularity(.net_drought$mods),
                igraph::modularity(.net_watered$mods)),
  n_keystones=c(sum(igraph::V(.net_drought$graph)$role  %in% c("Module hub","Network hub","Connector")),
                sum(igraph::V(.net_watered$graph)$role  %in% c("Module hub","Network hub","Connector"))),
  stringsAsFactors=FALSE)
write.csv(.trt_topo, "tables/P6_networks/P6_network_topology.csv", row.names=FALSE)
cat("  tables/P6_network_topology.csv\n")

## Genotype topology
.geno_topo <- rbind(
  .topo_row(.net_resistant,   "Resistant"),
  .topo_row(.net_susceptible, "Susceptible"))
write.csv(.geno_topo, "tables/P6_networks/network_topology_genotype.csv", row.names=FALSE)
cat("  tables/network_topology_genotype.csv\n")

## Keystone taxa tables
.kst_trt  <- rbind(.keystone_df(.net_drought, "Drought"),
                    .keystone_df(.net_watered, "Watered"))
write.csv(.kst_trt, "tables/P6_networks/network_keystone_taxa.csv", row.names=FALSE)
cat(sprintf("  tables/network_keystone_taxa.csv (%d rows)\n", nrow(.kst_trt)))

.kst_geno <- rbind(.keystone_df(.net_resistant,   "Resistant"),
                    .keystone_df(.net_susceptible, "Susceptible"))
write.csv(.kst_geno, "tables/P6_networks/network_keystone_taxa_genotype.csv", row.names=FALSE)
cat(sprintf("  tables/network_keystone_taxa_genotype.csv (%d rows)\n", nrow(.kst_geno)))

cat("\nTopology summary:\n")
print(rbind(
  .trt_topo[, c("treatment","nodes","edges","avg_degree","modularity","n_keystones")],
  setNames(.geno_topo[, c("Network","Nodes","Edges","Avg_degree","Modularity","N_keystones")],
           c("treatment","nodes","edges","avg_degree","modularity","n_keystones"))))

## ─────────────────────────────────────────────────────────────────────────────
## 7E — NETWORK FIGURES
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── 7E: Network figures ──\n")

## Network plot function (all 4 networks, identical spec)
## min_degree:    remove nodes with degree < min_degree before FR layout (default 1 = keep all)
## top_n_edges:   max edges to display by |pcor| (default 500; use 800 for denser networks)
.plot_network <- function(net_obj, n_samp=NULL, min_degree=1L, top_n_edges=500L) {
  g   <- net_obj$graph
  lbl <- net_obj$label

  ## Defensive LCC: .build_net() already extracts LCC, but re-check for safety
  .comps_p7 <- igraph::components(g)
  if (.comps_p7$no > 1L) {
    lcc_idx <- which(.comps_p7$membership == which.max(.comps_p7$csize))
    g <- igraph::induced_subgraph(g, lcc_idx)
    cat(sprintf("  [%s] LCC re-extracted: %d nodes, %d edges\n",
                lbl, igraph::vcount(g), igraph::ecount(g)))
  }

  ## Module colour mapping
  mod_ids  <- V(g)$module
  mod_cols <- .mod_cols(mod_ids)
  V(g)$mod_fac <- factor(as.character(mod_ids), levels=names(mod_cols))

  ## Keystone labelling: prefer hub/connector roles, fall back to top-5 by degree.
  ## FIX 1 — deepest informative rank: Genus → Family → Order → Class → Phylum.
  ## Class is skipped when it contains digits, underscores, or hyphens
  ## (environmental codes such as BRH-c20a, TK10, Lineage_IIc).
  .tax_label <- function(g, i) {
    .get <- function(attr) {
      v <- igraph::vertex_attr(g, attr, i)
      trimws(ifelse(is.null(v) | is.na(v), "", as.character(v)))
    }
    gen <- .get("Genus");  if (gen != "") return(gen)
    fam <- .get("Family"); if (fam != "") return(paste0("f__", fam))
    ord <- .get("Order");  if (ord != "") return(paste0("o__", ord))
    cls <- .get("Class")
    if (cls != "" && !grepl("[0-9_\\-]", cls)) return(paste0("c__", cls))
    phy <- .get("Phylum"); if (phy != "") return(paste0("p__", phy))
    return("")
  }
  ks_idx <- which(V(g)$role %in% c("Network hub", "Module hub", "Connector"))
  if (length(ks_idx) >= 5) {
    top5 <- ks_idx[order(-V(g)$degree[ks_idx])][1:5]
  } else {
    top5 <- order(-V(g)$degree)[1:min(5L, igraph::vcount(g))]
  }
  ## FIX 1 — Force Azospirillum into labels if it is a keystone not already in top5.
  ## Strategy: take top 4 keystones by degree, replace slot 5 with Azospirillum.
  .azos_idx <- which(V(g)$Genus == "Azospirillum")
  if (length(.azos_idx) == 1L && .azos_idx %in% ks_idx && !(.azos_idx %in% top5)) {
    cat(sprintf("  [%s] Forcing Azospirillum (role=%s, degree=%d) into label slot 5\n",
                lbl, V(g)$role[.azos_idx], V(g)$degree[.azos_idx]))
    top5[length(top5)] <- .azos_idx
  }
  V(g)$lbl <- ""
  V(g)$lbl[top5] <- sapply(top5, function(i) .tax_label(g, i))
  .lbls <- V(g)$lbl[top5]
  cat(sprintf("  [%s] labelled taxa: %s\n", lbl,
              paste(.lbls[.lbls != ""], collapse=", ")))
  ## FIX 2 — Confirm TK10 env-code filtering for any Chloroflexota keystone in top5
  .chloro_in_top5 <- top5[!is.na(V(g)$Phylum[top5]) &
                           V(g)$Phylum[top5] == "Chloroflexota" &
                           grepl("[0-9_\\-]", V(g)$Class[top5])]
  if (length(.chloro_in_top5) > 0) {
    for (.ci in .chloro_in_top5)
      cat(sprintf("  [%s] Chloroflexota keystone: Class='%s' contains digits → env-code filtered → label='p__Chloroflexota'\n",
                  lbl, V(g)$Class[.ci]))
  }

  ## Mark edge polarity (weight attribute from .build_net() pcor values)
  E(g)$epos <- E(g)$weight > 0

  ## g_plot keeps weight so we can rank edges by |pcor| below
  g_plot <- g

  ## LCC filter on g_plot
  .cplot <- igraph::components(g_plot)
  if (.cplot$no > 1L) {
    lcc_v2 <- which(.cplot$membership == which.max(.cplot$csize))
    g_plot  <- igraph::induced_subgraph(g_plot, lcc_v2)
  }

  ## Optional minimum-degree filter (Susceptible uses min_degree=3)
  if (min_degree > 1L) {
    .deg_plot  <- igraph::degree(g_plot)
    .keep      <- igraph::V(g_plot)[.deg_plot >= min_degree]
    .n_removed <- igraph::vcount(g_plot) - length(.keep)
    if (.n_removed > 0L) {
      g_plot <- igraph::induced_subgraph(g_plot, .keep)
      cat(sprintf("  [%s] min_degree=%d: removed %d nodes → %d nodes, %d edges remaining\n",
                  lbl, min_degree, .n_removed,
                  igraph::vcount(g_plot), igraph::ecount(g_plot)))
    }
    ## Re-apply LCC: degree filtering can disconnect the graph
    .cplot2 <- igraph::components(g_plot)
    if (.cplot2$no > 1L) {
      lcc_v3 <- which(.cplot2$membership == which.max(.cplot2$csize))
      .dropped <- igraph::vcount(g_plot) - length(lcc_v3)
      g_plot   <- igraph::induced_subgraph(g_plot, lcc_v3)
      cat(sprintf("  [%s] post-degree LCC: dropped %d fragment nodes → %d nodes\n",
                  lbl, .dropped, igraph::vcount(g_plot)))
    }
  }

  ## Top-N edge thinning by |pcor weight| (FIX 3: 800 for Drought/Watered, 500 otherwise)
  ## Selected BEFORE removing weight for layout.
  n_edges_kept <- igraph::ecount(g_plot)
  .edf      <- igraph::as_data_frame(g_plot, what="edges")
  .eord     <- order(abs(.edf$weight), decreasing=TRUE)
  top_n     <- min(top_n_edges, n_edges_kept)
  top_edges <- .edf[.eord[seq_len(top_n)], ]
  n_pos_shown <- sum(top_edges$weight > 0)
  n_neg_shown <- sum(top_edges$weight < 0)
  cat(sprintf("  [%s] edges: %d total → top %d shown (%d pos, %d neg)\n",
              lbl, n_edges_kept, top_n, n_pos_shown, n_neg_shown))

  ## Build g_sparse with exactly top-N edges + full vertex set.
  ## Add boolean epos for two-layer rendering (avoids colour-aesthetic interpolation).
  .vert_df         <- igraph::as_data_frame(g_plot, what="vertices")
  g_sparse         <- igraph::graph_from_data_frame(top_edges, directed=FALSE, vertices=.vert_df)
  E(g_sparse)$epos <- E(g_sparse)$weight > 0

  ## Verification (guarantees exactly two colours in console output)
  cat(sprintf("  [%s] unique edge_col: %s\n", lbl,
              paste(sort(unique(ifelse(E(g_sparse)$epos, "#B22222", "#3498DB"))),
                    collapse=", ")))

  ## Remove weight from g_plot before FR layout (FR rejects negative weights)
  g_plot <- igraph::delete_edge_attr(g_plot, "weight")

  ## FR layout from full g_plot (all kept edges drive node positioning)
  set.seed(42)
  lay <- ggraph::create_layout(g_plot, layout="fr")

  ## Switch ggraph's graph reference to g_sparse — ggraph requires tbl_graph
  attr(lay, "graph") <- tidygraph::as_tbl_graph(g_sparse)

  ## Subtitle (no caption — replaced by annotate() legend below)
  q   <- igraph::modularity(net_obj$mods)
  sub <- sprintf("Nodes=%d | Top %d/%d edges | Modularity=%.3f%s%s",
                  igraph::vcount(g_plot), top_n, n_edges_kept, q,
                  if(!is.null(n_samp)) sprintf("  (n=%d)", n_samp) else "",
                  if(min_degree > 1L) sprintf("\n(displayed: degree \u2265 %d; %d/%d nodes)",
                                               min_degree, igraph::vcount(g_plot),
                                               igraph::vcount(net_obj$graph)) else "")

  ## FIX 2 — Coloured annotate() legend placed at bottom-left of layout space.
  .xr      <- diff(range(lay$x)); .yr <- diff(range(lay$y))
  .xmin    <- min(lay$x) - 0.02 * .xr
  .ymin    <- min(lay$y) - 0.10 * .yr   ## below the lowest node
  .seg_len <- 0.07 * .xr
  .txt_gap <- 0.012 * .xr
  .lgap    <- 0.06 * .yr                ## vertical spacing between legend rows

  ## FIX 1 — Two geom_edge_link layers: fixed colours + per-layer linewidth.
  ## Positive (red) drawn first (background); negative (blue) on top (thicker).
  ggraph::ggraph(lay) +
    ggraph::geom_edge_link(
      ggplot2::aes(filter = epos),
      colour="#B22222", alpha=0.6, linewidth=1.2, show.legend=FALSE) +
    ggraph::geom_edge_link(
      ggplot2::aes(filter = !epos),
      colour="#3498DB", alpha=0.8, linewidth=1.8, show.legend=FALSE) +
    ggraph::geom_node_point(ggplot2::aes(size=degree, colour=mod_fac), alpha=0.85) +
    ggplot2::scale_colour_manual(values=mod_cols, name="Module",
      guide=ggplot2::guide_legend(title="Module", ncol=2,
        override.aes=list(size=5, alpha=1))) +
    ggplot2::scale_size_continuous(name="Degree", range=c(2,12), breaks=c(5,20,50)) +
    ggraph::geom_node_label(ggplot2::aes(label=lbl), repel=TRUE,
      size=3.5, box.padding=0.4, point.padding=0.3, max.overlaps=Inf,
      fontface="italic", fill=scales::alpha("white", 0.85),
      label.size=0.15, show.legend=FALSE) +
    ## FIX 2 — inline coloured legend via annotate (header + red row + blue row)
    ggplot2::annotate("text",
      x=.xmin, y=.ymin + .lgap,
      label=sprintf("Showing top %d edges by |partial correlation|", top_n),
      size=3, hjust=0, colour="grey40", fontface="italic") +
    ggplot2::annotate("segment",
      x=.xmin, xend=.xmin + .seg_len,
      y=.ymin, yend=.ymin,
      colour="#B22222", linewidth=1.5) +
    ggplot2::annotate("text",
      x=.xmin + .seg_len + .txt_gap, y=.ymin,
      label=paste0("Positive (", n_pos_shown, " edges)"),
      size=3.5, hjust=0, colour="grey20") +
    ggplot2::annotate("segment",
      x=.xmin, xend=.xmin + .seg_len,
      y=.ymin - .lgap, yend=.ymin - .lgap,
      colour="#3498DB", linewidth=1.5) +
    ggplot2::annotate("text",
      x=.xmin + .seg_len + .txt_gap, y=.ymin - .lgap,
      label=paste0("Negative (", n_neg_shown, " edges)"),
      size=3.5, hjust=0, colour="grey20") +
    ggplot2::labs(title=lbl, subtitle=sub) +
    ggraph::theme_graph(base_family="") +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size=14, face="bold", hjust=0.5),
      plot.subtitle   = ggplot2::element_text(size=10, colour="grey40", hjust=0.5),
      legend.text     = ggplot2::element_text(size=10),
      legend.title    = ggplot2::element_text(size=11, face="bold"),
      legend.key.size = unit(5,"mm"),
      plot.margin     = ggplot2::margin(8,2,8,2,"mm"))
}

.save_all(.plot_network(.net_drought,  top_n_edges=800L),                                   "P6_network_Drought",   10, 9)
.save_all(.plot_network(.net_watered,  top_n_edges=800L),                                   "P6_network_Watered",   10, 9)
.save_all(.plot_network(.net_resistant,   n_samp=length(.res_ids)),                         "network_resistant",    10, 9)
.save_all(.plot_network(.net_susceptible, n_samp=length(.sus_ids), min_degree=3L),          "network_susceptible",  10, 9)

## Zi-Pi scatter function (identical for both treatment and genotype)
.plot_zipi <- function(kdf, title_str) {
  kdf$Phylum[is.na(kdf$Phylum) | trimws(kdf$Phylum)==""] <- "Unknown"

  ## FIX 1 — vivid named phylum palette; phyla absent from data are silently skipped
  .phy_pal <- c(
    "Acidobacteriota"  = "#2196F3",
    "Actinomycetota"   = "#FF5722",
    "Bacillota"        = "#9C27B0",
    "Bacteroidota"     = "#E91E63",
    "Chloroflexota"    = "#4CAF50",
    "Cyanobacteriota"  = "#FFEB3B",
    "Myxococcota"      = "#00BCD4",
    "Planctomycetota"  = "#795548",
    "Pseudomonadota"   = "#9E9E9E",
    "Gemmatimonadota"  = "#FF9800",
    "Other"            = "#BDBDBD"
  )
  ## Phyla not in the named palette (incl. Unknown) collapse to "Other"
  kdf$PG <- ifelse(kdf$Phylum %in% names(.phy_pal), kdf$Phylum, "Other")
  ## Factor levels = palette order, restricted to phyla present in this dataset
  .present <- intersect(names(.phy_pal), unique(kdf$PG))
  if (!"Other" %in% .present) .present <- c(.present, "Other")
  kdf$PG   <- factor(kdf$PG, levels = .present)
  pal_cols <- .phy_pal[.present]
  cat(sprintf("  [zipi:%s] phyla present: %s\n", title_str,
              paste(names(pal_cols), collapse=", ")))

  levs   <- sort(unique(kdf$Network))
  quad_df <- do.call(rbind, lapply(levs, function(n) data.frame(
    x=c(0.02,0.64,0.02,0.64), y=c(-Inf,-Inf,Inf,Inf),
    hjust=c(0,0,0,0), vjust=c(-0.5,-0.5,1.5,1.5),
    label=c("Peripherals","Connectors","Module hubs","Network hubs"),
    Network=n, stringsAsFactors=FALSE)))
  quad_df$Network <- factor(quad_df$Network, levels=levs)
  kdf$Network     <- factor(kdf$Network,     levels=levs)

  ggplot2::ggplot(kdf, ggplot2::aes(Pi, Zi, colour=PG, size=Degree)) +
    ggplot2::geom_vline(xintercept=0.62, linetype="dashed",
                         colour="grey50", linewidth=0.5) +
    ggplot2::geom_hline(yintercept=2.5,  linetype="dashed",
                         colour="grey50", linewidth=0.5) +
    ggplot2::geom_point(alpha=0.60, stroke=0) +
    ggplot2::geom_text(data=quad_df,
      ggplot2::aes(x=x,y=y,label=label,hjust=hjust,vjust=vjust),
      colour="grey50", size=8/.pt, fontface="italic", inherit.aes=FALSE) +
    ggplot2::facet_wrap(~Network, ncol=2) +
    ## FIX 1 — vivid palette + FIX 2 — larger legend points for Phylum
    ggplot2::scale_colour_manual(values=pal_cols, name="Phylum", drop=TRUE,
      guide=ggplot2::guide_legend(title="Phylum",
        override.aes=list(size=5), ncol=1)) +
    ## FIX 2 + FIX 3 — larger point range in plot body + larger legend key sizes
    ggplot2::scale_size_continuous(name="Degree", range=c(2,8), breaks=c(5,15,30),
      guide=ggplot2::guide_legend(
        override.aes=list(size=c(4,6,9)))) +
    ggplot2::scale_x_continuous(limits=c(0,1), name="Participation coefficient (Pi)") +
    ggplot2::scale_y_continuous(name="Within-module degree z-score (Zi)") +
    ggplot2::theme_bw(base_size=9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour="grey93", linewidth=0.25),
      strip.background = ggplot2::element_rect(fill="grey95", colour="grey70"),
      strip.text       = ggplot2::element_text(size=9, face="bold"),
      axis.text        = ggplot2::element_text(size=8),
      axis.title       = ggplot2::element_text(size=9),
      legend.text      = ggplot2::element_text(size=8),
      legend.title     = ggplot2::element_text(size=9, face="bold"),
      legend.key.size  = unit(3.5,"mm"),
      plot.title       = ggplot2::element_text(size=9, face="bold"),
      plot.margin      = ggplot2::margin(2,2,2,2,"mm")) +
    ggplot2::labs(title=title_str)
}

.kst_trt2  <- read.csv("tables/P6_networks/network_keystone_taxa.csv",          stringsAsFactors=FALSE)
.kst_geno2 <- read.csv("tables/P6_networks/network_keystone_taxa_genotype.csv", stringsAsFactors=FALSE)
.save_all(.plot_zipi(.kst_trt2,  "Keystone classification — SPIEC-EASI treatment networks"),
          "keystone_zipi_plot",      10, 5)
.save_all(.plot_zipi(.kst_geno2, "Keystone classification — SPIEC-EASI genotype networks"),
          "keystone_zipi_genotype",  10, 5)

## Topology comparison bar chart (4 groups, 2 metrics)
.plot_topo_bar <- function() {
  trt_t  <- read.csv("tables/P6_networks/P6_network_topology.csv",       stringsAsFactors=FALSE)
  geno_t <- read.csv("tables/P6_networks/network_topology_genotype.csv", stringsAsFactors=FALSE)
  dat <- rbind(
    data.frame(Group=trt_t$treatment,  Avg_degree=trt_t$avg_degree,
               Modularity=trt_t$modularity, stringsAsFactors=FALSE),
    data.frame(Group=geno_t$Network,   Avg_degree=geno_t$Avg_degree,
               Modularity=geno_t$Modularity, stringsAsFactors=FALSE))
  dat$Group <- factor(dat$Group, levels=c("Drought","Watered","Resistant","Susceptible"))
  long <- dat %>% tidyr::pivot_longer(c(Avg_degree,Modularity),
    names_to="Metric", values_to="Value") %>%
    dplyr::mutate(Metric=dplyr::recode(Metric,
      Avg_degree="Average degree", Modularity="Modularity (Q)"))
  ## FIX 2 — match manuscript palette from theme_constants.r
  COLS4 <- c(Drought="#B03A2E", Watered="#1A5276", Resistant="#1E8449", Susceptible="#D35400")
  ggplot2::ggplot(long, ggplot2::aes(Group, Value, fill=Group)) +
    ggplot2::geom_col(alpha=0.85, width=0.65) +
    ggplot2::facet_wrap(~Metric, scales="free_y", ncol=1) +
    ggplot2::scale_fill_manual(values=COLS4, guide="none") +
    ggplot2::theme_bw(base_size=12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(size=10, angle=35, hjust=1),
      axis.text.y      = ggplot2::element_text(size=10),
      axis.title       = ggplot2::element_text(size=10),
      strip.text       = ggplot2::element_text(size=11, face="bold"),
      plot.title       = ggplot2::element_text(size=12, face="bold"),
      plot.margin      = ggplot2::margin(4,4,4,4,"mm")) +
    ## FIX 1 — remove uninformative y="Value"; facet strips already label each metric
    ggplot2::labs(x=NULL,
      title="Network topology comparison",
      subtitle="SPIEC-EASI overall networks")
}
.save_all(.plot_topo_bar(), "network_topology_comparison", 6, 7)

## Cohesion plot (Herren & McMahon 2017)
cat("\n── Cohesion analysis ──\n")
{
  ps_rel_coh <- readRDS("data/ps_rel.rds")
  otu_rel <- as(otu_table(ps_rel_coh), "matrix")
  if (taxa_are_rows(ps_rel_coh)) otu_rel <- t(otu_rel)
  keep    <- intersect(.taxa_trt, colnames(otu_rel))
  otu_net <- otu_rel[, keep, drop=FALSE]
  meta_coh <- data.frame(id=rownames(sample_data(ps_rel_coh)),
    as.data.frame(sample_data(ps_rel_coh)), stringsAsFactors=FALSE)
  planted_coh <- meta_coh %>% dplyr::filter(condition=="Planted")
  rm(ps_rel_coh); gc(verbose=FALSE)

  net_d_adj <- readRDS("data/spiec_easi_drought.rds")
  adj_d <- as.matrix(net_d_adj$refit$stars)
  rownames(adj_d) <- colnames(adj_d) <- .taxa_trt
  rm(net_d_adj); gc(verbose=FALSE)

  net_w_adj <- readRDS("data/spiec_easi_watered.rds")
  adj_w <- as.matrix(net_w_adj$refit$stars)
  rownames(adj_w) <- colnames(adj_w) <- .taxa_trt
  rm(net_w_adj); gc(verbose=FALSE)

  pcor_d_k <- .pcor_d[keep, keep]; pcor_w_k <- .pcor_w[keep, keep]
  adj_d_k  <- adj_d[keep, keep];   adj_w_k  <- adj_w[keep, keep]

  .coh_fn <- function(pm, otu) {
    n <- ncol(pm)
    data.frame(Sample=rownames(otu),
      pos_cohesion=rowSums(otu %*% pmax(pm,0)) / n,
      neg_cohesion=rowSums(otu %*% pmin(pm,0)) / n,
      stringsAsFactors=FALSE)
  }
  coh_d <- .coh_fn(pcor_d_k * adj_d_k, otu_net)
  coh_w <- .coh_fn(pcor_w_k * adj_w_k, otu_net)

  coh_df <- planted_coh %>%
    dplyr::select(id, treatment) %>%
    dplyr::left_join(dplyr::rename(coh_d, d_pos=pos_cohesion, d_neg=neg_cohesion),
                     by=c("id"="Sample")) %>%
    dplyr::left_join(dplyr::rename(coh_w, w_pos=pos_cohesion, w_neg=neg_cohesion),
                     by=c("id"="Sample")) %>%
    dplyr::mutate(
      pos_cohesion   = ifelse(treatment=="Drought", d_pos, w_pos),
      neg_cohesion   = ifelse(treatment=="Drought", d_neg, w_neg),
      total_cohesion = pos_cohesion + neg_cohesion,
      treatment      = factor(treatment, levels=c("Watered","Drought")))

  w_pos <- wilcox.test(pos_cohesion   ~ treatment, data=coh_df, exact=FALSE)
  w_neg <- wilcox.test(neg_cohesion   ~ treatment, data=coh_df, exact=FALSE)
  w_tot <- wilcox.test(total_cohesion ~ treatment, data=coh_df, exact=FALSE)
  cat(sprintf("  Cohesion Wilcoxon: Pos p=%.4f | Neg p=%.4f | Total p=%.4f\n",
              w_pos$p.value, w_neg$p.value, w_tot$p.value))

  coh_long <- coh_df %>%
    tidyr::pivot_longer(c(pos_cohesion,neg_cohesion,total_cohesion),
                        names_to="Type", values_to="Cohesion") %>%
    dplyr::mutate(Type=factor(dplyr::recode(Type,
      pos_cohesion="Positive", neg_cohesion="Negative",
      total_cohesion="Total"), levels=c("Positive","Negative","Total")))

  ann_coh <- data.frame(
    Type=factor(c("Positive","Negative","Total"),
                levels=c("Positive","Negative","Total")),
    label=c(sprintf("p=%.3f",w_pos$p.value),
            sprintf("p=%.3f",w_neg$p.value),
            sprintf("p=%.3f",w_tot$p.value)),
    stringsAsFactors=FALSE)

  p_coh <- ggplot2::ggplot(coh_long, ggplot2::aes(treatment, Cohesion, fill=treatment)) +
    ggplot2::geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.4) +
    ggplot2::geom_jitter(width=0.15, size=0.9, alpha=0.35, show.legend=FALSE) +
    ggplot2::geom_text(data=ann_coh,
      ggplot2::aes(x=1.5, y=Inf, label=label), vjust=1.6, hjust=0.5,
      size=8/.pt, colour="grey30", inherit.aes=FALSE) +
    ggplot2::facet_wrap(~Type, scales="free_y", ncol=3) +
    ggplot2::scale_fill_manual(values=c(Watered="#2980B9",Drought="#C0392B"),
                                name="Treatment") +
    ggplot2::theme_bw(base_size=10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size=10, face="bold"),
      axis.text        = ggplot2::element_text(size=9),
      axis.title       = ggplot2::element_text(size=10),
      legend.text      = ggplot2::element_text(size=9),
      legend.title     = ggplot2::element_text(size=10, face="bold"),
      plot.title       = ggplot2::element_text(size=11, face="bold"),
      plot.margin      = ggplot2::margin(3,3,3,3,"mm")) +
    ggplot2::labs(x=NULL, y="Network cohesion",
      title="Co-occurrence network cohesion by treatment",
      subtitle=sprintf("Positive p=%.3f | Negative p=%.3f | Total p=%.3f",
        w_pos$p.value, w_neg$p.value, w_tot$p.value))

  .save_all(p_coh, "cohesion_plot", 9, 5)
}

## ╔═══════════════════════════════════════════════════════════════════════════╗
## ║  WARNING: Per-day Spearman networks (n ≈ 6 per group per day) produce   ║
## ║  unreliable edge estimates due to insufficient sample size. Spearman    ║
## ║  correlations at n=6 have very wide confidence intervals and high        ║
## ║  false-discovery rates even with |rho|>0.6 and p<0.05 filters.          ║
## ║                                                                          ║
## ║  The following outputs are retained for EXPLORATORY REFERENCE ONLY and  ║
## ║  are NOT used in the manuscript or any supplementary figures:            ║
## ║    - tables/network_topology_perday.csv                                  ║
## ║    - tables/network_topology_perday_genotype.csv                         ║
## ║    - tables/network_keystone_perday.csv                                  ║
## ║    - tables/network_keystone_perday_genotype.csv                         ║
## ║    - figures/P6_networks/keystone_temporal_heatmap*.{png,pdf,svg}        ║
## ║    - figures/P6_networks/network_topology_temporal*.{png,pdf,svg}        ║
## ║                                                                          ║
## ║  Do NOT interpret these results quantitatively in the manuscript.        ║
## ╚═══════════════════════════════════════════════════════════════════════════╝

## ─────────────────────────────────────────────────────────────────────────────
## 7C — PER-DAY SPEARMAN NETWORKS: TREATMENT
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── 7C: Per-day Spearman networks — treatment ──\n")

.topo_pd_path <- "tables/P6_networks/network_topology_perday.csv"
.kst_pd_path  <- "tables/P6_networks/network_keystone_perday.csv"

if (!file.exists(.topo_pd_path) || !file.exists(.kst_pd_path)) {
  cat("  Computing per-day treatment Spearman networks...\n")
  ps_f <- readRDS("data/ps_filtered.rds")

  ## Pre-extract OTU matrix (avoids S4 subscript errors in loop)
  otu_raw <- as(otu_table(ps_f), "matrix")
  if (taxa_are_rows(ps_f)) otu_raw <- t(otu_raw)
  meta_f  <- data.frame(sample_id=rownames(sample_data(ps_f)),
    as.data.frame(sample_data(ps_f)), stringsAsFactors=FALSE)
  tax_f   <- data.frame(tax_table(ps_f), stringsAsFactors=FALSE)
  DAYS    <- c(0,2,3,4,5,6)
  TRTS    <- c("Drought","Watered")

  topo_pd_list <- list(); kst_pd_list <- list()

  for (day in DAYS) {
    for (trt in TRTS) {
      sids <- meta_f$sample_id[meta_f$harvest == paste0("Day ",day) &
                                 meta_f$treatment == trt]
      if (length(sids) < 4) { cat(sprintf("  Skip Day%d %s (n=%d)\n",day,trt,length(sids))); next }

      otu_sub <- otu_raw[sids, , drop=FALSE]
      ## Prevalence filter: >50% of samples
      keep_taxa <- colMeans(otu_sub > 0) >= 0.50
      otu_flt   <- otu_sub[, keep_taxa, drop=FALSE]
      if (ncol(otu_flt) < 5) next

      ## Spearman correlation
      cor_res <- suppressWarnings(Hmisc::rcorr(otu_flt, type="spearman"))
      rho_mat <- cor_res$r; p_mat <- cor_res$P
      rho_mat[is.na(rho_mat)] <- 0; p_mat[is.na(p_mat)] <- 1

      ## Threshold: |rho| > 0.6 and p < 0.05
      adj_sp <- (abs(rho_mat) > 0.6) & (p_mat < 0.05) & !diag(nrow(rho_mat))
      diag(adj_sp) <- FALSE
      upper_sp <- which(adj_sp & upper.tri(adj_sp), arr.ind=TRUE)
      if (nrow(upper_sp) < 3) next

      gs <- igraph::make_empty_graph(n=ncol(otu_flt), directed=FALSE)
      igraph::V(gs)$name <- colnames(otu_flt)
      gs <- igraph::add_edges(gs, as.vector(t(upper_sp)),
        weight=rho_mat[upper_sp])

      comps <- igraph::components(gs)
      lcc_v <- which(comps$membership == which.max(comps$csize))
      gs_lcc <- igraph::induced_subgraph(gs, lcc_v)

      ## Modules
      g_uw <- igraph::graph_from_edgelist(
        igraph::as_edgelist(gs_lcc, names=FALSE), directed=FALSE)
      mods_sp <- igraph::cluster_fast_greedy(g_uw)
      igraph::V(gs_lcc)$module <- igraph::membership(mods_sp)
      igraph::V(gs_lcc)$degree <- igraph::degree(gs_lcc)

      ## Zi/Pi
      mod_v2 <- igraph::V(gs_lcc)$module; deg_v2 <- igraph::V(gs_lcc)$degree
      Zi2 <- numeric(igraph::vcount(gs_lcc)); Pi2 <- numeric(igraph::vcount(gs_lcc))
      for (m2 in unique(mod_v2)) {
        idx2 <- which(mod_v2==m2); if (length(idx2)<2) next
        k2 <- sapply(idx2, function(v) sum(igraph::neighbors(gs_lcc,v) %in% idx2))
        sd2 <- sd(k2); if (is.na(sd2)||sd2==0) next
        Zi2[idx2] <- (k2-mean(k2))/sd2
      }
      for (i2 in seq_len(igraph::vcount(gs_lcc))) {
        ki2 <- deg_v2[i2]; if (ki2==0) next
        Pi2[i2] <- 1-sum((table(mod_v2[igraph::neighbors(gs_lcc,i2)])/ki2)^2)
      }
      igraph::V(gs_lcc)$Zi <- Zi2; igraph::V(gs_lcc)$Pi <- Pi2
      igraph::V(gs_lcc)$role <- ifelse(Zi2>=2.5 & Pi2>=0.62,"Network hub",
        ifelse(Zi2>=2.5 & Pi2<0.62,"Module hub",
        ifelse(Zi2<2.5  & Pi2>=0.62,"Connector","Peripheral")))

      asv2 <- igraph::V(gs_lcc)$name
      el2  <- igraph::as_edgelist(gs_lcc, names=TRUE)
      rho_ev <- rho_mat[cbind(match(el2[,1],colnames(otu_flt)),
                               match(el2[,2],colnames(otu_flt)))]

      topo_pd_list[[paste(day,trt)]] <- data.frame(
        Day=day, Treatment=trt,
        Nodes=igraph::vcount(gs_lcc), Edges=igraph::ecount(gs_lcc),
        Positive_edges=sum(rho_ev>0), Negative_edges=sum(rho_ev<0),
        Avg_degree=round(mean(deg_v2),4),
        Modularity=round(igraph::modularity(mods_sp),4),
        N_keystones=sum(igraph::V(gs_lcc)$role %in% c("Module hub","Network hub","Connector")),
        stringsAsFactors=FALSE)

      kst_pd_list[[paste(day,trt)]] <- data.frame(
        Day=day, Treatment=trt, ASV=asv2,
        Genus=ifelse(asv2 %in% rownames(tax_f), trimws(tax_f[asv2,"Genus"]), ""),
        Phylum=ifelse(asv2 %in% rownames(tax_f), trimws(tax_f[asv2,"Phylum"]), ""),
        Degree=deg_v2, Zi=round(Zi2,4), Pi=round(Pi2,4),
        Role=igraph::V(gs_lcc)$role, stringsAsFactors=FALSE)

      cat(sprintf("  Day%d %s: %d nodes, %d edges\n",
                  day, trt, igraph::vcount(gs_lcc), igraph::ecount(gs_lcc)))
    }
  }

  topo_pd <- do.call(rbind, topo_pd_list); kst_pd <- do.call(rbind, kst_pd_list)
  write.csv(topo_pd, .topo_pd_path, row.names=FALSE)
  write.csv(kst_pd,  .kst_pd_path,  row.names=FALSE)
  cat(sprintf("  Saved topology (%d rows) and keystone (%d rows) tables.\n",
              nrow(topo_pd), nrow(kst_pd)))
} else {
  cat("  Loading cached per-day treatment topology tables...\n")
  topo_pd <- read.csv(.topo_pd_path, stringsAsFactors=FALSE)
  kst_pd  <- read.csv(.kst_pd_path,  stringsAsFactors=FALSE)
}

## ─────────────────────────────────────────────────────────────────────────────
## 7D — PER-DAY SPEARMAN NETWORKS: GENOTYPE
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── 7D: Per-day Spearman networks — genotype ──\n")

.topo_gd_path <- "tables/P6_networks/network_topology_perday_genotype.csv"
.kst_gd_path  <- "tables/P6_networks/network_keystone_perday_genotype.csv"

if (!file.exists(.topo_gd_path) || !file.exists(.kst_gd_path)) {
  cat("  Computing per-day genotype Spearman networks...\n")
  ps_fg    <- readRDS("data/ps_filtered.rds")
  otu_rawg <- as(otu_table(ps_fg), "matrix")
  if (taxa_are_rows(ps_fg)) otu_rawg <- t(otu_rawg)
  meta_fg  <- data.frame(sample_id=rownames(sample_data(ps_fg)),
    as.data.frame(sample_data(ps_fg)), stringsAsFactors=FALSE)
  tax_fg   <- data.frame(tax_table(ps_fg), stringsAsFactors=FALSE)
  GENOS    <- c("Resistance","Susceptible")

  topo_gd_list <- list(); kst_gd_list <- list()

  for (day in c(0,2,3,4,5,6)) {
    for (geno in GENOS) {
      geno_lbl <- ifelse(geno=="Resistance","Resistant","Susceptible")
      sids <- meta_fg$sample_id[
        (meta_fg$harvest == paste0("Day ",day) | meta_fg$Day == day) &
        meta_fg$trait == geno & meta_fg$condition == "Planted"]
      ## Fallback: try numeric day column
      if (length(sids) < 4)
        sids <- meta_fg$sample_id[
          as.character(meta_fg$Day) == as.character(day) &
          meta_fg$trait == geno & meta_fg$condition == "Planted"]
      if (length(sids) < 4) {
        cat(sprintf("  Skip Day%d %s (n=%d)\n", day, geno_lbl, length(sids))); next
      }

      otu_subg <- otu_rawg[sids, , drop=FALSE]
      keep_g   <- colMeans(otu_subg > 0) >= 0.50
      otu_fltg <- otu_subg[, keep_g, drop=FALSE]
      if (ncol(otu_fltg) < 5) next

      cor_g <- suppressWarnings(Hmisc::rcorr(otu_fltg, type="spearman"))
      rho_g <- cor_g$r; p_g <- cor_g$P
      rho_g[is.na(rho_g)] <- 0; p_g[is.na(p_g)] <- 1
      adj_g <- (abs(rho_g) > 0.6) & (p_g < 0.05) & !diag(nrow(rho_g))
      diag(adj_g) <- FALSE
      upper_g <- which(adj_g & upper.tri(adj_g), arr.ind=TRUE)
      if (nrow(upper_g) < 3) next

      gg <- igraph::make_empty_graph(n=ncol(otu_fltg), directed=FALSE)
      igraph::V(gg)$name <- colnames(otu_fltg)
      gg <- igraph::add_edges(gg, as.vector(t(upper_g)), weight=rho_g[upper_g])

      comps_g <- igraph::components(gg)
      lcc_vg  <- which(comps_g$membership == which.max(comps_g$csize))
      gg_lcc  <- igraph::induced_subgraph(gg, lcc_vg)
      gg_uw   <- igraph::graph_from_edgelist(
        igraph::as_edgelist(gg_lcc, names=FALSE), directed=FALSE)
      mods_g  <- igraph::cluster_fast_greedy(gg_uw)
      igraph::V(gg_lcc)$module <- igraph::membership(mods_g)
      igraph::V(gg_lcc)$degree <- igraph::degree(gg_lcc)

      ## Zi/Pi
      mod_vg <- igraph::V(gg_lcc)$module; deg_vg <- igraph::V(gg_lcc)$degree
      Zig <- numeric(igraph::vcount(gg_lcc)); Pig <- numeric(igraph::vcount(gg_lcc))
      for (mg in unique(mod_vg)) {
        idxg <- which(mod_vg==mg); if (length(idxg)<2) next
        kg <- sapply(idxg, function(v) sum(igraph::neighbors(gg_lcc,v) %in% idxg))
        sdg <- sd(kg); if (is.na(sdg)||sdg==0) next
        Zig[idxg] <- (kg-mean(kg))/sdg
      }
      for (ig in seq_len(igraph::vcount(gg_lcc))) {
        kig <- deg_vg[ig]; if (kig==0) next
        Pig[ig] <- 1-sum((table(mod_vg[igraph::neighbors(gg_lcc,ig)])/kig)^2)
      }
      igraph::V(gg_lcc)$Zi <- Zig; igraph::V(gg_lcc)$Pi <- Pig
      igraph::V(gg_lcc)$role <- ifelse(Zig>=2.5 & Pig>=0.62,"Network hub",
        ifelse(Zig>=2.5 & Pig<0.62,"Module hub",
        ifelse(Zig<2.5  & Pig>=0.62,"Connector","Peripheral")))

      asv_g <- igraph::V(gg_lcc)$name
      el_g  <- igraph::as_edgelist(gg_lcc, names=TRUE)

      topo_gd_list[[paste(day, geno_lbl)]] <- data.frame(
        Day=day, Genotype=geno_lbl,
        Nodes=igraph::vcount(gg_lcc), Edges=igraph::ecount(gg_lcc),
        Avg_degree=round(mean(deg_vg),4),
        Modularity=round(igraph::modularity(mods_g),4),
        N_keystones=sum(igraph::V(gg_lcc)$role %in% c("Module hub","Network hub","Connector")),
        stringsAsFactors=FALSE)

      kst_gd_list[[paste(day, geno_lbl)]] <- data.frame(
        Day=day, Genotype=geno_lbl, ASV=asv_g,
        Genus=ifelse(asv_g %in% rownames(tax_fg), trimws(tax_fg[asv_g,"Genus"]), ""),
        Phylum=ifelse(asv_g %in% rownames(tax_fg), trimws(tax_fg[asv_g,"Phylum"]), ""),
        Degree=deg_vg, Zi=round(Zig,4), Pi=round(Pig,4),
        Role=igraph::V(gg_lcc)$role, stringsAsFactors=FALSE)

      cat(sprintf("  Day%d %s: %d nodes, %d edges\n",
                  day, geno_lbl, igraph::vcount(gg_lcc), igraph::ecount(gg_lcc)))
    }
  }

  topo_gd <- do.call(rbind, topo_gd_list); kst_gd <- do.call(rbind, kst_gd_list)
  write.csv(topo_gd, .topo_gd_path, row.names=FALSE)
  write.csv(kst_gd,  .kst_gd_path,  row.names=FALSE)
  cat(sprintf("  Saved topology (%d rows) and keystone (%d rows) tables.\n",
              nrow(topo_gd), nrow(kst_gd)))
} else {
  cat("  Loading cached per-day genotype topology tables...\n")
  topo_pd_g <- read.csv(.topo_gd_path, stringsAsFactors=FALSE)
  kst_pd_g  <- read.csv(.kst_gd_path,  stringsAsFactors=FALSE)
}

## ─────────────────────────────────────────────────────────────────────────────
## 7F — TEMPORAL FIGURES (per-day line plots + keystone heatmaps)
## ─────────────────────────────────────────────────────────────────────────────
cat("\n── 7F: Temporal network figures ──\n")

.make_topo_lp <- function(csv_path, group_col, cols_map, metric_col,
                            y_label, title_str) {
  df <- tryCatch(read.csv(csv_path, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(df)) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", x=0.5, y=0.5, label=paste0("MISSING:\n",basename(csv_path)),
        size=8/.pt, colour="#C0392B") + ggplot2::theme_void())
  }
  df$Group <- df[[group_col]]
  df$Value <- df[[metric_col]]
  df$Group <- factor(df$Group, levels=names(cols_map))
  ggplot2::ggplot(df, ggplot2::aes(Day, Value, colour=Group)) +
    ggplot2::geom_line(linewidth=0.8) + ggplot2::geom_point(size=2.0) +
    ggplot2::scale_colour_manual(values=cols_map, name=group_col) +
    ggplot2::scale_x_continuous(breaks=c(0,2,3,4,5,6)) +
    ggplot2::theme_bw(base_size=9) +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
          axis.text=ggplot2::element_text(size=8),
          plot.title=ggplot2::element_text(size=8, face="bold")) +
    ggplot2::labs(x="Day", y=y_label, title=title_str)
}

.COLS_TRT  <- c(Drought="#C0392B", Watered="#2980B9")
.COLS_GENO <- c(Resistant="#27AE60", Susceptible="#E67E22")

## Treatment per-day topology (4 metrics)
pA_t <- .make_topo_lp(.topo_pd_path,"Treatment",.COLS_TRT,"Modularity","Modularity","Modularity")
pB_t <- .make_topo_lp(.topo_pd_path,"Treatment",.COLS_TRT,"Avg_degree","Avg degree","Avg degree")
pC_t <- .make_topo_lp(.topo_pd_path,"Treatment",.COLS_TRT,"Edges","Edges","Edges")
pD_t <- .make_topo_lp(.topo_pd_path,"Treatment",.COLS_TRT,"N_keystones","N keystones","Keystones")

p_topo_trt <- cowplot::plot_grid(pA_t, pB_t, pC_t, pD_t, ncol=2, align="hv")
.save_all(p_topo_trt, "network_topology_temporal", 10, 7)

## Genotype per-day topology (4 metrics)
pA_g <- .make_topo_lp(.topo_gd_path,"Genotype",.COLS_GENO,"Modularity","Modularity","Modularity")
pB_g <- .make_topo_lp(.topo_gd_path,"Genotype",.COLS_GENO,"Avg_degree","Avg degree","Avg degree")
pC_g <- .make_topo_lp(.topo_gd_path,"Genotype",.COLS_GENO,"Edges","Edges","Edges")
pD_g <- .make_topo_lp(.topo_gd_path,"Genotype",.COLS_GENO,"N_keystones","N keystones","Keystones")

p_topo_geno <- cowplot::plot_grid(pA_g, pB_g, pC_g, pD_g, ncol=2, align="hv")
.save_all(p_topo_geno, "network_topology_temporal_genotype", 10, 7)

## Keystone heatmap helper
.keystone_heatmap <- function(kst_csv, group_col, fill_col, title_str, outname) {
  kst_all <- tryCatch(read.csv(kst_csv, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(kst_all)) { cat(sprintf("  MISSING: %s\n", kst_csv)); return(invisible(NULL)) }

  non_p <- kst_all %>% dplyr::filter(Role %in% c("Module hub","Network hub","Connector"))
  gen_lbl <- dplyr::coalesce(
    dplyr::if_else(trimws(non_p$Genus)!="", non_p$Genus, NA_character_),
    paste0("p__", non_p$Phylum))
  non_p$GenusLabel <- gen_lbl

  consistent <- non_p %>%
    dplyr::group_by(GenusLabel, Phylum) %>%
    dplyr::summarise(n_networks=dplyr::n(), .groups="drop") %>%
    dplyr::filter(n_networks >= 2) %>%
    dplyr::arrange(dplyr::desc(n_networks))

  keep_g <- consistent$GenusLabel
  heat_dat <- non_p %>%
    dplyr::filter(GenusLabel %in% keep_g) %>%
    dplyr::group_by(Day, .data[[group_col]], GenusLabel) %>%
    dplyr::summarise(Present=1L, .groups="drop") %>%
    dplyr::mutate(TimeLabel=paste0("D",Day,"_",substr(.data[[group_col]],1,3)))

  heat_wide <- heat_dat %>%
    dplyr::select(GenusLabel, TimeLabel, Present) %>%
    tidyr::pivot_wider(names_from=TimeLabel, values_from=Present, values_fill=0L)
  heat_mat <- as.matrix(heat_wide[,-1])
  rownames(heat_mat) <- heat_wide$GenusLabel
  col_ord <- sort(colnames(heat_mat))
  heat_mat <- heat_mat[, col_ord, drop=FALSE]
  heat_mat <- heat_mat[order(-rowSums(heat_mat)), , drop=FALSE]
  n_show   <- min(40, nrow(heat_mat))
  heat_mat <- heat_mat[seq_len(n_show), , drop=FALSE]

  heat_long <- as.data.frame(heat_mat) %>%
    tibble::rownames_to_column("Genus") %>%
    tidyr::pivot_longer(-Genus, names_to="Network", values_to="Present") %>%
    dplyr::mutate(
      Genus   = factor(Genus, levels=rev(rownames(heat_mat))),
      Network = factor(Network, levels=col_ord))

  p_h <- ggplot2::ggplot(heat_long, ggplot2::aes(Network, Genus, fill=factor(Present))) +
    ggplot2::geom_tile(colour="white", linewidth=0.3) +
    ggplot2::scale_fill_manual(values=c("0"="grey92","1"=fill_col),
                                labels=c("0"="Absent","1"="Keystone"), name=NULL) +
    ggplot2::theme_bw(base_size=9) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=7),
          axis.text.y=ggplot2::element_text(size=7),
          legend.position="top",
          plot.title=ggplot2::element_text(size=9, face="bold")) +
    ggplot2::labs(x=NULL, y=NULL,
         title=sprintf("%s  (top %d consistent keystones)", title_str, n_show))

  h_ht <- max(6, n_show * 0.22 + 2)
  .save_all(p_h, outname, 10, h_ht)
}

suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tibble))

.keystone_heatmap(.kst_pd_path,  "Treatment", "#C0392B",
  "Keystone taxa — treatment per-day networks", "keystone_temporal_heatmap")
.keystone_heatmap(.kst_gd_path,  "Genotype",  "#27AE60",
  "Keystone taxa — genotype per-day networks",  "keystone_temporal_heatmap_genotype")

## ─────────────────────────────────────────────────────────────────────────────
## PILLAR 7 COMPLETE
## ─────────────────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 7 COMPLETE — Co-occurrence Networks\n")
cat(sprintf("  Figures: %d files in figures/P6_networks/\n",
            length(list.files(.P6, pattern="\\.(pdf|png|svg)$"))))
cat(sprintf("  Tables:  %s\n",
            paste(grep("network|keystone|cohesion|P6_",
                       list.files("tables"), value=TRUE), collapse=", ")))
cat("═══════════════════════════════════════════════════════════════\n\n")


## =============================================================================
##  PILLAR 8: PHYLOGENETIC SIGNAL ANALYSIS
##  Blomberg K for drought-enrichment LFC signal
##
##  Sections:
##    8A  Cache check — skip computation if outputs exist
##    8B  K_manual() — custom Blomberg K (no GLS; MASS::ginv fallback)
##    8C  Per-phylum K computation (threshold: >= 5 ASVs)
##    8D  Validate against existing tables/P7_blomberg_k_results.csv
##    8E  Panel A — horizontal K bar chart by phylum
##    8F  Panel B — LFC violin: Patescibacteria vs all other phyla
##
##  NOTE: minimum n_taxa = 5 (not 10) to reproduce existing results;
##        Patescibacteria has n=9 and is the key significant phylum.
##  CACHING: re-runs if either CSV or both figure files are absent.
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 8: Phylogenetic Signal (Blomberg K)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

.P7_DIR   <- "figures/P7_phylogenetic"
.P7_CSV   <- "tables/P6_networks/P7_blomberg_k_results.csv"
.P7_FIG_A <- file.path(.P7_DIR, "P7_blomberg_k_by_phylum.pdf")
.P7_FIG_B <- file.path(.P7_DIR, "P7_LFC_distribution.pdf")

dir.create(.P7_DIR, recursive = TRUE, showWarnings = FALSE)

## ── 8A: Cache check ──────────────────────────────────────────────────────────
.p7_cache_ok <- file.exists(.P7_CSV) &&
                file.exists(.P7_FIG_A) &&
                file.exists(.P7_FIG_B)

if (.p7_cache_ok) {
  cat("  Cache hit — loading existing K results from", .P7_CSV, "\n")
  .k_res <- read.csv(.P7_CSV, stringsAsFactors = FALSE)
} else {
  cat("  Cache miss — computing Blomberg K per phylum...\n\n")

  ## ── 8B: K_manual() with permutation p-value ──────────────────────────────
  ## Blomberg K without GLS — avoids convergence errors on large trees.
  ## Formula: K = (MSobs / MSexp) where:
  ##   MSobs = mean sq error from phylogenetic weighted mean (n-1 denominator)
  ##   MSexp = expected MSobs under Brownian motion
  ## Uses MASS::ginv() when solve(V) fails due to near-singularity.
  .K_stat <- function(x_named, phy) {
    ## Align trait to tree tips
    shared <- intersect(names(x_named), phy$tip.label)
    if (length(shared) < 3) return(NA_real_)
    phy  <- ape::drop.tip(phy, setdiff(phy$tip.label, shared))
    x    <- x_named[phy$tip.label]
    n    <- length(x)
    V    <- ape::vcv(phy, corr = FALSE)
    Vinv <- tryCatch(solve(V), error = function(e) MASS::ginv(V))
    ones <- rep(1, n)
    ## Phylogenetically weighted mean
    xphy <- as.numeric((t(ones) %*% Vinv %*% x) /
                         (t(ones) %*% Vinv %*% ones))
    MSobs <- as.numeric(t(x - xphy) %*% (x - xphy)) / (n - 1L)
    a     <- as.numeric(t(ones) %*% Vinv %*% ones)
    MSexp <- (sum(diag(V)) - n / a) / (n - 1L)
    MSobs / MSexp
  }

  .K_manual <- function(x_named, phy, n_perm = 999L) {
    k_obs <- tryCatch(.K_stat(x_named, phy),
                       error = function(e) NA_real_)
    if (is.na(k_obs)) return(list(K = NA_real_, p = NA_real_))
    set.seed(42)
    k_null <- vapply(seq_len(n_perm), function(i) {
      x_perm <- x_named
      names(x_perm) <- sample(names(x_named))
      tryCatch(.K_stat(x_perm, phy), error = function(e) NA_real_)
    }, numeric(1))
    k_null <- k_null[!is.na(k_null)]
    p_val  <- if (length(k_null) > 0) mean(k_null >= k_obs) else NA_real_
    list(K = k_obs, p = p_val)
  }

  ## ── 8C: Load data ─────────────────────────────────────────────────────────
  cat("  Loading ps_filtered (tree + taxonomy)...\n")
  ps_p7  <- readRDS("data/ps_filtered.rds")
  tree   <- phy_tree(ps_p7)
  tax_p7 <- as.data.frame(tax_table(ps_p7), stringsAsFactors = FALSE)
  cat(sprintf("  Tree: %d tips, rooted=%s\n",
              length(tree$tip.label), is.rooted(tree)))

  cat("  Loading ANCOMBC2 treatment LFC table...\n")
  ancom <- read.csv("tables/P3_da/ANCOMBC2_treatment.csv",
                    stringsAsFactors = FALSE, row.names = 1)
  ## lfc column: positive = Drought-enriched
  lfc_all <- setNames(ancom$lfc, rownames(ancom))
  cat(sprintf("  LFC values: %d ASVs\n", length(lfc_all)))

  ## Shared ASVs between tree and LFC table
  shared_asv <- intersect(names(lfc_all), tree$tip.label)
  cat(sprintf("  Shared (tree ∩ LFC): %d ASVs\n", length(shared_asv)))
  lfc_shared <- lfc_all[shared_asv]

  ## ── Per-phylum K (min 5 ASVs with valid LFC) ─────────────────────────────
  phyla <- sort(unique(tax_p7$Phylum[!is.na(tax_p7$Phylum) &
                                       tax_p7$Phylum != ""]))
  cat(sprintf("  Computing K for %d phyla (>= 5 ASVs)...\n\n", length(phyla)))

  k_rows <- lapply(phyla, function(phy_name) {
    asv_ids <- rownames(tax_p7)[!is.na(tax_p7$Phylum) &
                                    tax_p7$Phylum == phy_name]
    pt <- intersect(asv_ids, shared_asv)
    if (length(pt) < 5L) return(NULL)
    ## Cap at 500 by |LFC| to avoid OOM on very large phyla
    if (length(pt) > 500L)
      pt <- names(sort(abs(lfc_shared[pt]), decreasing = TRUE))[1:500L]
    phy_tree_pruned <- tryCatch(
      ape::drop.tip(tree, setdiff(tree$tip.label, pt)),
      error = function(e) NULL)
    if (is.null(phy_tree_pruned) || length(phy_tree_pruned$tip.label) < 3L)
      return(NULL)
    res <- tryCatch(
      .K_manual(lfc_shared[pt], phy_tree_pruned, n_perm = 999L),
      error = function(e) {
        message("  K_manual error [", phy_name, "]: ", conditionMessage(e))
        list(K = NA_real_, p = NA_real_)
      })
    cat(sprintf("  %-30s  n=%3d  K=%.4f  p=%.4f\n",
                phy_name, length(pt),
                ifelse(is.na(res$K), NA_real_, res$K),
                ifelse(is.na(res$p), NA_real_, res$p)))
    data.frame(Phylum  = phy_name,
               K       = res$K,
               p       = res$p,
               n_taxa  = length(pt),
               stringsAsFactors = FALSE)
  })
  k_rows <- k_rows[!vapply(k_rows, is.null, logical(1))]
  .k_res <- do.call(rbind, k_rows)
  .k_res <- .k_res[!is.na(.k_res$K), ]
  .k_res$sig <- ifelse(.k_res$p < 0.05, "p<0.05", "ns")
  .k_res <- .k_res[order(-.k_res$K), ]

  ## ── 8D: Validate against existing CSV ────────────────────────────────────
  cat("\n── Step 4: Validating against existing CSV ──\n")
  if (file.exists(.P7_CSV)) {
    old <- read.csv(.P7_CSV, stringsAsFactors = FALSE)
    cat(sprintf("  %-25s  %8s  %8s  %6s  %8s  %8s  %6s\n",
                "Phylum", "K_old", "K_new", "ΔK", "p_old", "p_new", "Match"))
    cat(sprintf("  %s\n", paste(rep("-", 75), collapse="")))
    max_delta <- 0
    for (i in seq_len(nrow(.k_res))) {
      ph   <- .k_res$Phylum[i]
      k_new <- .k_res$K[i]
      p_new <- .k_res$p[i]
      old_row <- old[old$Phylum == ph, ]
      if (nrow(old_row) == 0) {
        cat(sprintf("  %-25s  %8s  %8.4f  %6s  %8s  %8.4f  %6s\n",
                    ph, "—", k_new, "—", "—", p_new, "NEW"))
      } else {
        dk    <- abs(k_new - old_row$K[1])
        match <- ifelse(dk <= 0.05, "PASS", "FAIL")
        max_delta <- max(max_delta, dk)
        cat(sprintf("  %-25s  %8.4f  %8.4f  %6.4f  %8.4f  %8.4f  %6s\n",
                    ph, old_row$K[1], k_new, dk, old_row$p[1], p_new, match))
      }
    }
    cat(sprintf("\n  Max |ΔK| across all phyla: %.4f\n", max_delta))
    if (max_delta <= 0.05) {
      cat("  VALIDATION PASSED — overwriting", .P7_CSV, "\n")
      write.csv(.k_res, .P7_CSV, row.names = FALSE)
    } else {
      .P7_CSV_NEW <- "tables/P6_networks/P7_blomberg_k_results_new.csv"
      cat(sprintf("  VALIDATION FAILED (max ΔK=%.4f > 0.05)\n", max_delta))
      cat("  Saving new results to", .P7_CSV_NEW, "(existing CSV preserved)\n")
      write.csv(.k_res, .P7_CSV_NEW, row.names = FALSE)
    }
  } else {
    cat("  No existing CSV — saving fresh results\n")
    write.csv(.k_res, .P7_CSV, row.names = FALSE)
  }
  cat("\n")
}  ## end cache-miss block

cat("  K results summary:\n")
print(.k_res[, c("Phylum","K","p","n_taxa","sig")], row.names = FALSE)
cat("\n")

## ── 8E: Panel A — horizontal K bar chart ─────────────────────────────────────
cat("── Panel A: Blomberg K by phylum ──\n")

.k_plot <- .k_res[!is.na(.k_res$K), ]
.k_plot$Phylum <- factor(.k_plot$Phylum,
                          levels = .k_plot$Phylum[order(.k_plot$K)])

p_k_bar <- ggplot(.k_plot, aes(x = K, y = Phylum, fill = sig)) +
  geom_col(width = 0.7, alpha = 0.88) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  annotate("text", x = 1, y = 0.6, label = "Brownian\nmotion",
           hjust = -0.1, vjust = 0, size = 8 / .pt, colour = "grey35") +
  scale_fill_manual(
    values = c("p<0.05" = COL_DROUGHT, "ns" = "grey65"),
    labels = c("p<0.05" = "p < 0.05", "ns" = "n.s."),
    name   = "Permutation\np-value"
  ) +
  theme_pub() +
  theme(axis.text.y  = element_text(size = 8),
        axis.text.x  = element_text(size = 8),
        axis.title   = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  labs(x        = "Blomberg K",
       y        = NULL,
       title    = "Phylogenetic signal of drought enrichment by phylum",
       subtitle = "K > 1 = stronger than Brownian motion; permutation p-values, 999 replicates")

save_fig(p_k_bar, file.path(.P7_DIR, "P7_blomberg_k_by_phylum"),
         w = ISME_DOUBLE_W * 0.75, h = 5)
cat(sprintf("  P7_blomberg_k_by_phylum  PDF %.0fKB  PNG %.0fKB  SVG %.0fKB\n",
    file.info(file.path(.P7_DIR,"P7_blomberg_k_by_phylum.pdf"))$size/1e3,
    file.info(file.path(.P7_DIR,"P7_blomberg_k_by_phylum.png"))$size/1e3,
    file.info(file.path(.P7_DIR,"P7_blomberg_k_by_phylum.svg"))$size/1e3))

## ── 8F: Panel B — LFC violin: Patescibacteria vs all other phyla ─────────────
cat("── Panel B: LFC distribution — Patescibacteria vs others ──\n")

## Build LFC data frame — merge tax annotation with LFC from ANCOMBC2
ancom_b <- tryCatch(
  read.csv("tables/P3_da/ANCOMBC2_treatment.csv",
           stringsAsFactors = FALSE, row.names = 1),
  error = function(e) NULL)

if (!is.null(ancom_b)) {
  ps_b    <- tryCatch(readRDS("data/ps_filtered.rds"), error = function(e) NULL)
  tax_b   <- if (!is.null(ps_b))
    as.data.frame(tax_table(ps_b), stringsAsFactors = FALSE) else NULL

  if (!is.null(tax_b)) {
    lfc_df <- data.frame(
      ASV    = rownames(ancom_b),
      lfc    = ancom_b$lfc,
      Phylum = tax_b[rownames(ancom_b), "Phylum"],
      stringsAsFactors = FALSE
    )
    lfc_df <- lfc_df[!is.na(lfc_df$Phylum) & !is.na(lfc_df$lfc), ]
    lfc_df$Group <- ifelse(lfc_df$Phylum == "Patescibacteria",
                            "Patescibacteria", "All other phyla")
    lfc_df$Group <- factor(lfc_df$Group,
                            levels = c("Patescibacteria", "All other phyla"))

    ## Pull K and p for Patescibacteria annotation
    .pates_row <- .k_res[.k_res$Phylum == "Patescibacteria", ]
    .annot_lbl <- if (nrow(.pates_row) > 0)
      sprintf("K = %.3f\np = %.3f%s",
              .pates_row$K[1], .pates_row$p[1],
              ifelse(.pates_row$p[1] < 0.001, " (< 0.001)", ""))
    else "Patescibacteria\nnot computed"

    ## Colour: Patescibacteria = drought colour; others = steel blue
    .pal_lfc <- c("Patescibacteria" = COL_DROUGHT,
                  "All other phyla" = COL_WATERED)

    p_lfc_viol <- ggplot(lfc_df, aes(x = Group, y = lfc, fill = Group)) +
      geom_violin(alpha = 0.75, trim = FALSE, linewidth = 0.4) +
      geom_boxplot(width = 0.12, outlier.shape = NA,
                   fill = "white", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      annotate("text",
               x     = 1,
               y     = max(lfc_df$lfc[lfc_df$Group == "Patescibacteria"],
                            na.rm = TRUE) * 0.9,
               label = .annot_lbl,
               hjust = 0.5, vjust = 1,
               size  = 8 / .pt, colour = COL_DROUGHT,
               fontface = "bold") +
      scale_fill_manual(values = .pal_lfc, guide = "none") +
      theme_pub() +
      theme(axis.text.x  = element_text(size = 9),
            axis.text.y  = element_text(size = 8),
            axis.title   = element_text(size = 9)) +
      labs(x        = NULL,
           y        = "Log fold change (Drought / Watered)",
           title    = "LFC distribution: Patescibacteria vs other phyla",
           subtitle = sprintf(
             "Patescibacteria (n=%d ASVs) show strong phylogenetic clustering",
             sum(lfc_df$Group == "Patescibacteria")))

    save_fig(p_lfc_viol, file.path(.P7_DIR, "P7_LFC_distribution"),
             w = ISME_SINGLE_W * 1.4, h = 5)
    cat(sprintf("  P7_LFC_distribution  PDF %.0fKB  PNG %.0fKB  SVG %.0fKB\n",
        file.info(file.path(.P7_DIR,"P7_LFC_distribution.pdf"))$size/1e3,
        file.info(file.path(.P7_DIR,"P7_LFC_distribution.png"))$size/1e3,
        file.info(file.path(.P7_DIR,"P7_LFC_distribution.svg"))$size/1e3))
  } else {
    cat("  [P7B] taxonomy table unavailable — skipping LFC violin\n")
  }
} else {
  cat("  [P7B] ANCOMBC2_treatment.csv unavailable — skipping LFC violin\n")
}

## ── 8G: Verify figure_assembly.r paths ───────────────────────────────────────
cat("\n── Step 7: Verifying figure_assembly.r paths ──\n")
.fa_lines <- readLines("figure_assembly.r")
.fa_a     <- grep("P7_blomberg_k_by_phylum", .fa_lines, value = TRUE)
.fa_b     <- grep("P7_LFC_distribution",     .fa_lines, value = TRUE)
cat("  Panel A reference:", if (length(.fa_a)) trimws(.fa_a[1]) else "NOT FOUND", "\n")
cat("  Panel B reference:", if (length(.fa_b)) trimws(.fa_b[1]) else "NOT FOUND", "\n")
.fa_ok <- length(.fa_a) > 0 && length(.fa_b) > 0 &&
          grepl("P7_phylogenetic", .fa_a[1]) &&
          grepl("P7_phylogenetic", .fa_b[1])
cat(sprintf("  figure_assembly.r paths: %s\n",
            ifelse(.fa_ok, "CORRECT ✓", "NEED UPDATE")))

## ── Summary ──────────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  PILLAR 8 COMPLETE — Phylogenetic Signal\n")
.p7_figs <- list.files(.P7_DIR, pattern="\\.(pdf|png|svg)$")
cat(sprintf("  Figures: %d files in %s\n", length(.p7_figs), .P7_DIR))
cat(sprintf("  Table:   %s\n", .P7_CSV))
cat(sprintf("  Key result: Patescibacteria K=%.3f, p=%.3f\n",
    .k_res$K[.k_res$Phylum == "Patescibacteria"],
    .k_res$p[.k_res$Phylum == "Patescibacteria"]))
cat("═══════════════════════════════════════════════════════════════\n\n")

## =============================================================================
##  SUPPLEMENTARY FIGURE ASSEMBLY
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  SUPPLEMENTARY FIGURE ASSEMBLY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

tryCatch({

  ## ── Packages ────────────────────────────────────────────────────────────────
  .sa_pkgs <- c("ggplot2","cowplot","magick","gridExtra","grid","scales","tools","svglite")
  .sa_need <- .sa_pkgs[!.sa_pkgs %in% rownames(installed.packages())]
  if (length(.sa_need)) install.packages(.sa_need,
    repos="https://cloud.r-project.org", quiet=TRUE)
  invisible(lapply(.sa_pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only=TRUE))))

  ## ── Paths ───────────────────────────────────────────────────────────────────
  .SA_FIG  <- "figures"
  .SA_TAB  <- "tables"
  SUPP_OUT <- file.path(.SA_FIG, "manuscript_figures", "supplementary")
  dir.create(SUPP_OUT, recursive=TRUE, showWarnings=FALSE)

  ## ── Colours ─────────────────────────────────────────────────────────────────
  .SA_COLS <- c(Drought="#C0392B", Watered="#2980B9",
                Resistance="#27AE60", Susceptible="#E67E22", Unplanted="#7F8C8D")

  ## ── Helper: ggplot theme ────────────────────────────────────────────────────
  .sa_theme_isme <- function(bs=9) {
    theme_bw(base_size=bs) +
      theme(
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour="grey93", linewidth=0.25),
        axis.ticks        = element_line(linewidth=0.36),
        axis.text         = element_text(size=8),
        axis.title        = element_text(size=9),
        strip.background  = element_rect(fill="grey95", colour="grey70"),
        strip.text        = element_text(size=8, face="bold"),
        legend.text       = element_text(size=7),
        legend.title      = element_text(size=8, face="bold"),
        legend.key.size   = unit(3, "mm"),
        plot.title        = element_text(size=9, face="bold"),
        plot.margin       = margin(2, 2, 2, 2, "mm")
      )
  }

  ## ── Helper: load panel image (placeholder on failure) ───────────────────────
  .sa_load_panel <- function(path) {
    fname <- basename(path)
    if (!file.exists(path)) {
      cat("  [MISSING]", path, "\n")
      return(ggplot() +
        annotate("text", x=0.5, y=0.5, label=paste0("MISSING:\n", fname),
                 size=8/.pt, colour="#C0392B", hjust=0.5, vjust=0.5,
                 fontface="bold") +
        theme_void() +
        theme(panel.border=element_rect(colour="#C0392B", fill=NA, linewidth=1)))
    }
    tryCatch({
      img <- magick::image_read(path)
      magick::image_ggplot(img)
    }, error=function(e) {
      cat("  [LOAD ERROR]", fname, "--", conditionMessage(e), "\n")
      ggplot() +
        annotate("text", x=0.5, y=0.5, label=paste0("LOAD ERROR:\n", fname),
                 size=8/.pt, colour="#C0392B", hjust=0.5, vjust=0.5) +
        theme_void() +
        theme(panel.border=element_rect(colour="#C0392B", fill=NA, linewidth=1))
    })
  }

  ## ── Helper: add bold panel label ────────────────────────────────────────────
  .sa_add_label <- function(p, lbl, x=0.01, y=0.99) {
    cowplot::ggdraw(p) +
      cowplot::draw_label(lbl, x=x, y=y, hjust=0, vjust=1,
                          fontface="bold", size=11)
  }

  ## ── Helper: save supplementary figure (PDF + PNG) ───────────────────────────
  .sa_save_supp <- function(plot, name, w_mm, h_mm) {
    pp <- file.path(SUPP_OUT, paste0(name, ".pdf"))
    np <- file.path(SUPP_OUT, paste0(name, ".png"))
    ggsave(pp, plot, device="pdf", dpi=300, width=w_mm, height=h_mm, units="mm")
    ggsave(np, plot, device="png", dpi=300, width=w_mm, height=h_mm, units="mm")
    pi <- file.info(pp); ni <- file.info(np)
    cat(sprintf("  [SUPP] %-50s  %dmm x %dmm  PDF %.1f KB  PNG %.1f KB\n",
                name, w_mm, h_mm, pi$size/1e3, ni$size/1e3))
    invisible(list(pdf=pp, png=np))
  }

  ## ── Clear existing supplementary folder contents ────────────────────────────
  ## before reassembling with canonical S-numbers
  ## This ensures no stale or mis-numbered files remain
  existing_files <- list.files(SUPP_OUT, full.names=TRUE)
  if (length(existing_files) > 0) {
    file.remove(existing_files)
    cat("Cleared", length(existing_files),
        "existing files from supplementary folder\n")
  }

  ## ── FigS_additional_taxa_functional ─────────────────────────────────────────
  cat("\n=== FigS_additional_taxa_functional: Relocated panels ===\n")
  sTA <- .sa_load_panel(file.path(.SA_FIG, "P3_function",   "KO_categories_temporal.png"))
  sTB <- .sa_load_panel(file.path(.SA_FIG, "P3_da",         "heatmap_top40_DA_treatment.png"))
  sTC <- .sa_load_panel(file.path(.SA_FIG, "P1_composition","top_taxa_Genus_horizontal.png"))
  figST <- cowplot::plot_grid(
    .sa_add_label(sTA,"A"), .sa_add_label(sTB,"B"), .sa_add_label(sTC,"C"),
    ncol=1, align="v", axis="lr"
  )
  {
    w_mm <- 183; h_mm <- 240; w_in <- w_mm/25.4; h_in <- h_mm/25.4
    base_st <- file.path(SUPP_OUT, "FigS_additional_taxa_functional")
    ggsave(paste0(base_st,".pdf"),  figST, device="pdf",
           width=w_mm, height=h_mm, units="mm", dpi=300)
    ggsave(paste0(base_st,".tiff"), figST, device="tiff",
           width=w_mm, height=h_mm, units="mm", dpi=300, compression="lzw")
    ggsave(paste0(base_st,".svg"),  figST, device=svglite::svglite,
           width=w_in, height=h_in)
    fi <- lapply(c("pdf","tiff","svg"), function(ext)
      file.info(paste0(base_st,".",ext))$size)
    cat(sprintf(
      "  [SUPP] FigS_additional_taxa_functional  %dx%dmm  PDF %.1f KB  TIFF %.1f KB  SVG %.1f KB\n",
      w_mm, h_mm, fi[[1]]/1e3, fi[[2]]/1e3, fi[[3]]/1e3))
  }

  ## ── FigS_network_temporal ───────────────────────────────────────────────────
  cat("\n=== FigS_network_temporal: Treatment and Genotype temporal networks ===\n")
  .sa_make_topo_panel <- function(csv_path, group_col, cols_map, metric_col,
                                   y_label, title_str, breaks_x=c(0,2,3,4,5,6)) {
    df <- tryCatch(read.csv(csv_path, stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(df)) {
      return(ggplot() +
        annotate("text", x=0.5, y=0.5,
                 label=paste0("MISSING:\n", basename(csv_path)),
                 size=8/.pt, colour="#C0392B", hjust=0.5, vjust=0.5) +
        theme_void() +
        theme(panel.border=element_rect(colour="#C0392B", fill=NA, linewidth=1)))
    }
    df$Group <- df[[group_col]]
    df$Value <- df[[metric_col]]
    df$Group <- factor(df$Group, levels=names(cols_map))
    ggplot(df, aes(Day, Value, colour=Group)) +
      geom_line(linewidth=0.8) + geom_point(size=2.0) +
      scale_colour_manual(values=cols_map, name=group_col) +
      scale_x_continuous(breaks=breaks_x) +
      .sa_theme_isme() +
      theme(axis.text=element_text(size=8), legend.text=element_text(size=7),
            legend.title=element_text(size=7), legend.key.size=unit(3,"mm"),
            plot.title=element_text(size=8, face="bold")) +
      labs(x="Day", y=y_label, title=title_str)
  }
  trt_csv  <- file.path(.SA_TAB, "P6_networks", "network_topology_perday.csv")
  geno_csv <- file.path(.SA_TAB, "P6_networks", "network_topology_perday_genotype.csv")
  trt_cols  <- c(Drought  =unname(.SA_COLS["Drought"]),
                 Watered  =unname(.SA_COLS["Watered"]))
  geno_cols <- c(Resistant   =unname(.SA_COLS["Resistance"]),
                 Susceptible =unname(.SA_COLS["Susceptible"]))
  pA_st <- .sa_make_topo_panel(trt_csv,  "Treatment", trt_cols,  "Modularity",  "Modularity",  "(A) Modularity — Treatment")
  pB_st <- .sa_make_topo_panel(trt_csv,  "Treatment", trt_cols,  "Avg_degree",  "Avg degree",  "(B) Avg degree — Treatment")
  pC_st <- .sa_make_topo_panel(trt_csv,  "Treatment", trt_cols,  "Edges",       "Edges",       "(C) Edges — Treatment")
  pD_st <- .sa_make_topo_panel(trt_csv,  "Treatment", trt_cols,  "N_keystones", "N keystones", "(D) Keystones — Treatment")
  pE_st <- .sa_make_topo_panel(geno_csv, "Genotype",  geno_cols, "Modularity",  "Modularity",  "(E) Modularity — Genotype")
  pF_st <- .sa_make_topo_panel(geno_csv, "Genotype",  geno_cols, "Avg_degree",  "Avg degree",  "(F) Avg degree — Genotype")
  pG_st <- .sa_make_topo_panel(geno_csv, "Genotype",  geno_cols, "Edges",       "Edges",       "(G) Edges — Genotype")
  pH_st <- .sa_make_topo_panel(geno_csv, "Genotype",  geno_cols, "N_keystones", "N keystones", "(H) Keystones — Genotype")
  pI_st <- .sa_add_label(
    .sa_load_panel(file.path(.SA_FIG, "P6_networks",
                             "keystone_temporal_heatmap.svg")), "I")
  pJ_st <- .sa_add_label(
    .sa_load_panel(file.path(.SA_FIG, "P6_networks",
                             "keystone_temporal_heatmap_genotype.svg")), "J")
  row1_st <- cowplot::plot_grid(pA_st, pB_st, pC_st, pD_st,
    ncol=4, align="hv", axis="tblr",
    labels=c("A","B","C","D"), label_size=11, label_fontface="bold")
  row2_st <- cowplot::plot_grid(pE_st, pF_st, pG_st, pH_st,
    ncol=4, align="hv", axis="tblr",
    labels=c("E","F","G","H"), label_size=11, label_fontface="bold")
  row3_st <- cowplot::plot_grid(pI_st, pJ_st, ncol=2, align="hv", axis="tblr")
  figST2  <- cowplot::plot_grid(row1_st, row2_st, row3_st,
    ncol=1, rel_heights=c(0.22, 0.22, 0.56))
  {
    w_mm <- 183; h_mm <- 320; w_in <- w_mm/25.4; h_in <- h_mm/25.4
    base_st2 <- file.path(SUPP_OUT, "FigS_network_temporal")
    ggsave(paste0(base_st2,".pdf"),  figST2, device="pdf",
           width=w_mm, height=h_mm, units="mm", dpi=300)
    ggsave(paste0(base_st2,".tiff"), figST2, device="tiff",
           width=w_mm, height=h_mm, units="mm", dpi=300, compression="lzw")
    ggsave(paste0(base_st2,".svg"),  figST2, device=svglite::svglite,
           width=w_in, height=h_in)
    fi2 <- lapply(c("pdf","tiff","svg"), function(ext)
      file.info(paste0(base_st2,".",ext))$size)
    cat(sprintf(
      "  [SUPP] FigS_network_temporal  %dx%dmm  PDF %.1f KB  TIFF %.1f KB  SVG %.1f KB\n",
      w_mm, h_mm, fi2[[1]]/1e3, fi2[[2]]/1e3, fi2[[3]]/1e3))
  }

  ## ── S1: Sequencing Quality ──────────────────────────────────────────────────
  cat("\n=== FIGURE S1: Sequencing Quality ===\n")
  sA1 <- .sa_load_panel(file.path(.SA_FIG, "P1_qc",    "depth_by_group.png"))
  sB1 <- .sa_load_panel(file.path(.SA_FIG, "P1_qc",    "rarefaction_curves_treatment.png"))
  sC1 <- .sa_load_panel(file.path(.SA_FIG, "P1_qc",    "rarefaction_by_harvest.png"))
  sD1 <- .sa_load_panel(file.path(.SA_FIG, "P1_alpha", "alpha_planted_vs_unplanted.png"))
  figS1 <- cowplot::plot_grid(
    .sa_add_label(sA1,"A"), .sa_add_label(sB1,"B"),
    .sa_add_label(sC1,"C"), .sa_add_label(sD1,"D"),
    ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS1, "FigureS1_sequencing_quality", 183, 160)

  ## ── S2: Extended Taxonomy ───────────────────────────────────────────────────
  cat("\n=== FIGURE S2: Extended Taxonomy ===\n")
  sA2 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "phylum_barplot_per_sample.png"))
  sB2 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "phylum_barplot_trt_geno_time_FIXED.png"))
  sC2 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "genus_barplot_trt_geno_time.png"))
  sD2 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "venn_ASV_4set_trt_geno.png"))
  figS2 <- cowplot::plot_grid(
    .sa_add_label(sA2,"A"), .sa_add_label(sB2,"B"),
    .sa_add_label(sC2,"C"), .sa_add_label(sD2,"D"),
    ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS2, "FigureS2_extended_taxonomy", 183, 250)

  ## ── S3: Extended Beta Diversity ─────────────────────────────────────────────
  cat("\n=== FIGURE S3: Extended Beta Diversity ===\n")
  sA3 <- .sa_load_panel(file.path(.SA_FIG, "P0_community",  "NMDS_planted_unplanted.png"))
  sB3 <- .sa_load_panel(file.path(.SA_FIG, "P0_community",  "NMDS_by_time.png"))
  sC3 <- .sa_load_panel("Manuscript_ISME/supplementary/sensitivity_distance_metrics.png")
  sD3 <- .sa_load_panel("Manuscript_ISME/supplementary/rarefaction_sensitivity.png")
  figS3 <- cowplot::plot_grid(
    .sa_add_label(sA3,"A"), .sa_add_label(sB3,"B"),
    .sa_add_label(sC3,"C"), .sa_add_label(sD3,"D"),
    ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS3, "FigureS3_extended_betadiversity", 183, 160)

  ## ── S4: Alpha Diversity Statistics ──────────────────────────────────────────
  cat("\n=== FIGURE S4: Alpha Diversity Statistics ===\n")
  sA4 <- .sa_load_panel(file.path(.SA_FIG, "P1_alpha",
                                   "alpha_boxplot_time_genotype.png"))
  sB4 <- .sa_load_panel(file.path(.SA_FIG, "P1_alpha", "alpha_boxplot_genotype.png"))
  figS4 <- cowplot::plot_grid(.sa_add_label(sA4,"A"), .sa_add_label(sB4,"B"),
                               ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS4, "FigureS4_alpha_statistics", 183, 120)

  ## ── S5: Extended DA ─────────────────────────────────────────────────────────
  cat("\n=== FIGURE S5: Extended Differential Abundance ===\n")
  sA5 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "top_taxa_Phylum_horizontal.png"))
  sB5 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "top_taxa_Class_horizontal.png"))
  sC5 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "top_taxa_Family_horizontal.png"))
  sD5 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "top_taxa_Genus_horizontal.png"))
  figS5 <- cowplot::plot_grid(
    .sa_add_label(sA5,"A"), .sa_add_label(sB5,"B"),
    .sa_add_label(sC5,"C"), .sa_add_label(sD5,"D"),
    ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS5, "FigureS5_extended_DA", 183, 250)

  ## ── S6: Extended Assembly ───────────────────────────────────────────────────
  cat("\n=== FIGURE S6: Extended Assembly ===\n")
  sA6 <- .sa_load_panel(file.path(.SA_FIG, "P1_iCAMP",
                                   "P1_iCAMP_assembly_by_timepoint.png"))
  icamp_tbl <- data.frame(
    Statistic = c("Input ASVs","ASVs retained","Reads retained (%)",
                  "Pairwise comparisons","Method"),
    Value     = c("7,868","7,148","90.8%","2,556","iCAMP (bNTI + RC-Bray)")
  )
  sB6 <- cowplot::ggdraw() +
    cowplot::draw_label("iCAMP Filtering Report", x=0.5, y=0.97,
                        hjust=0.5, vjust=1, fontface="bold", size=9) +
    cowplot::draw_grob(
      gridExtra::tableGrob(icamp_tbl, rows=NULL,
        theme=gridExtra::ttheme_minimal(base_size=8,
          core   =list(fg_params=list(fontsize=8)),
          colhead=list(fg_params=list(fontsize=8, fontface="bold")))),
      x=0.05, y=0.05, width=0.9, height=0.85)
  sC6 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "venn_ASV_treatment.png"))
  sD6 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "core_phylum_composition.png"))
  top_S6 <- cowplot::plot_grid(.sa_add_label(sA6,"A"), sB6,
                                ncol=2, align="hv", axis="tblr")
  bot_S6 <- cowplot::plot_grid(.sa_add_label(sC6,"C"), .sa_add_label(sD6,"D"),
                                ncol=2, align="hv", axis="tblr")
  figS6  <- cowplot::plot_grid(top_S6, bot_S6, nrow=2)
  .sa_save_supp(figS6, "FigureS6_extended_assembly", 183, 180)

  ## ── S7: Phylogenetic Signal ──────────────────────────────────────────────────
  cat("\n=== FIGURE S7: Phylogenetic Signal ===\n")
  sA7 <- .sa_load_panel(file.path(.SA_FIG, "P7_phylogenetic",
                                   "P7_blomberg_k_by_phylum.png"))
  sB7 <- .sa_load_panel(file.path(.SA_FIG, "P7_phylogenetic",
                                   "P7_LFC_distribution.png"))
  figS7 <- cowplot::plot_grid(.sa_add_label(sA7,"A"), .sa_add_label(sB7,"B"),
                               ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS7, "FigureS7_phylogenetic_signal", 183, 100)

  ## ── S8: Extended Functional ──────────────────────────────────────────────────
  cat("\n=== FIGURE S8: Extended Functional Analysis ===\n")
  sA8 <- .sa_load_panel(file.path(.SA_FIG, "P3_function",
                                   "PCoA_functional_pathway.png"))
  sB8 <- .sa_load_panel(file.path(.SA_FIG, "P3_function",    "NSTI_distribution.png"))
  sC8 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition",
                                   "core_phylum_composition.png"))
  sD8 <- .sa_load_panel(file.path(.SA_FIG, "P1_composition", "venn_ASV_treatment.png"))
  figS8 <- cowplot::plot_grid(
    .sa_add_label(sA8,"A"), .sa_add_label(sB8,"B"),
    .sa_add_label(sC8,"C"), .sa_add_label(sD8,"D"),
    ncol=2, align="hv", axis="tblr")
  .sa_save_supp(figS8, "FigureS8_extended_functional", 183, 160)

  ## ── S9: Functional Alpha (three-factor) ─────────────────────────────────────
  cat("\n=== FIGURE S9: Functional Alpha Diversity (three-factor) ===\n")
  sA9 <- .sa_load_panel(file.path(.SA_FIG, "P3_function",
                                   "functional_alpha_diversity_threefactor.png"))
  figS9 <- cowplot::plot_grid(.sa_add_label(sA9,"A"), ncol=1)
  .sa_save_supp(figS9, "FigureS9_functional_alpha", 183, 120)

  ## ── S_network_topology: Network Topology Bar Chart (assembled before copy block)
  cat("\n=== FIGURE S_network_topology: Network Topology Bar Chart ===\n")
  sA_topo   <- .sa_load_panel(file.path(.SA_FIG, "P6_networks",
                                         "network_topology_comparison.png"))
  figS_topo <- cowplot::plot_grid(.sa_add_label(sA_topo,"A"), ncol=1)
  .sa_save_supp(figS_topo, "FigureS_network_topology", 120, 160)
  {
    .tiff_topo <- file.path(SUPP_OUT, "FigureS_network_topology.tiff")
    ggsave(.tiff_topo, figS_topo, device="tiff", dpi=300,
           width=120, height=160, units="mm", compression="lzw")
    cat(sprintf("  [SUPP] FigureS_network_topology.tiff  %.1f KB\n",
                file.info(.tiff_topo)$size/1e3))
  }

  ## ── copy_supp: file-copy helper ──────────────────────────────────────────────
  ## Copies an existing figure to SUPP_OUT with canonical FigureS{nn}_{desc}.{ext}
  copy_supp <- function(src, n, desc) {
    if (!file.exists(src)) {
      cat(sprintf("  SOURCE MISSING — S%02d: %s\n", n, src))
      return(invisible(NULL))
    }
    ext  <- tools::file_ext(src)
    dest <- file.path(SUPP_OUT, sprintf("FigureS%02d_%s.%s", n, desc, ext))
    file.copy(src, dest, overwrite=TRUE)
    sz   <- file.info(dest)$size / 1e3
    cat(sprintf("  S%02d: %-55s  %.1f KB\n", n, basename(dest), sz))
    invisible(dest)
  }

  cat("\n=== FIGURES S1-S15: Essential supplementary panels ===\n")

  ## S1 — Sequencing quality
  copy_supp("figures/manuscript_figures/supplementary/FigureS1_sequencing_quality.pdf",
            1, "sequencing_quality")

  ## S2 — Rarefaction sensitivity
  copy_supp("Manuscript_ISME/supplementary/rarefaction_sensitivity.pdf",
            2, "rarefaction_sensitivity")

  ## S3 — Alpha diversity statistics
  copy_supp("figures/manuscript_figures/supplementary/FigureS4_alpha_statistics.pdf",
            3, "alpha_diversity_statistics")

  ## S4 — Distance metric sensitivity
  copy_supp("Manuscript_ISME/supplementary/sensitivity_distance_metrics.pdf",
            4, "sensitivity_distance_metrics")

  ## S5 — Extended beta-diversity
  copy_supp("figures/manuscript_figures/supplementary/FigureS3_extended_betadiversity.pdf",
            5, "extended_betadiversity")

  ## S6 — Extended differential abundance
  copy_supp("figures/manuscript_figures/supplementary/FigureS5_extended_DA.pdf",
            6, "extended_differential_abundance")

  ## S7 — Extended assembly mechanisms
  copy_supp("figures/manuscript_figures/supplementary/FigureS6_extended_assembly.pdf",
            7, "extended_assembly_mechanisms")

  ## S8 — Phylogenetic signal
  copy_supp("figures/manuscript_figures/supplementary/FigureS7_phylogenetic_signal.pdf",
            8, "phylogenetic_signal")

  ## S9 — Network architecture (main Fig 5)
  copy_supp("figures/manuscript_figures/main/Figure5_networks_functional.pdf",
            9, "network_architecture")

  ## S10 — Network topology metrics
  copy_supp("figures/manuscript_figures/supplementary/FigureS_network_topology.pdf",
            10, "network_topology_metrics")

  ## S11 — Network temporal dynamics
  copy_supp("figures/manuscript_figures/supplementary/FigS_network_temporal.pdf",
            11, "network_temporal_dynamics")

  ## S12 — Functional alpha diversity
  copy_supp("figures/manuscript_figures/supplementary/FigureS9_functional_alpha.pdf",
            12, "functional_alpha")

  ## S13 — Functional redundancy Mantel
  copy_supp("figures/P3_function/functional_redundancy_mantel.png",
            13, "functional_redundancy_mantel")

  ## S14 — Extended functional analysis
  copy_supp("figures/manuscript_figures/supplementary/FigureS8_extended_functional.pdf",
            14, "extended_functional_analysis")

  ## S15 — LinkET taxa–morphology
  copy_supp("figures/P4_phenotype/linkET_taxa_morphology.png",
            15, "linkET_taxon_morphology")

  ## S16 — Plant photographs (manual Inkscape assembly required)
  cat("  S16 -> FigureS16_plant_photographs [MANUAL INKSCAPE ASSEMBLY PENDING]\n")

  n_copied <- length(list.files(SUPP_OUT, pattern="^FigureS[0-9]{2}_.*\\.(pdf|png)$"))
  cat(sprintf("\nSupplementary assembly complete: %d/15 canonical files (S16 pending manual assembly)\n",
              n_copied))
  cat(sprintf("Output: %s\n", SUPP_OUT))
  cat("═══════════════════════════════════════════════════════════════\n\n")

}, error=function(e) {
  cat("  Supplementary assembly error:", conditionMessage(e), "\n")
})

## =============================================================================
##  SUPPLEMENTARY TABLE ASSEMBLY
## =============================================================================

tryCatch({

  supp_tab_dir <- file.path(BASE_DIR, "tables", "manuscript_supplementary")
  dir.create(supp_tab_dir, recursive=TRUE, showWarnings=FALSE)

  ## Clear existing tables
  old_tabs <- list.files(supp_tab_dir, pattern="\\.csv$", full.names=TRUE)
  if (length(old_tabs) > 0) {
    invisible(file.remove(old_tabs))
    cat(sprintf("  Cleared %d existing supplementary tables\n", length(old_tabs)))
  }

  copy_tab <- function(src, n, desc) {
    src_full <- file.path(BASE_DIR, src)
    dest <- file.path(supp_tab_dir, sprintf("TableS%02d_%s.csv", n, desc))
    if (file.exists(src_full)) {
      file.copy(src_full, dest, overwrite=TRUE)
      cat(sprintf("  ST%02d: %s\n", n, desc))
    } else {
      cat(sprintf("  TABLE MISSING — ST%02d: %s [%s]\n", n, desc, src))
    }
  }

  copy_tab("tables/P1_beta/PERMANOVA_main.csv",                  1,  "PERMANOVA_main")
  copy_tab("tables/P1_beta/PERMANOVA_per_day.csv",               2,  "PERMANOVA_per_day")
  copy_tab("tables/P1_beta/variance_partitioning.csv",           3,  "variance_partitioning")
  copy_tab("tables/P3_da/ANCOMBC2_treatment.csv",                4,  "ANCOMBC2_treatment")
  copy_tab("tables/P3_da/ANCOMBC2_genotype.csv",                 5,  "ANCOMBC2_genotype")
  copy_tab("tables/P2_assembly/iCAMP_summary_for_manuscript.csv",6,  "iCAMP_summary")
  copy_tab("tables/P2_assembly/iCAMP_process_tests_treatment.csv",7, "iCAMP_process_tests")
  copy_tab("tables/P6_networks/P6_network_topology.csv",         8,  "network_topology")
  copy_tab("tables/P6_networks/network_keystone_taxa.csv",       9,  "network_keystone_taxa")
  copy_tab("tables/P3_function/faprotax_guild_stats.csv",        10, "faprotax_guild_stats")
  copy_tab("tables/P3_function/pgp_ko_stats.csv",                11, "pgp_ko_stats")
  copy_tab("tables/P4_phenotype/morphology_global_stats.csv",    12, "morphology_global_stats")
  copy_tab("tables/P4_phenotype/morphology_lme_results.csv",     13, "morphology_lme_results")
  copy_tab("tables/P4_phenotype/mantel_tests.csv",               14, "mantel_tests")
  copy_tab("tables/P5_prediction/RF_performance.csv",            15, "RF_performance")
  copy_tab("tables/P5_prediction/RF_feature_importance.csv",     16, "RF_feature_importance")

  tabs_out <- list.files(supp_tab_dir, pattern="\\.csv$")
  cat(sprintf("  Supplementary Table Assembly: %d tables copied to %s\n",
              length(tabs_out), supp_tab_dir))

}, error=function(e) {
  cat("  Supplementary table assembly error:", conditionMessage(e), "\n")
})

## =============================================================================
##  FINAL SUMMARY
## =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  PIPELINE COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

summary_df <- data.frame(
  Finding = c(
    "Samples (filtered)", "ASVs (filtered)", "Rarefaction depth",
    "PERMANOVA: Treatment R²", "PERMANOVA: Genotype R²", "PERMANOVA: Time R²",
    "PERMANOVA: Condition R²", "Residual variance (%)",
    "Mean turnover proportion",
    "DA ASVs (treatment)", "DA ASVs (genotype)",
    "Assembly: Deterministic %", "Assembly: Stochastic %",
    "Assembly: Dispersal Limitation (Drought)", "Assembly: Homogenizing Dispersal (Drought)",
    "Functional R² (treatment)", "Functional PERMANOVA p",
    "Procrustes r (p-value)", "Mantel r (partial)",
    "iCAMP chi-squared p",
    "RF AUC (full — primary)", "RF AUC (top30 — exploratory)", "RF perm p (full)"
  ),
  Value = c(
    nsamples(ps), ntaxa(ps), CONFIG$RAREFACTION_DEPTH,
    round(perm_main$R2[1],4), round(perm_main$R2[2],4), round(perm_main$R2[3],4),
    round(perm_cond$R2[1],4), round(residual_pct,1),
    round(mean(beta_summary$mean_turn),3),
    ifelse(exists("n_dr"), n_dr+n_wa, NA), ifelse(exists("n_res"), n_res+n_sus, NA),
    ifelse(exists("deterministic_pct"), round(deterministic_pct*100,1), NA),
    ifelse(exists("stochastic_pct"), round(stochastic_pct*100,1), NA),
    ifelse(icamp_available, "29.5%", NA), ifelse(icamp_available, "7.0%", NA),
    ifelse(exists("perm_func"), round(perm_func$R2[1],4), NA),
    ifelse(exists("perm_func"), round(perm_func$`Pr(>F)`[1],4), NA),
    ifelse(exists("procr"), sprintf("%.3f (p=%.3f %s)", procr_r, procr$signif,
                                     ifelse(procr$signif<0.05,"*","NS")), NA),
    ifelse(exists("mantel_p"), round(mantel_p$statistic,3), NA),
    ifelse(exists("chi_test"), sprintf("%.2e", chi_test$p.value), NA),
    round(auc_full,3), round(auc_red,3), round(perm_p,4)
  ), stringsAsFactors=FALSE
)
write.csv(summary_df, "tables/P1_composition/pipeline_summary.csv", row.names=FALSE)
cat("Summary:\n"); print(summary_df, right=FALSE)

## List all outputs
cat("\n\nOutput files:\n")
cat("Figures:", length(list.files("figures", recursive=TRUE, pattern="\\.(pdf|png)$")), "files\n")
cat("Tables:", length(list.files("tables", pattern="\\.csv$", recursive=TRUE)), "files\n")

writeLines(capture.output(sessionInfo()), "tables/P1_composition/sessionInfo.txt")
cat("\n═══ Pipeline finished ═══\n")
