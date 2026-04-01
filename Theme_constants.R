## =============================================================================
## Global typography, colors and themes
## Source this file at the top of all analysis and assembly scripts
## =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
  library(showtext)
  library(RColorBrewer)
  library(svglite)
})

## ── Font setup ───────────────────────────────────────────────────────────────
## Load Arial via showtext for consistent rendering across platforms
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf",
         bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
         bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
showtext_auto()
showtext_opts(dpi = 300)
FONT_FAMILY <- "Arial"

## ──  print dimensions ────────────────────────────────────────────────────
## Single column = 89mm = 3.50 inches
## Double column = 183mm = 7.20 inches
## Maximum height = 220mm = 8.66 inches
SINGLE_W  <- 3.50
DOUBLE_W  <- 7.20
MAX_H     <- 8.66

## ── Base font size calibration ───────────────────────────────────────────────
## At 300 DPI, base_size=14 in R renders at approximately 8-9pt at  print size
## This is the minimum acceptable size for  figures
BASE_SIZE <- 16

## ── Global color grammar ─────────────────────────────────────────────────────
## These colors must be used consistently across ALL figures

## Primary treatment colors
COL_DROUGHT  <- "#B03A2E"   # Deep brick red
COL_WATERED  <- "#1A5276"   # Deep navy blue

## Genotype colors
COL_RESIST   <- "#1E8449"   # Deep forest green
COL_SUSCEPT  <- "#D35400"   # Deep burnt orange
COL_UNPLANT  <- "#626567"   # Mid grey

## Assembly process colors — single canonical source of truth
## Canonical names (user-specified order): Variable Selection, Homogeneous Selection,
## Dispersal Limitation, Homogenizing Dispersal, Drift
ASSEMBLY_COLORS <- c(
  "Variable Selection"     = "#E63946",   # Clear red
  "Homogeneous Selection"  = "#2E4057",   # Dark slate blue
  "Dispersal Limitation"   = "#F4A261",   # Warm amber
  "Homogenizing Dispersal" = "#457B9D",   # Steel blue
  "Drift"                  = "#A8DADC"    # Pale teal
)

## Clean-name alias (backward compat) — identical to ASSEMBLY_COLORS
COL_ASSEMBLY <- ASSEMBLY_COLORS

## Dot-notation alias for iCAMP output column names — derived from ASSEMBLY_COLORS
## so colours are guaranteed identical for the same process
COL_ASSEMBLY_DOT <- setNames(
  ASSEMBLY_COLORS[c("Variable Selection", "Homogeneous Selection",
                    "Dispersal Limitation", "Homogenizing Dispersal",
                    "Drift", "Drift")],
  c("Variable.Selection", "Homogeneous.Selection",
    "Dispersal.Limitation", "Homogenizing.Dispersal",
    "Drift.and.Others", "Undominated..Drift.")
)

## Phylum colors — 8 primary phyla + Other, anchored consistently
COL_PHYLUM <- c(
  "Pseudomonadota"     = "#E41A1C",   # Red
  "Actinomycetota"     = "#377EB8",   # Blue
  "Acidobacteriota"    = "#4DAF4A",   # Green
  "Gemmatimonadota"    = "#984EA3",   # Purple
  "Bacteroidota"       = "#FF7F00",   # Orange
  "Chloroflexota"      = "#A65628",   # Brown
  "Myxococcota"        = "#F781BF",   # Pink
  "Planctomycetota"    = "#999999",   # Grey
  "Bacillota"          = "#FFFF33",   # Yellow
  "Cyanobacteriota"    = "#66C2A5",   # Teal
  "Other"              = "#CCCCCC"    # Light grey
)

## Time sequential scale — single-hue purple, never blue or red (reserved for treatment)
COL_TIME <- setNames(
  colorRampPalette(c("#F2F0F7", "#54278F"))(6),
  paste("Day", c(0, 2, 3, 4, 5, 6))
)

## Heatmap diverging scale — orange-red to blue, NOT red-blue (avoids treatment confusion)
COL_HEATMAP_LOW  <- "#2166AC"   # Blue = low abundance
COL_HEATMAP_MID  <- "#FFFFFF"   # White = midpoint
COL_HEATMAP_HIGH <- "#D6604D"   # Orange-red = high abundance

## Significance colors
COL_SIG_UP   <- COL_DROUGHT    # Enriched in drought
COL_SIG_DOWN <- COL_WATERED    # Enriched in watered
COL_SIG_NS   <- "#CCCCCC"      # Not significant

## ── Build PAL list (drop-in replacement for existing PAL) ────────────────────
PAL <- list(
  treatment = c(Watered = COL_WATERED, Drought = COL_DROUGHT),
  genotype  = c(Resistance = COL_RESIST, Susceptible = COL_SUSCEPT,
                Unplanted = COL_UNPLANT),
  condition = c(Planted = COL_RESIST, Unplanted = COL_UNPLANT),
  assembly       = COL_ASSEMBLY_DOT,
  assembly_clean = ASSEMBLY_COLORS,
  phylum    = COL_PHYLUM,
  time      = COL_TIME
)

## ── Taxon name italics helper ─────────────────────────────────────────────────
## Wraps genus/species names in markdown italics for use with ggtext
## Usage: italicise_taxa(vector_of_names)
## Returns character vector with *Name* markdown formatting
## Non-genus names (phyla, classes etc.) are returned unchanged
italicise_taxa <- function(x, rank = "genus") {
  ## Only italicise genus level and below
  ## Phylum, Class, Order, Family — NOT italicised by convention
  if (rank %in% c("genus", "species", "asv")) {
    ## Skip hash IDs, NAs, and names with underscores that are not binomials
    is_name <- !is.na(x) &
               !grepl("^[a-f0-9]{8,}$", x, perl = TRUE) &  # not hash IDs
               !grepl("^[A-Z]{2}[0-9]", x)                  # not accession numbers
    x[is_name] <- paste0("*", x[is_name], "*")
  }
  x
}

## ── Publication theme ─────────────────────────────────────────────────────────
## All figures must use this theme
theme_pub <- function(base_size = BASE_SIZE, legend_pos = "bottom") {
  theme_bw(base_size = base_size) %+replace%
    theme(
      ## Typography — Arial throughout
      text              = element_text(family = "Arial", size = base_size),
      axis.title        = element_text(family = "Arial", size = base_size,
                                       colour = "black"),
      axis.text         = element_text(family = "Arial", size = max(8L, base_size - 2L),
                                       colour = "black"),
      strip.text        = element_text(family = "Arial", size = base_size - 1,
                                       face = "bold"),
      legend.title      = element_text(family = "Arial", size = base_size - 1,
                                       face = "bold"),
      legend.text       = element_text(family = "Arial", size = max(8L, base_size - 2L)),
      plot.title        = element_text(family = "Arial", size = base_size,
                                       face = "bold", hjust = 0),
      plot.subtitle     = element_text(family = "Arial", size = max(8L, base_size - 2L),
                                       colour = "grey40", hjust = 0),
      ## Grid and background
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      panel.border      = element_rect(colour = "grey70", fill = NA,
                                       linewidth = 0.5),
      ## Legend
      legend.position   = legend_pos,
      legend.key.size   = unit(0.4, "cm"),
      legend.background = element_blank(),
      ## Margins
      plot.margin       = margin(4, 4, 4, 4, "mm")
    )
}

## ── Save figure helper ────────────────────────────────────────────────────────
## Saves PDF + PNG at 300 DPI using -calibrated dimensions
## w and h are in inches (use SINGLE_W, DOUBLE_W constants)
save_fig <- function(p, stem, w = DOUBLE_W, h = 5) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h,
         device = cairo_pdf, units = "in")
  ggsave(paste0(stem, ".png"), p, width = w, height = h,
         dpi = 300, units = "in", type = "cairo")
  ggsave(paste0(stem, ".svg"), p, width = w, height = h,
         device = svglite::svglite, units = "in")
  invisible(p)
}

## ── Significance star helper ──────────────────────────────────────────────────
sig_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

cat("theme_constants.r loaded: Arial font, BASE_SIZE=", BASE_SIZE,
    ",  dimensions calibrated\n")

## End of theme_constants.r
