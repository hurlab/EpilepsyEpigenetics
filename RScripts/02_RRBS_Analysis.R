################################################################################
# RRBS (Methylation) Analysis: Epileptogenesis Models (Kindling & Kainic Acid)
################################################################################
# This script performs reduced representation bisulfite sequencing (RRBS)
# analysis to identify DNA methylation changes in two rat models of
# epileptogenesis: intrahippocampal electrical kindling (KI) and systemic
# kainic acid-induced status epilepticus (KA).
#
# Analysis workflow:
#   1. Load methylation beta values
#   2. Perform PCA visualization
#   3. Differential methylation analysis (DMC calling) with methylKit
#   4. Differentially methylated region (DMR) analysis with 1-kb tiling
#   5. Gene annotation with genomation and Ensembl
#   6. DMG (differentially methylated genes) aggregation
#   7. Genomic context analysis (promoter/exon/intron/intergenic)
#   8. Cross-model comparisons and overlap analysis
#   9. Functional enrichment (GO/KEGG)
#  10. Permutation testing for overlap significance
#
# Author: [Your Name]
# Date: December 2025
################################################################################

################################################################################
# 1. Initialize Workspace
################################################################################

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  if (nchar(script_path) > 0) {
    setwd(dirname(script_path))
  }
}

cat("Working directory:", getwd(), "\n")

################################################################################
# 2. Load Required Packages
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  # Methylation analysis
  methylKit, genomation,

  # Genomic annotations
  GenomicRanges, GenomicFeatures, biomaRt,
  TxDb.Rnorvegicus.UCSC.rn6.refGene,
  org.Rn.eg.db, AnnotationDbi,

  # Visualization
  ggplot2, ggpubr, pheatmap, RColorBrewer,
  gridExtra, grid,

  # Data manipulation
  dplyr, tidyr, data.table, stringr,

  # Functional enrichment
  richR, VennDetail, UpSetR, ComplexHeatmap,

  # Multivariate analysis
  factoextra,

  # Parallel processing
  parallel
)

################################################################################
# 3. Define Global Variables and Directories
################################################################################

# Analysis parameters
analysis_date <- '250502'
meth_diff_threshold <- 25  # Minimum methylation difference (%)
qvalue_threshold <- 0.01   # Q-value cutoff
padj_threshold <- 0.05     # Adjusted p-value for enrichment
thread <- min(60, parallel::detectCores() - 1)

# Directory structure
output_dir <- 'output_RRBS/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
sample_file <- '../GEOData/samples_nooutlier.txt'
beta_file <- '../GEOData/beta_values.txt.gz'
comparison_file <- 'comparison_250502.txt'

# Annotation files
gene_obj_file <- 'gene.obj.rdata'
ensembl_gene_file <- 'ensembl.gene.rdata'

################################################################################
# 4. Load or Create Gene Annotation Objects
################################################################################

cat("\n=== Loading gene annotation objects ===\n")

# Load or create gene.obj (TxDb-based gene annotations)
if (file.exists(gene_obj_file)) {
  cat("Loading existing gene.obj from", gene_obj_file, "\n")
  load(gene_obj_file)
} else {
  cat("Creating gene.obj from TxDb.Rnorvegicus.UCSC.rn6.refGene...\n")
  txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
  gene.obj <- genomation::readTranscriptFeatures(txdb)
  save(gene.obj, file = gene_obj_file)
  cat("Saved gene.obj to", gene_obj_file, "\n")
}

# Load or create ensembl.gene (gene name annotations from Ensembl/biomaRt)
if (file.exists(ensembl_gene_file)) {
  cat("Loading existing ensembl.gene from", ensembl_gene_file, "\n")
  load(ensembl_gene_file)
} else {
  cat("Creating ensembl.gene from biomaRt...\n")
  ensembl <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  ensembl.gene <- biomaRt::getBM(
    attributes = c('ensembl_transcript_id', 'ensembl_gene_id',
                   'external_gene_name', 'description'),
    mart = ensembl
  )
  save(ensembl.gene, file = ensembl_gene_file)
  cat("Saved ensembl.gene to", ensembl_gene_file, "\n")
}

cat("Gene annotation objects ready\n")

################################################################################
# 5. Build or Load GO and KEGG Annotation Objects
################################################################################

cat("\n=== Building/Loading GO and KEGG annotations ===\n")

annot_file <- paste0('GO-KEGG-annotation_', analysis_date, '.RData')

if (file.exists(annot_file)) {
  cat("Loading existing annotation file:", annot_file, "\n")
  load(annot_file)
} else {
  cat("Building new annotation objects...\n")
  go_anno <- richR::buildAnnot(species = "rat", keytype = "SYMBOL", anntype = "GO")
  ko_anno <- richR::buildAnnot(species = "rat", keytype = "SYMBOL", anntype = "KEGG", builtin = FALSE)
  save(go_anno, ko_anno, file = annot_file)
  cat("Annotation objects saved to:", annot_file, "\n")
}

################################################################################
# 6. Load Sample Information and Comparisons
################################################################################

cat("\n=== Loading sample metadata ===\n")

samples <- read.delim(sample_file, sep = "\t", header = TRUE)
samples <- samples[!grepl("5ADs|1stlimbic", samples$group), ]

cat("Sample groups:\n")
print(table(samples$group))

comparisons <- read.delim(comparison_file, sep = "\t", header = TRUE)
cat("\nComparisons to be performed:\n")
print(comparisons)

################################################################################
# 7. Load Beta Values
################################################################################

cat("\n=== Loading methylation beta values ===\n")

# Load beta values (pre-processed from Bismark)
beta_data <- read.delim(beta_file, sep = "\t", header = TRUE)

cat("Beta values dimensions:", nrow(beta_data), "CpGs ×", ncol(beta_data), "columns\n")

# Create CpG site identifier
if (!"cpg.site" %in% colnames(beta_data)) {
  beta_data$cpg.site <- paste0(beta_data$chr, "_", beta_data$start)
}

# Reorder columns: genomic info first, then sample beta values
genomic_cols <- c("cpg.site", "chr", "start", "end", "strand")
sample_cols <- setdiff(colnames(beta_data), genomic_cols)
beta_data <- beta_data[, c(genomic_cols, sample_cols)]

cat("Sample columns:", length(sample_cols), "\n")

################################################################################
# 8. Principal Component Analysis (PCA)
################################################################################

cat("\n=== Performing PCA on methylation data ===\n")

# Select top 10% most variable CpGs for PCA
beta_matrix <- as.matrix(beta_data[, sample_cols])
rownames(beta_matrix) <- beta_data$cpg.site

# Remove CpGs with missing values
beta_matrix <- beta_matrix[complete.cases(beta_matrix), ]

# Calculate variance for each CpG
cpg_var <- apply(beta_matrix, 1, var, na.rm = TRUE)
top_var_cpgs <- names(sort(cpg_var, decreasing = TRUE)[1:ceiling(0.1 * length(cpg_var))])

# Perform PCA on top variable CpGs
beta_pca <- beta_matrix[top_var_cpgs, ]
pca_result <- prcomp(t(beta_pca), scale. = TRUE)

# Create PCA data frame
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  sample = rownames(pca_result$x)
)

# Match with sample metadata
sample_match <- match(pca_data$sample, samples$sample)
pca_data$group <- samples$group[sample_match]
pca_data$model <- ifelse(grepl("3ClassV1|Sham$", pca_data$group), "KI", "KA")
pca_data$condition <- ifelse(grepl("ClassV", pca_data$group), "Seizure", "Control")

# Calculate variance explained
var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Plot PCA
pdf(file.path(output_dir, "01_PCA_Methylation.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))

# Kindling model
ki_data <- pca_data[pca_data$model == "KI", ]
plot(ki_data$PC1, ki_data$PC2,
     col = ifelse(ki_data$condition == "Seizure", "red", "black"),
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", var_explained[1], "%)"),
     ylab = paste0("PC2 (", var_explained[2], "%)"),
     main = "Kindling Model - Methylation PCA")
legend("topright", legend = c("Sham", "3ClassV"),
       col = c("black", "red"), pch = 19)

# Kainic acid model
ka_data <- pca_data[pca_data$model == "KA", ]
plot(ka_data$PC1, ka_data$PC2,
     col = ifelse(ka_data$condition == "Seizure", "green", "black"),
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", var_explained[1], "%)"),
     ylab = paste0("PC2 (", var_explained[2], "%)"),
     main = "Kainic Acid Model - Methylation PCA")
legend("topright", legend = c("Sham", "3ClassV"),
       col = c("black", "green"), pch = 19)

dev.off()

cat("PCA plot saved\n")

################################################################################
# 9. Load Pre-computed DMC Results
################################################################################

cat("\n=== Loading differential methylation results ===\n")

# Note: This script assumes DMCs have been called previously using methylKit
# For a complete pipeline from raw data, see the methylKit documentation
# Here we load pre-computed results

# Containers for DMC data
all_dmc <- list()
all_dmc_hyper <- list()
all_dmc_hypo <- list()

# Load DMC results for each comparison
dmc_dir <- "output_250805/"  # Directory with pre-computed DMC results

for (i in 1:nrow(comparisons)) {
  experiment <- as.character(comparisons[i, 1])
  reference <- as.character(comparisons[i, 2])
  comparison_name <- paste0(experiment, "_vs_", reference)

  dmc_file <- file.path(dmc_dir, paste0(comparison_name, "_DECpG-diff25q001.txt"))

  if (file.exists(dmc_file)) {
    cat("Loading DMCs for", comparison_name, "\n")
    dmc_data <- read.delim(dmc_file, sep = "\t", header = TRUE)

    all_dmc[[comparison_name]] <- dmc_data
    all_dmc_hyper[[comparison_name]] <- dmc_data[dmc_data$meth.diff > 0, ]
    all_dmc_hypo[[comparison_name]] <- dmc_data[dmc_data$meth.diff < 0, ]

    cat("  Total DMCs:", nrow(dmc_data),
        "(", nrow(all_dmc_hyper[[comparison_name]]), "hyper,",
        nrow(all_dmc_hypo[[comparison_name]]), "hypo)\n")
  } else {
    cat("Warning: DMC file not found:", dmc_file, "\n")
    cat("  To generate DMCs from scratch, use methylKit::calculateDiffMeth()\n")
  }
}

################################################################################
# 10. Annotate DMCs with Genes
################################################################################

cat("\n=== Annotating DMCs with genes ===\n")

# Function to annotate CpGs with genes
annotate_cpgs_with_genes <- function(cpg_data, gene.obj, ensembl.gene) {
  cat("  Annotating", nrow(cpg_data), "DMCpGs with genes...\n")

  # Create GRanges object from CpG coordinates
  cpg_gr <- GRanges(
    seqnames = cpg_data$chr,
    ranges = IRanges(start = cpg_data$start, end = cpg_data$end),
    strand = ifelse("strand" %in% colnames(cpg_data), cpg_data$strand, "*")
  )

  # Annotate with gene parts (promoter, exon, intron, intergenic)
  suppressWarnings({
    gene_annot <- genomation::annotateWithGeneParts(cpg_gr, gene.obj)
    tss_annot <- genomation::getAssociationWithTSS(gene_annot)
  })

  # Extract annotations
  annot_df <- as.data.frame(gene_annot@members)
  annot_df$cpg_index <- 1:nrow(annot_df)

  # Add genomic context
  cpg_data$prom <- as.integer(annot_df$prom)
  cpg_data$exon <- as.integer(annot_df$exon)
  cpg_data$intron <- as.integer(annot_df$intron)

  # Get nearest gene information
  tss_df <- as.data.frame(tss_annot)
  cpg_data$feature.name <- tss_df$feature.name
  cpg_data$dist.to.feature <- tss_df$dist.to.feature

  # Map NCBI gene IDs to gene symbols using Ensembl
  gene_id_to_symbol <- setNames(
    ensembl.gene$external_gene_name,
    ensembl.gene$ensembl_gene_id
  )

  # Try to map feature names to symbols
  cpg_data$external_gene_name <- gene_id_to_symbol[cpg_data$feature.name]

  # If no match, keep original feature name
  cpg_data$external_gene_name <- ifelse(
    is.na(cpg_data$external_gene_name),
    cpg_data$feature.name,
    cpg_data$external_gene_name
  )

  cat("  Annotated", sum(!is.na(cpg_data$external_gene_name)), "CpG-gene associations\n")

  return(cpg_data)
}

# Annotate all DMC sets
all_dmc_annotated <- list()

for (comp_name in names(all_dmc)) {
  cat("\nAnnotating", comp_name, "\n")
  all_dmc_annotated[[comp_name]] <- annotate_cpgs_with_genes(
    all_dmc[[comp_name]],
    gene.obj,
    ensembl.gene
  )

  # Save annotated DMCs
  write.table(
    all_dmc_annotated[[comp_name]],
    file = file.path(output_dir, paste0(comp_name, "_DMCs_annotated.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

################################################################################
# 11. Genomic Context Analysis
################################################################################

cat("\n=== Analyzing genomic context of DMCs ===\n")

# Function to create pie chart of genomic features
create_genomic_context_pie <- function(cpg_data, title) {
  # Count CpGs in each genomic context
  counts <- c(
    Promoter = sum(cpg_data$prom == 1, na.rm = TRUE),
    Exon = sum(cpg_data$exon == 1 & cpg_data$prom == 0, na.rm = TRUE),
    Intron = sum(cpg_data$intron == 1 & cpg_data$exon == 0 & cpg_data$prom == 0, na.rm = TRUE)
  )

  # Intergenic = rest
  counts["Intergenic"] <- nrow(cpg_data) - sum(counts)

  # Create data frame
  pie_data <- data.frame(
    Context = names(counts),
    Count = as.numeric(counts),
    Percentage = round(100 * as.numeric(counts) / nrow(cpg_data), 1)
  )

  return(pie_data)
}

# Create pie charts for each comparison
pdf(file.path(output_dir, "02_Genomic_Context_PieCharts.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2))

genomic_context_summary <- list()

for (comp_name in names(all_dmc_annotated)) {
  context_data <- create_genomic_context_pie(all_dmc_annotated[[comp_name]], comp_name)
  genomic_context_summary[[comp_name]] <- context_data

  # Create pie chart
  pie(context_data$Count,
      labels = paste0(context_data$Context, "\n",
                     context_data$Count, " (",
                     context_data$Percentage, "%)"),
      main = comp_name,
      col = c("red", "blue", "green", "gray"))
}

dev.off()

# Save genomic context summary
for (comp_name in names(genomic_context_summary)) {
  write.table(
    genomic_context_summary[[comp_name]],
    file = file.path(output_dir, paste0(comp_name, "_genomic_context.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

################################################################################
# 12. Aggregate DMCs to DMGs (Differentially Methylated Genes)
################################################################################

cat("\n=== Aggregating DMCs to DMGs ===\n")

# Function to create DMG list from annotated DMCs
create_dmg_list <- function(cpg_data) {
  # Filter for CpGs with valid gene annotations
  cpg_with_genes <- cpg_data[!is.na(cpg_data$external_gene_name) &
                              cpg_data$external_gene_name != "", ]

  # Get unique genes
  dmg_list <- unique(cpg_with_genes$external_gene_name)

  # Calculate average methylation difference per gene
  gene_meth_summary <- cpg_with_genes %>%
    group_by(external_gene_name) %>%
    summarise(
      mean_meth_diff = mean(meth.diff, na.rm = TRUE),
      n_cpgs = n(),
      min_qvalue = min(qvalue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(mean_meth_diff)))

  return(list(
    genes = dmg_list,
    summary = gene_meth_summary
  ))
}

# Create DMG lists
all_dmg <- list()
dmg_summaries <- list()

for (comp_name in names(all_dmc_annotated)) {
  dmg_result <- create_dmg_list(all_dmc_annotated[[comp_name]])
  all_dmg[[comp_name]] <- dmg_result$genes
  dmg_summaries[[comp_name]] <- dmg_result$summary

  cat(comp_name, ":", length(dmg_result$genes), "DMGs\n")

  # Save DMG summary
  write.table(
    dmg_result$summary,
    file = file.path(output_dir, paste0(comp_name, "_DMG_summary.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

################################################################################
# 13. DMR Analysis (1-kb Tiling)
################################################################################

cat("\n=== Loading DMR results ===\n")

# Load pre-computed DMR results
# DMRs are identified using methylKit's tileMethylCounts() with 1-kb windows

all_dmr <- list()
all_dmr_annotated <- list()

dmr_dir <- "output_250805/"  # Directory with pre-computed DMR results

for (i in 1:nrow(comparisons)) {
  experiment <- as.character(comparisons[i, 1])
  reference <- as.character(comparisons[i, 2])
  comparison_name <- paste0(experiment, "_vs_", reference)

  dmr_file <- file.path(dmr_dir, paste0(comparison_name, "_DMR-diff25q001.txt"))

  if (file.exists(dmr_file)) {
    cat("Loading DMRs for", comparison_name, "\n")
    dmr_data <- read.delim(dmr_file, sep = "\t", header = TRUE)

    all_dmr[[comparison_name]] <- dmr_data

    # Annotate DMRs with genes
    dmr_annotated <- annotate_cpgs_with_genes(dmr_data, gene.obj, ensembl.gene)
    all_dmr_annotated[[comparison_name]] <- dmr_annotated

    cat("  Total DMRs:", nrow(dmr_data), "\n")
    cat("  Unique genes:", length(unique(dmr_annotated$external_gene_name)), "\n")

    # Save annotated DMRs
    write.table(
      dmr_annotated,
      file = file.path(output_dir, paste0(comparison_name, "_DMRs_annotated.txt")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

################################################################################
# 14. DMG Overlap Analysis
################################################################################

cat("\n=== Analyzing DMG overlaps ===\n")

# Simplify comparison names
rename_map <- c(
  "KI_3ClassV_vs_KI_Sham" = "KI",
  "KA_3ClassV_vs_KA_ShamSaline" = "KA",
  "KI_3ClassV_vs_KA_3ClassV" = "KIKA",
  "KI_Sham_vs_KA_ShamSaline" = "Sham"
)

# Create simplified DMG lists
dmg_lists <- lapply(names(rename_map), function(comp_name) {
  if (comp_name %in% names(all_dmg)) {
    return(all_dmg[[comp_name]])
  } else {
    return(character(0))
  }
})
names(dmg_lists) <- rename_map

# Print DMG counts
dmg_summary <- data.frame(
  Comparison = names(dmg_lists),
  DMG_Count = sapply(dmg_lists, length)
)
print(dmg_summary)

write.table(
  dmg_summary,
  file = file.path(output_dir, "Summary_DMG_counts.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 15. Venn Diagrams for DMGs
################################################################################

cat("\n=== Creating Venn diagrams for DMGs ===\n")

library(VennDetail)

# Four-way Venn
venn_four_dmg <- venndetail(dmg_lists)

pdf(file.path(output_dir, "03_Venn_DMG_FourWay.pdf"), width = 10, height = 10)
plot(venn_four_dmg, type = "venn")
dev.off()

# Two-way Venn: KI vs KA
venn_two_dmg <- venndetail(list(
  KI = dmg_lists$KI,
  KA = dmg_lists$KA
))

pdf(file.path(output_dir, "03_Venn_DMG_KI_vs_KA.pdf"), width = 8, height = 8)
plot(venn_two_dmg, type = "venn")
dev.off()

# Extract shared DMGs
shared_dmg_ki_ka <- intersect(dmg_lists$KI, dmg_lists$KA)
cat("Shared DMGs between KI and KA:", length(shared_dmg_ki_ka), "\n")

# Get methylation difference values for shared genes
shared_dmg_data <- data.frame(
  SYMBOL = shared_dmg_ki_ka
)

# Add mean methylation differences from gene summaries
ki_summary <- dmg_summaries$KI_3ClassV_vs_KI_Sham
ka_summary <- dmg_summaries$KA_3ClassV_vs_KA_ShamSaline

shared_dmg_data <- shared_dmg_data %>%
  left_join(
    ki_summary %>% dplyr::select(external_gene_name, mean_meth_diff, n_cpgs),
    by = c("SYMBOL" = "external_gene_name")
  ) %>%
  dplyr::rename(KI_meth_diff = mean_meth_diff, KI_n_cpgs = n_cpgs) %>%
  left_join(
    ka_summary %>% dplyr::select(external_gene_name, mean_meth_diff, n_cpgs),
    by = c("SYMBOL" = "external_gene_name")
  ) %>%
  dplyr::rename(KA_meth_diff = mean_meth_diff, KA_n_cpgs = n_cpgs)

# Classify as concordant or discordant
shared_dmg_data$direction <- ifelse(
  sign(shared_dmg_data$KI_meth_diff) == sign(shared_dmg_data$KA_meth_diff),
  "Concordant",
  "Discordant"
)

cat("Concordant vs Discordant DMGs:\n")
print(table(shared_dmg_data$direction))

write.table(
  shared_dmg_data,
  file = file.path(output_dir, "Shared_DMGs_KI_KA.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 16. Scatter Plot of Shared DMGs
################################################################################

cat("\n=== Creating scatter plot of shared DMGs ===\n")

p_scatter_dmg <- ggplot(shared_dmg_data,
                        aes(x = KI_meth_diff, y = KA_meth_diff, color = direction)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Concordant" = "blue", "Discordant" = "orange")) +
  labs(
    x = "Kindling Mean Methylation Difference (%)",
    y = "Kainic Acid Mean Methylation Difference (%)",
    title = "Shared DMGs between Kindling and Kainic Acid Models",
    color = "Direction"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "04_Scatter_Shared_DMGs.pdf"), p_scatter_dmg,
       width = 8, height = 7)

################################################################################
# 17. Functional Enrichment Analysis (GO and KEGG)
################################################################################

cat("\n=== Performing functional enrichment analysis for DMGs ===\n")

go_results_dmg <- list()
kegg_results_dmg <- list()

for (comp in names(dmg_lists)) {
  if (length(dmg_lists[[comp]]) >= 5) {
    cat("  Enrichment for", comp, ":", length(dmg_lists[[comp]]), "genes\n")

    # GO enrichment
    go_res <- richR::enrichment(
      x = dmg_lists[[comp]],
      annotation = go_anno,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      ont = "ALL"
    )

    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      go_results_dmg[[comp]] <- go_res
      write.table(
        go_res@result,
        file = file.path(output_dir, paste0("GO_enrichment_DMG_", comp, ".txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }

    # KEGG enrichment
    kegg_res <- richR::enrichment(
      x = dmg_lists[[comp]],
      annotation = ko_anno,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )

    if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
      kegg_results_dmg[[comp]] <- kegg_res
      write.table(
        kegg_res@result,
        file = file.path(output_dir, paste0("KEGG_enrichment_DMG_", comp, ".txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
}

# Create comparison plots
if (length(go_results_dmg) > 1) {
  go_compare_dmg <- richR::compareResult(go_results_dmg)

  pdf(file.path(output_dir, "05_GO_Enrichment_DMG_Comparison.pdf"),
      width = 12, height = 10)
  print(richR::ggdot(go_compare_dmg, topN = 10))
  dev.off()

  # UpSet plot
  go_sets_dmg <- lapply(go_results_dmg, function(x) unique(x@result$Term))
  pdf(file.path(output_dir, "05_GO_DMG_UpSet.pdf"), width = 10, height = 6)
  upset(fromList(go_sets_dmg), nsets = length(go_sets_dmg), order.by = "freq")
  dev.off()
}

if (length(kegg_results_dmg) > 1) {
  kegg_compare_dmg <- richR::compareResult(kegg_results_dmg)

  pdf(file.path(output_dir, "05_KEGG_Enrichment_DMG_Comparison.pdf"),
      width = 12, height = 10)
  print(richR::ggdot(kegg_compare_dmg, topN = 10))
  dev.off()
}

################################################################################
# 18. Permutation Testing for DMG Overlap Significance
################################################################################

cat("\n=== Performing permutation testing for DMG overlap significance ===\n")

permutation_test_overlap <- function(set1, set2, universe, n_perm = 10000) {
  observed_overlap <- length(intersect(set1, set2))
  size1 <- length(set1)
  size2 <- length(set2)

  null_overlaps <- replicate(n_perm, {
    random_set1 <- sample(universe, size1)
    random_set2 <- sample(universe, size2)
    length(intersect(random_set1, random_set2))
  })

  mean_null <- mean(null_overlaps)
  sd_null <- sd(null_overlaps)
  z_score <- (observed_overlap - mean_null) / sd_null
  p_value <- sum(null_overlaps >= observed_overlap) / n_perm

  return(list(
    observed = observed_overlap,
    expected = mean_null,
    sd = sd_null,
    z_score = z_score,
    p_value = p_value
  ))
}

# Define universe (all genes with methylation measurements)
universe_genes <- unique(unlist(lapply(all_dmc_annotated, function(x) {
  unique(x$external_gene_name[!is.na(x$external_gene_name)])
})))

# Test KI vs KA DMG overlap
perm_dmg_ki_ka <- permutation_test_overlap(
  dmg_lists$KI,
  dmg_lists$KA,
  universe_genes,
  n_perm = 10000
)

cat("\nKI vs KA DMG Overlap Permutation Test:\n")
cat("  Observed overlap:", perm_dmg_ki_ka$observed, "\n")
cat("  Expected overlap:", round(perm_dmg_ki_ka$expected, 2), "\n")
cat("  Z-score:", round(perm_dmg_ki_ka$z_score, 2), "\n")
cat("  P-value:", perm_dmg_ki_ka$p_value, "\n")

perm_results_dmg <- data.frame(
  Comparison = "KI_vs_KA_DMG",
  Observed = perm_dmg_ki_ka$observed,
  Expected = perm_dmg_ki_ka$expected,
  SD = perm_dmg_ki_ka$sd,
  Z_score = perm_dmg_ki_ka$z_score,
  P_value = perm_dmg_ki_ka$p_value
)

write.table(
  perm_results_dmg,
  file = file.path(output_dir, "Permutation_Test_DMG_Results.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 19. Save Workspace
################################################################################

cat("\n=== Saving workspace ===\n")

save.image(file = file.path(output_dir, "RRBS_Analysis_Workspace.RData"))

cat("\n=== RRBS analysis complete ===\n")
cat("Output directory:", output_dir, "\n")

################################################################################
# End of Script
################################################################################
