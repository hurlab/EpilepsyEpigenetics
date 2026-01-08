################################################################################
# RNA-Seq Analysis: Epileptogenesis Models (Kindling & Kainic Acid)
################################################################################
# This script performs differential gene expression analysis comparing two rat
# models of epileptogenesis: intrahippocampal electrical kindling (KI) and
# systemic kainic acid-induced status epilepticus (KA).
#
# Analysis workflow:
#   1. Load and prepare RNA-Seq data (FPKM values)
#   2. Perform PCA visualization
#   3. Differential expression analysis with DESeq2
#   4. Cross-model comparisons (pairwise)
#   5. Overlap analysis (Venn diagrams, alluvial plots)
#   6. Functional enrichment (GO/KEGG)
#   7. Cell-type deconvolution (CIBERSORT)
#   8. Permutation testing for overlap significance
#
# Author: Junguk Hur
# Date: December 2025
################################################################################

################################################################################
# 1. Initialize Workspace
################################################################################

# Set working directory (adjust as needed)
# If using RStudio, this will automatically set to script location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  if (nchar(script_path) > 0) {
    setwd(dirname(script_path))
  }
}

# Print working directory for confirmation
cat("Working directory:", getwd(), "\n")

################################################################################
# 2. Load Required Packages
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  # Core analysis
  Rsubread, DESeq2, edgeR,

  # Genomic annotations
  GenomicFeatures, biomaRt, org.Rn.eg.db,

  # Visualization
  ggplot2, ggpubr, pheatmap, RColorBrewer,
  EnhancedVolcano, ggalluvial, ggtern,

  # Data manipulation
  dplyr, tidyr, data.table,

  # Functional enrichment
  richR, VennDetail, UpSetR, ComplexHeatmap,

  # Cell-type deconvolution
  IOBR,

  # Multivariate analysis
  factoextra, NMF,

  # Parallel processing
  parallel
)

################################################################################
# 3. Define Global Variables and Directories
################################################################################

# Analysis parameters
analysis_date <- '250502'
padj_threshold <- 0.05
thread <- min(32, parallel::detectCores() - 1)

# Directory structure
output_dir <- 'output_RNASeq/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
sample_file <- '../GEOData/samples_nooutlier.txt'
count_file <- '../GEOData/subread_counts.csv'  # Raw count matrix from featureCounts
fpkm_file <- '../GEOData/fpkmDatawithMean_220205.txt'  # Alternative: FPKM values
bam_dir <- 'bam_dir/'  # Directory containing aligned BAM files (if regenerating counts)
genes_gtf_file <- './rn6.refGene.gtf'  # For featureCounts if regenerating counts

# Check for GTF file in alternative locations
if (!file.exists(genes_gtf_file)) {
  genes_gtf_file <- '/data/genomes/rn6/rn6.refGene.gtf'
}

# Comparison file defining experimental contrasts
comparison_file <- 'comparison_250502.txt'

################################################################################
# 4. Build or Load GO and KEGG Annotation Objects
################################################################################

cat("\n=== Building/Loading GO and KEGG annotations ===\n")

annot_file <- paste0('GO-KEGG-annotation_', analysis_date, '.RData')

if (file.exists(annot_file)) {
  cat("Loading existing annotation file:", annot_file, "\n")
  load(annot_file)
} else {
  cat("Building new annotation objects...\n")
  # Gene Ontology annotation
  go_anno <- richR::buildAnnot(
    species = "rat",
    keytype = "SYMBOL",
    anntype = "GO"
  )

  # KEGG pathway annotation (from KEGG website, not builtin)
  ko_anno <- richR::buildAnnot(
    species = "rat",
    keytype = "SYMBOL",
    anntype = "KEGG",
    builtin = FALSE
  )

  # Save for future use
  save(go_anno, ko_anno, file = annot_file)
  cat("Annotation objects saved to:", annot_file, "\n")
}

################################################################################
# 5. Load Sample Information and Define Comparisons
################################################################################

cat("\n=== Loading sample metadata ===\n")

# Load sample information
samples <- read.delim(sample_file, sep = "\t", header = TRUE)
rownames(samples) <- samples$sample

# Remove unwanted groups if present
samples <- samples[!grepl("5ADs|1stlimbic", samples$group), ]

# Display sample summary
cat("Sample groups:\n")
print(table(samples$group))

# Load comparison definitions
comparisons <- read.delim(comparison_file, sep = "\t", header = TRUE)
cat("\nComparisons to be performed:\n")
print(comparisons)

################################################################################
# 6. Load Count Data
################################################################################

cat("\n=== Loading count data ===\n")

# Load raw count matrix from featureCounts output
if (file.exists(count_file)) {
  cat("Loading raw counts from:", count_file, "\n")
  count_data <- read.csv(count_file, header = TRUE, row.names = 1)

  # Check if sample names need to be formatted (e.g., add 'S' prefix if needed)
  # Assuming columns are already named S1, S13, S14, etc.
  cat("Count data dimensions:", nrow(count_data), "genes ×", ncol(count_data), "samples\n")

} else {
  cat("Count file not found. Loading FPKM data as alternative...\n")

  # Fallback to FPKM data if counts not available
  fpkm_data <- read.delim(fpkm_file, sep = "\t", header = TRUE, row.names = 1)
  fpkm_data <- fpkm_data[, !grepl("Mean", colnames(fpkm_data))]

  # Approximate counts from FPKM (not ideal, but workable)
  # Convert FPKM back to approximate counts (for DESeq2)
  count_data <- round(fpkm_data)  # Simple approximation

  cat("Using FPKM-derived counts (not ideal):", nrow(count_data), "genes ×", ncol(count_data), "samples\n")
}

# Ensure sample order matches metadata
# Match sample names (handle potential differences in naming conventions)
sample_cols <- colnames(count_data)
sample_match <- samples$sample %in% sample_cols

if (sum(sample_match) < nrow(samples)) {
  cat("Warning: Not all samples in metadata found in count data\n")
  cat("  Samples in count data:", paste(sample_cols, collapse = ", "), "\n")
  cat("  Samples in metadata:", paste(samples$sample, collapse = ", "), "\n")

  # Try to match by removing/adding 'S' prefix
  if (all(grepl("^S", sample_cols))) {
    samples$count_name <- samples$sample
  } else {
    samples$count_name <- paste0("S", samples$sample)
  }

  common_samples <- intersect(samples$count_name, sample_cols)
  count_data <- count_data[, common_samples]
  samples <- samples[samples$count_name %in% common_samples, ]
  cat("  Matched", length(common_samples), "samples\n")
} else {
  count_data <- count_data[, samples$sample]
}

cat("Final count matrix:", nrow(count_data), "genes ×", ncol(count_data), "samples\n")

################################################################################
# 7. Principal Component Analysis (PCA)
################################################################################

cat("\n=== Performing PCA analysis ===\n")

# PCA on log-transformed count values (add pseudocount to avoid log(0))
count_log <- log2(count_data + 1)

# Perform PCA
pca_result <- prcomp(t(count_log), scale. = TRUE)

# Create PCA data frame
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  sample = rownames(pca_result$x)
)

# Add group information
pca_data <- merge(pca_data, samples[, c("sample", "group")], by = "sample")

# Define model grouping
pca_data$model <- ifelse(grepl("3ClassV1|Sham$", pca_data$group), "KI", "KA")
pca_data$condition <- ifelse(grepl("ClassV", pca_data$group), "Seizure", "Control")

# Calculate variance explained
var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Plot PCA for each model separately
pdf(file.path(output_dir, "01_PCA_ByModel.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))

# Kindling model
ki_data <- pca_data[pca_data$model == "KI", ]
plot(ki_data$PC1, ki_data$PC2,
     col = ifelse(ki_data$condition == "Seizure", "red", "black"),
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", var_explained[1], "%)"),
     ylab = paste0("PC2 (", var_explained[2], "%)"),
     main = "Kindling Model PCA")
legend("topright", legend = c("Sham", "3ClassV"),
       col = c("black", "red"), pch = 19)

# Kainic acid model
ka_data <- pca_data[pca_data$model == "KA", ]
plot(ka_data$PC1, ka_data$PC2,
     col = ifelse(ka_data$condition == "Seizure", "green", "black"),
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", var_explained[1], "%)"),
     ylab = paste0("PC2 (", var_explained[2], "%)"),
     main = "Kainic Acid Model PCA")
legend("topright", legend = c("Sham", "3ClassV"),
       col = c("black", "green"), pch = 19)

dev.off()

# ggplot2 version with ellipses
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = model)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    title = "PCA of RNA-Seq Data"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(file.path(output_dir, "01_PCA_Combined.pdf"), p, width = 10, height = 6)

cat("PCA plots saved\n")

################################################################################
# 8. Differential Expression Analysis with DESeq2
################################################################################

cat("\n=== Running DESeq2 differential expression analysis ===\n")

# Use the loaded count data
# Note: count_data was loaded from subread_counts.csv in section 6

# OPTIONAL: Generate counts from BAM files if regenerating from scratch
# Uncomment this section if you want to regenerate counts from BAM files
# cat("Generating count matrix from BAM files...\n")
# bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
# fc_result <- Rsubread::featureCounts(
#   files = bam_files,
#   annot.ext = genes_gtf_file,
#   isGTFAnnotationFile = TRUE,
#   GTF.featureType = "exon",
#   GTF.attrType = "gene_id",
#   useMetaFeatures = TRUE,
#   allowMultiOverlap = FALSE,
#   countMultiMappingReads = FALSE,
#   nthreads = thread
# )
# count_data <- fc_result$counts
# colnames(count_data) <- gsub(".bam", "", basename(colnames(count_data)))

# Store all DESeq2 results
deseq_results <- list()
deg_tables <- list()

# Perform differential expression for each comparison
for (i in 1:nrow(comparisons)) {
  experiment <- as.character(comparisons[i, 1])
  reference <- as.character(comparisons[i, 2])
  comparison_name <- paste0(experiment, "_vs_", reference)

  cat("\nProcessing:", comparison_name, "\n")

  # Subset samples for this comparison
  comparison_samples <- samples[samples$group %in% c(experiment, reference), ]
  counts_subset <- count_data[, rownames(comparison_samples)]

  # Create DESeq2 dataset
  coldata <- data.frame(
    condition = factor(comparison_samples$group, levels = c(reference, experiment))
  )
  rownames(coldata) <- rownames(comparison_samples)

  # Filter low-count genes (at least 10 reads in at least 3 samples)
  keep <- rowSums(counts_subset >= 10) >= 3
  counts_subset <- counts_subset[keep, ]

  cat("  Genes after filtering:", nrow(counts_subset), "\n")

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = coldata,
    design = ~ condition
  )

  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)

  # Extract results
  res <- results(dds, contrast = c("condition", experiment, reference))
  res <- as.data.frame(res)
  res$SYMBOL <- rownames(res)

  # Filter for significant genes
  deg <- res[!is.na(res$padj) & res$padj < padj_threshold, ]
  deg <- deg[order(deg$padj), ]

  cat("  DEGs identified:", nrow(deg), "\n")

  # Store results
  deseq_results[[comparison_name]] <- res
  deg_tables[[comparison_name]] <- deg

  # Save to file
  write.table(
    res,
    file = file.path(output_dir, paste0(comparison_name, "_AllGenes_", analysis_date, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  write.table(
    deg,
    file = file.path(output_dir, paste0(comparison_name, "_DEtable_padj0.05_", analysis_date, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

################################################################################
# 9. DEG Summary and Overlap Analysis
################################################################################

cat("\n=== Analyzing DEG overlaps ===\n")

# Create simplified names for key comparisons
rename_map <- c(
  "KI_3ClassV_vs_KI_Sham" = "KI",
  "KA_3ClassV_vs_KA_ShamSaline" = "KA",
  "KI_3ClassV_vs_KA_3ClassV" = "KIKA",
  "KI_Sham_vs_KA_ShamSaline" = "Sham"
)

# Create gene lists for overlaps
deg_lists <- lapply(names(rename_map), function(comp_name) {
  if (comp_name %in% names(deg_tables)) {
    return(deg_tables[[comp_name]]$SYMBOL)
  } else {
    return(character(0))
  }
})
names(deg_lists) <- rename_map

# Print DEG counts
deg_summary <- data.frame(
  Comparison = names(deg_lists),
  DEG_Count = sapply(deg_lists, length)
)
print(deg_summary)

write.table(
  deg_summary,
  file = file.path(output_dir, "Summary_DEG_counts.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 10. Venn Diagrams
################################################################################

cat("\n=== Creating Venn diagrams ===\n")

# Four-way Venn diagram
library(VennDetail)

venn_four <- venndetail(deg_lists)

pdf(file.path(output_dir, "02_Venn_FourWay.pdf"), width = 10, height = 10)
plot(venn_four, type = "venn")
dev.off()

# Two-way Venn: KI vs KA (main models)
venn_two <- venndetail(list(
  KI = deg_lists$KI,
  KA = deg_lists$KA
))

pdf(file.path(output_dir, "02_Venn_KI_vs_KA.pdf"), width = 8, height = 8)
plot(venn_two, type = "venn")
dev.off()

# Extract shared genes between KI and KA
shared_deg_ki_ka <- intersect(deg_lists$KI, deg_lists$KA)
cat("Shared DEGs between KI and KA:", length(shared_deg_ki_ka), "\n")

# Get log2FC values for shared genes
shared_deg_data <- data.frame(
  SYMBOL = shared_deg_ki_ka,
  KI_log2FC = deseq_results$KI_3ClassV_vs_KI_Sham[shared_deg_ki_ka, "log2FoldChange"],
  KA_log2FC = deseq_results$KA_3ClassV_vs_KA_ShamSaline[shared_deg_ki_ka, "log2FoldChange"]
)

# Classify as concordant or discordant
shared_deg_data$direction <- ifelse(
  sign(shared_deg_data$KI_log2FC) == sign(shared_deg_data$KA_log2FC),
  "Concordant",
  "Discordant"
)

# Count concordant vs discordant
table(shared_deg_data$direction)

write.table(
  shared_deg_data,
  file = file.path(output_dir, "Shared_DEGs_KI_KA.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 11. Scatter Plot of Shared DEGs
################################################################################

cat("\n=== Creating scatter plot of shared DEGs ===\n")

p_scatter <- ggplot(shared_deg_data, aes(x = KI_log2FC, y = KA_log2FC, color = direction)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Concordant" = "blue", "Discordant" = "orange")) +
  labs(
    x = "Kindling log2 Fold Change",
    y = "Kainic Acid log2 Fold Change",
    title = "Shared DEGs between Kindling and Kainic Acid Models",
    color = "Direction"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "03_Scatter_Shared_DEGs.pdf"), p_scatter, width = 8, height = 7)

################################################################################
# 12. Alluvial Plot
################################################################################

cat("\n=== Creating alluvial plot ===\n")

# Get three-way comparison data
alluvial_data <- data.frame(
  SYMBOL = unique(c(deg_lists$KI, deg_lists$KA, deg_lists$KIKA))
)

# Add directionality for each comparison
for (comp in c("KI", "KA", "KIKA")) {
  full_name <- names(rename_map)[rename_map == comp]
  genes_in_comp <- deg_lists[[comp]]

  alluvial_data[[comp]] <- ifelse(
    alluvial_data$SYMBOL %in% genes_in_comp,
    ifelse(
      deseq_results[[full_name]][alluvial_data$SYMBOL, "log2FoldChange"] > 0,
      "Up",
      "Down"
    ),
    "None"
  )
}

# Remove genes not in any comparison
alluvial_data <- alluvial_data[
  rowSums(alluvial_data[, c("KI", "KA", "KIKA")] != "None") > 0,
]

# Count frequencies
alluvial_freq <- alluvial_data %>%
  group_by(KI, KA, KIKA) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  filter(Freq > 0)

# Create alluvial plot
p_alluvial <- ggplot(
  alluvial_freq,
  aes(axis1 = KI, axis2 = KA, axis3 = KIKA, y = Freq)
) +
  geom_alluvium(aes(fill = KI), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("KI", "KA", "KIKA"), expand = c(.05, .05)) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "None" = "gray")) +
  labs(
    title = "Directional Changes in DEGs Across Comparisons",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "04_Alluvial_DEG_Direction.pdf"), p_alluvial, width = 10, height = 8)

################################################################################
# 13. Functional Enrichment Analysis (GO and KEGG)
################################################################################

cat("\n=== Performing functional enrichment analysis ===\n")

# Containers for enrichment results
go_results <- list()
kegg_results <- list()

# Perform enrichment for each DEG set
for (comp in names(deg_lists)) {
  if (length(deg_lists[[comp]]) >= 5) {  # Minimum 5 genes for enrichment
    cat("  Enrichment for", comp, ":", length(deg_lists[[comp]]), "genes\n")

    # GO enrichment
    go_res <- richR::enrichment(
      x = deg_lists[[comp]],
      annotation = go_anno,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      ont = "ALL"  # Biological Process, Molecular Function, Cellular Component
    )

    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      go_results[[comp]] <- go_res

      # Save GO results
      write.table(
        go_res@result,
        file = file.path(output_dir, paste0("GO_enrichment_", comp, ".txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }

    # KEGG enrichment
    kegg_res <- richR::enrichment(
      x = deg_lists[[comp]],
      annotation = ko_anno,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )

    if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
      kegg_results[[comp]] <- kegg_res

      # Save KEGG results
      write.table(
        kegg_res@result,
        file = file.path(output_dir, paste0("KEGG_enrichment_", comp, ".txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
}

# Create comparison plots for GO enrichment
if (length(go_results) > 1) {
  go_compare <- richR::compareResult(go_results)

  # Dot plot comparison
  pdf(file.path(output_dir, "05_GO_Enrichment_Comparison.pdf"), width = 12, height = 10)
  print(richR::ggdot(go_compare, topN = 10))
  dev.off()

  # UpSet plot for overlapping terms
  go_sets <- lapply(go_results, function(x) {
    unique(x@result$Term)
  })

  pdf(file.path(output_dir, "05_GO_UpSet.pdf"), width = 10, height = 6)
  upset(fromList(go_sets), nsets = length(go_sets), order.by = "freq")
  dev.off()
}

# Create comparison plots for KEGG enrichment
if (length(kegg_results) > 1) {
  kegg_compare <- richR::compareResult(kegg_results)

  # Dot plot comparison
  pdf(file.path(output_dir, "05_KEGG_Enrichment_Comparison.pdf"), width = 12, height = 10)
  print(richR::ggdot(kegg_compare, topN = 10))
  dev.off()
}

################################################################################
# 14. Cell-Type Deconvolution (CIBERSORT)
################################################################################

cat("\n=== Performing cell-type deconvolution ===\n")

# Note: This requires a cell-type signature matrix
# Using human brain cell-type signatures mapped to rat

# Load or create signature matrix
# For this example, we'll use a pre-existing human brain signature
# and map to rat genes using biomaRt

# Map human genes to rat orthologs
cat("Mapping human genes to rat orthologs...\n")

# This is a placeholder - in practice, you would:
# 1. Load a human brain cell-type signature (e.g., from Darmanis et al.)
# 2. Use biomaRt to map human genes to rat orthologs
# 3. Run CIBERSORT using IOBR package

# Example workflow (uncomment and modify as needed):
# library(biomaRt)
# human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# rat_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
#
# # Get human-rat ortholog mapping
# orthologs <- getLDS(
#   attributes = "hgnc_symbol",
#   filters = "hgnc_symbol",
#   values = rownames(signature_matrix),
#   mart = human_mart,
#   attributesL = "rgd_symbol",
#   martL = rat_mart
# )
#
# # Run CIBERSORT
# cibersort_results <- IOBR::CIBERSORT(
#   sig_matrix = signature_matrix_rat,
#   mixture_file = fpkm_data,
#   perm = 100,
#   QN = TRUE
# )
#
# # Analyze cell-type differences
# cibersort_data <- as.data.frame(cibersort_results)
# cibersort_data$group <- samples[rownames(cibersort_data), "group"]
#
# # Statistical testing (Kruskal-Wallis + Dunn's test)
# library(FSA)
# for (cell_type in colnames(cibersort_data)[1:(ncol(cibersort_data)-1)]) {
#   kw_test <- kruskal.test(
#     as.formula(paste(cell_type, "~ group")),
#     data = cibersort_data
#   )
#
#   if (kw_test$p.value < 0.05) {
#     dunn_test <- dunnTest(
#       as.formula(paste(cell_type, "~ group")),
#       data = cibersort_data,
#       method = "bh"
#     )
#     print(dunn_test)
#   }
# }

cat("Cell-type deconvolution section requires cell-type signature matrix\n")
cat("See commented code for implementation details\n")

################################################################################
# 15. Permutation Testing for Overlap Significance
################################################################################

cat("\n=== Performing permutation testing for overlap significance ===\n")

# Function to perform permutation test
permutation_test_overlap <- function(set1, set2, universe, n_perm = 10000) {
  observed_overlap <- length(intersect(set1, set2))
  size1 <- length(set1)
  size2 <- length(set2)

  # Generate null distribution
  null_overlaps <- replicate(n_perm, {
    random_set1 <- sample(universe, size1)
    random_set2 <- sample(universe, size2)
    length(intersect(random_set1, random_set2))
  })

  # Calculate statistics
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

# Define universe (all genes tested)
universe <- rownames(count_data)

# Test KI vs KA overlap
perm_ki_ka <- permutation_test_overlap(
  deg_lists$KI,
  deg_lists$KA,
  universe,
  n_perm = 10000
)

cat("\nKI vs KA DEG Overlap Permutation Test:\n")
cat("  Observed overlap:", perm_ki_ka$observed, "\n")
cat("  Expected overlap:", round(perm_ki_ka$expected, 2), "\n")
cat("  Z-score:", round(perm_ki_ka$z_score, 2), "\n")
cat("  P-value:", perm_ki_ka$p_value, "\n")

# Save permutation test results
perm_results <- data.frame(
  Comparison = "KI_vs_KA",
  Observed = perm_ki_ka$observed,
  Expected = perm_ki_ka$expected,
  SD = perm_ki_ka$sd,
  Z_score = perm_ki_ka$z_score,
  P_value = perm_ki_ka$p_value
)

write.table(
  perm_results,
  file = file.path(output_dir, "Permutation_Test_Results.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 16. Save Workspace
################################################################################

cat("\n=== Saving workspace ===\n")

save.image(file = file.path(output_dir, "RNASeq_Analysis_Workspace.RData"))

cat("\n=== RNA-Seq analysis complete ===\n")
cat("Output directory:", output_dir, "\n")

################################################################################
# End of Script
################################################################################
