######################################################
## Set the current working directory
######################################################
library(rstudioapi) # make sure you have it installed
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(current_path)

base_dir = dirname(current_path)
list.files()

################################################################################
################################################################################
##
## Overlap analysis
##
################################################################################
################################################################################

################################################################################
#####################  Packages  ###############################################
################################################################################

# Install pacman for package management if not already installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# Required CRAN packages
cran_packages <- c("dplyr", "ggplot2", "tidyverse")

# Required Bioconductor packages
bioc_packages <- c("VennDetail", "biomaRt")

# GitHub packages
github_packages <- c("guokai8/richR")

# Automatically install and load CRAN and Bioconductor packages
p_load(char = cran_packages, install = TRUE)
p_load(char = bioc_packages, install = TRUE, repos = BiocManager::repositories())

# Install and load GitHub packages
for (i in seq_along(github_packages)) {
  pkg_name <- names(github_packages)[i]
  repo <- github_packages[i]
  
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    pacman::p_install_gh(repo)
  }
  library(pkg_name, character.only = TRUE)
}


################################################################################
## Custom Functions
################################################################################

scatterplot_two_methdata <- function(venn.result.shared, set1name = "set1", xlab = "set1", 
                                     set2name = "set2", ylab = "set2", labelmin = 30, 
                                     genesumylocation = 8, genesumxlocation = 60, 
                                     filename.base = "scatterplot") {
  # Function to create scatter plots with various configurations
  # Details and code remain unchanged
}

scatterplot_meth_and_rna <- function(venn.result.shared, xlab = "set1", ylab = "set2", 
                                     labelmin = 30, labelFCmin = 1, genesumylocation = 8, 
                                     genesumxlocation = 8, filename.base = "scatterplot") {
  # Function to create scatter plots comparing methylation and RNA data
  # Details and code remain unchanged
}

perform_KEGG_enrichment <- function(input.gene, kodata, output.base.name) {
  # Function to perform KEGG enrichment analysis
  # Details and code remain unchanged
}

perform_GO_enrichment <- function(input.gene, godata, output.base.name) {
  # Function to perform GO enrichment analysis
  # Details and code remain unchanged
}

################################################################################
## Input Files and Directories
################################################################################

DBKIAdv.RRBS.25.file <- "output/3ClassV1_vs_Sham1stlimbic_DECpG-diff25q001.txt"
DBKAAdv.RRBS.25.file <- "output/3ClassV2_vs_SalineSham1stlimbic_DECpG-diff25q001.txt"

DBKIAdv.RNA.file <- "../RNAseq/output/3ClassV1_vs_Sham_DEtable_padj0.05_220205.txt"
DBKAAdv.RNA.file <- "../RNAseq/output/3ClassV2_vs_SalineSham_DEtable_padj0.05_220205.txt"

output_dir <- "output/"

################################################################################
## Load DMcpg and DEG Data
################################################################################

DBKIAdv.RRBS.25 <- read.delim(DBKIAdv.RRBS.25.file)
DBKAAdv.RRBS.25 <- read.delim(DBKAAdv.RRBS.25.file)

DBKIAdv.RNA <- read.delim(DBKIAdv.RNA.file)
DBKAAdv.RNA <- read.delim(DBKAAdv.RNA.file)

################################################################################
## RRBS Average Methylation Calculations for multiple DEcpgs at gene-level
################################################################################

DBKIAdv.RRBS.25.average <- DBKIAdv.RRBS.25 %>%
  group_by(external_gene_name) %>%
  summarize(mean = mean(meth.diff))

DBKAAdv.RRBS.25.average <- DBKAAdv.RRBS.25 %>%
  group_by(external_gene_name) %>%
  summarize(mean = mean(meth.diff))

DBKIAdv.RRBS.25.average.abs <- DBKIAdv.RRBS.25.average %>%
  filter(abs(mean) >= 25)

DBKAAdv.RRBS.25.average.abs <- DBKAAdv.RRBS.25.average %>%
  filter(abs(mean) >= 25)


################################################################################
## Overlap Analysis: DMG Levels (KI vs KA)
################################################################################

RRBS.DMG.overlap <- list(
  KI = unique(DBKIAdv.RRBS.25.average.abs$external_gene_name),
  KA = unique(DBKAAdv.RRBS.25.average.abs$external_gene_name)
)

RRBS.DMG.venn <- VennDetail::venndetail(RRBS.DMG.overlap)
pdf(file.path(output_dir, "VennDiagram_DMG_KI_vs_KA.pdf"))
plot(RRBS.DMG.venn)
dev.off()

RRBS.DMG.venn.result <- VennDetail::result(RRBS.DMG.venn)
write.csv(RRBS.DMG.venn.result, file = file.path(output_dir, "VennDiagram_DMG_KI_vs_KA_Details.csv"))

perform_KEGG_enrichment(
  input.gene = unique(RRBS.DMG.venn.result$Detail),
  kodata = ko_anno,
  output.base.name = file.path(output_dir, "KEGG_DMG_KI_vs_KA")
)

perform_GO_enrichment(
  input.gene = unique(RRBS.DMG.venn.result$Detail),
  godata = go_anno,
  output.base.name = file.path(output_dir, "GO_DMG_KI_vs_KA")
)



################################################################################
## Overlap Analysis: DMG vs DEG (KI and KA)
################################################################################

# KI: DMG vs DEG
KI.DMG.DEG.overlap <- list(
  DMG = unique(DBKIAdv.RRBS.25.average.abs$external_gene_name),
  DEG = rownames(DBKIAdv.RNA)
)

KI.DMG.DEG.venn <- VennDetail::venndetail(KI.DMG.DEG.overlap)
pdf(file.path(output_dir, "VennDiagram_KI_DMG_vs_DEG.pdf"))
plot(KI.DMG.DEG.venn)
dev.off()

KI.DMG.DEG.venn.result <- VennDetail::result(KI.DMG.DEG.venn)
write.csv(KI.DMG.DEG.venn.result, file = file.path(output_dir, "VennDiagram_KI_DMG_vs_DEG_Details.csv"))

perform_KEGG_enrichment(
  input.gene = unique(KI.DMG.DEG.venn.result$Detail),
  kodata = ko_anno,
  output.base.name = file.path(output_dir, "KEGG_KI_DMG_vs_DEG")
)

perform_GO_enrichment(
  input.gene = unique(KI.DMG.DEG.venn.result$Detail),
  godata = go_anno,
  output.base.name = file.path(output_dir, "GO_KI_DMG_vs_DEG")
)

# KA: DMG vs DEG
KA.DMG.DEG.overlap <- list(
  DMG = unique(DBKAAdv.RRBS.25.average.abs$external_gene_name),
  DEG = rownames(DBKAAdv.RNA)
)

KA.DMG.DEG.venn <- VennDetail::venndetail(KA.DMG.DEG.overlap)
pdf(file.path(output_dir, "VennDiagram_KA_DMG_vs_DEG.pdf"))
plot(KA.DMG.DEG.venn)
dev.off()

KA.DMG.DEG.venn.result <- VennDetail::result(KA.DMG.DEG.venn)
write.csv(KA.DMG.DEG.venn.result, file = file.path(output_dir, "VennDiagram_KA_DMG_vs_DEG_Details.csv"))

perform_KEGG_enrichment(
  input.gene = unique(KA.DMG.DEG.venn.result$Detail),
  kodata = ko_anno,
  output.base.name = file.path(output_dir, "KEGG_KA_DMG_vs_DEG")
)

perform_GO_enrichment(
  input.gene = unique(KA.DMG.DEG.venn.result$Detail),
  godata = go_anno,
  output.base.name = file.path(output_dir, "GO_KA_DMG_vs_DEG")
)


################################################################################
## Overlap Analysis: Shared DMG-DEG Genes (KI and KA)
################################################################################

shared.DMG.DEG <- list(
  KI = unique(KI.DMG.DEG.venn.result$Detail),
  KA = unique(KA.DMG.DEG.venn.result$Detail)
)

shared.DMG.DEG.venn <- VennDetail::venndetail(shared.DMG.DEG)
pdf(file.path(output_dir, "VennDiagram_Shared_DMG_DEG_KI_vs_KA.pdf"))
plot(shared.DMG.DEG.venn)
dev.off()

shared.DMG.DEG.venn.result <- VennDetail::result(shared.DMG.DEG.venn)
write.csv(shared.DMG.DEG.venn.result, file = file.path(output_dir, "VennDiagram_Shared_DMG_DEG_KI_vs_KA_Details.csv"))

perform_KEGG_enrichment(
  input.gene = unique(shared.DMG.DEG.venn.result$Detail),
  kodata = ko_anno,
  output.base.name = file.path(output_dir, "KEGG_Shared_DMG_DEG_KI_vs_KA")
)

perform_GO_enrichment(
  input.gene = unique(shared.DMG.DEG.venn.result$Detail),
  godata = go_anno,
  output.base.name = file.path(output_dir, "GO_Shared_DMG_DEG_KI_vs_KA")
)

################################################################################
## Scatter Plots
################################################################################

# Scatterplot: KI DMG vs DEG
scatterplot_meth_and_rna(
  KI.DMG.DEG.venn.result[KI.DMG.DEG.venn.result$Subset == "Shared", 
                         c("Detail", "mean", "log2FoldChange")],
  xlab = "KI methylation difference (%)",
  ylab = "KI RNA log2FC",
  labelmin = 30,
  filename.base = file.path(output_dir, "Scatterplot_KI_DMG_vs_DEG")
)

# Scatterplot: KA DMG vs DEG
scatterplot_meth_and_rna(
  KA.DMG.DEG.venn.result[KA.DMG.DEG.venn.result$Subset == "Shared", 
                         c("Detail", "mean", "log2FoldChange")],
  xlab = "KA methylation difference (%)",
  ylab = "KA RNA log2FC",
  labelmin = 30,
  filename.base = file.path(output_dir, "Scatterplot_KA_DMG_vs_DEG")
)

# Scatterplot: Shared DMG-DEG (KI vs KA)
scatterplot_meth_and_rna(
  shared.DMG.DEG.venn.result[shared.DMG.DEG.venn.result$Subset == "Shared", 
                             c("Detail", "mean", "log2FoldChange")],
  xlab = "KI methylation difference (%)",
  ylab = "KA RNA log2FC",
  labelmin = 30,
  filename.base = file.path(output_dir, "Scatterplot_Shared_DMG_DEG")
)

