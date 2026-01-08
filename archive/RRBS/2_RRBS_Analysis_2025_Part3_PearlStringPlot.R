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
## Pearl String Plot - Methylation
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
cran_packages <- c(
  "dplyr", "tidyr", "ggplot2", "reshape2", 
  "rlang", "RColorBrewer", "svglite", "scales"
)

# Required Bioconductor packages
bioc_packages <- c(
  "GenomicRanges", "GenomicFeatures", "methylKit", "Gviz", "genomation"
)

# Automatically install and load CRAN and Bioconductor packages
p_load(char = cran_packages, install = TRUE)
p_load(char = bioc_packages, install = TRUE, repos = BiocManager::repositories())



################################################################################
## functions - begin
################################################################################

################################################################################
## function: calculate_methylation_levels
################################################################################
## this function will calculate percent methylation from a methylBase object
calculate_methylation_levels <- function(data, sample_ids = NULL, coverage_prefix = "coverage", numCs_prefix = "numCs") {
  
  # Find columns with coverage information
  cov_cols <- grep(paste0("^", coverage_prefix), colnames(data))
  
  # Determine the number of samples and CpGs
  if (is.null(sample_ids)) {
    # chec if 
        num_samples <- length(unique(sub(paste0("^", coverage_prefix), "", colnames(data)[cov_cols])))
    sample_ids <- unique(sub(paste0("^", coverage_prefix), "", colnames(data)[cov_cols]))
  } else {
    num_samples <- length(sample_ids)
  }
  num_CpGs <- nrow(data)
  
  # Initialize a data frame to store the methylation levels
  meth_levels <- data.frame(matrix(ncol = num_samples + 4, nrow = num_CpGs))
  colnames(meth_levels)[1:4] <- colnames(data)[1:4]
  colnames(meth_levels)[5:(num_samples + 4)] <- sample_ids
  rownames(meth_levels) <- rownames(data)
  meth_levels[, 1:4] <- data[, 1:4]
  
  # Loop over the samples and calculate the methylation level for each CpG
  for (i in 1:num_samples) {
    # Create a data frame with the coverage and numCs columns for the current sample
    cov_col <- paste0(coverage_prefix, i)
    numCs_col <- paste0(numCs_prefix, i)
    repDF <- data[, c(cov_col, numCs_col)]
    
    # Calculate the methylation level for the current sample
    meth_levels[, i+4] <- repDF[, 2] / repDF[, 1]
  }
  
  # Update the column names
  colnames(meth_levels)[5:(num_samples + 4)] <- paste0(colnames(meth_levels)[5:(num_samples + 4)], "_pmeth")
  
  # Return the methylation levels
  return(meth_levels)
}

################################################################################
create_epitargets <- function (methsubset.obj, gene.obj) {
  
  ## create annotated objects
  methsubset.obj.annotated <- annotateWithGeneParts(as(methsubset.obj, "GRanges"), gene.obj)
  
  ## merge with Ensembl gene annotation object
  ## this is to be able to work using symbols
  methsubset.obj.disttoTSS.epitargets <- methsubset.obj.annotated@dist.to.TSS %>% 
    merge (ensembl.gene, by.x="feature.name", by.y="ensembl_transcript_id", 
           all.x=TRUE) %>% filter(external_gene_name %in% epitargets) 
  
  ## extract meth data and calculate % methylation
  methsubset.obj.data <- getData(methsubset.obj)
  methsubset.obj.data.pmeth <- calculate_methylation_levels(methsubset.obj.data,
                                                            sample_ids=methsubset.obj@sample.ids)
  
  ## merge two data frames (annotation and methylation levels)
  methsubset.obj.epitargets <- cbind(methsubset.obj.disttoTSS.epitargets,
                                     methsubset.obj.data.pmeth[methsubset.obj.disttoTSS.epitargets$target.row,])

  ## identify the longest transcript for each gene
  methsubset.obj.epitargets.transcript <- methsubset.obj.epitargets %>%
    group_by(feature.name, ensembl_gene_id) %>%
    summarise(transcript_length = max(start) - min(start) + 1) %>%
    ungroup()
  
  ## Keep only the longest transcript per unique ensembl_gene_id
  methsubset.obj.epitargets.longest <- methsubset.obj.epitargets.transcript %>%
    group_by(ensembl_gene_id) %>%
    slice(which.max(transcript_length)) %>%
    ungroup()
  
  # Subset methsubset.obj.epitargets based on the longest transcripts
  methsubset.obj.epitargets_subset <- subset(methsubset.obj.epitargets,
                                             feature.name %in% methsubset.obj.epitargets.longest$feature.name)
  
  
  ## Prepare output
  methsubset.obj.epitargets_subset$coor <- paste(methsubset.obj.epitargets_subset$chr, methsubset.obj.epitargets_subset$start, sep = "_")
  methsubset.obj.epitargets_subset <- methsubset.obj.epitargets_subset[,c(ncol(methsubset.obj.epitargets_subset),
                                                            1:(ncol(methsubset.obj.epitargets_subset)-1))]
  
  return(methsubset.obj.epitargets_subset)
}


################################################################################
## Load gene annotation
################################################################################

################################################################################
# Note: I am evaluating the differences from these two gtf files
#       I think the latter one "Rattus_norvegicus.Rnor_6.0.104_chr.gtf" 
#       was used in methylseq analysis; however, 
#       'gene.obj' is from Kai, which seems to be based on
#       different annotation. 
################################################################################

# Specify the path to the GTF file containing gene annotations
gtfFile <- "rn6.refGene.gtf"


## define the target
#epitargets <- c("Stc2","Mpped1","Adcy2","Nedd9","Ptpre","Dpf3","Thra","Etv4","Sirt2","Tspan2")
epitargets <- c("Fcgbpl1","Loxhd1","Nedd9","Ptpre","Adcy4")



################################################################################
## functions - end
################################################################################

## create epi-target ojbect
KIadvmethsubset.epitargets <- create_epitargets(methsubset.obj=KIadvmethsubset, 
                                                gene.obj=gene.obj)

KAadvmethsubset.epitargets <- create_epitargets(methsubset.obj=KAadvmethsubset, 
                                                gene.obj=gene.obj)

create_individual_methylation_gene_plots_v11 <- function(methylation_data, 
                                                         P25.sig.obj = NULL,
                                                         P10.sig.obj = NULL,
                                                         basename = NULL,
                                                         cordType = "original",
                                                         plotMode = "individual",
                                                         use_binned_scale = FALSE,
                                                         show_grid = TRUE,
                                                         flat_y = FALSE) {
  ## Subset the methylation data for this group
  methylation_data <- methylation_data %>%
    dplyr::select(starts_with("coor") | ends_with("_pmeth")) %>%
    cbind(external_gene_name = methylation_data$external_gene_name) %>% 
    dplyr::select(1, last_col(), everything()) %>% 
    separate(coor, into = c("chromosome", "position"), sep = "_") %>%
    mutate(position = as.numeric(position))
  
  ## Calculate evenly spaced values for 'position'
  methylation_data <- methylation_data %>%
    group_by(external_gene_name) %>%
    arrange(chromosome, position) %>%
    mutate(evenly_spaced_position = seq(min(position), max(position), length.out = n())) %>%
    dplyr::select(1:3, last_col(), everything())
  
  ## Calculate bin size - for significant CpGs
  length_out <- methylation_data %>%
    group_by(external_gene_name) %>%
    summarize(length_out = as.integer(unique(diff(evenly_spaced_position)))) %>% unique()
  
  ## merge into methylation_data
  methylation_data <- merge(methylation_data, length_out, 
                            by = c("external_gene_name")) %>%
    dplyr::select(1:4, last_col(), everything())
  
  ## Reshape the data to long format using gather()
  methylation_data_long <- methylation_data %>% as.data.frame %>%
    gather(sample, methylation, 6:ncol(methylation_data)) %>% 
    group_by(external_gene_name)
  
  ##############################################################################
  ## create color palette and updated sample names
  ##############################################################################
  color.df <- data.frame(sample = unique(methylation_data_long$sample),
                         sample_short = unique(sub("_pmeth", "", methylation_data_long$sample)))
  color.df$group <- group[color.df$sample_short, "group"]
  color.df.unique.groups <- unique(group[color.df$sample_short, "group"])
  color.df$sample_updated <- gsub("_pmeth", "", paste0(color.df$sample_short, "_", color.df$group))
  rownames(color.df) <- color.df$sample_short
  
  ## Create a set of colors similar to red
  red_colors <- colorRampPalette(c("#FFEDA0", "#F03B20"))(nrow(subset(color.df, group == color.df.unique.groups[1])))
  
  ## Create a set of colors similar to blue
  blue_colors <- colorRampPalette(c("#BFD3E6", "#08306B"))(nrow(subset(color.df, group == color.df.unique.groups[2])))
  
  ## merge two palettes
  my_colors <- c(red_colors, blue_colors)
  
  ## convert sample to factor
  methylation_data_long$sample <- factor(methylation_data_long$sample,
                                         level = unique(methylation_data_long$sample))
  
  ## add updated_sample column
  methylation_data_long$sample_updated <- color.df[gsub("_pmeth", "", methylation_data_long$sample), "sample_updated"]
  methylation_data_long$sample_updated <- factor(methylation_data_long$sample_updated,
                                                 level = unique(methylation_data_long$sample_updated))
  
  ## highlight those with significantly different
  methylation_data_long$significant <- rep("No", length(methylation_data_long[, 1]))
  methylation_data_long$significant[which(paste0(methylation_data_long$chromosome, "_", methylation_data_long$position) %in% P25.sig.obj$cpg.site)] <- "Yes"
  
  ## highlight those with significantly different
  methylation_data_long$significant_10p <- rep("No", length(methylation_data_long[, 1]))
  
  # Check if "cpg.site" column exists
  if (!"cpg.site" %in% names(P10.sig.obj)) {
    # If not, concatenate "start" and "end" columns and create a new "cpg.site" column
    P10.sig.obj$cpg.site <- paste(P10.sig.obj$chr, P10.sig.obj$start, sep = "_")
  }
  
  methylation_data_long$significant_10p[which(paste0(methylation_data_long$chromosome, "_", methylation_data_long$position) %in% P10.sig.obj$cpg.site)] <- "Yes"
  
  methylation_data_long$significant_color <- ifelse(methylation_data_long$significant == "Yes", "red", 
                                                    ifelse(methylation_data_long$significant_10p == "Yes", "yellow", "white"))
  
  ## Split the data by `external_gene_name`
  gene_data <- split(methylation_data_long, methylation_data_long$external_gene_name)
  
  # Create the subdirectory if it doesn't exist
  if (!dir.exists("output_gene")) {
    dir.create("output_gene")
  }
  
  ###################################################################################
  # Logic to decide which x-axis to use based on cordType
  x_cord <- "position"
  if (!is.na(cordType)) {
    cordType <- tolower(cordType)
    x_cord <- ifelse(cordType == "original", "position", 
                     ifelse(cordType == "evenlydistributed", "evenly_spaced_position", NA))
  }
  
  # Create a separate plot for each gene
  for (i in seq_along(gene_data)) {
    
    ## Add a rectangle for the subset of x-axis significant ones
    significant_data <- unique(gene_data[[i]][gene_data[[i]]$significant_10p == "Yes",
                                              c("position", "evenly_spaced_position", "length_out", "significant_color")])
    
    ## Create fill_scale to handle the continuous or binned gradients
    fill_scale <- if (use_binned_scale) {
      scale_fill_binned(
        low = "white",
        high = "black",
        n.breaks = 11,
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.1),
        labels = seq(0, 1, by = 0.1)
      )
    } else {
      scale_fill_gradient(
        low = "white",
        high = "black",
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      )
    }
    
    ## Define the grid theme based on the show_grid parameter
    grid_theme <- if (show_grid) {
      theme_minimal()
    } else {
      theme_minimal() + theme(panel.grid = element_blank())
    }
    
    if (plotMode == "individual") {
      p <- ggplot(gene_data[[i]]) 
      if (flat_y) {
        p <- p + geom_point(aes(x = .data[[x_cord]], 
                                y = as.numeric(as.factor(sample_updated)), 
                                fill = methylation, 
                                group = sample_updated), 
                            color = "black", shape = 21, stroke = 0.25, size = 5) +
          geom_line(aes(x = .data[[x_cord]], 
                        y = as.numeric(as.factor(sample_updated)), 
                        group = sample_updated),
                    linewidth = 0.5, alpha = 0.40) +
          xlab(paste0("CpG on ", unique(gene_data[[i]]$chromosome))) +
          ylab("Samples") +
          ggtitle(paste0("Methylation levels for gene ", names(gene_data)[i])) +
          scale_y_continuous(breaks = 1:length(unique(gene_data[[i]]$sample_updated)), 
                             labels = unique(gene_data[[i]]$sample_updated),
                             limits = c(0.5, length(unique(gene_data[[i]]$sample_updated)) + 0.5))
      } else {
        p <- p + geom_line(aes(x = .data[[x_cord]], 
                               y = as.numeric(methylation), 
                               group = sample_updated,
                               color = sample_updated),
                           linewidth = 0.2, alpha = 0.40) +
          geom_point(aes(x = .data[[x_cord]], 
                         y = as.numeric(methylation), 
                         fill = methylation, 
                         color = sample_updated, 
                         group = sample_updated),
                     shape = 21, stroke = 0.25, size = 3) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))
      }
      
      p <- p + scale_color_manual(values = my_colors) +
        
        # The rest of the plot details
        xlab(paste0("CpG on ", unique(gene_data[[i]]$chromosome))) +
        ylab(ifelse(flat_y, "Samples", "Methylation level")) +
        ggtitle(paste0("Methylation levels for gene ", names(gene_data)[i])) +
        grid_theme +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_blank()
        ) +
        
        # Drawing the rectangles
        geom_rect(data = significant_data,
                  aes(xmin = .data[[x_cord]] - length_out/2, 
                      xmax = .data[[x_cord]] + length_out/2, 
                      ymin = -Inf, ymax = Inf),
                  fill = significant_data$significant_color,
                  alpha = 0.15, inherit.aes = FALSE) +
        fill_scale
    } else if (plotMode == "average") {
      
      # Calculate the average and standard deviation for each CpG site and sample
      gene_summary <- gene_data[[i]] %>% 
        separate(sample_updated, into = c("sample", "group"), sep = "_") %>%
        group_by(!!sym(x_cord), group) %>%
        summarize(mean_methylation = mean(as.numeric(methylation)),
                  sd_methylation = sd(as.numeric(methylation)))
      
      # Define the colors dynamically based on the unique groups
      unique_groups <- unique(gene_summary$group)
      if (length(unique_groups) == 2) {
        my_colors <- setNames(c("#FF0000", "#0000FF"), unique_groups)
      } else {
        my_colors <- setNames(colorRampPalette(c("#FF0000", "#0000FF"))(length(unique_groups)), unique_groups)
      }
      
      p <- ggplot(gene_summary)
      if (flat_y) {
        p <- p + geom_point(aes(x = .data[[x_cord]], 
                                y = as.numeric(as.factor(group)), 
                                fill = mean_methylation, 
                                group = group), 
                            color = "black", shape = 21, stroke = 0.25, size = 5) +
          geom_line(aes(x = .data[[x_cord]], 
                        y = as.numeric(as.factor(group)), 
                        group = group),
                    linewidth = 0.5, alpha = 0.40) +
          xlab(paste0("CpG on ", unique(gene_data[[i]]$chromosome))) +
          ylab("Groups") +
          ggtitle(paste0("Average Methylation levels for gene ", names(gene_data)[i])) +
          grid_theme +
          theme(plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_blank()) +
          scale_y_continuous(breaks = 1:length(unique(gene_summary$group)), 
                             labels = unique(gene_summary$group),
                             limits = c(0.5, length(unique(gene_summary$group)) + 0.5))
      } else {
        p <- p + geom_line(aes(x = .data[[x_cord]], 
                               y = mean_methylation, 
                               group = group,
                               color = group),
                           linewidth = 0.2) +
          geom_point(aes(x = .data[[x_cord]], 
                         y = mean_methylation, 
                         fill = mean_methylation, 
                         color = group,
                         group = group),
                     shape = 21, stroke = 0.25, size = 3) +
          geom_errorbar(aes(x = .data[[x_cord]], 
                            ymin = mean_methylation - sd_methylation, 
                            ymax = mean_methylation + sd_methylation,
                            group = group,
                            color = group),
                        width = 0.1,
                        size = 0.1) +
          scale_color_manual(values = my_colors) +
          xlab(paste0("CpG on ", unique(gene_data[[i]]$chromosome))) +
          ylab("Average Methylation level") +
          ggtitle(paste0("Average Methylation levels for gene ", names(gene_data)[i])) +
          grid_theme +
          theme(plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_blank()) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
          scale_x_continuous(breaks = seq(min(gene_summary[[x_cord]]), max(gene_summary[[x_cord]]), length.out = 5))
        
      }
      
      p <- p + fill_scale
      
      if (!is.null(significant_data) && nrow(significant_data) > 0) {
        p <- p + geom_rect(data = significant_data,
                           aes(xmin = .data[[x_cord]] - length_out/2, 
                               xmax = .data[[x_cord]] + length_out/2, 
                               ymin = -Inf, ymax = Inf),
                           fill = significant_data$significant_color,
                           alpha = 0.15, inherit.aes = FALSE)
      }
    }
    
    # Define the height based on plot mode to adjust the plot size
    height <- if (plotMode == "individual") 5 else 3
    
    # Define the file name components
    file_name_component <- paste0("_", substr(plotMode, 1, 3), 
                                  "_", substr(cordType, 1, 3), "_", 
                                  ifelse(use_binned_scale, "binned", "continuous"),
                                  ifelse(show_grid, "", "_noGrid"),
                                  ifelse(flat_y, "_flatY", ""))
    
    ## Save the plot as a PDF file with the appropriate file name
    ggsave(file = paste0("output_gene/", basename, names(gene_data)[i], file_name_component, ".pdf"), plot = p, width = 12, height = height)
    ggsave(file = paste0("output_gene/", basename, names(gene_data)[i], file_name_component, ".svg"), plot = p, width = 12, height = height)
    
  }
}




## Different combinations
## KI
create_individual_methylation_gene_plots_v11(methylation_data=KIadvmethsubset.epitargets,
                                             P25.sig.obj=DBKIAdv.RRBS.25,
                                             P10.sig.obj=DBKIAdv.RRBS.10,
                                             basename="KI_",
                                             cordType="original",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KIadvmethsubset.epitargets,
                                             P25.sig.obj=DBKIAdv.RRBS.25,
                                             P10.sig.obj=DBKIAdv.RRBS.10,
                                             basename="KI_",
                                             cordType="original",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE,
                                             flat_y = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KIadvmethsubset.epitargets,
                                             P25.sig.obj=DBKIAdv.RRBS.25,
                                             P10.sig.obj=DBKIAdv.RRBS.10,
                                             basename="KI_",
                                             cordType="evenlyDistributed",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KIadvmethsubset.epitargets,
                                             P25.sig.obj=DBKIAdv.RRBS.25,
                                             P10.sig.obj=DBKIAdv.RRBS.10,
                                             basename="KI_",
                                             cordType="evenlyDistributed",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE,
                                             flat_y = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KIadvmethsubset.epitargets,
                                             P25.sig.obj=DBKIAdv.RRBS.25,
                                             P10.sig.obj=DBKIAdv.RRBS.10,
                                             basename="KI_",
                                             cordType="evenlyDistributed",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = FALSE)





## KA
create_individual_methylation_gene_plots_v11(methylation_data=KAadvmethsubset.epitargets,
                                             P25.sig.obj=DBKAAdv.RRBS.25,
                                             P10.sig.obj=DBKAAdv.RRBS.10,
                                             basename="KA_",
                                             cordType="original",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KAadvmethsubset.epitargets,
                                             P25.sig.obj=DBKAAdv.RRBS.25,
                                             P10.sig.obj=DBKAAdv.RRBS.10,
                                             basename="KA_",
                                             cordType="original",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE,
                                             flat_y = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KAadvmethsubset.epitargets,
                                             P25.sig.obj=DBKAAdv.RRBS.25,
                                             P10.sig.obj=DBKAAdv.RRBS.10,
                                             basename="KA_",
                                             cordType="evenlyDistributed",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KAadvmethsubset.epitargets,
                                             P25.sig.obj=DBKAAdv.RRBS.25,
                                             P10.sig.obj=DBKAAdv.RRBS.10,
                                             basename="KA_",
                                             cordType="evenlyDistributed",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = TRUE,
                                             flat_y = TRUE)

create_individual_methylation_gene_plots_v11(methylation_data=KAadvmethsubset.epitargets,
                                             P25.sig.obj=DBKAAdv.RRBS.25,
                                             P10.sig.obj=DBKAAdv.RRBS.10,
                                             basename="KA_",
                                             cordType="evenlyDistributed",
                                             plotMode="average",
                                             use_binned_scale=TRUE,
                                             show_grid = FALSE)

