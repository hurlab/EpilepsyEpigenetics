######################################################
## Set the current working directory
######################################################
library(rstudioapi) # make sure you have it installed
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
base_dir = dirname(current_path)

######################################################
## Define some variables
######################################################
analysis_date 		= '050124'
output_dir		  	= 'output/'
bam_dir			    = 'bam_dir/'
genes.gtf.file		= 'rn6.refGene.gtf'
sample_file		  	= 'selected_samples.txt'
basename			= 'epiRNA'
additional_tag  	= ''
thread			    = 32
dir.create(file.path("./", output_dir), showWarnings = FALSE)


################################################################################
#####################  Packages  ###############################################
################################################################################

# Install pacman for package management if not already installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# Required packages
# CRAN packages
cran_packages <- c(
  "gplots", "reshape2", "ggplot2", "NMF", "dplyr", 
  "RColorBrewer", "pheatmap", "ggpubr", "dendextend", "pca3d"
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2", "Rsubread", "GenomicFeatures", "org.Hs.eg.db",
  "org.Rn.eg.db", "EnsDb.Hsapiens.v86", "GO.db", "VennDetail"
)

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
####################  Functions  ###############################################
################################################################################
## This Function will allow us to make a MA plot from a list of results objects

maPlot.lists <- function(x,i) {
  pdf(paste(output_dir, i, '_maPlot_', analysis_date, '.pdf', sep=''))
  plotMA(x, main=paste(i, 'alpha=0.01', sep=' '), 
         alpha=0.01, ylim=c(-6,6))
  abline(h=c(2,-2), col='red')
  dev.off()
}


perform_DE_analysis <- function (dds, experiment, reference) {
  # run DESeq2
  res <- results(dds,contrast = c("condition",experiment,reference))
  
  # add symbol to res
  res$SYMBOL <- rownames(geneAnno)

  # check summary
  #summary(res)
  
  # base.name
  base.name = paste0(experiment, '_vs_', reference)
  
  # create ggmaplot
  pdf(file=paste0(output_dir, base.name, "_ggmaplot_",analysis_date,".pdf"))
  ggmaplot(as.data.frame(res), main = paste(experiment, " vs ", reference, sep=""),
           fdr = 0.05, fc = 2, size = 0.4,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           genenames = as.vector(res$SYMBOL),
           legend = "top", top = 50,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal())
  dev.off()

  # write out all results 
  write.table(res[order(res$padj), ], quote=F, sep="\t",
              file=paste0(output_dir, base.name, '_allResults_',analysis_date,'.txt'))
  
  # write out sign results 
  sigRes <- subset(res, padj<=0.05)
  
  write.table(sigRes, quote=F,sep="\t",
              file=paste0(output_dir, base.name,'_DEtable_padj0.05_',analysis_date,'.txt'))
  
  # Another MA plot
  # let's take a look at the distribution of the DE genes with a MA plot
  pdf(paste(output_dir, base.name, "_MAPlot_",analysis_date,".pdf",sep=""))
  plotMA(res, main=paste(experiment, " vs ", reference, sep=""),
         alpha=0.05, ylim=c(-4,4))
  abline(h=c(1,-1), col='red')
  dev.off()
  
  ## Heatmap of DEGs
  sig_rld <- rld_mat_ext[rownames(sigRes), grep('Mean',colnames(rld_mat_ext))]
  if (length(rownames(sigRes))>2) {
    aheatmap(sig_rld, color = "-RdBu:50", scale = "row", breaks = 0,
             distfun = "spearman",treeheight=c(15,10),
             fontsize=10, cexCol=.7, verbose=T, 
             labRow = geneAnno[rownames(sig_rld),"SYMBOL"],
             filename=paste0(output_dir, base.name, "_heatmap_",analysis_date,".pdf"))  
    
    ## GO enrichment
    input.gene<-sub("\\.\\d+$","",rownames(sigRes))
    GO.res <- richGO(input.gene, godata=go_anno, ontology ="BP", 
                     pvalue = 0.05, padj = 0.10, minSize = 2, maxSize = 500, 
                     keepRich = TRUE, filename=paste0(output_dir,base.name,"_GO.out"), 
                     padj.method="BH", sep=", ")
    GO.sig <- filter(GO.res, Padj<0.05)
    
    ## generate plots
    if (length(rownames(GO.sig@result))>0) {
      pdf(paste0(output_dir,base.name,"_GO-bar_Top50.pdf"), width=12)
      ggbar(GO.sig, top = 50, usePadj = T, horiz=T) + ylab("% in genome") + ggtitle("Enrichment analysis")
      dev.off()
      pdf(paste0(output_dir,base.name,"_GO-dot_Top50.pdf"), width=12)
      ggdot(GO.sig, top = 50, usePadj = T) + ylab("% in genome") + ggtitle("Enrichment analysis")
      dev.off()  
    }

    ## KEGG enrichment
    KEGG.res <- richKEGG(input.gene, kodata=ko_anno, pvalue=0.05, builtin = FALSE)
    KEGG.sig <- filter(KEGG.res, Padj < 0.05)
    write.table (KEGG.res, file=paste0(output_dir,base.name,"_KEGG.out.txt"),
                 sep="\t", quote=F)
    
    ## generate plots
    if (length(rownames(KEGG.sig@result))>0) {
      pdf(paste0(output_dir,base.name,"_KEGG-bar_Top50.pdf"), width=12)
      ggbar(KEGG.sig, top = 50, usePadj = T, horiz=T) + ylab("% in genome") + ggtitle("Enrichment analysis")
      dev.off()
      pdf(paste0(output_dir,base.name,"_KEGG-dot_Top50.pdf"), width=12)
      ggdot(KEGG.sig,top=10,usePadj = F)
      dev.off()
    }
  }
}




################################################################################
####################  Counting using Rsubread  #################################
################################################################################
## Get the files and run featureCounts (Rsubread) using a linked gtf file
## You will only need to adjust the nthreads option to accomodate your hardware
## Assuming that you are running this script in the directory containing all bam
## files along with a linked gtf file

## working directory changed to my home directory
## still using the files in Atrayee's directory
setwd(paste(base_dir,bam_dir,sep="/"))
file_list <- list.files(pattern='*.bam$')
counts <- featureCounts(files=file_list, annot.ext=genes.gtf.file, 
                        isGTFAnnotationFile=T, GTF.featureType='exon', 
                        GTF.attrType='gene_id', nthreads=thread,
                        isPairedEnd=F)
setwd(base_dir)
counts.bak<-counts

## process file names 
## remove '.sort.bam'
counts$targets<-colnames(counts$counts)<-
  sub('^\\d+\\_','',sub('\\.sorted\\.bam','',colnames(counts$counts)))
colnames(counts$stat)<-tmp.names<-
  sub('^\\d+\\_','',sub('\\.sorted\\.bam','',colnames(counts$stat)))

## Write the counting counts/statistics to files 
write.csv(counts$stat, file=paste(output_dir, 'subread_countStats_', 
                                    analysis_date, additional_tag, '.csv', sep=""), quote=F, row.names=F)
write.csv(counts$counts, file=paste(output_dir, 'subread_counts_', 
                                    analysis_date, additional_tag, '.csv', sep=""), quote=F)


################################################################################
####################  Update counts data     ###################################
################################################################################
## many of the ensembl gene IDs are obsolete
## let's remove them before we do any analysis

## create a geneAnno data frame
geneAnno <- data.frame (SYMBOL = rownames(counts$counts))

## make the symbol unique
geneAnno$SYMBOL_UNIQ<-make.unique(geneAnno$SYMBOL, sep="^")
rownames(geneAnno)<-geneAnno$SYMBOL


################################################################################
####################  Build GO and KEGG annotation objects   ###################
################################################################################

## Gene ontology using GO.db
go_anno <- buildAnnot(species="rat", keytype="SYMBOL", anntype = "GO")

## KEGG pathway (up-to-date information from KEGG website)
ko_anno <- buildAnnot(species="rat", keytype="SYMBOL", anntype = "KEGG", builtin = FALSE)



################################################################################
####################  Group and comparison   ###################################
################################################################################
group<-read.delim("selected_samples.txt",sep="\t",header=T)
group.bak<-group
rownames(group)<-group$sample

## Reorder the count data based on the current sample file
## Note that the sample order in count object is not the same as in the sample
## list file
counts$counts<-counts$counts[,group$sample]
counts$targets<-group$sample

## Load comparison.txt file
comparison<-read.delim("comparison.txt",sep="\t")

## define condition and replicate object
condition<-group[counts$targets,"group"]
replicates<-group[counts$targets,"replicate"]



################################################################################
####################  DDS object initialization   ##############################
################################################################################

## create dds object
dds<-DESeqDataSetFromMatrix(counts$counts,DataFrame(condition),~condition)

## Add metadata from the gtf.  This will allow us to retrieve things like 
## gene identity from the assembled data set.
mcols(dds) <- DataFrame(mcols(dds), counts$annotation)

## Call algorithm on the assembled data set
dds <- DESeq(dds)

## Plot the dispersion of the experiment to verify that the Algorithm's 
## assumptions are valid for this dataset.  This will also show us if 
## the variance is too LOW in the samples indicating an error in replication
pdf(paste(output_dir, 'dispModel_', analysis_date, '.pdf', sep=""))
plotDispEsts(dds)
dev.off()



################################################################################
####################  Count and FPKM data   ####################################
################################################################################
## Retrieve count data and clean up the data frame 
## gives basemean
count_data <- counts(dds, normalized=T)

## this may not be necessary if colnames are already replicate names
colnames(count_data) <- replicates 

## We need the gene lengths to calculate fpkm.  
## the following code will retrieve and calculate the gene lengths from the gtf
## file then add them to the dds object.  This may take some time...
txdb <- makeTxDbFromGFF(genes.gtf.file, format='gtf')
exons.list.per.gene <- exonsBy(txdb, by='gene')
exonic.gene.sizes <- lapply(exons.list.per.gene,
                            function(x){sum(width(reduce(x)))})
gene_size <- data.frame(unlist(exonic.gene.sizes))
gene_size <- data.frame(cbind(gene_size, 
                              kb=apply(gene_size, 1, 
                                       function(x){x/1000})))
colnames(gene_size)[1] <- 'length'

## add to dds object
mcols(dds)$basepairs <- gene_size[,'length']

## Get fpkm data
fpkm_data <- fpkm(dds)
colnames(fpkm_data) <- colnames(count_data)

## save outfiles
write.table(count_data, file=paste(output_dir,'depthNormCountFiltered_',
                                   analysis_date,'.txt', sep=""), quote=F, sep="\t")
write.table(fpkm_data, file=paste(output_dir,'fpkmDataFiltered_',
                                  analysis_date,'.txt', sep=""), quote=F, sep="\t")


## Here we will add mean counts for each sample to the count and fpkm data frames.
## The is a purely convenience operation to make the output tables more useful.
## Then write those tables to file.
## You will need to match these column names to your condition variable.
## Additional lines will need to be added if you have more than two conditions
count_data <- data.frame(cbind(count_data, 
    ShamMean=apply(count_data[,(condition == "Sham")],1,mean),
    ThreeClassV1Mean=apply(count_data[,(condition == "3Class V")],1,mean),
    SalineShamMean=apply(count_data[,(condition == "Saline Sham")],1,mean),
    ThreeClassV2Mean=apply(count_data[,(condition == "3 Class V")],1,mean))
fpkm_data <- data.frame(cbind(fpkm_data, 
    ShamMean=apply(fpkm_data[,(condition == "Sham")],1,mean),
    ThreeClassV1Mean=apply(fpkm_data[,(condition == "3Class V")],1,mean),
    SalineShamMean=apply(fpkm_data[,(condition == "Saline Sham")],1,mean),
    ThreeClassV2Mean=apply(fpkm_data[,(condition == "3 Class V")],1,mean))

write.table(count_data, file=paste(output_dir,'depthNormCountWithMean_',
            analysis_date,'.txt', sep=""), quote=F, sep="\t")
write.table(fpkm_data, file=paste(output_dir,'fpkmDatawithMean_',
            analysis_date,'.txt', sep=""), quote=F, sep="\t")



###############################################################################
####################  Conditional QC   ########################################
###############################################################################
                                           
## Perform rLog transformation on the data 
rld <- rlog(dds) 
rld.bak <- rld

## First we'll make a distance clustering HM to look at replicate grouping
## Set the colors to be used. help(brewer.pal) will show some options
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

## Here we will use the base functions for calculating euclidian distance and 
## then cluster that data.  Take a look at the intermediate variables with 
## str() to see what we are doing with these fxns
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- colnames(mat) <- as.character(replicates)
hc <- hclust(dists)

## create a plot
par(mar = c(7.1, 4.1, 4.1, 2.1))
plot(hc)
dev.copy2pdf(file=paste0(output_dir,'hierarchicalClustering_v1_',
                         analysis_date,'.pdf', sep=""))

## create a temporary vector and handle the naming issue
## this step really depends on your file names
## so it needs to be adjusted

## alternative hierarchical clustering and label changes
hc2 <- as.dendrogram(hc)
plot(hc2)
dev.copy2pdf(file=paste0(output_dir,'hierarchicalClustering_v2_',
                         analysis_date,'.pdf', sep=""))

## Use replicate names
## alternative hierarchical clustering and label changes
colnames(rld)<-replicates
dists <- dist(t(assay(rld)))
hc <- hclust(dists)
hc2 <- as.dendrogram(hc)
plot(hc2)
dev.copy2pdf(file=paste0(output_dir,'hierarchicalClustering_v3_',
                         analysis_date,'.pdf', sep=""))






#################################################################################
## Plot HM and save to file.  Run without the pdf() and dev.off() fxns to plot
## in the R window before saving to file.  There are several other file type 
## options besides pdf  
#################################################################################
pdf(paste(output_dir, 'distClustering_', analysis_date, '.pdf', sep=""),
    width=12,height=10)
heatmap.2(mat, Rowv=as.dendrogram(hc2), symm=T, trace='none',
          col=rev(hmcol), margin=c(13,13), main='Distance heatmap')
dev.off()


#################################################################################
## Check correlation
#################################################################################
rld_mat <- assay(rld)

## Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
#head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

pdf(paste0(output_dir, 'correlationHeatmap_', analysis_date, '.pdf'))
pheatmap(rld_cor)
dev.off()

## create extended rld_mat_ext
rld_mat_ext <- data.frame(cbind(rld_mat, 
    ShamMean=apply(rld_mat[,(condition == "Sham")],1,mean),
    ThreeClassV1Mean=apply(rld_mat[,(condition == "3ClassV1")],1,mean),
    SalineShamMean=apply(rld_mat[,(condition == "SalineSham")],1,mean),
    ThreeClassV2Mean=apply(rld_mat[,(condition == "3ClassV2")],1,mean))


#################################################################################
## a PCA plot is another way to look at the data in a similar way.  Here we 
## plot one and save it to file.
#################################################################################
pdf(paste(output_dir, 'PCA_',analysis_date,'.pdf', sep=''))
plotPCA(rld, intgroup=c('condition')) +  stat_ellipse()
dev.off()

## Let's use different package to plot different PC combinations
library(pca3d)
pca <- prcomp(t(assay(rld)), center=T, scale.=F)
pca2d(pca, component=c(1,3), group=condition, labels=T, legend=T)
pca2d(pca, component=c(1,2), group=condition, labels=T, legend=T)
pca2d(pca, component=c(2,3), group=condition, labels=T, legend=T)

## Different versions of plots
scores <- data.frame(condition, pca$x[,1:3])

## PC1 vs PC2
pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=factor(condition),size=I(3),
               xlab=paste0("PC1 (",round(summary(pca)$importance[2,]*100,1)[1],"%)"),
               ylab=paste0("PC2 (",round(summary(pca)$importance[2,]*100,1)[2],"%)"))+ 
          theme(legend.position="right", text=element_text(size=10),
            axis.text=element_text(size=15), axis.title=element_text(size=20)) +
          geom_text(aes(label=rownames(scores)), size=3 ,vjust="inward", hjust="inward")
pdf(file=paste0(output_dir,'PCAplot_PC1-PC2_',analysis_date,'.pdf'))
print(pc1.2)
dev.off()

## PC1 vs PC3
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=factor(condition),size=I(3),
               xlab=paste0("PC1 (",round(summary(pca)$importance[2,]*100,1)[1],"%)"),
               ylab=paste0("PC3 (",round(summary(pca)$importance[2,]*100,1)[3],"%)"))+ 
  theme(legend.position="right", text=element_text(size=10),
        axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  geom_text(aes(label=rownames(scores)), size=3 ,vjust="inward", hjust="inward")
pdf(file=paste0(output_dir,'PCAplot_PC1-PC3_',analysis_date,'.pdf'))
print(pc1.3)
dev.off()

## PC2 vs PC3
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=factor(condition),size=I(3),
               xlab=paste0("PC2 (",round(summary(pca)$importance[2,]*100,1)[2],"%)"),
               ylab=paste0("PC3 (",round(summary(pca)$importance[2,]*100,1)[3],"%)"))+ 
  theme(legend.position="right", text=element_text(size=10),
        axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  geom_text(aes(label=rownames(scores)), size=3 ,vjust="inward", hjust="inward")
pdf(file=paste0(output_dir,'PCAplot_PC2-PC3_',analysis_date,'.pdf'))
print(pc2.3)
dev.off()



################################################################################
#################### DE Tests
################################################################################
## Set up the results comparisons.  ORDER MATTERS!!  The second label will be 
## the denominator (reference) in the comparison, 
## FC will be read relative to this group.
## Many of these results comparisons objects will need to be created for more
## complicated experiments.

## this for loop will run differential analyses automatically as defined
## in the 'comparison.txt' file
for (i in c(1:length(comparison[,1]))) {
  perform_DE_analysis (dds, experiment=toString(comparison[i,1]), reference=toString(comparison[i,2]))
}



################################################################################
#################### Overlap
################################################################################

setwd(output_dir)
DEGfiles<-list.files(pattern="*DEtable*")

## create DEG list summary
DEGSummary<-as.data.frame(t(matrix(unlist(lapply(DEGfiles, function(x){
  shortname<-sub(paste0(output_dir, "DEtable_padj0.05_",analysis_date,".txt"),'',x)
  df<-read.delim(x)
  return(c(shortname, dim(df)[1]))})),nrow=2)))
colnames(DEGSummary)<-c("Comparison","#ofDEGs")

## save the DEG list summary
sink(paste0(output_dir, "Summary_DE_count_",analysis_date,".txt"))
DEGSummary
sink()

## get the DEG lists
DEGMaster<-lapply(DEGfiles, function(x){
  shortname<-sub(paste0("DEtable_padj0.05_",analysis_date,".txt"),'',x)
  return(read.delim(x))})
names(DEGMaster)<-DEGSummary$Comparison


## Some select overlap analysis
DE.output.KIKA<-list(ThreeClassV1_vs_Sham = rownames(DEGMaster[["3ClassV1_vs_Sham"]]),
                     ThreeClassV2_vs_SalineSham = rownames(DEGMaster[["3ClassV2_vs_SalineSham"]]))

## venn diagram overlap analysis
venn <- VennDetail::venndetail(DEoutputs)
pdf(paste0(output_dir, "VennDiagram_two-DEGsets.pdf"))
plot(venn)
dev.off()

VennDetail::result(venn)

