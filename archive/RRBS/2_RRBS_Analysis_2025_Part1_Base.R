######################################################
## Set the current working directory
######################################################
library(rstudioapi) # make sure you have it installed
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(current_path)

base_dir = dirname(current_path)
list.files()

######################################################
## Define some variables
######################################################
analysis_date 		= '050124'
output_dir		  	= 'output/'
bam_dir			    = 'bam_dir/'
genes_gtf_file		= 'Rattus_norvegicus.Rnor_6.0.104_chr.gtf'
sample_file		  	= 'selected_samples.txt'
basename			= 'epiRRBS'
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

# Required CRAN packages
cran_packages <- c("dplyr", "tidyverse", "rstudioapi", "stringr", "data.table")

# Required Bioconductor packages
bioc_packages <- c(
  "methylKit", "genomation", "GenomicFeatures", 
  "TxDb.Rnorvegicus.UCSC.rn6.refGene", "biomaRt"
)

# Automatically install and load CRAN and Bioconductor packages
p_load(char = cran_packages, install = TRUE)
p_load(char = bioc_packages, install = TRUE, repos = BiocManager::repositories())



################################################################################
####################  Load Rn6 gene annotation - biomart   #####################
################################################################################

ensembl = useMart ("ensembl", dataset="rnorvegicus_gene_ensembl")
ensembl.gene = getBM(attributes=c('ensembl_transcript_id',
                                  'ensembl_gene_id', 'external_gene_name',
                                  'description'), mart = ensembl)


################################################################################
####################  Group and comparison   ###################################
################################################################################
setwd(base_dir)

# group backup
group<-read.delim(sample_file,sep="\t",header=T)
rownames(group)<-group$sample
group$sampleFileShort<-paste0(sub("S","",group$sample),"_",group$sample)

## Load comparison.txt file
comparison<-read.delim("comparison.txt",sep="\t")

## define condition and replicate object
condition<-group$group
replicates<-group$replicate
treatment<-group$groupNum

###prepare all sample with sam format and remove the first row with name
setwd(file.path(base_dir, bam_dir))
filenames<-list.files(pattern="*.sorted.bam$")
filenums<-as.numeric(str_split(sub("\\.sorted\\.bam","",filenames),"_",simplify=T)[,1])
filenames.df<-data.frame(filename=filenames,filenum=filenums)
filenames.df<-filenames.df[order(filenames.df$filenum),]

## prepare sample id
file_list<-as.list(filenames.df$filename)

methobj<-processBismarkAln(location=file_list, 
                         sample.id=as.list(group$sampleFileShort),
                         treatment=treatment,
                         save.folder="./methyl/", save.context="CpG", 
                         read.context="CpG", assembly="rn6",
                         save.db=FALSE, nolap=FALSE, phred64=FALSE)

setwd(base_dir)

## Get all methylation stats
for(i in 1:length(methobj)){
  getMethylationStats(methobj[[i]],plot=TRUE,both.strands=FALSE)
  dev.print(pdf,file=paste0(output_dir,filenames.df$filename[i],"_percentage_meth.pdf"))
}

for(i in 1:length(methobj)){
  getCoverageStats(methobj[[i]],plot=TRUE,both.strands=FALSE)
  dev.print(pdf,file=paste0(output_dir,filenames.df$filename[i],"_coverage_meth.pdf"))
}

## filter can be used based on coverage or minimum count
## no huge difference in general
#filterobj=filterByCoverage(methobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)

## unite all matrix into one big matrix
## minimum coverage could be set below 
## make huge different if you used  
meth<-methylKit::unite(methobj, destrand=FALSE)
meth.bak<-meth

## simply the sample name
meth@sample.ids<-str_split(meth@sample.ids,"_",simplify=T)[,2]

## make clustering plot
clusterSamples(meth,dist="correlation", method="ward", plot=TRUE)
dev.print(pdf,file=paste0(output_dir,"meth_clustering_v1.pdf"))

## make PCA plot
PCASamples(meth)
dev.print(pdf,file=paste0(output_dir,"meth_PCA_v1.pdf"))

## screeplot for variance
PCASamples(meth,screeplot=T)
dev.print(pdf,file=paste0(output_dir,"meth_PCA_screeplot.pdf"))

## use replicate names instead of sample names
meth@sample.ids<-group$replicate

## make clustering plot
clusterSamples(meth,dist="correlation", method="ward", plot=TRUE)
dev.print(pdf,file=paste0(output_dir,"meth_clustering_v2.pdf"))

## make PCA plot
PCASamples(meth)
dev.print(pdf,file=paste0(output_dir,"meth_PCA_v2.pdf"))

## checking only the most variable ones
meth.tmp<-meth
methTopMean<-matrixMean(percMethylation(meth))
meth.tmp$Variance<-matrixVar(percMethylation(meth),methTopMean)
methTop <- arrange(meth.tmp, -Variance)

## Check top 1000
methTopF1000 <- methTop[1:1000]
PCASamples(methTopF1000)


###########################################################
## filter sex chromosome
###########################################################
methNoSexChr<-meth.bak[!meth.bak$chr%in%c("chrY","chrX"),]

## simply the sample name
methNoSexChr@sample.ids<-str_split(methNoSexChr@sample.ids,"_",simplify=T)[,2]

## make clustering plot
clusterSamples(methNoSexChr,dist="correlation", method="ward", plot=TRUE)
dev.print(pdf,file=paste0(output_dir,"methNoSexChr_clustering_v1.pdf"))

## make PCA plot
PCASamples(methNoSexChr)
dev.print(pdf,file=paste0(output_dir,"methNoSexChr_PCA_v1.pdf"))

## screeplot for variance
PCASamples(methNoSexChr,screeplot=T)
dev.print(pdf,file=paste0(output_dir,"methNoSexChr_PCA_screeplot.pdf"))

## use replicate names instead of sample names
methNoSexChr@sample.ids<-group$replicate

## make clustering plot
clusterSamples(methNoSexChr,dist="correlation", method="ward", plot=TRUE)
dev.print(pdf,file=paste0(output_dir,"methNoSexChr_clustering_v2.pdf"))

## make PCA plot
PCASamples(methNoSexChr)
dev.print(pdf,file=paste0(output_dir,"methNoSexChr_PCA_v2.pdf"))

## save the image
save.image(file="completeData_08272021_0119.RData", compress="gzip")


###########################################################
## Differential analysis
###########################################################
## this for loop will run differential analyses automatically as defined
## in the 'comparison.txt' file
for (i in c(1:length(comparison[,1]))) {
  experiment=toString(comparison[i,1])
  reference=toString(comparison[i,2])
  
  ## create a subset object
  selectIndex<-which(group$group %in% c(experiment,reference))
  methobjsubset<-reorganize(methobj, 
                            sample.ids = group$sampleFileShort[selectIndex], 
                            treatment = methobj@treatment[selectIndex])
  
  ## unite
  methsubset<-methylKit::unite(methobjsubset, destrand=FALSE)
  
  ## simply the sample name
  methsubset@sample.ids<-str_split(methsubset@sample.ids,"_",simplify=T)[,2]
  
  subsetDiff<-calculateDiffMeth(methsubset,num.cores=thread)
  subsetDiff25p<-getMethylDiff(subsetDiff,difference=25,qvalue=0.01)
  
  write.table(subsetDiff25p,
              file=paste0(output_dir,experiment,"_vs_",reference,"_DECpG-diff25q001.txt"),
              sep="\t",quote=FALSE, row.names = FALSE)
}

## save the image
save.image(file=paste0("RRBS_Temp_",analysis_date,".RData"), compress="gzip")

