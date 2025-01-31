library("knitr")
library("BiocStyle")
#install.packages("ggplot2")
library(ggplot2)
library(dada2)
library(DECIPHER)
library(readr)
library("DECIPHER")

library(phangorn)
library(DiagrammeR)
library(DiagrammeRsvg)
library(Rgraphviz)
library(phyloseq)
library(dendextend)
library(vegan)
library("Matrix")
library("ShortRead")

getwd()
setwd("Desktop/Latest/")
ill_path <- "../Data_1/"
list.files(ill_path)
#Sort  ensures f/r reads are in same order
sortF <- sort(list.files(ill_path,pattern='_R1_001.fastq.gz'))
sortR <- sort(list.files(ill_path,pattern='_R2_001.fastq.gz'))

sampleName <- sapply(strsplit(sortF,"_"),'[',1)
fnF <- file.path(ill_path,sortF)

fnR <- file.path(ill_path,sortR)

filt_path <- file.path(ill_path,"filtered")
if (!file_test('-d', filt_path)) dir.create(filt_path)
filtF <- file.path(filt_path,paste0(sampleName,"_F_filt.fastq.gz"))
print(filtF)
filtR <- file.path(filt_path,paste0(sampleName,"_R_filt.fastq.gz"))
print(filtR)
#Filter f/r reads
out <- filterAndTrim(fnF,filtF,fnR,filtR, truncLen=c(280,280), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE)

#Dereplication reads
derepF <- derepFastq(filtF,verbose=TRUE)
derepR <- derepFastq(filtR,verbose=TRUE)
names(derepF) <- sampleName
names(derepR) <- sampleName
errF <- learnErrors(filtF,multithread=TRUE)
errR <- learnErrors(filtR,multithreads=TRUE)

dadaFs <- dada(derepF, err=errF, multithreads=TRUE)
dadaRs <- dada(derepR, err=errR, multithreads=TRUE)
mergers <- mergePairs(dadaFs,derepF, dadaRs, derepR)
print(mergers)
seqtabAll <- makeSequenceTable(mergers)
seqtabAll.nochim <- removeBimeraDenovo(seqtabAll, verbose = TRUE, multithread = TRUE)
sum(seqtabAll.nochim)/sum(seqtabAll)

getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=sampleName, dada2_input = out[,1], filtered=out[,2], dada_f=sapply(dadaFs,getN), dada_r=sapply(dadaRs, getN), merged=sapply(mergers,getN),
                          nonchim=rowSums(seqtabAll.nochim), final_perc_reads_retained = round(rowSums(seqtabAll.nochim)/out[,1]*100,1))
write.table(summary_tab,"read_count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

## loading reference taxonomy object

asv_seqs <- colnames(seqtabAll.nochim)
asv_headers <- vector(dim(seqtabAll.nochim)[2], mode="character")
sequence <- getSequences(seqtabAll.nochim)
names(sequence) <- gsub(pattern=">", replacement="", x=asv_headers)
alignment <- AlignSeqs(DNAStringSet(sequence), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
treeNJ_rooted=midpoint(treeNJ)
for (i in 1:dim(seqtabAll.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
load("../SILVA_SSU_r138_2019.RData")
dna <- DNAStringSet(getSequences(seqtabAll.nochim))
seqtabAll.nochim<- as.data.frame(seqtabAll.nochim)
seqtabAll.nochim$new<- row.names(seqtabAll.nochim)
seqtabAll.nochim$new <- sapply(strsplit(seqtabAll.nochim$new, "-"), `[`, 3)
seqtabAll.nochim$type <- sapply(strsplit(row.names(seqtabAll.nochim), "-"), `[`, 1)
seqtabAll.nochim<-seqtabAll.nochim[!seqtabAll.nochim$new=="76",]
seqtabAll.nochim$new <- NULL
row.names(seqtabAll.nochim)<-as.character(seqtabAll.nochim$new)
seqtabAll.nochim[68,"new"]="157m"
asvtab <- t(seqtabAll.nochim)
row.names(asvtab) <- sub(">", "", asv_headers)
colnames(asvtab)
write.table(asvtab, "ASVs_counts_new.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(asv_tax, "ASVs_taxonomy_new.tsv", sep = "\t", quote=F, col.names=NA)


sample_info_tab <- read.table(file = "../metadata_p.tsv", row.names = 1, sep="\t", header = T)
sample_info_tab

relative_asv <- asvtab_1/apply(asvtab_1,2,sum)
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(asvtab, taxa_are_rows=T)
vst_count_phy
sample_info_tab_phy <- sample_data(sample_info_tab)
sample_info_tab_phy
sample_names(vst_count_phy)
sample_names(sample_info_tab_phy)
taxa <-tax_table(asv_tax)
vst_physeq <- phyloseq(otu_table(asvtab, taxa_are_rows=T), sample_names(sample_info_tab_phy), tax_table(asv_tax))
rel_phyl <-phyloseq(otu_table(relative_asv, taxa_are_rows=T), sample_info_tab_phy, taxa)
##alpha-richness
richness = estimate_richness(vst_physeq, measures = c("Shannon", "Chao1", "Observed"))
write.table(richness, "alpha_richness.tsv", sep="\t", quote=F, col.names=NA)
pdf(file="alpha_diversity_Shannon.pdf")

#Оставила без изменений
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
plot_ordination(vst_physeq, vst_pcoa, shape="Treatment") + geom_point(size=2, aes(color = Time)) ++scale_color_manual(values=unique(as.vector(sample_info_tab$conditions)))+coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1]))+ ggtitle("PCoA (VST Euc. dist.; N=48)")+
  theme_bw() + theme(legend.position="bottom")


#Alpha diversity
asvtab

pdf(file="alpha_diversity_ordination.pdf")
#ggrare(vst_physeq, step=1000, color='SampleType', label='Sample')
rarecurve(t(asvtab), step=100, col=sample_info_tab$color, lwd=2, ylab='ASVs', label=T)+theme(legend.position="bottom")
#abline(v=(min(rowSums(t(asvtab)))))
dev.off()



vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples


pdf(file="Reads_counts.pdf")
hist(sample_sums(vst_physeq), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=12)
dev.off()

#Phyla
phyla_asv<-asvtab
row.names(phyla_asv)<-asv_tax[,2]
phyla_samp<- t(rowsum(phyla_asv, row.names(phyla_asv)))
write.table(phyla_samp, "Phyla_count.tsv", sep = "\t", quote=F, col.names=NA)
#Class
class_asv<- asvtab
row.names(class_asv)<- asv_tax[,3]
class_asv<- t(rowsum(class_asv, row.names(class_asv)))
write.table(class_asv,'Class_count.tsv',sep='\t', quote=F, col.names = NA)
#Order
order_asv<- asvtab
row.names(order_asv)<-asv_tax[,4]
order_samp<- t(rowsum(order_asv, row.names(order_asv)))
write.table(order_samp, "Order_count.tsv", sep = "\t", quote=F, col.names=NA)
#Family
family_asv<- asvtab
row.names(family_asv)<- asv_tax[,5]
family_samp<- t(rowsum(family_asv, row.names(family_asv)))
write.table(family_samp, "Family_count.tsv", sep = "\t", quote=F, col.names=NA)
#Genus

genus_asv<- asvtab
row.names(genus_asv)<- asv_tax[,6]
genus_samp<- t(rowsum(genus_asv, row.names(genus_asv)))
write.table(genus_samp, "Genus_count.tsv", sep = "\t", quote=F, col.names=NA)

#Heatmap
rarecurve(asvtab, step=100, col=color_c, lwd=2, ylab='ASVs', label=F,cex.main=4,cex.names=4, cex =1, cex.lab=3,cex.axis=2.5,main=sprintf("Rarefaction curve colored by %s",cond))+theme(legend.position="bottom",axis.text.y=element_text(size =4, hjust=1))
rarecurve(vst_physeq, step=100, lwd=2, ylab='ASVs', label=F,cex.main=4,cex.names=4, cex =1, cex.lab=3,cex.axis=2.5,main=sprintf("Rarefaction curve"))+theme(legend.position="bottom",axis.text.y=element_text(size =4, hjust=1))

#ps<-phyloseq(vst_count_phy,  sample_info_tab_phy, tax_tab_phy)
rarecurve(asvtab, step=100, col=color_c, lwd=2, ylab='ASVs', label=F,cex.main=4,cex.names=4, cex =1, cex.lab=3,cex.axis=2.5,main=sprintf("Rarefaction curve colored by %s",cond))+theme(legend.position="bottom",axis.text.y=element_text(size =4, hjust=1))