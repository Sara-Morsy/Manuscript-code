#the data was checked first on GEO database ... study tables was downloaded, cleaned and all control group was deleted. 
library(GEOquery)
library(oligo)
library(readr)
library(ggplot2)
library(dplyr)
library(data.table)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")
getGEOSuppFiles("GSE18123")
untar("GSE18123/GSE18123_RAW.tar", exdir = "GSE18123/CEL")
celfilesGSE18123 <- list.files("GSE18123/CEL/GPL570", full = TRUE)


#download characeteristics from the GEO and delete the control group files (did it manually and was reviewed by two members of the team)
GPL570pheno <- read_csv("GPL570pheno.csv")
SDRF <- AnnotatedDataFrame(GPL570pheno)
GPL570 <- read.celfiles(celfilesGSE18123,phenoData = SDRF)

#GPL570 is Affymetrix Human Genome U133 Plus 2.0 Array
#check the data through principle components analysis and boxplots
exp_raw <- log2(Biobase::exprs(GPL570))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = GPL570@phenoData@data[["Dx"]]
)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Disease)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) 

oligo::boxplot(GPL570, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")

#normalization then checking the data again
GPL570_norm <- oligo::rma(GPL570)

exp_GPL570 <- Biobase::exprs(GPL570_norm)
PCA <- prcomp(t(exp_GPL570), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease = 
                       Biobase::pData(GPL570_norm)$Dx)


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Disease)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) 

oligo::boxplot(GPL570_norm, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")

#filtering low intensity genes
GPL570_f <- rowMedians(Biobase::exprs(GPL570_norm))

hist_res <- hist(GPL570_f, 100, col="#e7efd8", freq = FALSE,
	     main = "Histogram of the median intensities",
	     xlab = "Median intensities")

emp_mu <- hist_res$breaks[which.max(hist_res$density)]
emp_sd <- mad(GPL570_f)/2
prop_cental <- 0.50

hist(GPL570_f, 100, col="#e7efd8", freq = FALSE,
	     main = "Histogram of the median intensities",
	     xlab = "Median intensities")+lines(sort(GPL570_f), prop_cental*dnorm(sort(GPL570_f),
				 mean = emp_mu , sd = emp_sd),
				 col = "grey10", lwd = 4)


cut_val <- 0.05 / prop_cental
thresh_median <- qnorm(0.05 / prop_cental, emp_mu, emp_sd)
no_of_samples <- table(paste0(pData(GPL570_norm)$Dx))
no_of_samples
samples_cutoff <- min(no_of_samples)

idx_thresh_median <- apply(exprs(GPL570_norm), 1, function(x){
				   sum(x > thresh_median) >= samples_cutoff})
table(idx_thresh_median)
#There is no needed filtration 
GPL570_filtered <- subset(GPL570_norm, idx_thresh_median)

#annotation to remove gene with NA volume
library(hgu133plus2.db)
anno_GPL570 <- AnnotationDbi::select(hgu133plus2.db,
                                       keys =(featureNames(GPL570_filtered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")
probe_stats <- anno_GPL570 %>%
    group_by(PROBEID) %>%
    summarize(no_of_matches = n_distinct(SYMBOL)) %>%
    filter(no_of_matches > 1)

ids_to_exlude <- ((featureNames(GPL570_filtered) %in% probe_stats$PROBEID) |
		 featureNames(GPL570_filtered)  %in% subset(anno_GPL570 ,
							         is.na(SYMBOL))$PROBEID)
table(ids_to_exlude)

GPL570_final <- subset(GPL570_filtered, !ids_to_exlude)

fData(GPL570_final)$PROBEID <- rownames(fData(GPL570_final))
fData(GPL570_final) <- left_join(fData(GPL570_final), anno_GPL570)
rownames(fData(GPL570_final)) <-fData(GPL570_final)$PROBEID
validObject(GPL570_final)

#batch corection using Combat
library(sva)
#Since the microarray for each sample was performed at different times, this was identified as one of the reasons for the batch effect
GPL570batch<-c("B3","B4","B3","B5","B6","B5","B5","B3","B8","B3","B9","B5","B1","B1","B2","B5","B6","B8","B10","B1","B6","B4","B10","B10","B12","B12","B13","B14","B14","B9","B13","B9","B10","B9","B10","B14","B11","B16","B12","B16")

GPL570_filteredb<-ComBat(dat=exprs(GPL570_final),batch=GPL570batch)
dataset1<-as.data.frame(GPL570_filteredb)

#averaging gene expression values for duplicate genes
dataset1$gene<-GPL570_final@featureData@data[["SYMBOL"]]
dataset1<- dataset1 %>%
           group_by(gene) %>%
           summarise_all(mean)

#ID2 Affymetrix Human Exon 1.0 ST Array (GPL6244)
celfilesg <- list.files("GSE18123/CEL/GPL6244", full = TRUE)
phenoGPL6244 <- read_csv("phenoGPL6244.csv")
SDRF <- AnnotatedDataFrame(phenoGPL6244)
GPL6244 <- read.celfiles(celfilesg,phenoData = SDRF)
#############################################
#Same steps as before
exp_raw <- log2(Biobase::exprs(GPL6244))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = GPL6244@phenoData@data[["ID"]]
)
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Disease)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

oligo::boxplot(GPL6244,  target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")

library(pd.hugene.1.0.st.v1)
GPL6244_norm <- oligo::rma(GPL6244)

exp_GPL6244 <- Biobase::exprs(GPL6244_norm)
PCA <- prcomp(t(exp_GPL6244), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease = 
                       Biobase::pData(GPL6244_norm)$Dx)


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Disease)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) 

oligo::boxplot(GPL6244_norm, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data")

#filtering low intensity genes
GPL6244_f <- rowMedians(Biobase::exprs(GPL6244_norm))

hist_res <- hist(GPL6244_f, 100, col = "cornsilk1", freq = FALSE, 
            main = "Histogram of the median intensities", 
            border = "antiquewhite4",
            xlab = "Median intensities")
emp_mu <- hist_res$breaks[which.max(hist_res$density)]
emp_sd <- BiocGenerics::mad(GPL6244_f)/2
prop_cental <- 0.50
hist(GPL6244_f, 100, col = "cornsilk1", freq = FALSE, 
            main = "Histogram of the median intensities", 
            border = "antiquewhite4",
            xlab = "Median intensities")+lines(sort(GPL6244_f), prop_cental*dnorm(sort(GPL6244_f),
				 mean = emp_mu , sd = emp_sd),
				 col = "grey10", lwd = 4)

#filtering 
cut_val <- 0.05 / prop_cental
thresh_median <- qnorm(0.05 / prop_cental, emp_mu, emp_sd)
no_of_samples <- table(paste0(pData(GPL6244_norm)$Dx))
no_of_samples
samples_cutoff <- min(no_of_samples)

idx_thresh_median <- apply(exprs(GPL6244_norm), 1, function(x){
				   sum(x > thresh_median) >= samples_cutoff})
table(idx_thresh_median)

GPL6244_filtered <- subset(GPL6244_norm, idx_thresh_median)
#annotation to remove genes with NA values
library(hugene10sttranscriptcluster.db)
GSE103965_anno <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                  keys=(featureNames(GPL6244_filtered)),
                                  columns = c("SYMBOL"),
                                  keytype = "PROBEID")
## 'select()' returned 1:many mapping between keys and columns
probe_stats <- GSE103965_anno %>%
    group_by(PROBEID) %>%
    summarize(no_of_matches = n_distinct(SYMBOL)) %>%
    filter(no_of_matches > 1)
ids_to_exlude <- ((featureNames(GPL6244_filtered) %in% probe_stats$PROBEID) |
         featureNames(GPL6244_filtered)  %in% subset(GSE103965_anno ,
                                     is.na(SYMBOL))$PROBEID)
table(ids_to_exlude)
GSE103965_final <- subset(GPL6244_filtered, !ids_to_exlude)

fData(GSE103965_final)$SYMBOL <- rownames(fData(GSE103965_final))
fData(GSE103965_final) <- left_join(fData(GSE103965_final), GSE103965_anno)
## Joining, by = "PROBEID"
rownames(fData(GSE103965_final)) <-fData(GSE103965_final)$SYMBOL

validObject(GSE103965_final)
#Batch correction using Combat
x <- c("B2","B3","B1","B1","B2","B2","B3","B2","B2","B1","B1","B4","B5","B5","B2","B2","B1","B2","B3","B1","B1","B1","B1","B1","B1","B2","B1","B2","B2","B2","B2","B2","B2","B2","B2","B2","B2","B3","B2","B2","B2","B5","B5","B5","B4","B5","B4","B5","B5","B5","B4","B2","B2","B2","B2","B2")
GPL6244_filteredb<-ComBat(dat=exprs(GSE103965_final), batch=x,mod=phenoGPL6244$D)
dataset2<-as.data.frame(GPL6244_filteredb)
dataset2<-setDT(dataset2, keep.rownames = TRUE)[]
dataset2$Gene<-AnnotationDbi::mapIds(hugene10sttranscriptcluster.db,
                                  keys=dataset2$rn,
                                  column = c("SYMBOL"),
                                  keytype = "PROBEID")


full<-merge(dataset2, dataset1, by.x = "Gene", by.y = "gene")
fulldf<-as.data.frame(full)
fulldf<-aggregate(fulldf,list(fulldf$Gene),mean)
fulldf$Gene<-NULL
fulldf$rn<-NULL
row.names(fulldf)<-fulldf$Group.1
fulldf$Group.1<-NULL
final<-t(fulldf)

fullpheno <- read_csv("fullpheno.csv")
rownames(final)<-fullpheno$ID

#This file is the one that will be used for the downstream analysis
#removing samples with structural genetic disorders
remove<-c("GSM650525","GSM650516","GSM650526","GSM650533","GSM650715","GSM650721","GSM650512","GSM650656","GSM650723")
datTraits<-fullpheno[!(fullpheno$ID %in% remove), ]
final2<-final[ !(rownames(final) %in% remove), ]

## Batch effect correction after combining data from two different sources to remove any noise
mod = matrix(1, nrow = dim(final2)[1], ncol = 1)
n.pc = num.sv(t(final2), mod, method="be")
datExpr<-sva_network(t(final2), n.pc)
