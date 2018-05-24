# Prioritize features picked using MZ Mine 
# P. Vervliet 09/2017
# Modified for upload to GitHub 05/2018

#' Script was developed based on the data obtained in the in vitro metabolism study assay as
#' described by Mortelé O., Vervliet P. et al (doi: 10.1016/j.jpba.2018.02.32).
#' An overview of the samples in the assay can be found in the supplied CSV file (HLM_samples.csv)
#' A visual representation can be found in the supplied image (HLMassay.png)
#' 
#' @author Philippe Vervliet (ORCID: 0000-0003-2644-6820)
#' 
#' @references 
#' Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” 
#' Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

##Load libraries
library("tidyverse")
library("DESeq2")
library("ggrepel")
library("gridExtra")

##Load data
PeaksOrig <- read.csv("Sample-InVitroMassList R.csv", sep = ",", header = T, row.names = 1)
ExpDataOrig <- read.csv("HLM_samples.csv", sep = ",", header = T, row.names = 1)


#Make a copy of the data to work with and perform necessary transformations
ExpData <- ExpDataOrig
ExpData <- ExpData %>% 
  mutate(rownames = rownames(ExpDataOrig)) %>% 
  mutate(Phase.I = as.factor(ExpData$Phase.I)) %>% 
  mutate(Phase.II = as.factor(ExpData$Phase.II)) %>%
  mutate(SolvBlank = as.factor(ExpData$SolvBlank)) %>%
  mutate(CompPresent = as.factor(ExpData$CompPresent)) %>%
  mutate(CoFPresent = as.factor(ExpData$CoFPresent)) %>%
  mutate(HLMPres = as.factor(ExpData$HLMPres)) %>% 
  column_to_rownames(var = "rownames")

Peaks <- PeaksOrig
Peaks <- Peaks %>% 
  select(3:30) %>% 
  mutate(Sample.0 = as.integer(Peaks$Sample.0)) %>% 
  mutate(Sample.1 = as.integer(Peaks$Sample.1)) %>% 
  mutate(Sample.2 = as.integer(Peaks$Sample.2)) %>% 
  mutate(Sample.3 = as.integer(Peaks$Sample.3)) %>% 
  mutate(Sample.4 = as.integer(Peaks$Sample.4)) %>% 
  mutate(Sample.5 = as.integer(Peaks$Sample.5)) %>% 
  mutate(Sample.6 = as.integer(Peaks$Sample.6)) %>% 
  mutate(Sample.7 = as.integer(Peaks$Sample.7)) %>% 
  mutate(Sample.8 = as.integer(Peaks$Sample.8)) %>% 
  mutate(Sample.9 = as.integer(Peaks$Sample.9)) %>% 
  mutate(Sample.10 = as.integer(Peaks$Sample.10)) %>% 
  mutate(Sample.11 = as.integer(Peaks$Sample.11)) %>% 
  mutate(Sample.12 = as.integer(Peaks$Sample.12)) %>% 
  mutate(Sample.13 = as.integer(Peaks$Sample.13)) %>% 
  mutate(Sample.14 = as.integer(Peaks$Sample.14)) %>% 
  mutate(Sample.15 = as.integer(Peaks$Sample.15)) %>% 
  mutate(Sample.15a = as.integer(Peaks$Sample.15a)) %>% 
  mutate(Sample.16 = as.integer(Peaks$Sample.16)) %>% 
  mutate(Sample.17 = as.integer(Peaks$Sample.17)) %>% 
  mutate(Sample.18 = as.integer(Peaks$Sample.18)) %>% 
  mutate(Sample.19 = as.integer(Peaks$Sample.19)) %>% 
  mutate(Sample.20 = as.integer(Peaks$Sample.20)) %>% 
  mutate(Sample.21 = as.integer(Peaks$Sample.21)) %>% 
  mutate(Sample.21a = as.integer(Peaks$Sample.21a)) %>% 
  mutate(Blank.1 = as.integer(Peaks$Blank.1)) %>% 
  mutate(Blank.2 = as.integer(Peaks$Blank.2)) %>% 
  mutate(Blank.3 = as.integer(Peaks$Blank.3)) %>% 
  mutate(Blank.4 = as.integer(Peaks$Blank.4)) %>% 
  mutate(Blank.5 = as.integer(Peaks$Blank.5)) %>% 
  mutate(mzvalue = rownames(PeaksOrig)) %>% 
  column_to_rownames(var = "mzvalue")

##PCA FULL DATASET
#Create DESeqDataSet object structure
PeaksDds <- DESeqDataSetFromMatrix(countData = Peaks, colData = ExpData,
                              design = ~ ColourCode)

#Execute DESeq 2. This will create a new data object, called dds here, that contains the fitted
#statistical model. This will have performed the following steps; normalizing
#the read counts for library sizes, estimating the variance for each gene included
#in the analysis, calculating the different log-fold change values for all included
#factors
PeaksDds <- DESeq(PeaksDds)

#PCA
PeaksSE1.2 <- SummarizedExperiment(log2(counts(PeaksDds, normalized=TRUE) + 1),
                           colData=colData(PeaksDds))

pcase<-plotPCA(DESeqTransform(PeaksSE1.2), intgroup=c("ColourCode")
               ,ntop=10000,returnData=TRUE)

ggplot(pcase,aes(x=PC1,y=PC2,label=rownames(ExpData))) +
  geom_point(aes(col=ExpData$ColourCode),alpha=1, size=2) +
  geom_text_repel(aes(colour=ExpData$ColourCode), nudge_y = 20, nudge_x = 10) +
  theme_minimal() +
  labs(title="PCA - All Samples" , x="PC1", y="PC2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")


##PCA ONLY SAMPLES
#Create subsets
PeaksSamples <- subset(Peaks, select=c(5:7,8:10,13:15,20:22))
ExpDataSamples <- subset(ExpData, ExpData$SampleType=="Sample")

#Create DESeqDataSet object structure
PeaksSamplesDds <- DESeqDataSetFromMatrix(countData = PeaksSamples, colData = ExpDataSamples,
                                   design = ~ ColourCode)

#Execute DESeq 2
PeaksSamplesDds <- DESeq(PeaksSamplesDds)

#PCA
PeaksSamplesSE1.2 <- SummarizedExperiment(log2(counts(PeaksSamplesDds, normalized=TRUE) + 1),
                                   colData=colData(PeaksSamplesDds))

pcasesamples<-plotPCA(DESeqTransform(PeaksSamplesSE1.2), intgroup=c("ColourCode")
               ,ntop=10000,returnData=TRUE)

ggplot(pcasesamples,aes(x=PC1,y=PC2,label=rownames(ExpDataSamples))) +
  geom_point(aes(col=ExpDataSamples$ColourCode),alpha=1, size=2) +
  geom_text_repel(aes(colour=ExpDataSamples$ColourCode), nudge_y = 20, nudge_x = 10) +
  theme_minimal() +
  labs(title="PCA - Replicate Samples" , x="PC1", y="PC2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")


###VOLCANO PLOTS

##Ph1.1h
Ph1.1h <- subset(PeaksOrig, select=c(4:5,7:9,1:2)) %>% 
  mutate(MeanCtrl = rowMeans(PeaksOrig[c("Sample.1", "Sample.2")])) %>% 
  mutate(MeanSmpl = rowMeans(PeaksOrig[c("Sample.4", "Sample.5", "Sample.6")])) %>% 
  mutate(FC = MeanSmpl/MeanCtrl)
#Calculate p-value, unpaired student T-test
t.result1 <- apply(Ph1.1h[,1:5], 1, function (x) t.test(x[1:2],x[3:5],paired=FALSE))
Ph1.1h$p_value <- unlist(lapply(t.result1, function(x) x$p.value))
#Add rownames
Ph1.1h$mzvalue <- rownames(Peaks)
#Log transformations fold change & p-value
Ph1.1h$negLog10P <- -log10(Ph1.1h$p_value)
Ph1.1h$Log10FC <- log10(Ph1.1h$FC)
#Signal significant p-value
Ph1.1h <- mutate(Ph1.1h, sigP=ifelse(Ph1.1h$p_value<0.05, 1, 0))
Ph1.1h <- mutate(Ph1.1h, sigFC= ifelse(abs(Ph1.1h$Log10FC)>1, 1, 0))
Ph1.1h <- mutate(Ph1.1h, Legend= ifelse(c(Ph1.1h$sigP>0.5 & Ph1.1h$sigFC>0.5), "p < 0.05, log10FC > 1", 
                                     ifelse(c(Ph1.1h$sigP>0.5 & Ph1.1h$sigFC<0.5),"p < 0.05, log10FC < 1",
                                            ifelse(c(Ph1.1h$sigP<0.5 & Ph1.1h$sigFC>0.5),"p > 0.05, log10FC > 1","p > 0.05, log10FC < 1"))))

#Subset of significant data
Ph1.1h.subset <- Ph1.1h %>% 
  filter(sigP > 0.5) %>% 
  select(12,7,6,1:5,10,11,17)

write.csv(Ph1.1h.subset, file = "Significant-Ph1-1hbis.csv",row.names=FALSE)
  
##Ph1.3h
Ph1.3h <- subset(PeaksOrig, select=c(4:5,10:12,1:2)) %>% 
  mutate(MeanCtrl = rowMeans(PeaksOrig[c("Sample.1", "Sample.2")])) %>% 
  mutate(MeanSmpl = rowMeans(PeaksOrig[c("Sample.7", "Sample.8", "Sample.9")])) %>% 
  mutate(FC = MeanSmpl/MeanCtrl)
#Calculate p-value, unpaired student T-test
t.result2 <- apply(Ph1.3h[,1:5], 1, function (x) t.test(x[1:2],x[3:5],paired=FALSE))
Ph1.3h$p_value <- unlist(lapply(t.result2, function(x) x$p.value))
#Add rownames
Ph1.3h$mzvalue <- rownames(Peaks)
#Log transformations fold change & p-value
Ph1.3h$negLog10P <- -log10(Ph1.3h$p_value)
Ph1.3h$Log10FC <- log10(Ph1.3h$FC)
#Signal significant p-value
Ph1.3h <- mutate(Ph1.3h, sigP=ifelse(Ph1.3h$p_value<0.05, 1, 0))
Ph1.3h <- mutate(Ph1.3h, sigFC= ifelse(abs(Ph1.3h$Log10FC)>1, 1, 0))
Ph1.3h <- mutate(Ph1.3h, Legend= ifelse(c(Ph1.3h$sigP>0.5 & Ph1.3h$sigFC>0.5), "p < 0.05, log10FC > 1", 
                                        ifelse(c(Ph1.3h$sigP>0.5 & Ph1.3h$sigFC<0.5),"p < 0.05, log10FC < 1",
                                               ifelse(c(Ph1.3h$sigP<0.5 & Ph1.3h$sigFC>0.5),"p > 0.05, log10FC > 1","p > 0.05, log10FC < 1"))))
#Subset of significant data
Ph1.3h.subset <- Ph1.3h %>% 
  filter(sigP > 0.5) %>% 
  select(12,7,6,1:5,10,11,17)

write.csv(Ph1.3h.subset, file = "Significant-Ph1-3hbis.csv",row.names=FALSE)

##Ph2.Gluc
Ph2.Gluc <- subset(PeaksOrig, select=c(13:17,1:2)) %>% 
  mutate(MeanCtrl = rowMeans(PeaksOrig[c("Sample.10", "Sample.11")])) %>% 
  mutate(MeanSmpl = rowMeans(PeaksOrig[c("Sample.12", "Sample.13", "Sample.14")])) %>% 
  mutate(FC = MeanSmpl/MeanCtrl)
#Calculate p-value, unpaired student T-test
t.result3 <- apply(Ph2.Gluc[,1:5], 1, function (x) t.test(x[1:2],x[3:5],paired=FALSE))
Ph2.Gluc$p_value <- unlist(lapply(t.result3, function(x) x$p.value))
#Add rownames
Ph2.Gluc$mzvalue <- rownames(Peaks)
#Log transformations fold change & p-value
Ph2.Gluc$negLog10P <- -log10(Ph2.Gluc$p_value)
Ph2.Gluc$Log10FC <- log10(Ph2.Gluc$FC)
#Signal significant p-value
Ph2.Gluc <- mutate(Ph2.Gluc, sigP=ifelse(Ph2.Gluc$p_value<0.05, 1, 0))
Ph2.Gluc <- mutate(Ph2.Gluc, sigFC= ifelse(abs(Ph2.Gluc$Log10FC)>1, 1, 0))
Ph2.Gluc <- mutate(Ph2.Gluc, Legend= ifelse(c(Ph2.Gluc$sigP>0.5 & Ph2.Gluc$sigFC>0.5), "p < 0.05, log10FC > 1", 
                                        ifelse(c(Ph2.Gluc$sigP>0.5 & Ph2.Gluc$sigFC<0.5),"p < 0.05, log10FC < 1",
                                               ifelse(c(Ph2.Gluc$sigP<0.5 & Ph2.Gluc$sigFC>0.5),"p > 0.05, log10FC > 1","p > 0.05, log10FC < 1"))))
#Subset of significant data
Ph2.Gluc.subset <- Ph2.Gluc %>% 
  filter(sigP > 0.5) %>% 
  select(12,7,6,1:5,10,11,17)

write.csv(Ph2.Gluc.subset, file = "Significant-Ph2-Glucbis.csv",row.names=FALSE)

##Ph2.Sulf
Ph2.Sulf <- subset(PeaksOrig, select=c(20:24,1:2)) %>% 
  mutate(MeanCtrl = rowMeans(PeaksOrig[c("Sample.16", "Sample.17")])) %>% 
  mutate(MeanSmpl = rowMeans(PeaksOrig[c("Sample.18", "Sample.19", "Sample.20")])) %>% 
  mutate(FC = MeanSmpl/MeanCtrl)
#Calculate p-value, unpaired student T-test
t.result4 <- apply(Ph2.Sulf[,1:5], 1, function (x) t.test(x[1:2],x[3:5],paired=FALSE))
Ph2.Sulf$p_value <- unlist(lapply(t.result4, function(x) x$p.value))
#Add rownames
Ph2.Sulf$mzvalue <- rownames(Peaks)
#Log transformations fold change & p-value
Ph2.Sulf$negLog10P <- -log10(Ph2.Sulf$p_value)
Ph2.Sulf$Log10FC <- log10(Ph2.Sulf$FC)
#Signal significant p-value
Ph2.Sulf <- mutate(Ph2.Sulf, sigP=ifelse(Ph2.Sulf$p_value<0.05, 1, 0))
Ph2.Sulf <- mutate(Ph2.Sulf, sigFC= ifelse(abs(Ph2.Sulf$Log10FC)>1, 1, 0))
Ph2.Sulf <- mutate(Ph2.Sulf, Legend= ifelse(c(Ph2.Sulf$sigP>0.5 & Ph2.Sulf$sigFC>0.5), "p < 0.05, log10FC > 1", 
                                            ifelse(c(Ph2.Sulf$sigP>0.5 & Ph2.Sulf$sigFC<0.5),"p < 0.05, log10FC < 1",
                                                   ifelse(c(Ph2.Sulf$sigP<0.5 & Ph2.Sulf$sigFC>0.5),"p > 0.05, log10FC > 1","p > 0.05, log10FC < 1"))))
#Subset of significant data
Ph2.Sulf.subset <- Ph2.Sulf %>% 
  filter(sigP > 0.5) %>% 
  select(12,7,6,1:5,10,11,17)

write.csv(Ph2.Sulf.subset, file = "Significant-Ph2-Sulfbis.csv",row.names=FALSE)

##Volcano Plot Colors
VolcColors <- c("p < 0.05, log10FC < 1", "p < 0.05, log10FC > 1", "p > 0.05, log10FC < 1", "p > 0.05, log10FC > 1")
VolcColors.col <- c("#BDB69C", "#67C5AB","#717D8C", "#41AAC4")
names(VolcColors.col) <- VolcColors

##Volcano Plot
V1 <- ggplot(Ph1.1h, aes(x=Ph1.1h$Log10FC, y=Ph1.1h$negLog10P)) +
  geom_point(aes(col=Legend),alpha=0.7, size=2) +
  scale_color_manual(values = VolcColors.col) +
  theme_minimal() +
  labs(title="Phase 1 - 1 hour" , x="Log10 Fold Change", y="-Log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "bottom")
V2 <- ggplot(Ph1.3h, aes(x=Ph1.3h$Log10FC, y=Ph1.3h$negLog10P)) +
  geom_point(aes(col=Legend),alpha=0.7, size=2) +
  scale_color_manual(values = VolcColors.col) +
  theme_minimal() +
  labs(title="Phase 1 - 3 hour" , x="Log10 Fold Change", y="-Log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "bottom")
V3 <- ggplot(Ph2.Gluc, aes(x=Ph2.Gluc$Log10FC, y=Ph2.Gluc$negLog10P)) +
  geom_point(aes(col=Legend),alpha=0.7, size=2) +
  scale_color_manual(values = VolcColors.col) +
  theme_minimal() +
  labs(title="Phase 2 - Gluc" , x="Log10 Fold Change", y="-Log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "bottom")
V4 <- ggplot(Ph2.Sulf, aes(x=Ph2.Sulf$Log10FC, y=Ph2.Sulf$negLog10P)) +
  geom_point(aes(col=Legend),alpha=0.7, size=2) +
  scale_color_manual(values = VolcColors.col) +
  theme_minimal() +
  labs(title="Phase 2 - Sulf" , x="Log10 Fold Change", y="-Log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "bottom")

#Plot all separately
V1

V2

V3

V4

#Plot all in grid
grid.arrange(V1, V2, V3, V4, ncol=2, padding=unit(10, "line"))


