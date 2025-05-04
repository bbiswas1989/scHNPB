
# This package has been developed to implement the scHNPB method.
#
# The main function of this package is scHNPB.
# It contains the object data matrix, condition (case/control) and condition_F for adding the column names for NODES method.


###########################################
################scHNPB#####################
###########################################

###### Import these Packages########
library("scDEA")
library(plyr)
library(dplyr)
library(ROSeq)
library(NODES)

#####Main function for scHNPB method#######

scHNPB<-function(scData,condition,condition_F){

#This looping is necessary when a number of genes of a dataset contains all zeros across cell#

  index<-c()
  for (i in 1:dim(scData)[1]){
    index[i]=length(which(as.numeric(scData[i,])!=0))/length(as.numeric(scData[i,]))
  }
  keep1<-which(index>=0.001)
  filterData<-scData[keep1,]

############NODES################

pred.NODES <- Sys.time()
colnames(scData) =condition_F

norm_data<-pQ(scData)
a<-table(colnames(norm_data))
condition1 <- c(rep("C1",as.numeric(a[1])),rep("C2",as.numeric(a[2])))
Res <- NODES(norm_data,condition1,r=50)
dat2<-data.frame(Res)
dat1<-data.frame(scData[,1:2])
library(tibble)
dat1 = rownames_to_column(dat1, "Plants")
dat2 = rownames_to_column(dat2, "Plants")
library(dplyr)
dat = full_join(dat1, dat2, )
dat = dat %>% replace(is.na(.), 0.25)

NODES<-dat$Fisher
#NODES.adj<-p.adjust(NODES,method="BH")
#NODES.adjF<-Pvals_full$NODES.adjF
end.NODES<- Sys.time()
time.NODES<- difftime(end.NODES, pred.NODES, units = "mins")
    cat("Run time for NODES: ", time.NODES, "min","\n")

##########wilcoxon##############
pred.Wilcoxon<- Sys.time()
Pvals_wilcoxon<-scDEA_individual_methods(
  raw.count=scData,
  cell.label=condition,
  is.normalized = FALSE,
  verbose = TRUE,
  BPSC = F,
  DEsingle = F,
  DESeq2 = F,
  edgeR = F,
  MAST = F,
  monocle = F,
  scDD = F,
  Ttest = F,
  Wilcoxon = TRUE,
  limma = F,
  Seurat = F,
  zingeR.edgeR = F,
  BPSC.coef = 2,
  BPSC.normalize = "CPM",
  BPSC.parallel = TRUE,
  DEsingle.parallel = TRUE,
  DEsingle.normalize = "CPM",
  DESeq2.test = "LRT",
  DESeq2.parallel = TRUE,
  DESeq2.beta.prior = TRUE,
  DESeq2.fitType = "parametric",
  DESeq2.normalize = "CPM",
  edgeR.Test = "QLFT",
  edgeR.normalize = "TMM",
  limma.method.fit = "ls",
  limma.trend = TRUE,
  limma.robust = TRUE,
  limma.normalize = "CPM",
  Seurat.normalize = "CPM",
  Seurat.method = "bimod",
  MAST.method = "bayesglm",
  MAST.normalize = "CPM",
  MAST.parallel = TRUE,
  monocle.cores = 1,
  monocle.normalize = "CPM",
  scDD.alpha1 = 0.01,
  scDD.mu0 = 0,
  scDD.s0 = 0.01,
  scDD.a0 = 0.01,
  scDD.b0 = 0.01,
  scDD.normalize = "CPM",
  scDD.permutation = 0,
  Ttest.normalize = "CPM",
  Wilcoxon.normalize = "TMM",
  zingeR.edgeR.normalize = "CPM",
  zingeR.edgeR.maxit.EM = 100
)

Wilcoxon<-Pvals_wilcoxon[,1]
end.Wilcoxon<- Sys.time()
time.Wilcoxon<- difftime(end.Wilcoxon, pred.Wilcoxon, units = "mins")
    cat("Run time for Wilcoxon: ", time.Wilcoxon, "min","\n")

#############ROSeq##################
pred.ROSeq <- Sys.time()
output_full<-ROSeq(countData=filterData, condition = condition, numCores=1)
data_ROSeq_full<-data.frame(output_full)
ROSeq<-data_ROSeq_full$pVals
end.ROSeq<- Sys.time()
time.ROSeq <- difftime(end.ROSeq, pred.ROSeq, units = "mins")
    cat("Run time for ROSeq: ", time.ROSeq, "min","\n")

#############scHNPB############
pred.HNPB<- Sys.time()
Pvals_HNPB_full<-as.matrix(data.frame(ROSeq,NODES,Wilcoxon))
Pvals_HNPB_full= Pvals_HNPB_full%>% replace(is.na(.), 0.20)

#############Lancasters#################
combination.Pvals_HNPB_full<- lancaster.combination(Pvals_HNPB_full, weight = TRUE, trimmed = 0.2)
adjusted.Pvals_HNPB_full<- scDEA.p.adjust(combination.Pvals_HNPB_full, adjusted.method = "BH")
end.HNPB<- Sys.time()
time.HNPB1<- difftime(end.HNPB, pred.HNPB, units = "mins")
time.HNPB<-time.ROSeq+time.NODES+time.Wilcoxon+time.HNPB1
    cat("Run time for HNPB: ", time.HNPB, "min","\n")

data_HNPB_full<-data.frame(Pvals_HNPB_full,combination.Pvals_HNPB_full,adjusted.Pvals_HNPB_full)
return(data_HNPB_full)
}
################Curated scRNA-seq data################
#scDataB<-read.csv("F:\\PhD Folder\\Analysis for P-values Combination Method\\Liver Cancer Data\\GSE77288_Top3000.csv",header=T,row.names=1)

##############Group####################################
#condition<- factor(c(rep("1",201),rep("2",221)))
#condition_F<-c(rep("C1",201),rep("C2",221))
#scData<-scDataB
#Pvals<-scHNPB(scData,condition,condition_F)
