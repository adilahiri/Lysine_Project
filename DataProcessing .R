setwd("~/Desktop/Research/Septi/GSE98455_RAW/TextFile")
source('~/Desktop/Research/Septi/GSE98455_RAW/TextFile/binarized_data.R')
source('~/Desktop/Research/Septi/GSE98455_RAW/TextFile/get_quartile.R')
source('~/Desktop/Research/Septi/GSE98455_RAW/TextFile/get_discrete.R')
source('~/Desktop/Research/Septi/GSE98455_RAW/TextFile/get_compare.R')

library(DESeq2)
library(plotly)
library(dplyr)

Count_df = read.table("Counts.csv",sep=",",check.names = TRUE) # Read the csv file

row.names(Count_df) <- Count_df$V1 # Extract the gene names
Gene_Name <-Count_df[,1, drop=FALSE]
Gene_Name<-t(Gene_Name)
Gene_Name<-Gene_Name[-1]

countdata <- as.matrix(Count_df)
countdata<-countdata[,-1]
countdata<-countdata[-1,]
countdata<-apply(countdata,2,as.numeric)
rownames(countdata)<- Gene_Name

arr_con <- c("con","con","tr","tr")
arr_con<-rep(arr_con,92)
con_tr = factor(arr_con)
coldata <- data.frame(row.names=colnames(countdata), con_tr) # Setup column data for deseq normalization

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~con_tr) # Deseq Analysis
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) # Normalized Counts


Illumina_Names<- c("13101.t06254","13103.t05563","13107.t01862",
                   "13108.t02313","13109.t01024","13103.t04797",
                   "13104.t01554","13104.t04365","13102.t02109",
                   "13103.t01222","13103.t00810","13103.t01671",
                   "13112.t03489","13102.t02141")

MSU_Names <- c("LOC_Os01g70300","LOC_Os03g63330",	"LOC_Os07g20544",	
  "LOC_Os08g25390",	"LOC_Os09g12290","LOC_Os03g55280",
  "LOC_Os04g18200",	"LOC_Os04g48540",	"LOC_Os02g24020",	
  "LOC_Os03g14120",	"LOC_Os03g09910",	"LOC_Os03g18810",	
  "LOC_Os12g37960",	"LOC_Os02g24354")

Index_List_Illumina<-as.array(numeric())
for(iter in 1:length(Illumina_Names)){
  Index_List_Illumina[iter] <- match(Illumina_Names[iter],Gene_Name)
}

PGM_Data <- matrix(nrow=14,ncol=368)
for (iter in 1:14){
  PGM_Data[iter,] <- normalized_counts[Index_List_Illumina[iter],]
}
rownames(PGM_Data) <-Gene_Name[Index_List_Illumina]

control_treatment_matrix <- as.character(coldata$con_tr)
treatment_index <- which(control_treatment_matrix=="tr")
control_index <- which(control_treatment_matrix=="con")

data_treatment <-PGM_Data[,treatment_index]
data_control <-PGM_Data[,control_index]

data_treatment<-t(data_treatment)
data_control<-t(data_control)
colnames(data_control)<- MSU_Names
colnames(data_treatment)<-MSU_Names

data_treatment_disc <-get_discrete(data_treatment)
data_control_disc <- get_discrete(data_control)

Activ_Inhib_Mat_Tr <- matrix(nrow=14,ncol=4)
colnames(Activ_Inhib_Mat_Tr)<-c("Inhibit","Dormant","Active","total")

Activ_Inhib_Mat_Cn <- matrix(nrow=14,ncol=4)
colnames(Activ_Inhib_Mat_Cn)<-c("Inhibit","Dormant","Active","total")

for (iter in 1:14){
  Activ_Inhib_Mat_Tr [iter,1] <- table(data_treatment_disc[,iter])[1]
  Activ_Inhib_Mat_Tr [iter,2] <- table(data_treatment_disc[,iter])[2]
  Activ_Inhib_Mat_Tr [iter,3] <- table(data_treatment_disc[,iter])[3]
  Activ_Inhib_Mat_Tr [iter,4] <- table(data_treatment_disc[,iter])[1]+
    table(data_treatment_disc[,iter])[2] +table(data_treatment_disc[,iter])[3]
  
  Activ_Inhib_Mat_Cn [iter,1] <- table(data_control_disc[,iter])[1]
  Activ_Inhib_Mat_Cn [iter,2] <- table(data_control_disc[,iter])[2]
  Activ_Inhib_Mat_Cn [iter,3] <- table(data_control_disc[,iter])[3]
  Activ_Inhib_Mat_Cn [iter,4] <- table(data_control_disc[,iter])[1]+
    table(data_control_disc[,iter])[2] +table(data_control_disc[,iter])[3]
  
  
  
}

##### Plots ########


Activ_Inhib_Mat_Tr<-as.data.frame(Activ_Inhib_Mat_Tr)
Activ_Inhib_Mat_Tr$GeneNames<-MSU_Names

Activ_Inhib_Mat_Cn<-as.data.frame(Activ_Inhib_Mat_Cn)
Activ_Inhib_Mat_Cn$GeneNames<-MSU_Names

p1<-plot_ly(Activ_Inhib_Mat_Tr,x=~GeneNames,y=~Inhibit,type='bar',
            name='Inhibited')%>%
  add_trace(y=~Dormant, name='Dormant')%>%
  add_trace(y=~Active, name='Activated')%>%
  layout(title="Activation/Inhibition of Genes under Saline Treatment",
         xaxis=list(title='Nodes'),yaxis=list(title='count'),barmode='stack')
show(p1)

p2<-plot_ly(Activ_Inhib_Mat_Cn,x=~GeneNames,y=~Inhibit,type='bar',
            name='Inhibited')%>%
  add_trace(y=~Dormant, name='Dormant')%>%
  add_trace(y=~Active, name='Activated')%>%
  layout(title="Activation/Inhibition of Genes under Control",
         xaxis=list(title='Nodes'),yaxis=list(title='count'),barmode='stack')
show(p2)

