setwd("~/Desktop/Research/Septi/Lysine_Project")
source('~/Desktop/Research/Septi/Lysine_Project/get_discrete_data.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_dirichlet.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_param.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_expected.R')

library(DESeq2)
library(plotly)
library(dplyr)
library(arules)

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


Illumina_Names<- c("13101.t06254","13103.t05563","13107.t01862","13108.t02313","13109.t01024",
                   "13103.t04797",
                   "13104.t01554","13104.t04365",
                   "13102.t02109","13103.t01222",
                   "13103.t00810","13103.t01671",
                   "13112.t03489",
                   "13102.t02141")

MSU_Names <- c("LOC_Os01g70300","LOC_Os03g63330",	"LOC_Os07g20544","LOC_Os08g25390","LOC_Os09g12290",
               "LOC_Os03g55280",
               "LOC_Os04g18200","LOC_Os04g48540",
               "LOC_Os02g24020","LOC_Os03g14120",
               "LOC_Os03g09910","LOC_Os03g18810",
               "LOC_Os12g37960",
               "LOC_Os02g24354")
Topological_Names <- LETTERS[1:14] # A-N Names

Index_List_Illumina<-as.array(numeric())
for(iter in 1:length(Illumina_Names)){
  Index_List_Illumina[iter] <- match(Illumina_Names[iter],Gene_Name)
}

Num_Genes = length(MSU_Names)
PGM_Data <- matrix(nrow=Num_Genes,ncol=368)
for (iter in 1:Num_Genes){
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

meth="cluster"
data_treatment_disc<- get_discrete_data(data_treatment,meth)
data_control_disc <- get_discrete_data(data_control,meth)

### Matrix for Plot P1 and P2
Activ_Inhib_Mat_Tr <- matrix(nrow=Num_Genes,ncol=4)
colnames(Activ_Inhib_Mat_Tr)<-c("Inhibit","Dormant","Active","total")

Activ_Inhib_Mat_Cn <- matrix(nrow=Num_Genes,ncol=4)
colnames(Activ_Inhib_Mat_Cn)<-c("Inhibit","Dormant","Active","total")

for (iter in 1:Num_Genes){
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

### Parameter Estimation Section ####
data_control_disc_topo <- data_control_disc
data_treatment_disc_topo <- data_treatment_disc
colnames(data_treatment_disc_topo)<-Topological_Names
colnames(data_control_disc_topo) <- Topological_Names
dirichlet_params_treatment<- get_dirichlet(data_treatment_disc_topo)
dirichlet_params_control <- get_dirichlet(data_control_disc_topo)

expected_val_treatment <- get_expected(dirichlet_params_treatment)
expected_val_control <- get_expected(dirichlet_params_control)

row.names(expected_val_treatment)<-NULL
row.names(expected_val_control)<-NULL

# Create the CPT Tables
status <- c("inhi","dorm", "act")

# Node A
Treatment_CPT_A <- matrix(c(expected_val_treatment[1,1], expected_val_treatment[1,2],
                            expected_val_treatment[1,3]), 
                          ncol=3, dimnames=list(NULL, status))
Control_CPT_A <- matrix(c(expected_val_control[1,1], expected_val_control[1,2],
                          expected_val_control[1,3]), 
                          ncol=3, dimnames=list(NULL, status))
# Node B
Treatment_CPT_B <- matrix(c(expected_val_treatment[2,1], expected_val_treatment[2,2],
                            expected_val_treatment[2,3]), 
                          ncol=3, dimnames=list(NULL, status))
Control_CPT_B <- matrix(c(expected_val_control[2,1], expected_val_control[2,2],
                          expected_val_control[2,3]), 
                        ncol=3, dimnames=list(NULL, status))
# Node C
Treatment_CPT_C <- matrix(c(expected_val_treatment[3,1], expected_val_treatment[3,2],
                            expected_val_treatment[3,3]), 
                          ncol=3, dimnames=list(NULL, status))
Control_CPT_C <- matrix(c(expected_val_control[3,1], expected_val_control[3,2],
                          expected_val_control[3,3]), 
                        ncol=3, dimnames=list(NULL, status))

# Node D
Treatment_CPT_D <- matrix(c(expected_val_treatment[4,1], expected_val_treatment[4,2],
                            expected_val_treatment[4,3]), 
                          ncol=3, dimnames=list(NULL, status))
Control_CPT_D <- matrix(c(expected_val_control[4,1], expected_val_control[4,2],
                          expected_val_control[4,3]), 
                        ncol=3, dimnames=list(NULL, status))


# Node E
Treatment_CPT_E <- matrix(c(expected_val_treatment[5,1], expected_val_treatment[5,2],
                            expected_val_treatment[5,3]), 
                          ncol=3, dimnames=list(NULL, status))
Control_CPT_E <- matrix(c(expected_val_control[5,1], expected_val_control[5,2],
                          expected_val_control[5,3]), 
                        ncol=3, dimnames=list(NULL, status))

# Node F
Treatment_CPT_F<-c(t(expected_val_treatment[6:248,1:3]))
dim(Treatment_CPT_F)<-c(3,3,3,3,3,3)
dimnames(Treatment_CPT_F)=list("F"=status,"A"=status,"B"=status,"C"=status,"D"=status,"E"=status)

Control_CPT_F<-c(t(expected_val_control[6:248,1:3]))
dim(Control_CPT_F)<-c(3,3,3,3,3,3)
dimnames(Control_CPT_F)=list("F"=status,"A"=status,"B"=status,"C"=status,"D"=status,"E"=status)

# Node G
Treatment_CPT_G<-c(t(expected_val_treatment[249:275,1:3]))
dim(Treatment_CPT_G)<-c(3,3,3,3)
dimnames(Treatment_CPT_G)=list("G"=status,"D"=status,"E"=status,"F"=status)

Control_CPT_G<-c(t(expected_val_control[249:275,1:3]))
dim(Control_CPT_G)<-c(3,3,3,3)
dimnames(Control_CPT_G)=list("G"=status,"D"=status,"E"=status,"F"=status)


# Node H
Treatment_CPT_H<-c(t(expected_val_treatment[276:302,1:3]))
dim(Treatment_CPT_H)<-c(3,3,3,3)
dimnames(Treatment_CPT_H)=list("H"=status,"D"=status,"E"=status,"F"=status)

Control_CPT_H<-c(t(expected_val_control[276:302,1:3]))
dim(Control_CPT_H)<-c(3,3,3,3)
dimnames(Control_CPT_H)=list("H"=status,"D"=status,"E"=status,"F"=status)


# Node I
Treatment_CPT_I<-c(t(expected_val_treatment[303:311,1:3]))
dim(Treatment_CPT_I)<-c(3,3,3)
dimnames(Treatment_CPT_I)=list("I"=status,"G"=status,"H"=status)

Control_CPT_I<-c(t(expected_val_control[303:311,1:3]))
dim(Control_CPT_I)<-c(3,3,3)
dimnames(Control_CPT_I)=list("I"=status,"G"=status,"H"=status)


# Node J
Treatment_CPT_J<-c(t(expected_val_treatment[312:320,1:3]))
dim(Treatment_CPT_J)<-c(3,3,3)
dimnames(Treatment_CPT_J)=list("J"=status,"G"=status,"H"=status)

Control_CPT_J<-c(t(expected_val_control[312:320,1:3]))
dim(Control_CPT_J)<-c(3,3,3)
dimnames(Control_CPT_J)=list("J"=status,"G"=status,"H"=status)

# Node K
Treatment_CPT_K<-c(t(expected_val_treatment[321:329,1:3]))
dim(Treatment_CPT_K)<-c(3,3,3)
dimnames(Treatment_CPT_K)=list("K"=status,"I"=status,"J"=status)

Control_CPT_K<-c(t(expected_val_control[321:329,1:3]))
dim(Control_CPT_K)<-c(3,3,3)
dimnames(Control_CPT_K)=list("K"=status,"I"=status,"J"=status)

#Node L
Treatment_CPT_L<-c(t(expected_val_treatment[330:338,1:3]))
dim(Treatment_CPT_L)<-c(3,3,3)
dimnames(Treatment_CPT_L)=list("L"=status,"I"=status,"J"=status)

Control_CPT_L<-c(t(expected_val_control[330:338,1:3]))
dim(Control_CPT_L)<-c(3,3,3)
dimnames(Control_CPT_L)=list("L"=status,"I"=status,"J"=status)

# Node M
Treatment_CPT_M<-c(t(expected_val_treatment[339:347,1:3]))
dim(Treatment_CPT_M)<-c(3,3,3)
dimnames(Treatment_CPT_M)=list("M"=status,"K"=status,"L"=status)

Control_CPT_M<-c(t(expected_val_control[339:347,1:3]))
dim(Control_CPT_M)<-c(3,3,3)
dimnames(Control_CPT_M)=list("M"=status,"K"=status,"L"=status)


# Node N
Treatment_CPT_N<-c(t(expected_val_treatment[348:350,1:3]))
dim(Treatment_CPT_N)<-c(3,3)
dimnames(Treatment_CPT_N)=list("N"=status,"M"=status)

Control_CPT_N<-c(t(expected_val_control[348:350,1:3]))
dim(Control_CPT_N)<-c(3,3)
dimnames(Control_CPT_N)=list("N"=status,"M"=status)


##### Plots ########
Activ_Inhib_Mat_Tr<- data.frame(Activ_Inhib_Mat_Tr)
Activ_Inhib_Mat_Cn<- data.frame(Activ_Inhib_Mat_Cn)

Activ_Inhib_Mat_Tr$GeneNames<-MSU_Names
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

