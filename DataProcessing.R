setwd("~/Desktop/Research/Septi/Lysine_Project")
source('~/Desktop/Research/Septi/Lysine_Project/get_discrete_data.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_dirichlet.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_param.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_expected.R')
source('~/Desktop/Research/Septi/Lysine_Project/get_nlargest.R')

library(DESeq2)
library(plotly)
library(dplyr)
library(arules)
library(bnlearn)

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

#detach("package:bnlearn", unload=TRUE) # as bnlearn has function with same name 

meth="cluster" ## interval or cluster
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


######### Define and plot the graph structure (DAG)
#library(bnlearn)
network_model<-model2network("[A][B][C][D][E][F|A:B:C:D:E][G|D:E:F][H|D:E:F][I|G:H][J|G:H][K|I:J][L|I:J][M|K:L][N|M]")





dfit_treatment <- custom.fit(network_model, dist=list(A=Treatment_CPT_A,B=Treatment_CPT_B,
                                                      C=Treatment_CPT_C,D=Treatment_CPT_D,
                                                      E=Treatment_CPT_E,F=Treatment_CPT_F,
                                                      G=Treatment_CPT_G,H= Treatment_CPT_H,
                                                      I=Treatment_CPT_I,J= Treatment_CPT_J,
                                                      K=Treatment_CPT_K, L=Treatment_CPT_L,
                                                      M=Treatment_CPT_M,N=Treatment_CPT_N))



dfit_control <- custom.fit(network_model, dist=list(A=Control_CPT_A,B=Control_CPT_B,
                                                      C=Control_CPT_C,D=Control_CPT_D,
                                                      E=Control_CPT_E,F=Control_CPT_F,
                                                      G=Control_CPT_G,H= Control_CPT_H,
                                                      I=Control_CPT_I,J= Control_CPT_J,
                                                      K=Control_CPT_K, L=Control_CPT_L,
                                                      M=Control_CPT_M,N=Control_CPT_N))
# Likelihood weighting for treatment
k=600000
Treatment_Matrix_Inhibit<-matrix(nrow=13,ncol=1)
Treatment_Matrix_Dormant<-matrix(nrow=13,ncol=1)
Treatment_Matrix_Active<-matrix(nrow=13,ncol=1)

colnames(Treatment_Matrix_Inhibit)<- c("N")
colnames(Treatment_Matrix_Dormant)<- c("N")
colnames(Treatment_Matrix_Active)<- c("N")

remove<-c("N")
Interevention_List<-Topological_Names
Interevention_List<-Interevention_List[! Interevention_List %in% remove]

rownames(Treatment_Matrix_Inhibit)<-Interevention_List
rownames(Treatment_Matrix_Dormant)<-Interevention_List
rownames(Treatment_Matrix_Active)<-Interevention_List

iter_N=3 # 1 for inhibition 2 for dormancy 3 for activation

set.seed(4)
Treatment_Matrix_Inhibit[1,1]=cpquery(dfit_treatment,
                                       (N==status[iter_N]),
                                       list(A="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[1,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(A="dorm"),method='lw',n=k)
Treatment_Matrix_Active[1,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(A="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[2,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(B="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[2,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(B="dorm"),method='lw',n=k)
Treatment_Matrix_Active[2,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(B="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[3,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(C="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[3,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(C="dorm"),method='lw',n=k)
Treatment_Matrix_Active[3,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(C="act"),method='lw',n=k)



Treatment_Matrix_Inhibit[4,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(D="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[4,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(D="dorm"),method='lw',n=k)
Treatment_Matrix_Active[4,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(D="act"),method='lw',n=k)



Treatment_Matrix_Inhibit[5,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(E="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[5,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(E="dorm"),method='lw',n=k)
Treatment_Matrix_Active[5,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(E="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[6,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(F="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[6,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(F="dorm"),method='lw',n=k)
Treatment_Matrix_Active[6,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(F="act"),method='lw',n=k)

Treatment_Matrix_Inhibit[7,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(G="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[7,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(G="dorm"),method='lw',n=k)
Treatment_Matrix_Active[7,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(G="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[8,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(H="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[8,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(H="dorm"),method='lw',n=k)
Treatment_Matrix_Active[8,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(H="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[9,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(I="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[9,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(I="dorm"),method='lw',n=k)
Treatment_Matrix_Active[9,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(I="act"),method='lw',n=k)

Treatment_Matrix_Inhibit[10,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(J="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[10,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(J="dorm"),method='lw',n=k)
Treatment_Matrix_Active[10,1]=cpquery(dfit_treatment,
                                     (N==status[iter_N]),
                                     list(J="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[11,1]=cpquery(dfit_treatment,
                                       (N==status[iter_N]),
                                       list(K="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[11,1]=cpquery(dfit_treatment,
                                       (N==status[iter_N]),
                                       list(K="dorm"),method='lw',n=k)
Treatment_Matrix_Active[11,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(K="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[12,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(L="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[12,1]=cpquery(dfit_treatment,
                                       (N==status[iter_N]),
                                       list(L="dorm"),method='lw',n=k)
Treatment_Matrix_Active[12,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(L="act"),method='lw',n=k)


Treatment_Matrix_Inhibit[13,1]=cpquery(dfit_treatment,
                                       (N==status[iter_N]),
                                       list(M="inhi"),method='lw',n=k)
Treatment_Matrix_Dormant[13,1]=cpquery(dfit_treatment,
                                       (N==status[iter_N]),
                                       list(M="dorm"),method='lw',n=k)
Treatment_Matrix_Active[13,1]=cpquery(dfit_treatment,
                                      (N==status[iter_N]),
                                      list(M="act"),method='lw',n=k)




# Likelihood weighting for control

Control_Matrix_Inhibit<-matrix(nrow=13,ncol=1)
Control_Matrix_Dormant<-matrix(nrow=13,ncol=1)
Control_Matrix_Active<-matrix(nrow=13,ncol=1)

colnames(Control_Matrix_Inhibit)<- c("N")
colnames(Control_Matrix_Dormant)<- c("N")
colnames(Control_Matrix_Active)<- c("N")


rownames(Control_Matrix_Inhibit)<-Interevention_List
rownames(Control_Matrix_Dormant)<-Interevention_List
rownames(Control_Matrix_Active)<-Interevention_List

iter_N=3 # 1 for inhibition 2 for dormancy 3 for activation

set.seed(4)
Control_Matrix_Inhibit[1,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(A="inhi"),method='lw',n=k)
Control_Matrix_Dormant[1,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(A="dorm"),method='lw',n=k)
Control_Matrix_Active[1,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(A="act"),method='lw',n=k)


Control_Matrix_Inhibit[2,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(B="inhi"),method='lw',n=k)
Control_Matrix_Dormant[2,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(B="dorm"),method='lw',n=k)
Control_Matrix_Active[2,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(B="act"),method='lw',n=k)


Control_Matrix_Inhibit[3,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(C="inhi"),method='lw',n=k)
Control_Matrix_Dormant[3,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(C="dorm"),method='lw',n=k)
Control_Matrix_Active[3,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(C="act"),method='lw',n=k)



Control_Matrix_Inhibit[4,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(D="inhi"),method='lw',n=k)
Control_Matrix_Dormant[4,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(D="dorm"),method='lw',n=k)
Control_Matrix_Active[4,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(D="act"),method='lw',n=k)



Control_Matrix_Inhibit[5,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(E="inhi"),method='lw',n=k)
Control_Matrix_Dormant[5,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(E="dorm"),method='lw',n=k)
Control_Matrix_Active[5,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(E="act"),method='lw',n=k)


Control_Matrix_Inhibit[6,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(F="inhi"),method='lw',n=k)
Control_Matrix_Dormant[6,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(F="dorm"),method='lw',n=k)
Control_Matrix_Active[6,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(F="act"),method='lw',n=k)

Control_Matrix_Inhibit[7,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(G="inhi"),method='lw',n=k)
Control_Matrix_Dormant[7,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(G="dorm"),method='lw',n=k)
Control_Matrix_Active[7,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(G="act"),method='lw',n=k)


Control_Matrix_Inhibit[8,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(H="inhi"),method='lw',n=k)
Control_Matrix_Dormant[8,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(H="dorm"),method='lw',n=k)
Control_Matrix_Active[8,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(H="act"),method='lw',n=k)


Control_Matrix_Inhibit[9,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(I="inhi"),method='lw',n=k)
Control_Matrix_Dormant[9,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(I="dorm"),method='lw',n=k)
Control_Matrix_Active[9,1]=cpquery(dfit_control,
                                     (N==status[iter_N]),
                                     list(I="act"),method='lw',n=k)

Control_Matrix_Inhibit[10,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(J="inhi"),method='lw',n=k)
Control_Matrix_Dormant[10,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(J="dorm"),method='lw',n=k)
Control_Matrix_Active[10,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(J="act"),method='lw',n=k)


Control_Matrix_Inhibit[11,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(K="inhi"),method='lw',n=k)
Control_Matrix_Dormant[11,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(K="dorm"),method='lw',n=k)
Control_Matrix_Active[11,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(K="act"),method='lw',n=k)


Control_Matrix_Inhibit[12,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(L="inhi"),method='lw',n=k)
Control_Matrix_Dormant[12,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(L="dorm"),method='lw',n=k)
Control_Matrix_Active[12,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(L="act"),method='lw',n=k)


Control_Matrix_Inhibit[13,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(M="inhi"),method='lw',n=k)
Control_Matrix_Dormant[13,1]=cpquery(dfit_control,
                                       (N==status[iter_N]),
                                       list(M="dorm"),method='lw',n=k)
Control_Matrix_Active[13,1]=cpquery(dfit_control,
                                      (N==status[iter_N]),
                                      list(M="act"),method='lw',n=k)
####################### TWO COMBINATIONS ##########################################################
Intercombo <-combn(Interevention_List, 2)

Treatment_Combo <- matrix (nrow=length(Intercombo)/2,ncol=3^2)
Control_Combo <- matrix (nrow=length(Intercombo)/2,ncol=3^2)

names_row<-NULL


for (iter in 1:78){
  names_row[iter]= paste(Intercombo[1,iter],"+",Intercombo[2,iter])
 }

rownames(Treatment_Combo)<-names_row
rownames(Control_Combo) <- names_row

col_assign <- expand.grid(c(-1,0,1),c(-1,0,1))

names_col <- c("-1-1","0-1","1-1","-10", "00","10","-11","01","11")

colnames(Treatment_Combo)<-names_col
colnames(Control_Combo)<-names_col

for (iter in 1: 78){
  set.seed(4)
  
  Treatment_Combo[iter,1] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="inhi")
                                   & eval(parse(text=Intercombo[2,iter]))=="inhi",n=k)
  
  Treatment_Combo[iter,2] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="dorm")
                                   & eval(parse(text=Intercombo[2,iter]))=="inhi",n=k)
  
  Treatment_Combo[iter,3] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="act")
                                   & eval(parse(text=Intercombo[2,iter]))=="inhi",n=k)
  
  Treatment_Combo[iter,4] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="inhi")
                                   & eval(parse(text=Intercombo[2,iter]))=="dorm",n=k)
  
  Treatment_Combo[iter,5] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="dorm")
                                   & eval(parse(text=Intercombo[2,iter]))=="dorm",n=k)
  
  Treatment_Combo[iter,6] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="act")
                                   & eval(parse(text=Intercombo[2,iter]))=="dorm",n=k)
  
  Treatment_Combo[iter,7] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="inhi")
                                   & eval(parse(text=Intercombo[2,iter]))=="act",n=k)
  
  Treatment_Combo[iter,8] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="dorm")
                                   & eval(parse(text=Intercombo[2,iter]))=="act",n=k)
  
  Treatment_Combo[iter,9] <-cpquery(dfit_treatment, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="act")
                                   & eval(parse(text=Intercombo[2,iter]))=="act",n=k)
  
  
  
  Control_Combo[iter,1] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="inhi")
                                   & eval(parse(text=Intercombo[2,iter]))=="inhi",n=k)
  
  Control_Combo[iter,2] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="dorm")
                                   & eval(parse(text=Intercombo[2,iter]))=="inhi",n=k)
  
  Control_Combo[iter,3] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="act")
                                   & eval(parse(text=Intercombo[2,iter]))=="inhi",n=k)
  
  Control_Combo[iter,4] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="inhi")
                                   & eval(parse(text=Intercombo[2,iter]))=="dorm",n=k)
  
  Control_Combo[iter,5] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="dorm")
                                   & eval(parse(text=Intercombo[2,iter]))=="dorm",n=k)
  
  Control_Combo[iter,6] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="act")
                                   & eval(parse(text=Intercombo[2,iter]))=="dorm",n=k)
  
  Control_Combo[iter,7] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="inhi")
                                   & eval(parse(text=Intercombo[2,iter]))=="act",n=k)
  
  Control_Combo[iter,8] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="dorm")
                                   & eval(parse(text=Intercombo[2,iter]))=="act",n=k)
  
  Control_Combo[iter,9] <-cpquery(dfit_control, (N=="act"),
                                   evidence = c(eval(parse(text=Intercombo[1,iter]))=="act")
                                   & eval(parse(text=Intercombo[2,iter]))=="act",n=k)

}

Treatment_Top_Five <- get_nlargest(Treatment_Combo,n=5)
Y_Treatment <- Treatment_Top_Five$values
X_Treatment <- names_col[Treatment_Top_Five$V2]
Treatment_Top_Five$X_Names<- paste(Treatment_Top_Five$Names,X_Treatment)

Control_Top_Five <- get_nlargest(Control_Combo,n=5)
Y_Control <- Control_Top_Five$values
X_Control <- names_col[Control_Top_Five$V2]
Control_Top_Five$X_Names<- paste(Control_Top_Five$Names,X_Control)


# Plot the results


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


par(cex=2)
show(graphviz.plot(network_model))


Treatment_Results <- cbind(Treatment_Matrix_Inhibit,Treatment_Matrix_Dormant,Treatment_Matrix_Active)
colnames(Treatment_Results)<-c("inhibition","dormancy","activation")
Treatment_Results<-as.data.frame(Treatment_Results)
Treatment_Results$Names <-Interevention_List

Control_Results <- cbind(Control_Matrix_Inhibit,Control_Matrix_Dormant,Control_Matrix_Active)
colnames(Control_Results)<-c("inhibition","dormancy","activation")
Control_Results<-as.data.frame(Control_Results)
Control_Results$Names <-Interevention_List

p3 <- Control_Results %>% plot_ly()
p3 <- p3 %>% add_trace(x = ~Names, y = ~inhibition, type = 'bar',name="inhibition",
                         text = Control_Results$inhibition, textposition = 'auto',
                         marker = list(color = '#1f77b4',
                                       line = list(color = '#1f77b4', width = 1.5)))
p3 <- p3 %>% add_trace(x = ~Names, y = ~dormancy, type = 'bar',name="dormancy",
                         text = Control_Results$dormancy, textposition = 'auto',
                         marker = list(color = '#ff7f0e',
                                       line = list(color = '#ff7f0e', width = 1.5)))
p3 <- p3 %>% add_trace(x = ~Names, y = ~activation, type = 'bar',name="activation",
                         text = Control_Results$activation, textposition = 'auto',
                         marker = list(color = '#2ca02c',
                                       line = list(color = '#2ca02c', width = 1.5)))

p3 <- p3 %>% layout(title = "Control LYSA activation probabilities ",
                      barmode = 'group',
                      xaxis = list(title = "Gene Names"),
                      yaxis = list(title = "Probability of Activation of N",range(0,0.7)))

show(p3)

p4 <- Treatment_Results %>% plot_ly()
p4 <- p4 %>% add_trace(x = ~Names, y = ~inhibition, type = 'bar',name="inhibition",
                       text = Treatment_Results$inhibition, textposition = 'auto',
                       marker = list(color = '#1f77b4',
                                     line = list(color = '#1f77b4', width = 1.5)))
p4 <- p4 %>% add_trace(x = ~Names, y = ~dormancy, type = 'bar',name="dormancy",
                       text = Treatment_Results$dormancy, textposition = 'auto',
                       marker = list(color = '#ff7f0e',
                                     line = list(color = '#ff7f0e', width = 1.5)))
p4 <- p4 %>% add_trace(x = ~Names, y = ~activation, type = 'bar',name="activation",
                       text = Treatment_Results$activation, textposition = 'auto',
                       marker = list(color = '#2ca02c',
                                     line = list(color = '#2ca02c', width = 1.5)))

p4 <- p4 %>% layout(title = "Saline Treatment LYSA activation probabilities ",
                    barmode = 'group',
                    xaxis = list(title = "Gene Names"),
                    yaxis = list(title = "Probability of Activation of N",range=c(0,0.6)))

show(p4)

p5 <- plot_ly(Treatment_Top_Five, x = ~X_Names, y = ~values,
              marker=list(color=c('#2ca02c','#ff7f0e','#1f77b4','#d62728','#9467bd')), 
              name='Saline Treatment Two Point Intervention Top Five',
              text=round(Treatment_Top_Five$values,4), textposition='auto',type = 'bar') %>% 
  layout(title="Saline Treatment Two Point Intervention Top Five",xaxis=list(title='Genes'),yaxis = list(title = 'Probabilitly ',range=c(0,0.7)))
show(p5)

p6 <- plot_ly(Control_Top_Five, x = ~X_Names, y = ~values,
              marker=list(color=c('#2ca02c','#ff7f0e','#1f77b4','#d62728','#9467bd')), 
              name='Control Two Point Intervention Top Five',
              text=round(Control_Top_Five$values,4), textposition='auto',type = 'bar') %>% 
  layout(title="Control Two Point Intervention Top Five",xaxis=list(title='Genes'),yaxis = list(title = 'Probabilitly ',range=c(0,0.6)))
show(p6)

