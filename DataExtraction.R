setwd("~/Desktop/Research/Septi/GSE98455_RAW/TextFile")
library(readr)
library(DESeq2)

filename=list.files(pattern = "*.txt")

df_raw<-data.frame(matrix(nrow=57845,ncol=1))
for(iter in 1:368){
  temp=read_delim(filename[iter], "\t", escape_double = FALSE, 
                 col_names = FALSE, 
                 trim_ws = TRUE)
  temp = subset(temp, select = -c(X1) )
  df_raw=cbind(df_raw,temp)
}
df_raw<-df_raw[,-1]

write.csv(df_raw,"Counts.csv")