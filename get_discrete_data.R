get_discrete_data <- function (data,meth)
{
  col_num = ncol(data)
  for (iter in 1:col_num){
    temp <- discretize(data[,iter],method= meth,
                              breaks=3,labels = c(-1,0,1))
    data[,iter] <- as.numeric(levels(temp))[temp]
  }
  
  return(data)
}