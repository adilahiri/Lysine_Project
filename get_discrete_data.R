get_discrete_data <- function (data,meth)
{
  for (iter in 1:14){
    temp <- discretize(data[,iter],method= meth,
                              breaks=3,labels = c(-1,0,1))
    data[,iter] <- as.numeric(levels(temp))[temp]
  }
  
  return(data)
}