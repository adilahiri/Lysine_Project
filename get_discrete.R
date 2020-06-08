get_discrete <- function (data)
{
  for(iter in 1:14){
    data[,iter] <- get_compare(data[,iter])
  }
  return(data)
}