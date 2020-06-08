get_quartile <- function(data) # Enter data columnwise
{
  Q= quantile(data)
  First_Third_Quart <- c(Q[2][[1]],Q[4][[1]])
  return(First_Third_Quart)
}