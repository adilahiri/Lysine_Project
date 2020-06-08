get_compare <- function (data){
  Q=get_quartile(data)
  q1=Q[1]
  q3=Q[2]
  index_q1 <- which(data < q1)
  index_q3 <- which(data > q3)
  index_q1_q3 <- which(data >= q1 & data <= q3)
  data[index_q1] <- (-1)
  data[index_q3] <- 1
  data[index_q1_q3] <- 0
  return(data)
}