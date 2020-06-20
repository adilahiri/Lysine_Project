get_param <- function(prior_param,node,parent){
  alpha_dirichlet_posterior <- matrix(ncol = 3)
  # Case when there are no parents
  if(is.null(parent)){
    counts_inhibit<- as.numeric(table(node)[1])
    counts_dormant<- as.numeric(table(node)[2])
    counts_active<- as.numeric(table(node)[3])
    
    alpha_inhibit <- prior_param[1] + counts_inhibit
    alpha_dormant <- prior_param[2] + counts_dormant
    alpha_active  <- prior_param[3] + counts_active
    
    
  }
  # Case when there is one parent
  
  # Case when there are two parents
  
  # Case when there are three parents
  
  # Case when there are five parents
  
  alpha_dirichlet_posterior <- c(alpha_inhibit,alpha_dormant,alpha_active)
  return(alpha_dirichlet_posterior)
}