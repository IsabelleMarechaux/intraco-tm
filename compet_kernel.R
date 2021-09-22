# competition functions

# assume c_i is the performance

compet_step <- function(x_i, x_j){
 res <- rep(0, length.out = length(x_i))
 if(length(x_j) == 1) {
  res[x_i > rep(x_j, length.out = length(x_i))] <- 1
 }else{
    stop("length of x_j greater than 1")
 }
 res
}


compet_lin <- function(x_i, x_j){
  if(length(x_j) == 1) {
    res <- max(x_i - x_j, 0)
  }else{
    stop("length of x_j greater than x_i")
  }
  res
}


compet_logit <- function(x_i, x_j, alpha = 2){
  if(length(x_j) == 1) {
    res <- exp(-alpha*x_j)/(exp(-alpha*x_i) + exp(-alpha*x_j))
  }else{
    stop("length of x_j greater than x_i")
  }
  res
}







