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
    res <- pmax(x_i - rep(x_j, length.out = length(x_i)), 0)
  }else{
    stop("length of x_j greater than x_i")
  }
  res
}


compet_logistic <- function(x_i, x_j, alpha = 20){
  if(length(x_j) == 1) {
    res <- exp(-alpha*x_j)/(exp(-alpha*x_i) + exp(-alpha*x_j))
  }else{
    stop("length of x_j greater than x_i")
  }
  res
}



# Plot of teh competition kernel
if(FALSE) { # example
 perf_vec <- (0:100)/100
 plot(perf_vec, compet_step(perf_vec, 0.5), type = "l")
 lines(perf_vec, compet_lin(perf_vec, 0.5), lty = 2)
 lines(perf_vec, compet_logistic(perf_vec, 0.5), lty = 2)
}



