att_to_perf_function = function(a, amin = -4, amax = 4, nonlin = 1) {
  if(!amin<amax) {
    return("error: amin must be < amax")
  } else if(nonlin<1 & nonlin>-1) {
    return("error: nonlin must be less than -1, or greater than +1")
  } else {
    result = numeric(length(a))
    
    result[a<=amin] = 0
    result[a>=amax] = 1
    
    ps = which(a<amax & a>amin)
    if(length(ps)>0) {
      tmp = ((a[ps]-amin)/(amax-amin))
      if(nonlin>0) {
        tmp = tmp^nonlin
      } else {
        tmp = 1-(1-tmp)^(-nonlin)
      }
      
      result[ps] = tmp
    }
    
    return(result)
  }
}


if(FALSE) { # example
  atest = sort(rnorm(1000, 0, 2))
  linrel = att_to_perf_function(a = atest, nonlin = 1)
  sublinrel = att_to_perf_function(a = atest, nonlin = -5)
  suplinrel = att_to_perf_function(a = atest, nonlin = 5)
  
  matplot(atest, cbind(linrel, sublinrel, suplinrel), type = "l", lty = c(1,2,4), col = 1, xlab = "a", ylab = "perf")
  abline(h=c(0,1),lty=3)
}

