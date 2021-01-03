# Negative log-likelhood function for binomial likelihood
seroprev_negative_log_likelihood <- function(pars,x){
  p <- seroprevalence(x,pars[1],pars[2])
  NLL <- -sum(x$k * log(p) + (x$n - x$k) * log(1-p))
}

# Age seroprevalence function
seroprevalence <- function(x,c,lambda){
  p <- c/lambda * (exp(-lambda * x$age_low)-exp(-lambda * x$age_upp))/(x$age_upp-x$age_low)
  return(p)
}