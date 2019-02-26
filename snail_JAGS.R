model{
  for (i in 1:nobs){
    #eta[i] <- b5*x1[i] + b12*x2[i] + b14*x3[i] + b15*x4[i] + b16*x5[i] + b18*x6[i] + b19*x7[i]
    #pi[i] <- 1/(1+exp(-eta[i]))
    #logit(pi[i]) <- -(b0+b5*x1[i] + b12*x2[i] + b14*x3[i] + b15*x4[i] + b16*x5[i] + 
                              #b18*x6[i] + b19*x7[i])
    y[i] ~ dbern(logit(-(b0+b5*x1[i] )))
    #+ b12*x2[i] + b14*x3[i] + b15*x4[i] + b16*x5[i] + b18*x6[i] + b19*x7[i]
  }
  
  #priors
  b0~dnorm(0,.01)
  b5~dnorm(0,.01)
  #b12~dnorm(0,.01)
  #b14~dnorm(0,.01)
  #b15~dnorm(0,.01)
  #b16~dnorm(0,.01)
  #b18~dnorm(0,.01)
  #b19~dnorm(0,.01)
}

