BSEM <- function(x, factors){
  #require(R2jags)
  #require(coda)
  ### Data====
  x       <- as.matrix(data.frame(scale(x)))
  factors <- as.vector(factors)
  Rmat    <- diag(length(unique(factors)))
  fixed   <- rep(which(duplicated(factors) == F), each=length(unique(factors)))
  
  n <- nrow(x)
  v <- ncol(x)
  l <- length(unique(factors))
    
  data <- list("x", "n", "v", "l", "fixed", "factors")
  params <- c("theta", "lambda", "sigma")
  inits <- function(){list(#"theta"=matrix(stats::rnorm(n*l, 0, 1e2), n, l) ,
                           "LL"=stats::rnorm(v, 0, 1e2)   ,
                           "sigma"=stats::rgamma(v, 1e-2, 1e-2))}

  ### Model====
  modelString = "
  data{
    for(i in 1:n){
      pseudo[i] ~ dnorm(0, 1)
    }
  }
  model{
    # Factor scores
    for(f in 1:l){
      corr[f] ~ dnorm(0, 1e-4)T(-1,1)
      #prec[f] ~ dgamma(1e-3, 1e-2)
    }
    for(i in 1:n){
      for(f in 1:l){
        #theta[i, f] ~ dnorm(pseudo[i]*corr[f], prec[f])
        theta[i, f] ~ dnorm(pseudo[i]*corr[f], 1)
      }
    }
    
    # Factor loadings
    for(j in 1:v){
      LL[j] ~ dnorm(0, 1e-4)
      sigma[j] ~ dgamma(1e-3, 1e-2)
    }
    for(j in 1:v){
      lambda[j] <- ifelse(j == fixed[j], 1, LL[j])
    }
    
    # Observed responses
    for(i in 1:n){
      for(j in 1:v) {
        x[i,j] ~ dnorm(lambda[j] * theta[i, factors[j]], sigma[j])
      }
    }
  }" # close quote for modelString
  model = textConnection(modelString)

  ### Fitting====
  nthin    = 1     # How Much Thinning?
  nchains  = 3     # How Many Chains?
  nburnin  = 1000  # How Many Burn-in Samples?
  nsamples = 11000 # How Many Recorded Samples?
  ## Calling JAGS to sample
  samples <- R2jags::jags(data, inits, params, model.file=model,
                          n.chains=nchains, n.iter=nsamples, n.burnin=nburnin,
                          n.thin=nthin, DIC=T, jags.seed=666)

  ### Results====
  # Factor score
  factorScore <- scale(sapply(1:l, function(g) {
    colMeans(samples$BUGSoutput$sims.list$theta[,,g])
  }))
  precision <- sapply(1:l, function(g) {
    apply(samples$BUGSoutput$sims.list$theta[,,g], 2, sd)
  })
  # Lambdas
  standardLambda <- as.vector(sapply(1:l, function(g) {
    cor(factorScore[,g], x[,which(factors == g)])
  }))
  lambda   <- data.frame("Variable"=colnames(x),
                         "Factor"=factors,
                         "Factor_Load"=colMeans(samples$BUGSoutput$sims.list$lambda),
                         "Std_Factor_Load"=standardLambda,
                         "Variance"= 1/colMeans(samples$BUGSoutput$sims.list$sigma))
  Rho      <- cor(factorScore)
  colnames(Rho) <- rownames(Rho) <- paste("F",1:l, sep="")
  colnames(factorScore) <- colnames(Rho)
  DIC     <- samples$BUGSoutput$DIC
  Results <- list("abil"=data.frame(factorScore),"output"=lambda,
                  "corr"=Rho,"dic"=DIC,"full"=samples)
  return(Results)
}
