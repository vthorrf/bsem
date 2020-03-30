BSEMIRT <- function(x, factors, k=max(myData)){
  #require(R2jags)
  #require(coda)
  ### Data====
  x <- x
  factors <- as.vector(factors)
  Rmat    <- diag(length(unique(factors)))
  
  n <- nrow(x)
  v <- ncol(x)
  l <- length(unique(factors))
  k <- k
    
  data <- list("x", "n", "v", "l", "k", "factors")
  params <- c("theta", "Disc", "Diff", "loads")
  #inits <- function(){list( #"theta"=matrix(stats::rnorm(n*l, 0, 1e2), n, l) ,
  #                          "Disc"=stats::rnorm(v, 0, 1e2)   ,
  #                          "Diff"=stats::rnorm(v, 0, 1e2) )}

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
    
    # Difficulty and Discrimination
    mu_a ~ dnorm(0, 1)T(0,)
    mu_b ~ dnorm(0, 1)
    var_a ~ dgamma(1e-2, 1e-2)
    var_b ~ dgamma(1e-2, 1e-2)
    
    for(j in 1:v){
      Disc[j] ~ dnorm(mu_a, var_a)
      Diff[j] ~ dnorm(mu_b, var_b)
      loads[j] <- Disc[j]/sqrt(1+(Disc[j]^2))
    }
    
    # Observed responses
    for(i in 1:n){
      for(j in 1:v) {
        Pr[i,j] <- ilogit( Disc[j] * (theta[i, factors[j]] - Diff[j]) )
        x[i,j] ~ dbin(Pr[i,j], k)
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
  samples <- R2jags::jags(data, inits=NULL, params, model.file=model,
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
  standardLambda <- colMeans(samples$BUGSoutput$sims.list$loads)
  Discrimination <- colMeans(samples$BUGSoutput$sims.list$Disc)
  Difficulty     <- colMeans(samples$BUGSoutput$sims.list$Diff)
  
  lambda   <- data.frame("Variable"=colnames(x),
                         "Factor"=factors,
                         "Discrimination"=Discrimination,
                         "Factor_Load"=standardLambda,
                         "Difficulty"=Difficulty)
  Rho      <- cor(factorScore)
  colnames(Rho) <- rownames(Rho) <- paste("F",1:l, sep="")
  colnames(factorScore) <- colnames(Rho)
  DIC     <- samples$BUGSoutput$DIC
  Results <- list("abil"=data.frame(factorScore),"output"=lambda,
                  "corr"=Rho,"dic"=DIC,"full"=samples)
  return(Results)
}
