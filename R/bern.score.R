bern.score <- function(x,what="person",method="bayes"){
  #require(fitdistrplus)
  #require(R2jags)
  #require(DPpackage)
  #require(coda)
  x <- as.matrix(x)
  if(ncol(x) == 1) x <- t(x)
  ### Item====
  if (what == "item"){
    ### Frequentist
    if (method == "freq"){
      v <- ncol(x)
      Easiness <- vector("numeric")
      sd <- vector("numeric")
      BIC <- vector("numeric")
      AIC <- vector("numeric")
      for (i in 1:v) {
        set.seed(i)
        A <- fitdistrplus::fitdist(data=x[,i], dist="binom",
                                   fix.arg=list(size=(l-1)),
                                   start=list(prob=.5))
        Easiness[i] <- A$estimate
        sd[i] <- A$sd
        BIC[i] <- A$bic
        AIC[i] <- A$aic
      }
      Results <- list("easi"=Easiness,"sd"=sd,"bic"=BIC,"aic"=AIC)
      ### Bayesian
    } else if (method == "bayes") {
      v <- ncol(x)
      Easiness <- vector("numeric")
      sd <- vector("numeric")
      DIC <- vector("numeric")
      for (i in 1:v) {
        y <- x[,i]
        k <- list()
        for (j in 1:length(y)) {k[[j]] <- c(rep(0,(l-1-y[j])),rep(1,y[j]))}
        k <- unlist(k)
        n <- length(k)
        k <- length(k[which( unlist(k) == 1)])

        data <- list("k", "n")
        params <- c("theta")
        inits <- function(){list("theta"=stats::rbeta(1,5,5))}

        modelString = "
        model{
        # Prior Distribution for Rate Theta
        theta ~ dbeta(1,1)
        # Observed Counts
        k ~ dbin(theta,n)
        }" # close quote for modelString
        model = textConnection(modelString)

        nthin    = 1    # How Much Thinning?
        nchains  = 3    # How Many Chains?
        nburnin  = 1000  # How Many Burn-in Samples?
        nsamples = 15000 # How Many Recorded Samples?
        ### Calling JAGS to sample
        samples <- R2jags::jags(data, inits, params, model.file=model,
                                n.chains=nchains, n.iter=nsamples, n.burnin=nburnin,
                                n.thin=nthin, DIC=T, jags.seed=666)
        Easiness[i] <- mean(samples$BUGSoutput$sims.list$theta)
        sd[i] <- sd(samples$BUGSoutput$sims.list$theta)
        DIC[i] <- samples$BUGSoutput$DIC
      }
      Results <- list("easi"=Easiness,"sd"=sd,"dic"=DIC)
      ### Bayesian Semiparametric
  } else if (method == "bsp") {
    DPest <- vector("numeric")
    DPsd <- vector("numeric")
    ess <- vector("numeric")
    v <- ncol(x)
    for (R in 1:v) {
      # Data
      y <- cbind(x[,R],rep((l-1),length(x[,R])))

      # Prior information
      prior<-list(a0=1, b0=1, a1=1, b1=1)

      # Initial state
      state <- NULL

      # MCMC parameters
      mcmc <- list(nburn=5000,
                   nsave=10000,
                   nskip=3,
                   ndisplay=1000)

      # Fitting the model
      fit <- DPpackage::DPbetabinom(y=y,ngrid=100,
                                    prior=prior,
                                    mcmc=mcmc,
                                    state=state,
                                    status=TRUE)

      DPest[R] <- mean(fit$save.state$randsave)
      DPsd[R] <- sd(fit$save.state$randsave)
      ess[R] <- coda::effectiveSize(fit$densp.m)
    }
    Easiness <- DPest
    sd <- DPsd
    ESS <- ess
    Results <- list("easi"=Easiness,"sd"=sd,"ess"=ESS)
  } else { warning("Unkown method!") }
    ### Person estimates====
} else if (what == "person") {
  ### Frequentist
  if (method == "freq"){
    n <- nrow(x)
    Ability <- vector("numeric")
    sd <- vector("numeric")
    BIC <- vector("numeric")
    AIC <- vector("numeric")
    for (i in 1:n) {
      set.seed(i)
      A <- fitdistrplus::fitdist(data=x[i,], dist="binom",
                                 fix.arg=list(size=(l-1)),
                                 start=list(prob=.5))
      Ability[i] <- A$estimate
      sd[i] <- A$sd
      BIC[i] <- A$bic
      AIC[i] <- A$aic
    }
    Results <- list("abil"=Ability,"sd"=sd,"bic"=BIC,"aic"=AIC)
    ### Bayesian
  } else if (method == "bayes") {
    n <- nrow(x)
    v <- ncol(x)
    
    data <- list("x","n", "v")
    params <- c("theta")
    inits <- function(){list("theta"=stats::rbeta(n,1,1))}

    modelString = "
    model{
      # Prior Distribution for Theta
      mu ~ dbeta(1,1)
      tau ~ dgamma(1e-2, 1e-2)
      for(i in 1:n){
        theta[i] ~ dbeta((mu * tau) + 1, ((1 - mu) * tau) + 1)
      }
      
      # Theta for all respondents
      for(i in 1:n){
        for(j in 1:v) {
          # Observed Counts
          x[i,j] ~ dbern(theta[i])
        }
      }
    }" # close quote for modelString
    model = textConnection(modelString)

    nthin    = 1    # How Much Thinning?
    nchains  = 3    # How Many Chains?
    nburnin  = 1000  # How Many Burn-in Samples?
    nsamples = 15000 # How Many Recorded Samples?
    ### Calling JAGS to sample
    samples <- R2jags::jags(data, inits, params, model.file=model,
                            n.chains=nchains, n.iter=nsamples, n.burnin=nburnin,
                            n.thin=nthin, DIC=T, jags.seed=666)
    Ability <- mean(samples$BUGSoutput$sims.list$theta)
    ABHDI   <- apply(samples$BUGSoutput$sims.list$theta, 2, quantile, c(.025,.974))
    fabil   <- samples$BUGSoutput$sims.list$theta
    sd      <- sd(samples$BUGSoutput$sims.list$theta)
    DIC     <- samples$BUGSoutput$DIC
      
    Results <- list("abil"=Ability,"abilHDI"=ABHDI,"abilFull"=fabil,"dic"=DIC)
    ### Bayesian Semiparametric
} else if (method == "bsp") {
  DPest <- vector("numeric")
  DPsd <- vector("numeric")
  ess <- vector("numeric")
  n <- nrow(x)
  for (R in 1:n) {
    # Data
    y <- cbind(x[R,],rep((l-1),length(x[R,])))

    # Prior information
    prior<-list(a0=1, b0=1, a1=1, b1=1)

    # Initial state
    state <- NULL

    # MCMC parameters
    mcmc <- list(nburn=5000,
                 nsave=10000,
                 nskip=3,
                 ndisplay=1000)

    # Fitting the model
    fit <- DPpackage::DPbetabinom(y=y,ngrid=100,
                                  prior=prior,
                                  mcmc=mcmc,
                                  state=state,
                                  status=TRUE)

    DPest[R] <- mean(fit$save.state$randsave)
    DPsd[R] <- sd(fit$save.state$randsave)
    ess[R] <- coda::effectiveSize(fit$densp.m)
  }
  Ability <- DPest
  sd <- DPsd
  ESS <- ess
  Results <- list("abil"=Ability,"sd"=sd,"ess"=ESS)
} else { warning("Unkown method!") }
  } else { warning("Unkown dimension!") }
  return(Results)
}
