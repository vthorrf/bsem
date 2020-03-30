\name{bern.score}
\alias{bern.score}
\title{Bernoulli Scores}
\description{
Estimate theta parameters of binomial distributions as scores for item easiness and person ability.
}
\usage{
bern.score(x,levels,what="item",method="freq")
}
\arguments{
   \item{x}{An data set.}
   \item{levels}{How many levels does the scale has.}
   \item{what}{Says if it estimates bias for "person" or for "item".}
   \item{method}{The statistical framework of the binomial test. "freq" if frequentist is desired. "bayes" if bayesian is desired. "bsp" if semiparametric bayesian is desired.}
}
\value{
A list that depends if Bayesian or frequentist method are used.
The frequentist contains the following components:
  \item{easi or abil}{A vector of estimates of easiness or ability, depending on the desired outcome.}
  \item{sd}{The standard deviation of the estimates.}
  \item{bic}{BIC fit index for each estimate.}
  \item{aic}{AIC fit index for each estimate.}

The parametric Bayesian contains the following components:
  \item{easi or abil}{A vector of estimates of easiness or ability, depending on the desired outcome.}
  \item{sd}{The standard deviation of the estimates.}
  \item{dic}{DIC fit index for each estimate.}

The semiparametric Bayesian contains the following components:
  \item{easi or abil}{A vector of estimates of easiness or ability, depending on the desired outcome.}
  \item{sd}{The standard deviation of the estimates.}
  \item{ess}{Effective sample size; the sample size adjusted for autocorrelation.}
}
\examples{
### Random Data with 100 observations and 10 variables with 5 levels of response
n = 100; v = 10; l = 5
Data <- matrix(NA,nrow=n,ncol=v)
psi <- seq(length=n,-3,3); theta <- seq(length=v,-3,3); eta <- matrix(NA,ncol=v,nrow=n)
seed <- matrix(c(1:(n*v)),nrow=n,ncol=v)
for(i in 1:n)
{
  for(j in 1:v)
  {
    set.seed(seed[i,j])
    eta[i,j] <- exp(psi[i]-theta[j])/(1+exp(psi[i]-theta[j]))
    Data[i,j] <- rbinom(1,l-1,eta[i,j])
  }
}
rm(n,v,eta,seed,i,j)

### Create Frequentist Scores
AbilF <- latent.vs(Data,l,what="person",method="freq")
Af <- AbilF$abil
EasiF <- latent.vs(Data,l,what="item",method="freq")
Ef <- EasiF$easi

### Create Bayesian Scores
AbilB <- latent.vs(Data,l,what="person",method="bayes")
Ab <- AbilB$abil
EasiB <- latent.vs(Data,l,what="item",method="bayes")
Eb <- EasiB$easi

### Create semiparametric Bayesian Scores
AbilS <- latent.vs(Data,l,what="person",method="bsp")
As <- AbilS$abil
EasiS <- latent.vs(Data,l,what="item",method="bsp")
Es <- EasiS$easi

### Compare
cor(data.frame(cbind(Af,Ab,As,psi)))
plot(data.frame(cbind(Af,Ab,As,psi)))

cor(data.frame(cbind(Ef,Eb,Es,theta)))
plot(data.frame(cbind(Ef,Eb,Es,theta)))
}
