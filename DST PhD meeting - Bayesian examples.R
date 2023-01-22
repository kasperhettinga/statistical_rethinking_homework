library(rethinking)
library(BayesFactor)
library(bayestestR)
library(dagitty)

# Load the data
grass_ketones <- read.csv("feed_acetone.csv")

#Log 10 transform the GCMS values
grass_ketones$log_acetone <- log(grass_ketones$Acetone, base=10)

####Contrast analysis####

#Build the MCMC model
model_acetone <- ulam(
  alist(
    log_acetone ~ dnorm(mu,sigma),
    mu <- a[Feed],
    a[Feed] ~ dnorm(7,2),
    sigma ~ dexp(1)
  ), data=grass_ketones, chains=4, cores=4
)

#Check the primary output of the MCMC model that has been build
precis(model_acetone,depth=2)
traceplot(model_acetone)
trankplot(model_acetone)

#Extract samples from the posterior for estimating the predicted mean
post <- extract.samples(model_acetone)
dens(post$a[,1], col=2, xlim=c(6.9,7.7))
dens(post$a[,2], col=3, add=TRUE)
dens(post$a[,3], col=4, add=TRUE)
dens(post$a[,4], col=5, add=TRUE)

#Calculate contrast in acetone mean relative to reference diet
contrast_mu_17_0 <- post$a[,2]-post$a[,1]
contrast_mu_33_0 <- post$a[,3]-post$a[,1]
contrast_mu_50_0 <- post$a[,4]-post$a[,1]
dens(contrast_mu_17_0, xlab="posterior mean acetone contrast", lwd=2, col=2,
     xlim=c(-0.6,0.3))
dens(contrast_mu_33_0, lwd=2, col=3, add=TRUE)
dens(contrast_mu_50_0, lwd=3, col=4, add=TRUE)

#Posterior simulated acetone data distribution
A_0 <- rnorm(1e4,post$a[,1],post$sigma)
A_17 <- rnorm(1e4,post$a[,2],post$sigma)
A_33 <- rnorm(1e4,post$a[,3],post$sigma)
A_50 <- rnorm(1e4,post$a[,4],post$sigma)
dens(A_0, col=2, xlim=c(6.5,8.1))
dens(A_17, col=3, add=TRUE)
dens(A_33, col=4, add=TRUE)
dens(A_50, col=5, add=TRUE)

#Contrast of simulated acetone data
contrast_0_17 <- A_17 - A_0
dens(contrast_0_17, xlab="posterior predicted acetone contrast", lwd=3)
sum(contrast_0_17<0)/1e4

contrast_0_33 <- A_33 - A_0
dens(contrast_0_33, xlab="posterior predicted acetone contrast", lwd=3)
sum(contrast_0_33<0)/1e4

contrast_0_50 <- A_50 - A_0
dens(contrast_0_50, xlab="posterior predicted acetone contrast", lwd=3)
sum(contrast_0_50<0)/1e4

####Two-sample Bayes Factor####
#Based on: https://bayesfactor.blogspot.com/2014/02/bayes-factor-t-tests-part-2-two-sample.html

#Take out 0 vs 33% as an example
grass_ketones_0_33 <- grass_ketones[c(0:9,19:27),]

#regular t-test
classical.test = t.test(log_acetone ~ Feed, data = grass_ketones_0_33, var.eq = TRUE)
classical.test

#two-sided BF test
bf2s = ttestBF(formula = log_acetone ~ Feed, data = grass_ketones_0_33)
bf2s
chains = posterior(bf2s, iterations = 1000)
plot(chains[,2])

bf1s = ttestBF(formula = log_acetone ~ Feed, data = grass_ketones_0_33, 
             nullInterval=c(0,Inf))
bf1s

grass_ketones$feed.factor <- as.factor(grass_ketones$Feed)

#Bayesian ANOVA alternative
bfA = anovaBF(formula = log_acetone ~ feed.factor, data=grass_ketones)
bfA

#Manual Bayes Factor calculation
# based on: https://easystats.github.io/bayestestR/

#Prior and posterior simulated acetone data distribution
prior_contrast <- rnorm(1e6,0,post$sigma)
posterior_0 <- rnorm(1e6,post$a[,1],post$sigma)
posterior_17 <- rnorm(1e6,post$a[,2],post$sigma)
posterior_33 <- rnorm(1e6,post$a[,3],post$sigma)
posterior_50 <- rnorm(1e6,post$a[,4],post$sigma)

posterior_contrast_0_33 <- posterior_0-posterior_33
dens(prior_contrast,col="blue")
dens(posterior_contrast_0_33,col="red",add=TRUE)

bayesfactor_parameters(posterior_contrast_0_33, prior_contrast, 
                       direction = "two-sided", null = 0)

####Regression example####
BDI <- read.csv("Abuse_BDI.csv")

#Build the MCMC model
model_bdi_1 <- ulam(
  alist(
    bdi ~ dnorm(mu,sigma),
    mu <- a + b*abuse,
    a ~ dnorm(0.5,0.5),
    b ~ dnorm(0.25,0.25),
    sigma ~ dexp(1)
  ), data=BDI, chains=4, cores=4
)

#Check the primary output of the MCMC model that has been build
precis(model_bdi_1)
traceplot(model_bdi_1)
trankplot(model_bdi_1)

#Check prior predictions
set.seed(1234)
prior1 <- extract.prior(model_bdi_1)
plot(NULL, xlim=c(0,5), ylim=c(-1,2), xlab="Abused milk added (mL/L)", 
     ylab="BDI value")
abuse_seq <- seq(from=0, to=5, length.out=30)
mu <- link(model_bdi_1, post=prior1, data=data.frame(abuse=abuse_seq))
for(i in 1:50) lines(abuse_seq,mu[i,],col=col.alpha("black",0.3))

#Can we do better? Cange priors to get more sensible prediction
model_bdi_2 <- ulam(
  alist(
    bdi ~ dnorm(mu,sigma),
    mu <- a + b*abuse,
    a ~ dnorm(0.5,0.2),
    b ~ dnorm(0.25,0.1),
    sigma ~ dexp(1)
  ), data=BDI, chains=4, cores=4
)
precis(model_bdi_2)

#Check prior predictions
set.seed(1234)
prior2 <- extract.prior(model_bdi_2)
plot(NULL, xlim=c(0,5), ylim=c(0,1.5), xlab="Abused milk added (mL/L)", 
     ylab="BDI value")
abuse_seq <- seq(from=0, to=5, length.out=30)
mu <- link(model_bdi_2, post=prior2, data=data.frame(abuse=abuse_seq))
for(i in 1:50) lines(abuse_seq,mu[i,],col=col.alpha("black",0.3))

#Plot posterior predictions & simulated data
post <- link(model_bdi, data=list(abuse=abuse_seq))
sim.bdi <- sim(model_bdi,data=list(abuse=abuse_seq))
bdi.ci <- apply(sim.bdi,2,PI,prob=0.89)
mu <- apply(post,2,mean)
ci <- apply(post,2,PI,prob=0.89)
plot(bdi~abuse,data=BDI, xlab="Abused milk added (mL/L)", ylab="BDI")
lines(abuse_seq,mu)
shade(ci,abuse_seq)
shade(bdi.ci,abuse_seq)

####DAG for explicit scientific model####

# The presence of mastitis (M) will impact proteolysis of caseins (P)
# Mastitis has three levels (health, subclinical, clinical; indexed 1, 2, 3)
# At any point in time on a farm, mastitis will be 0.8, 0.15, 0.05
# Mastitis will lead to higher somatic cell count (SCC)
# clinical mastitis will lead to more proteolysis & higher SSC
# Higher cell count will in general already lead to proteolysis
# Unobserved parameters: 
# days in lactation (Ud), with higher Ud leading to higher SSC & M
# b-casein ratio (Ub), with higher Ub leading to more proteolysis
dag <- dagitty("dag{
                  Ub [unobserved]
                  Ud [unobserved]
                  M -> P <- S
                  M <- Ud -> S <- M
                  Ub -> P
                  }")
coordinates(dag) <- list(x=c(M=0, S=2, P=1, Ud=1, Ub=2),
                            y=c(M=1, S=1, P=2, Ud=0, Ub=2))
drawdag(dag)
adjustmentSets(dag6H6a, exposure="S", outcome="P")
impliedConditionalIndependencies(dag)