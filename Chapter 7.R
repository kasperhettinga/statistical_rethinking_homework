library(rethinking)

sppnames <- c("afarensis","africanus","habilis","boisei","rudolfensis",
              "ergaster","sapiens")
brainvolcc <- c(438,452,612,521,752,871,1350)
masskg <- c(37.0,35.5,34.5,41.5,55.5,61.0,53.5)
d <- data.frame(species=sppnames,brain=brainvolcc,mass=masskg)

plot(d$brain~d$mass, xlab="approx body mass", ylab="brain volume")
d$mass_std <- (d$mass-mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain/max(d$brain)

model_brain1 <- quap(
  alist(
    brain_std ~ dnorm(mu,exp(log_sigma)),
    mu <- a+b*mass_std,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ), data=d
)

set.seed(12)
s <- sim(model_brain1)
r <- apply(s,2,mean)-d$brain_std
resid_var <- var2(r)
outcome_var <- var2(d$brain_std)
1-resid_var/outcome_var

r2_is_bad <- function(quap_fit){
  s <- sim(quap_fit,refresh=0)
  r <- apply(s,2,mean)-d$brain_std
  1-var2(r)/var2(d$brain_std)
}

model_brain2 <- quap(
  alist(
    brain_std ~ dnorm(mu,exp(log_sigma)),
    mu <- a+b[1]*mass_std+b[2]*mass_std^2,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ), data=d, start=list(b=rep(0,2))
)

model_brain3 <- quap(
  alist(
    brain_std ~ dnorm(mu,exp(log_sigma)),
    mu <- a+b[1]*mass_std+b[2]*mass_std^2+b[3]*mass_std^3,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ), data=d, start=list(b=rep(0,3))
)

model_brain4 <- quap(
  alist(
    brain_std ~ dnorm(mu,exp(log_sigma)),
    mu <- a+b[1]*mass_std+b[2]*mass_std^2+b[3]*mass_std^3+b[4]*mass_std^4,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ), data=d, start=list(b=rep(0,4))
)

model_brain5 <- quap(
  alist(
    brain_std ~ dnorm(mu,exp(log_sigma)),
    mu <- a+b[1]*mass_std+b[2]*mass_std^2+b[3]*mass_std^3+b[4]*mass_std^4+
    b[5]*mass_std^5,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ), data=d, start=list(b=rep(0,5))
)

model_brain6 <- quap(
  alist(
    brain_std ~ dnorm(mu,0.001),
    mu <- a+b[1]*mass_std+b[2]*mass_std^2+b[3]*mass_std^3+b[4]*mass_std^4+
      b[5]*mass_std^5+b[6]*mass_std^6,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10),
    log_sigma ~ dnorm(0,1)
  ), data=d, start=list(b=rep(0,6))
)

post <- extract.samples(model_brain1)
mass_seq <- seq(from=min(d$mass_std),to=max(d$mass_std),length.out=100)
l <- link(model_brain1, data=list(mass_std=mass_seq))
mu <- apply(l,2,mean)
ci <- apply(l,2,PI)
plot(brain_std~mass_std,data=d)
lines(mass_seq,mu)
shade(ci,mass_seq)

post <- extract.samples(model_brain6)
mass_seq <- seq(from=min(d$mass_std),to=max(d$mass_std),length.out=100)
l <- link(model_brain6, data=list(mass_std=mass_seq))
mu <- apply(l,2,mean)
ci <- apply(l,2,PI)
plot(brain_std~mass_std,data=d)
lines(mass_seq,mu)
shade(ci,mass_seq)

set.seed(1234)
lppd(model_brain1,n=1e4)

sapply(list(model_brain1,model_brain2,model_brain3,model_brain4,model_brain5,
            model_brain6), function(m) sum(lppd(m)))

#Manually calculate WAIC based on example in overthinking
data(cars)
model_waictest <- quap(
  alist(
    dist ~ dnorm(mu,sigma),
    mu <- a+b*speed,
    a ~ dnorm(0,100),
    b ~ dnorm(0,10),
    sigma ~ dexp(1)
  ), data=cars
)
plot(cars$dist~cars$speed)
post <- extract.samples(model_waictest,n=1e3)
WAIC(model_waictest)

n_samples <- 1e3
logprob <- sapply(1:n_samples,
                  function(s) {
                    mu <- post$a[s] + post$b[s]*cars$speed
                    dnorm(cars$dist, mu, post$sigma[s], log=TRUE)
                  })

n_cases <- nrow(cars)
lppd <- sapply(1:n_cases, function(i) log_sum_exp(logprob[i,])-log(n_samples))
sum(lppd)

pWAIC <- sapply(1:n_cases, function(i) var(logprob[i,]))
-2*(sum(lppd)-sum(pWAIC)) # Calculate WAIC
WAIC_vec <- -2*(lppd-pWAIC)
sqrt(n_cases*var(WAIC_vec)) # Calculate SE of WAIC

#Outlier detection in marriage example
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

model5.1 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
model5.2 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
model5.3 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

set.seed(1234)
compare(model5.1,model5.2,model5.3, func=PSIS)

set.seed(1234)
PSIS_5.3 <- PSIS(model5.3,pointwise=TRUE)
set.seed(1234)
WAIC_5.3 <- WAIC(model5.3,pointwise=TRUE)
plot(PSIS_5.3$k,WAIC_5.3$penalty,xlab="PSIS k",ylab="WAIC penalty")

#Excersises from the book

#7E1
#Uncertainty should be on a continuous scale to compare models
#Uncertainty should increase with the number of events
#Uncertainty should be additive

#7E2
entropy <- function(x){
  -sum(x*log(x))
}
entropy(c(0.7, 0.3))

#7E3
entropy(c(0.2,0.25,0.25,0.3))

#7E4
entropy(c(0.33,0.33,0.33))

#7M1
# WAIC & AIC give similar results if:
# the priors are flat 
# and/or posterior is overwhelmed by the likelihood
# and posterior is Gaussian.

#7M2
# With model comparison you look at markers for difference, which can also be
# used for detection of outliers
# With model selection, you want to find the specific best model to work with.

#7M3
# The calculation of the information criteria is summing up probabilities
# resulting in higher values for models with more observations

#7M4
# Re-use model5.3
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

model7M4a <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,1),
    bM ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)
model7M4b <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.1),
    bM ~ dnorm(0,0.1),
    sigma ~ dexp(1)
  ), data=d
)
model7M4c <- quap(
  alist(
    D ~ dstudent(2,mu,sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,1),
    bM ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)

set.seed(1234)
compare(model7M4a,model7M4b,model7M4c, func=PSIS)

set.seed(1234)
PSIS_7M4a <- PSIS(model7M4a,pointwise=TRUE)
set.seed(1234)
PSIS_7M4b <- PSIS(model7M4b,pointwise=TRUE)
plot(PSIS_7M4a$k,PSIS_7M4b$k,xlab="PSIS k model_a",ylab="PSIS k model_b")

#The penalty goes down with a model with narrower priors

#7M5
# Narrower priors leave less flexibility to the model. It will therefore tend to
# be less likely to come up with impossible models that overfit the data.
# In other words, the lower chance that very extreme values will drive the model
# is a benefit preventing overfitting.

#7M6
# Too narrow priors do not allow the model the capture the actual variation
# that is real, therefore leading to insufficient flexibility to capture this.

#7H1
data(Laffer)
head(Laffer)
d <- Laffer
d$Ra <- standardize(d$tax_rate)
d$Re <- standardize(d$tax_revenue)
head(d)

model7H1 <- quap(
  alist(
    Re ~ dnorm(mu,sigma),
    mu <- a + b*Ra,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)
precis(model7H1)
# There seems to be a linear increase?

model7H1q <- quap(
  alist(
    Re ~ dnorm(mu,sigma),
    mu <- a + b[1]*Ra + b[2]*Ra^2,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d, start=list(b=rep(0,2))
)
precis(model7H1q, depth=2)

model7H1c <- quap(
  alist(
    Re ~ dnorm(mu,sigma),
    mu <- a + b[1]*Ra + b[2]*Ra^2 + b[3]*Ra^3,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d, start=list(b=rep(0,3))
)
precis(model7H1c, depth=2)

set.seed(1234)
compare(model7H1,model7H1q,model7H1c, func=WAIC)
# The quadratic model seems marginally better, but not by much
# SE is much larger than the difference

Ra_seq <- seq(from=-3, to=3, length.out=30)
mu1 <- link(model7H1,data=list(Ra=Ra_seq))
mu.mean1 <- apply(mu1,2,mean)
mu.PI1 <- apply(mu1,2,PI,prob=0.91)
mu2 <- link(model7H1q,data=list(Ra=Ra_seq))
mu.mean2 <- apply(mu2,2,mean)
mu.PI2 <- apply(mu2,2,PI,prob=0.91)
mu3 <- link(model7H1c,data=list(Ra=Ra_seq))
mu.mean3 <- apply(mu3,2,mean)
mu.PI3 <- apply(mu3,2,PI,prob=0.91)
plot(d$Re~d$Ra)
lines(Ra_seq,mu.mean1,col="blue")
lines(Ra_seq,mu.mean2,col="red")
lines(Ra_seq,mu.mean3,col="green")
shade(mu.PI1,Ra_seq,col=col.alpha("blue",0.25))
shade(mu.PI2,Ra_seq,col=col.alpha("red",0.25))
shade(mu.PI3,Ra_seq,col=col.alpha("green",0.25))
#The one point on top is an outlier. No model comes close to newspapers.

#7H2
d2 <- d[d$tax_revenue != max(d$tax_revenue), ]
model7H2 <- quap(
  alist(
    Re ~ dnorm(mu,sigma),
    mu <- a + b*Ra,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d2
)
precis(model7H2)

model7H2q <- quap(
  alist(
    Re ~ dnorm(mu,sigma),
    mu <- a + b[1]*Ra + b[2]*Ra^2,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d2, start=list(b=rep(0,2))
)
precis(model7H2q, depth=2)

set.seed(1234)
compare(model7H1,model7H1q,model7H1c,func=WAIC)
set.seed(1234)
compare(model7H2,model7H2q,func=WAIC) 
# You get a much reduced SE 7 penalty term.

model7H2_r <- quap(
  alist(
    Re ~ dstudent(2,mu,sigma),
    mu <- a + b*Ra,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)
model7H2q_r <- quap(
  alist(
    Re ~ dstudent(2,mu,sigma),
    mu <- a + b[1]*Ra + b[2]*Ra^2,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d, start=list(b=rep(0,2))
)

set.seed(1234)
compare(model7H2_r,model7H2q_r, func=WAIC)
set.seed(1234)
compare(model7H2_r,model7H2q_r, func=PSIS)
# Differences between models is still small. No k-value errors anymore and 
# smaller SE than the dnorm models

#7H3
entropy <- function(x){
  -sum(x*log(x))
}

prob_1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
entropy(prob_1)
prob_2 <- c(0.8, 0.1, 0.05, 0.025, 0.025)
entropy(prob_2)
prob_3 <- c(0.05, 0.15, 0.7, 0.05, 0.05)
entropy(prob_3)
# Island 1 has the highest entropy, which is to be expected based on the equal 
# chances on each bird specie, giving a lot of uncertainty about the specie

D_12 <- sum(prob_1 * log(prob_1 / prob_2))
D_21 <- sum(prob_2 * log(prob_2 / prob_1))
D_13 <- sum(prob_1 * log(prob_1 / prob_3))
D_31 <- sum(prob_3 * log(prob_3 / prob_1))
D_23 <- sum(prob_2 * log(prob_2 / prob_3))
D_32 <- sum(prob_3 * log(prob_3 / prob_2))
# D23 has the largest distance: predicting island 3 based on island 2 gives
# the largest KL divergence: these two island are the most dissimilar
# D32 is smaller, as the deviation is smaller going from island 3->2

#7H4
# Collider bias
d <- sim_happiness(seed=1977, N_years=1000)
precis(d)
d2 <- d[d$age>17,]
d2$A <- (d2$age-18)/(65-18)
d2$mid <- d2$married+1

model_happymarried1 <- quap(
  alist(
    happiness ~ dnorm(mu,sigma),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm(0,1),
    bA ~ dnorm (0,2),
    sigma ~ dexp(1)
  ), data=d2
)
precis(model_happymarried1,depth=2)

model_happymarried2 <- quap(
  alist(
    happiness ~ dnorm(mu,sigma),
    mu <- a + bA*A,
    a ~ dnorm(0,1),
    bA ~ dnorm (0,2),
    sigma ~ dexp(1)
  ), data=d2
)
precis(model_happymarried2,depth=2)

compare(model_happymarried1,model_happymarried2,func=WAIC)
# Model 1 gives a better score, with the difference between models much 
# larger than the SE, even though there is a larger penalty for model 1
# However, this model is not causally correct, but it fits the data better
# still resulting in a lower WAIC.

data(foxes)
d <- list()
d$W <- standardize(foxes$weight)
d$T <- standardize(foxes$area)
d$S <- standardize(foxes$groupsize)
d$F <- standardize(foxes$avgfood)

model_F <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bF*F,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

model_T <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bT*T,
    a ~ dnorm(0,0.2),
    bT ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

model_ST <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bS*S + bT*T,
    a ~ dnorm(0,0.2),
    bS ~ dnorm(0,0.5),
    bT ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

model_FS <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bF*F + bS*S,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    bS ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

model_FST <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bF*F + bS*S + bT*T,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    bS ~ dnorm(0,0.5),
    bT ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

compare(model_F,model_T,model_ST,model_FS,model_FST,func=WAIC)
# The models include S (ST, FS and FST) are all performing similarly well
# whereas the models with only F or T are not as well performing
# with S both having an effect directly and part of a pipe, it contains the
# most data from the other parameters as well, being the most complete predictor
# with T->F being a pipe as well, also here the WAIC will be similar