library(rethinking)

# YT lecture examples
data("WaffleDivorce")
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

dd <- d[,14:16]
  
data_waffle <- alist(
  D ~ dnorm(mu,sigma),
  mu <- a + bA*A + bM*M,
  a ~ dnorm(0,0.2),
  bA ~ dnorm(0,0.5),
  bM ~ dnorm(0,0.5),
  sigma ~ dexp(1)
)

model_waffle_quap <- quap(data_waffle, data=d)
precis(model_waffle_quap)

model_waffle_mcmc <- ulam(
  data_waffle,
  data=dd,
  cores=4,
  chains=4
)
precis(model_waffle_mcmc)

# Book code examples

set.seed(1234)
num_weeks <- 1e5
positions <- rep(0, num_weeks) 
current   <- 10
for (i in 1:num_weeks) {
  # record current position 
  positions[i] <- current
  # flip coin to generate proposal
  proposal <- current + sample(c(-1, 1), size = 1)
  # now make sure he loops around the archipelago 
  if (proposal < 1) proposal <- 10
  if (proposal > 10) proposal <- 1
  # move?
  prob_move <- proposal / current
  current <- ifelse(runif(1) < prob_move, proposal, current)
}  

plot(1:100,positions[1:100])
plot(table(positions))

library(rethinking)

data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1,1,2)

model_rugged_gdp3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd
)

dd_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)
str(dd_slim)

model_rugged_gdp_mcmc <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1
)
precis(model_rugged_gdp3,depth=2)
precis(model_rugged_gdp_mcmc,depth=2)

#Diagnosing MCMC problems
y <- c(-1,1)
set.seed(1234)
wildchain <- ulam(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- alpha,
    alpha ~ dnorm(0,1000),
    sigma ~ dexp(0.0001)
  ), data=list(y=y), chains=3
)
precis(wildchain)
trankplot(wildchain)
traceplot(wildchain)
pairs(wildchain@stanfit)

set.seed(1234)
tamechain <- ulam(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- alpha,
    alpha ~ dnorm(1,10),
    sigma ~ dexp(1)
  ), data=list(y=y), chains=3
)
precis(tamechain)
trankplot(tamechain)

set.seed(1234)
y <- rnorm(100,mean=0,sd=1)

set.seed(1234)
unpredictable <- ulam(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- a1 + a2,
    c(a1,a2) ~ dnorm(0,1000),
    sigma ~ dexp(1)
  ), data=list(y=y), chains=3, cores=3
)
precis(unpredictable)
trankplot(unpredictable)
traceplot(unpredictable)

# Book practices
#9E1: the proposal distribution needs to be symmetric
#9E2: Gibbs is smarter by making adaptive proposals that make it more efficient
# The disadvantage is that in high dimensional situations, it can get stuck as
# well as when there is a correlation among parameters.
#9E3: HMC can not handle discrete parameters, because it has a physics model
# that moves around the distribution, requiring continuous parameters
#9E4: If the samples from HMC are not correlated, the n_eff will be equal
# to the actual number of samples, if they are strongly correlated, you get much
# less information per samples and n_eff decreases.
#9E5: Rhat should approach 1.00 for a proper chain
#9E6: A good HMC chain makes random moves, efficiently scanning the probability
# distribution. A bad chain will make loops, gets stuck in the wrong area, etc.
#9E7: A good chain will show different chains swapping positions all the time.
# A bad chain will show more extreme with specific chains showing extreme value
# for prolonged periods of time.
#9M1: 
library(rethinking)

data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1,1,2)
dd_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)

model_9M1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1
)
model_9M1_flat <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dunif(0,1)
  ), data=dd_slim, chains=1
)
precis(model_9M1,depth=2)
precis(model_9M1_flat,depth=2)
traceplot(model_9M1)
traceplot(model_9M1_flat)

#Plot the posterior sigma
sigma.exp <- extract.samples(model_9M1, n = 1e4)$sigma
sigma.flat <- extract.samples(model_9M1_flat, n = 1e4)$sigma
dens(sigma.exp, col="red")
dens(sigma.flat, col="blue", add=TRUE)

# In both cases, sigma moves very quickly during warm-up to the correct area
# There is not a major issues therefore, as the prior is quickly overruled by
# the likelihood of the data, giving a similar posterior, as also shown by the
# plots.

#9M2
model_9M2 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dexp(0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1
)
# The dexp(0.3) will have a broad distribution with a long tail, but will also
# make negative values for b impossible. This last one may be difficult as
# in the earlier model b[2] was negative.
precis(model_9M2,depth=2)
traceplot(model_9M2)
# The traceplot shows odd behaviour for b[2], as it can't go to negative values
# resulting in a hairy catterpillar, but only the positive half

#9M3
model_9M3a <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1, warmup=10, iter=990
)
model_9M3b <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1 # Effective warm-up/iter both 500
)
precis(model_9M3a,depth=2)
precis(model_9M3b,depth=2)
# With 500 samples, the warmup is actually too long, as only 10 samples for 
# warmup result in a higher n_eff and similar predicted parameter values. But
# can we do better by choosing something more extreme or in between?
model_9M3c <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1, warmup=1, iter=999
)
# This is no good, a bit less extreme
model_9M3d <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1, warmup=5, iter=995
)
# Still some issues, with much divergence. What about a bit higher than 10?
model_9M3e <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd_slim, chains=1, warmup=100, iter=900
)
precis(model_9M3c,depth=2)
precis(model_9M3d,depth=2)
precis(model_9M3e,depth=2)
# It seems 10 is actually not too bad

#9H1
mp <- ulam(
  alist(
    a ~ dnorm(0,1),
    b ~ dcauchy(0,1)
  ), data=list(y=1), chains=1
)
precis(mp)
# Very imprecise, with low n_eff (especially for cauchy)
traceplot(mp)
#a is a rough caterpillar while b has some odd spikes due to its long tails
#from the cauchy distribution.

samples <- extract.samples(mp)
dens(samples$a)
curve(dnorm(x, 0, 1), from = -4, to = 4, add = T, lty = 2)
# A bit rough, but still manages to capture the normal distribution rather well
dens(samples$b, ylim=c(0,0.4))
curve(dcauchy(x, 0, 1), from = -10, to = 10, add = T, lty = 2)
#Not captured as well

#9H2
data("WaffleDivorce")
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

#Make a model for divorce rate vs age at marriage
model5.1 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
#Make a model for divorce rate vs marriage rate
model5.2 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
#Make a multivariate model for divorce rate vs marriage rate & age
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
#Compare parameters of models 5.1-5.3
plot(coeftab(model5.1,model5.2,model5.3),par=c("bA","bM"))

d_slim <- list(
  A = d$A,
  D = d$D,
  M = d$M
)

#Make a model for divorce rate vs age at marriage
model5.1u <- ulam(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d_slim, chains=4, cores=4, log_lik=TRUE
)
#Make a model for divorce rate vs marriage rate
model5.2u <- ulam(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d_slim, chains=4, cores=4, log_lik=TRUE
)
#Make a multivariate model for divorce rate vs marriage rate & age
model5.3u <- ulam(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d_slim, chains=4, cores=4, log_lik=TRUE
)
#Compare parameters of models 5.1-5.3
plot(coeftab(model5.1,model5.2,model5.3,model5.1u,model5.2u,model5.3u),
     par=c("bA","bM"))
#Gives same predictions

compare(model5.1,model5.2,model5.3,func=WAIC)
compare(model5.1u,model5.2u,model5.3u,func=WAIC)

#The models with quap and stan give a similar picture, with model 1 being the 
#best predictor. Absolute values differ, but differences and ranking remain

#9H3
N <- 100
height <- rnorm(N,10,2)
leg_prop <- runif(N,0.4,0.5)
leg_left <- leg_prop*height + rnorm(N,2,0.02)
leg_right <- leg_prop*height + rnorm(N,2,0.02)
d <- data.frame(height,leg_left,leg_right)

model9H3a <- ulam(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)
  ), data=d, chains=4, cores=4, start=list(a=10,bl=0,br=0.1,sigma=1),
  log_lik=TRUE
)
model9H3b <- ulam(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)
  ), data=d, chains=4, cores=4, constraints=list(br="lower=0"),
  start=list(a=10,bl=0,br=0.1,sigma=1), log_lik=TRUE
)
posta <- extract.samples(model9H3a)
plot(bl~br,posta)
postb <- extract.samples(model9H3b)
plot(bl~br,postb)
#Model b can't have negative br values, and therefore has more density towards
# zero, giving a skewed result

#9H4
compare(model9H3a,model9H3b,func=WAIC)
compare(model9H3a,model9H3b,func=PSIS)
#Model performance is equal, model a has more effective parameters because the
# prior is slightly more informative (br>0)

#9H5-9H7
# Out of my comfort zone to make these assignments, also not that relevant?!
