#### From Book ####

## Grid approxination example
#define grid
p_grid = seq (from =0, to=1, length.out=20)
plot(p_grid)

# define prior
prior = rep (1, 20)

# likehood of each discrete p
likelihood = dbinom(6, size=9, prob=p_grid)
plot(likelihood)

# calculate posterior
unstd.posterior = prior * likelihood
plot(unstd.posterior)

posterior = unstd.posterior/(sum(unstd.posterior))
plot(posterior)

# plot posterior versus p
plot(p_grid, posterior, type="b")

# what about a prior saying only p>0.5 is possible
prior = ifelse(p_grid<0.5, 0, 1)
plot(prior)

unstd.posterior = prior * likelihood
posterior = unstd.posterior/(sum(unstd.posterior))
plot(p_grid, posterior, type="b")

# what about a prior peak at p=0.5
prior = exp(-5*abs(p_grid-0.5))
plot(prior)

unstd.posterior = prior * likelihood
posterior = unstd.posterior/(sum(unstd.posterior))
plot(p_grid, posterior, type="b")

##QUAP example
library(rethinking)
globe.qa <- quap(
  alist(
    W ~ dbinom(W+L, p),
    p ~ dunif(0,1)
  ),
  data=list(W=6,L=3))

#display main outcomes
precis(globe.qa)

# analytical calculation
W=6
L=3
curve(dbeta(x,W+1,L+1),from=0,to=1)
curve(dnorm(x, 0.67, 0.16), lty=2, add=TRUE)

# MCMC calculation of the posterior
n=100000
p=rep(NA,n)
p[1]=0.5
W=6
L=3
for (i in 2:n){
  p_new = rnorm(1, p[i-1], 0.1)
  if (p_new<0) p_new=abs(p_new)
  if (p_new>1) p_new=2-p_new
  q0=dbinom(W, W+L, p[i-1])
  q1=dbinom(W, W+L, p_new)
  p[i]=ifelse(runif(1)<q1/q0, p_new, p[i-1])
}

dens(p,xlim=c(0,1))
curve(dbeta(x,W+1,L+1), lty=2, add=TRUE)

#### Book Practices ####

# 2M1
p_grid <- seq(from=0 , to=1 , length.out=100)
likelihood <- dbinom(5 , size=7 , prob=p_grid)
prior <- rep(1,100) # uniform prior
posterior <- likelihood * prior
posterior <- posterior / sum(posterior) # standardize
plot(posterior)

# 2M2
p_grid <- seq(from=0 , to=1 , length.out=100)
likelihood <- dbinom(5 , size=7 , prob=p_grid)
prior <- ifelse(p_grid<0.5,0,1)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior) # standardize
plot(posterior)

# 2M3
p_grid <- seq(from=0 , to=1 , length.out=100)
likelihood <- dbinom(5 , size=7 , prob=p_grid)
prior <- ifelse(p_grid<0.5,0,1)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior) # standardize
plot(posterior)

# 2H1
# Species 1: 10% twins
# Species 2: 20% twins
# You have a twin: change that mom is species 2, is twice the chance (2/3)
# Specie 1: 1/3*1/10=1/30
# Specie 2: 2/3*2/10=4/30
# Total 5/30=1/6
# Specie 1: 1/3*
# singleton: 9vs8 : 9/17 vs 8/17 => 1/3/
