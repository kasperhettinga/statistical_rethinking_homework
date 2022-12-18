#### From Book ####

library(rethinking)

p_grid = seq (from =0, to=1, length.out=1000)
prob_p = rep(1,1000)
prob_data=dbinom(6, size=9, prob=p_grid)
posterior = prob_p*prob_data
posterior = posterior/sum(posterior)
plot(posterior)

samples=sample(p_grid,prob=posterior,size=1e5,replace=TRUE)
plot(samples)
dens(samples)

w=rbinom(1e7,size=9,prob=samples)
dens(w)

##### From lecture 2 ####
library(rethinking)

p_grid = seq(from=0, to=1, length.out=1000)
prob_p = rep(1, 1000)
prob_data = dbinom(6, size=9, prob=p_grid)
posterior = prop_data*prob_p
posterior = posterior/sum(posterior)

p_grid = seq(from=0, to=1, length.out=1000)
prob_p = dbeta(p_grid, 3, 1)
prob_data = dbinom(6, size=9, prob=p_grid)
posterior = prop_data*prob_p
posterior = posterior/sum(posterior)


#### Practice from book ####

p_grid = seq (from =0, to=1, length.out=1000)
prior = rep(1,1000)
likelihood = dbinom(6, size=9, prob=p_grid)
posterior = likelihood*prior
posterior = posterior/sum(posterior)

set.seed(100)
samples=sample(p_grid,prob=posterior,size=1e9,replace=TRUE)

#EASY

# 1-3
sum(samples<0.2)/1e4
sum(samples>0.8)/1e4
sum(samples>0.2&samples<0.8)/1e4

#4-5
dens(samples)
quantile(samples,0.2)
quantile(samples,0.8)

#6-7
HPDI(samples, prob=0.66)
PI(samples, prob=0.66)

#MEDIUM
p_grid = seq (from =0, to=1, length.out=1000)
prior = rep(1,1000)
likelihood = dbinom(8, size=15, prob=p_grid)
posterior = likelihood*prior
posterior = posterior/sum(posterior)

#1
plot(posterior~p_grid, type="l")

#2
set.seed(100)
samples=sample(p_grid,prob=posterior,size=1e4,replace=TRUE)
dens(samples)
HPDI(samples, prob=0.9)

#3
w=rbinom(1e4, size=15, prob=samples)
simplehist(w)
sum(w==8)/1e4

#4
w=rbinom(1e4, size=9, prob=samples)
simplehist(w)
sum(w==6)/1e4

#5
p_grid = seq (from =0, to=1, length.out=1000)
prior_new = c(rep(0,500),rep(1,500))
plot(prior_new)
likelihood = dbinom(8, size=15, prob=p_grid)
posterior = likelihood*prior_new
posterior = posterior/sum(posterior)
plot(posterior~p_grid, type="l")
set.seed(100)
samples=sample(p_grid,prob=posterior,size=1e4,replace=TRUE)
dens(samples)
HPDI(samples, prob=0.9)
w=rbinom(1e4, size=15, prob=samples)
simplehist(w)
sum(w==8)/1e4
w=rbinom(1e4, size=9, prob=samples)
simplehist(w)
sum(w==6)/1e4

#HARD
birth1 <- c(1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,0,
            0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,
            1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,
            1,0,1,1,1,0,1,1,1,1)
birth2 <- c(0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,
            1,1,1,0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,1,
            0,0,0,1,1,1,0,0,0,0)

#1
p_grid = seq (from =0, to=1, length.out=1000)
prior = rep(1,length(p_grid))
boys=sum(birth1)+sum(birth2)
likelihood = dbinom(boys, size=200, prob=p_grid)
posterior = likelihood*prior
posterior = posterior/sum(posterior)
plot(posterior~p_grid, type="l")
abline( v=0.5 , lty=2 )

p_grid[which.max(posterior)]

#2
samples=sample(p_grid,prob=posterior,size=1e4,replace=TRUE)
dens(samples)
HPDI(samples, prob=0.5)
HPDI(samples, prob=0.89)
HPDI(samples, prob=0.97)

#3
birth_sim=rbinom(1e4, size=200, prob=samples)
simplehist(birth_sim)
abline(v=111 , col="red")

#4
birth_sim=rbinom(1e4, size=100, prob=samples)
simplehist(birth_sim)
abline(v=sum(birth1) , col="red")

#5
birth_01 = birth2[birth1==0]
birth_01_sim = rbinom(1e4, size=length(birth_01), prob=samples)
simplehist(birth_01_sim)
abline(v=sum(birth_01) , col="red")

#### Homework from github ####
p_grid = seq(from=0, to=1, length.out=1000)
prob_p = c(rep(0, 500), rep(1, 500))
prob_data = dbinom(4, size=6, prob=p_grid)
posterior = prob_data*prob_p
posterior = posterior/sum(posterior)
samples = sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
w = rbinom(1e4, size=9, prob=samples)

plot(posterior)

HPDI(samples)
PI(samples)

W = rbinom (1, size=20, prob=0.7*0.8)
grid_p = seq(from=0, to=1, len=100)
pr_p = dbeta(grid_p,1,1)
prW = dbinom(W,20,grid_p*0.8)
post=pr_p*prW

plot(post)

post_bad = dbinom(W,20,grid_p)
plot(grid_p, post, type="l", lwd=4)
lines(grid_p, post_bad, col=2, lwd=4)