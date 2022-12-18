#### From book ####
library(rethinking)
set.seed(100)
pos = replicate(1000, sum(runif(16,-1,1)))

plot(density(pos))

prod(1+runif(12,0,0.1))
growth = replicate(1000,prod(1+runif(12,0,0.1)))
dens(growth, norm.comp=TRUE)

big = replicate(1000,prod(1+runif(12,0,0.5)))
dens(big, norm.comp=TRUE)

small = replicate(1000,prod(1+runif(12,0,0.01)))
dens(small, norm.comp=TRUE)

log.big = replicate(1000,log(prod(1+runif(12,0,0.5))))
dens(log.big, norm.comp=TRUE)

library(rethinking)
data(Howell1)
d=Howell1
str(d)
precis(d)
d2=d[d$age>=18,]
dens(d2$height)
dens(d2$height[d2$male==1])
dens(d2$height[d2$male==0])
curve(dnorm(x,178,20),from=100, to=250)
curve(dunif(x,0,50),from=-10, to=60)

sample_mu=rnorm(1e4,178,20)
sample_sigma=runif(1e4,0,50)
prior_pred_h=rnorm(1e4,sample_mu,sample_sigma)
dens(prior_pred_h)
abline(v=178,col="red")

## R code 4.16-18
mu.list <- seq( from=150, to=160 , length.out=100 )
sigma.list <- seq( from=7 , to=9 , length.out=100 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum(
  dnorm( d2$height , post$mu[i] , post$sigma[i] , log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )
contour_xyz( post$mu , post$sigma , post$prob )
image_xyz( post$mu , post$sigma , post$prob )

# Ch3: samples=sample(p_grid,prob=posterior,size=1e4,replace=TRUE)
sample.rows=sample(1:nrow(post),size=1e6,replace=TRUE,prob=post$prob)
sample.mu=post$mu[sample.rows]
sample.sigma=post$sigma[sample.rows]
plot(sample.mu,sample.sigma,cex=0.5,pch=16,col=col.alpha(rangi2,0.1))
dens(sample.mu)
PI(sample.mu)

# R code 4.23-25
d3 <- sample( d2$height , size=20 )
mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
  sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
              log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )
dens( sample2.sigma , norm.comp=TRUE )

#QUAP
flist = alist(
  height ~ dnorm(mu,sigma),
  mu ~ dnorm(178,20),
  sigma ~ dunif(0,50)
)

model4.1=quap(flist,data=d2)
precis(model4.1)

start = list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)

model4.1b=quap(flist,data=d2,start=start)
precis(model4.1b)

model4.2=quap(
  alist(
    height ~ dnorm(mu,sigma),
    mu ~ dnorm(178,0.1),
    sigma ~ dunif(0,50)
  ), data=d2
)
precis(model4.2)

#LINEAR PREDICTIONS
plot (d2$height ~ d2$weight)

## R code 4.38-41
set.seed(1234)
N <- 100                   # 100 lines
a <- rnorm(N,178,20)
b <- rnorm(N,0,10)

#prior simulation using normal distribution
plot(NULL, xlim=range(d2$weight), ylim=c(-100,400), 
     xlab="weight" , ylab="height" )
abline(h=0 , lty=2)
abline(h=272, lty=1 ,lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for (i in 1:N) curve(a[i] + b[i]*(x - xbar), from=min(d2$weight), 
                     to=max(d2$weight), add=TRUE, col=col.alpha("black",0.2))

#see a lognormal distribution
b <- rlnorm(1e4, 0, 1)
dens(b, xlim=c(0,5), adj=0.1)

#prior simulation using lognormal distribution
b <- rlnorm(N, 0, 1)
plot(NULL, xlim=range(d2$weight), ylim=c(-100,400), 
     xlab="weight" , ylab="height" )
abline(h=0 , lty=2)
abline(h=272, lty=1 ,lwd=0.5)
mtext("b ~ dnorm(0,10)")
for (i in 1:N) curve(a[i] + b[i]*(x - xbar), from=min(d2$weight), 
                     to=max(d2$weight), add=TRUE, col=col.alpha("black",0.2))

#create a posterior predictive model based on adult data
model4.3 <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a +b*(weight-xbar),
    a ~ dnorm(178,20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d2
)
#check model predictions of parameters
precis(model4.3)
round(vcov(model4.3),3)
pairs(model4.3)

#plot data against posterior parameter means
plot(height~weight, data=d2, col=rangi2)
post <- extract.samples(model4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map*(x-xbar), add=TRUE)

#make models for only part of the data to get indication of variation
N <- 10
dN <- d2[1:N,]
modelN <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a +b*(weight-mean(weight)),
    a ~ dnorm(178,20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=dN
)
#check model predictions of parameters
precis(modelN)
round(vcov(modelN),3)
pairs(modelN)
#extract samples from posterior
post <- extract.samples(modelN,n=20)
plot(dN$weight, dN$height, xlim=range(d2$weight),
     ylim=range(d2$height), col=rangi2)
mtext(concat("N = ",N))
for (i in 1:20) curve(post$a[i] + post$b[i]*(x - mean(dN$weight)),
                      col=col.alpha("black",0.3), add=TRUE)

#back to the complete model calculate mu 10,000 times for a weight of 50kg
post <- extract.samples(model4.3)
mu_at_50 <- post$a + post$b*(50 - xbar)
dens(mu_at_50)
PI(mu_at_50,prob=0.95)

#calculate (by default 1000 times) the mu values for all weights in the dataset
mu <- link(model4.3)

#calculate (by default 1000 times) the mu values for predefined weight range
weight.seq <- seq(from=25,to=70,by=1)
mu <- link(model4.3,data=data.frame(weight=weight.seq))
plot(height~weight, data=d2, type="n")
for (i in 1:100) points(weight.seq, mu[i,], pch=16, col=col.alpha(rangi2,0.1))

#plot the data with a shaded PI interval
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.95)
plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq,mu.mean)
shade(mu.PI, weight.seq)

#calculate (by default 1000 times) the mu values for predefined weight range
weight.seq <- seq(from=25,to=70,by=1)
sim.heights <- sim(model4.3,data=list(weight=weight.seq),n=1e4)
height.PI <- apply(sim.heights,2,PI,prob=0.95)

#plt weights, correlation line, simulated mu and simulated height as CI
plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq,mu.mean)
shade(mu.PI,weight.seq)
shade(height.PI,weight.seq)

#ch4.5 curved lines
plot(height~weight,data=d)

#normalize data
weight_s = (d$weight-mean(d$weight))/sd(d$weight)
weight_s2 = weight_s^2

#create a posterior predictive model based on adult & child data
model4.5 <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178,20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d
)
precis(model4.5)

#calculate the mu & simulate height for the standardized weight range
weight.seq <- seq(from=-2.2,to=1.9,length.out=60)
pred_data <- list(weight_s=weight.seq, weight_s2=weight.seq^2)
mu <- link(model4.5, data=pred_data)
mu.mean = apply(mu,2,mean)
mu.PI = apply(mu,2,PI,prob=0.95)
sim.heights <- sim(model4.5, data=pred_data)
heights.PI = apply(sim.heights,2,PI,prob=0.95)

#plot lines, and PI for mu as well as height
plot(height~weight_s,data=d)
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(heights.PI, weight.seq)

#REPEAT the above for a third order polynomial
weight_s = (d$weight-mean(d$weight))/sd(d$weight)
weight_s2 = weight_s^2
weight_s3 = weight_s^3
model4.6 <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3,
    a ~ dnorm(178,20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~ dnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d
)
precis(model4.6)

#calculate the mu & simulate height for the standardized weight range
weight.seq <- seq(from=-2.2,to=1.9,length.out=60)
pred_data <- list(weight_s=weight.seq, weight_s2=weight.seq^2, 
                  weight_s3=weight.seq^3)
mu <- link(model4.6, data=pred_data)
mu.mean = apply(mu,2,mean)
mu.PI = apply(mu,2,PI,prob=0.95)
sim.heights <- sim(model4.6, data=pred_data)
heights.PI = apply(sim.heights,2,PI,prob=0.95)

#plot lines, and PI for mu as well as height
plot(height~weight_s,data=d)
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(heights.PI, weight.seq)

#plot again, but now on normal x-axis with weights
plot(height~weight_s,data=d, xaxt="n")
at <- c(-2, -1, 0, 1, 2)
labels <- at*sd(d$weight) + mean(d$weight)
axis(side=1, at=at, labels=round(labels,1))

#CHERRY BLOSOM DATA - start using B-splines
library(rethinking)
data("cherry_blossoms")
d=cherry_blossoms
precis(d)
plot(doy~year,data=d)
d2 <- d[complete.cases(d$doy),]
num_knots <- 15
knotlist <- quantile(d2$year, probs=seq(0,1,length.out=num_knots))

library(splines)
B <- bs(d2$year, knots=knotlist[-c(1,num_knots)], degree=3, intercept=TRUE)
plot(NULL,xlim=range(d2$year), ylim=c(0,1), xlab="year", ylab="basis")
for (i in 1:ncol(B)) lines(d2$year, B[,i])

model4.7 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(100,10),
    w ~ dnorm(0,10),
    sigma ~ dexp(1)
  ), data=list(D=d2$doy, B=B),
  start=list(w=rep(0,ncol(B)))
)
precis(model4.7,depth=2)

post <- extract.samples(model4.7)
w <- apply(post$w,2,mean)
plot(NULL,xlim=range(d2$year), ylim=c(-6,6), xlab="year", ylab="basis*weight")
for (i in 1:ncol(B)) lines(d2$year, w[i]*B[,i])

# Predict mu and its 97% confidence interval
mu <- link(model4.7)
mu.PI = apply(mu,2,PI,prob=0.97)
plot(d2$year, d2$doy, col=col.alpha("blue",0.3),pch=16)
shade(mu.PI, d2$year, col=col.alpha("black",0.5))

#### Lecture examples ####

# From lecture 3
library(rethinking)
data("Howell1")

plot(Howell1$height, Howell1$weight)

adult = Howell1[Howell1$age>=18,]

alpha = 0
beta = 0.5
sigma = 5
n_individuals = 100

H = runif(n_individuals,130,170)

mu = alpha + beta*H
W = rnorm(n_individuals,mu,sigma)

plot(H,W)

n_samples = 10
alpha = rnorm(n_samples,0,1)
beta = rnorm(n_samples,0,1)
plot(NULL,xlim=c(-2,2),ylim=c(-2,2),xlab="x",ylab="y")
for(i in 1:n_samples)
  abline(alpha[i],beta[i],lwd=2,col=2)

n_samples = 10
alpha = rnorm(n_samples,60,10)
beta = rlnorm(n_samples,0,1)
H_bar = 150
Hseq = seq(from=130,to=170,len=30)
plot(NULL,xlim=c(130,170),ylim=c(10,100),xlab="length(cm)",ylab="weight(kg")
for(i in 1:n_samples)
  lines(Hseq,alpha[i]+beta[i]*(Hseq-H_bar),lwd=3,col=2)

# Model
# Prior-based simulation
alpha = 70
beta = 0.5
sigma = 5
n_individuals = 100
H = runif(n_individuals,130,170)
mu = alpha + beta*(H-mean(H))
W = rnorm(n_individuals,mu,sigma)
dat=list(H=H,W=W,Hbar=mean(H))

mvalidate = quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a+b*(H-Hbar),
    a ~ dnorm(60,10),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,10)
  ), data=dat)

# Actual data
data("Howell1")
adult = Howell1[Howell1$age>=18,]
dat=list(H=adult$height,W=adult$weight,Hbar=mean(adult$height))

m_adult = quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a+b*(H-Hbar),
    a ~ dnorm(60,10),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,10)
  ), data=dat)

## Plot data
# Existing data
col2=col.alpha(2,0.8)
plot(adult$height, adult$weight, col=col2, lwd=3, cex=1.2, xlab="length(cm)",ylab="weight(kg")

# Mean and expectation for a&b
xseq=seq(from=130,to=190,len=50)
mu=link(m_adult,data=list(H=xseq,Hbar=mean(adult$height)))
lines(xseq, apply(mu,2,mean), lwd=4)
shade(apply(mu,2,PI,prob=0.99),xseq,col=col.alpha(2,0.5))

# Prediction of 89% interval (simulates observations)
W_sim=sim(m_adult,data=list(H=xseq,Hbar=mean(adult$height)))
shade(apply(W_sim,2,PI,prob=0.89),xseq,col=col.alpha(1,0.3))

# Lecture 4

# Contrast analysis
library(rethinking)
data(Howell1)
d <- Howell1
str(d)
head(d)
d <- d[d$age>=18,]
dat <- list(W=d$weight, S=d$male+1)

model_sex.1 <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a[S],
    a[S] ~ dnorm(60,10),
    sigma ~ dexp(1)
  ), data=dat
)
precis(model_sex.1, depth=2)

#Posterior mean weight
post <- extract.samples(model_sex.1)
dens(post$a[,1], col=2, xlim=c(40,50))
dens(post$a[,2], col=3, add=TRUE)

#Posterior predicted weight distribution
W_F <- rnorm(1e4,post$a[,1],post$sigma)
W_M <- rnorm(1e4,post$a[,2],post$sigma)
dens(W_F, col=2, xlim=c(20,70))
dens(W_M, col=3, add=TRUE)

str(post)

#Calculate contrast in average weight
contrast_mu <- post$a[,2]-post$a[,1]
dens(contrast_mu, xlab="posterior mean weight contrast", lwd=3)

#Contrast of predicted weight
contrast_W <- W_M - W_F
dens(contrast_W, xlab="posterior predicted weight contrast", lwd=3)
sum(contrast_W>0)/1e4

#Regression with effect of sex
dat <- list(W=d$weight, 
            H=d$height,
            Hbar=mean(d$height),
            S=d$male+1)

model_sex.2 <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a[S] + b[S]*(H-Hbar),
    a[S] ~ dnorm(60,10),
    b[S] ~ dlnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat
)

xseq <- seq(from=100,to=190,len=50)
mu_F <- link(model_sex.2, data=list(S=rep(1,50),H=xseq,Hbar=mean(d$height)))
mu_F_mean <- apply(mu_F,2,mean)
lines(xseq,mu_F_mean,col=2)
mu_M <- link(model_sex.2, data=list(S=rep(2,50),H=xseq,Hbar=mean(d$height)))
mu_M_mean <- apply(mu_M,2,mean)
lines(xseq,mu_M_mean,col=3)

mu_contrast <- mu_M - mu_F
plot(NULL, xlim=range(xseq), ylim=c(-6,8))
for(p in c(0.5,0.6,0.7,0.8,0.9,0.99))
  shade(apply(mu_contrast,2,PI,prob=p),xseq)
abline(h=0,lty=2)

#### Practices from book ####

#4M1
#prior simulation
mu_prior=rnorm(1e4,0,10)
sigma_prior=rexp(1e4,1)
#sigma_prior=runif(1e4,0,10)
prior_pred_y=rnorm(1e4,mu_prior,sigma_prior)
dens(prior_pred_y)

#4M2
model4M2=quap(
  alist(
    y ~ dnorm(mu,sigma),
    mu ~ dnorm(10,1),
    sigma ~ dexp(1)
  ), data=
)

#4M3
# y ~ Normal(mu,sigma)
# mu = alpha+beta*x,i
# alpha ~ Normal(0,10)
# beta ~ Normal(0,1)
# sigma ~ Exponential(1)

#4M4
# height ~ Normal(mu,sigma)
# mu = alpha+beta*year,i
# alpha ~ Normal(0,50)
# beta ~ Normal(0,5)
# sigma ~ Exponential(1)

#4M5/4M6 -> done

#4M7
d=Howell1
d2=d[d$age>=18,]

model4.3b <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a + b*(weight),
    a ~ dnorm(178,20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d2
)
#check model predictions of parameters
precis(model4.3b)
round(vcov(model4.3b),3)
pairs(model4.3b)
precis(model4.3)
round(vcov(model4.3),3)
pairs(model4.3)
#calculate (by default 1000 times) the mu values for all weights in the dataset
mu <- link(model4.3b)
#calculate (by default 1000 times) the mu values for predefined weight range
weight.seq <- seq(from=25,to=70,by=1)
mu <- link(model4.3b,data=data.frame(weight=weight.seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.97)
plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq,mu.mean)
shade(mu.PI, weight.seq)
#calculate (by default 1000 times) the mu values for predefined weight range
weight.seq <- seq(from=25,to=70,by=1)
sim.heights <- sim(model4.3b,data=list(weight=weight.seq))
height.PI <- apply(sim.heights,2,PI,prob=0.97)
#plt weights, correlation line, simulated mu and simulated height as CI
plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq,mu.mean)
shade(mu.PI,weight.seq)
shade(height.PI,weight.seq)

#4M8
d=cherry_blossoms
d2 <- d[complete.cases(d$doy),]
set.seed(1234)
num_knots <- 25
knotlist <- quantile(d2$year, probs=seq(0,1,length.out=num_knots))
B <- bs(d2$year, knots=knotlist[-c(1,num_knots)], degree=3, intercept=TRUE)
model4.7b <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(100,10),
    w ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=list(D=d2$doy, B=B),
  start=list(w=rep(0,ncol(B)))
)
post <- extract.samples(model4.7b)
w <- apply(post$w,2,mean)
mu <- link(model4.7b)
mu.PI = apply(mu,2,PI,prob=0.97)
plot(d2$year, d2$doy, col=col.alpha("blue",0.3),pch=16)
shade(mu.PI, d2$year, col=col.alpha("black",0.5))

#4H1
d=Howell1
d2=d[d$age>=18,]
xbar <- mean(d2$weight)
model4.3c <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a +b*(weight-xbar),
    a ~ dnorm(150,30),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d2
)
new_weight <- c(46.95, 43.72, 64.78, 32.59, 54.63)
pred_height <- link(model4.3c, data = data.frame(weight = new_weight))
expected <- apply(pred_height, 2, mean)
interval <- apply(pred_height, 2, PI, prob = 0.89)

data.frame(
  individual = 1:5,
  weight = new_weight,
  expected = expected,
  lower = interval[1, ],
  upper = interval[2, ]
)

#4H2
d=Howell1
d2=d[d$age<18,]
xbar <- mean(d2$weight)
model4.3d <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a +b*(weight-xbar),
    a ~ dnorm(100,50),
    b ~ dlnorm(0,10),
    sigma ~ dunif(0,50)
  ), data=d2
)
precis(model4.3d)
mu <- link(model4.3d)
weight.seq <- seq(from=4,to=45,by=1)
mu <- link(model4.3d,data=data.frame(weight=weight.seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.97)
sim.heights <- sim(model4.3d,data=list(weight=weight.seq))
height.PI <- apply(sim.heights,2,PI,prob=0.97)
plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq,mu.mean)
shade(mu.PI,weight.seq)
shade(height.PI,weight.seq)

#4H3
d=Howell1
model4.3e <- quap (
  alist(
    height ~ dnorm (mu, sigma),
    mu <- a + b*(log(weight)),
    a ~ dnorm(150,30),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d
)
precis(model4.3e)


mu <- link(model4.3e)
weight.seq <- seq(from=4,to=75,by=1)
mu <- link(model4.3e,data=data.frame(weight=weight.seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.97)
sim.heights <- sim(model4.3e,data=list(weight=weight.seq))
height.PI <- apply(sim.heights,2,PI,prob=0.97)
plot(height~weight, data=d, col=col.alpha(rangi2,0.5))
lines(weight.seq,mu.mean)
shade(mu.PI,weight.seq)
shade(height.PI,weight.seq)

#4H4
d=Howell1
plot(height~weight,data=d)
weight_s = (d$weight-mean(d$weight))/sd(d$weight)
weight_s2 = weight_s^2

#prior simulation
set.seed(1234)
N=40
a_prior=rnorm(N,130,35)
b1_prior=rlnorm(N,1,1)
b2_prior=rlnorm(N,0,0.5)
plot(NULL, xlim=range(d2$weight), ylim=c(-100,400), 
     xlab="weight" , ylab="height" )
abline(h=0 , lty=2)
abline(h=272, lty=1 ,lwd=0.5)
mtext("b1/b2 ~ dlnorm(0,1) & dnorm(0,1")
for (i in 1:N) curve(a_prior[i] + b1_prior[i]*x + b2_prior[1]*x^2, from=min(d2$weight), 
                     to=max(d2$weight), add=TRUE, col=col.alpha("black",0.2))

#4H5
d=cherry_blossoms
d2 <- na.omit(d)
plot(d2$doy~d2$temp)
xbar <- mean(d2$temp)
model4.8 <- quap (
  alist(
    doy ~ dnorm (mu, sigma),
    mu <- a + b*(temp-xbar),
    a ~ dnorm(100,30),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ), data=d2
)
precis(model4.8)

mu <- link(model4.8)
temp.seq <- seq(from=4,to=9,by=0.1)
mu <- link(model4.8,data=data.frame(temp=temp.seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.97)
sim.doy <- sim(model4.8,data=list(temp=temp.seq))
doy.PI <- apply(sim.doy,2,PI,prob=0.97)
plot(doy~temp, data=d, col=col.alpha(rangi2,0.5))
lines(temp.seq,mu.mean)
shade(mu.PI,temp.seq)
shade(doy.PI,temp.seq)

#4H6
#TOO DIFFICULT

#4H7
d=cherry_blossoms
d2 <- na.omit(d)
set.seed(1234)
num_knots <- 15
knotlist <- quantile(d2$year, probs=seq(0,1,length.out=num_knots))
B <- bs(d2$year, knots=knotlist[-c(1,num_knots)], degree=3, intercept=FALSE)
model4.9 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- B %*% w,
    w ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=list(D=d2$doy, B=B),
  start=list(w=rep(0,ncol(B)))
)
post <- extract.samples(model4.9)
w <- apply(post$w,2,mean)
mu <- link(model4.9)
mu.PI = apply(mu,2,PI,prob=0.97)
plot(d2$year, d2$doy, col=col.alpha("blue",0.3),pch=16)
shade(mu.PI, d2$year, col=col.alpha("black",0.5))

#### Online homework ####
# 1
# Actual data
data("Howell1")
adult = Howell1[Howell1$age>=18,]
dat=list(H=adult$height,W=adult$weight,Hbar=mean(adult$height))

m_adult = quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a+b*(H-Hbar),
    a ~ dnorm(60,10),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,10)
  ), data=dat)

precis(m_adult)

dat2=list(H=c(140,160,175),Hbar=mean(adult$height))
w_sim=sim(m_adult,data=dat2)
ew=apply(w_sim,2,mean)
h_ci=apply(w_sim,2,PI,prob=0.89)

# 2
# Actual data
data("Howell1")
children = Howell1[Howell1$age<13,]
dat=list(A=children$age, W=children$weight)

# Simulate data
n_samples = 10
alpha = rnorm(n_samples,5,1)
beta = rlnorm(n_samples,0,1)
plot(NULL,xlim=range(children$age),ylim=range(children$weight),xlab="age(years)",ylab="weight(kg)")
for(i in 1:n_samples)
  abline(alpha[i], beta[i], lwd=3, col=2)

m_child_AW = quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a+b*A,
    a ~ dnorm(5,1),
    b ~ dlnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat)

precis(m_child_AW)

# 3
# Actual data
data("Howell1")
children = Howell1[Howell1$age<13,]
dat=list(A=children$age, W=children$weight, S=children$male+1)

# Model
m_child_AWS = quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a[S]+b[S]*A,
    a[S] ~ dnorm(5,1),
    b[S] ~ dlnorm(0,1),
    sigma ~ dexp(1)
  ), data=dat)

# Plot age vs weight
plot(children$age, children$weight, lwd=3, col=ifelse(children$male==1,4,2) , xlab="age (years)" , ylab="weight (kg)" )

# Separate plots for sexes
Aseq <- 0:12

# girls
muF <- link(m_child_AWS,data=list(A=Aseq,S=rep(1,13))) # Female mu prediction using sex-weighted model, predicting only for girls <13y
lines(Aseq, apply(muF,2,mean), lwd=3 , col=2 ) # Add a line for age <13, using the female prediction in color 2=red
shade(apply(muF,2,PI,0.99), Aseq, col=col.alpha(2,0.5) ) # Add a shade for 99% interval, using female prediction with 50% alpha

# boys
muM <- link(m_child_AWS,data=list(A=Aseq,S=rep(2,13)))
lines(Aseq, apply(muM,2,mean), lwd=3, col=4 )
shade(apply(muM,2,PI,0.99), Aseq, col=col.alpha(4,0.5) )

# contrast at each age
Aseq <- 0:12 # Only age <13y
mu1 <- sim(m_child_AWS,data=list(A=Aseq,S=rep(1,13))) # Simulate a mu for girls
mu2 <- sim(m_child_AWS,data=list(A=Aseq,S=rep(2,13))) # Simulate a mu for boys
mu_contrast <- mu1 # Perform a contrast analysis for girls by starting from the girl mu values
for (i in 1:13 ) mu_contrast[,i] <- mu2[,i] - mu1[,i] # Calculate for each age the difference in mu for girls and boys

# plot the mu contrast
plot(NULL , xlim=c(0,13), ylim=c(-15,15), xlab="age", ylab="weight difference (boys-girls)") # create axes
for (p in c(0.5,0.67,0.89,0.99)) 
  shade(apply(mu_contrast,2,PI,prob=p), Aseq) # Perform the analysis for different intervals, applying a shade for each
abline(h=0,lty=2,lwd=2) # Add a horizontal line for absence of an effect
for (i in 1:13) points(mu_contrast[1:1000,i], col=ifelse(mu_contrast[1:1000,i]>0,4,2), lwd=3) # plot the 1000 simulated datapoints, colour by sex

# 4
data(Oxboys)
d <- Oxboys

d$delta <- NA
for ( i in 1:nrow(d) ) {
  if ( d$Occasion[i] > 1 ) d$delta[i] <- d$height[i] - d$height[i-1]
}
d <- d[ !is.na(d$delta) , ]

dat=list(HD=d$delta)

# simulation from priors
n <- 1e3
alpha <- rnorm(n,0,0.1)
sigma <- rexp(n,3)
delta_sim <- rlnorm(n,alpha,sigma)
dens(delta_sim)

# Model
m_height_delta = quap(
  alist(
    HD ~ dlnorm(a,sigma),
    a ~ dnorm(0,0.1),
    sigma ~ dexp(3)
  ), data=dat)

# compute posterior sum of 8 increments
post <- extract.samples(m_height_delta)

dsim <- rlnorm(1e3,post$a,post$sigma)
dens(dsim)

inc_sum <- sapply(1:1000, function(s) sum(rlnorm(8,post$a[s],post$sigma[s])))
dens(inc_sum)