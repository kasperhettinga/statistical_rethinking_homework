# Book examples
setwd("C:/Users/kaspe/OneDrive - Wageningen University & Research/Courses/Statistical rethinking")
library(rethinking)
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

precis(model5.1)

#Plot the priors of the model
set.seed(1234)
prior <- extract.prior(model5.1)
mu <- link(model5.1, post=prior, data=list(A=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-2,2))
for(i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha("black",0.4))

#Calculate mu (mean and PI) based on the posterior
A_seq <- seq(from=-3, to=3.2, length.out=30)
mu <- link(model5.1,data=list(A=A_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(D~A,data=d,col=rangi2)
lines(A_seq,mu.mean,lwd=2)
shade(mu.PI,A_seq)

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

#Calculate mu (mean and PI) based on the posterior
M_seq <- seq(from=-3, to=3, length.out=30)
mu <- link(model5.2,data=list(M=M_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(D~M,data=d,col=rangi2)
lines(M_seq,mu.mean,lwd=2)
shade(mu.PI,M_seq)

#Draw the DAG for model 5.1
library(dagitty)
dag5.1 <- dagitty("dag{A->D;A->M;M->D}")
coordinates(dag5.1) <- list(x=c(A=0,D=1,M=2), y=c(A=0,D=1,M=0))
drawdag(dag5.1)
DMA_dag5.1 <- dagitty("dag{D<-A->M->D}")
impliedConditionalIndependencies(DMA_dag5.1)
DMA_dag5.2 <- dagitty("dag{D<-A->M}")
impliedConditionalIndependencies(DMA_dag5.2)

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
precis(model5.3)

#Compare parameters of models 5.1-5.3
plot(coeftab(model5.1,model5.2,model5.3),par=c("bA","bM"))

#Model age on mariage rate and calculate the residuals of the model
model5.4 <- quap(
  alist(
    M ~ dnorm(mu,sigma),
    mu <- a + bA*A,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

mu <- link(model5.4)
mu.mean <- apply(mu,2,mean)
mu.res <- d$M - mu.mean

A_seq <- seq(from=-3, to=3.2, length.out=30)
mu <- link(model5.4,data=list(A=A_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(M~A,data=d,col=rangi2)
lines(A_seq,mu.mean,lwd=2)

d$mu.res <- mu.res

#Model divorce rate on residuals of M vs A model
model5.5 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + br*mu.res,
    a ~ dnorm(0,0.2),
    br ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

res_seq <- seq(from=-2, to=2, length.out=40)
mu <- link(model5.5,data=list(mu.res=res_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(D~mu.res,data=d,col=rangi2)
lines(res_seq,mu.mean,lwd=2)
shade(mu.PI,res_seq)
abline(h=0 , lty=1)
abline(v=0 , lty=2)

#Model mariate rate on age and calculate the residuals of the model
model5.6 <- quap(
  alist(
    A ~ dnorm(mu,sigma),
    mu <- a + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

mu <- link(model5.4)
mu.mean <- apply(mu,2,mean)
mu.res <- d$A - mu.mean

M_seq <- seq(from=-3, to=3, length.out=30)
mu <- link(model5.6,data=list(M=M_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(A~M,data=d,col=rangi2)
lines(M_seq,mu.mean,lwd=2)

d$mu.res <- mu.res

#Model divorce rate on residuals of M vs A model
model5.7 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + br*mu.res,
    a ~ dnorm(0,0.2),
    br ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

res_seq <- seq(from=-5, to=5, length.out=40)
mu <- link(model5.7,data=list(mu.res=res_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(D~mu.res,data=d,col=rangi2)
lines(res_seq,mu.mean,lwd=2)
shade(mu.PI,res_seq)
abline(h=0 , lty=1)
abline(v=0 , lty=2)

#5152 Posterior prediction plots
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

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

mu <- link(model5.3)
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)

D_sim <- sim(model5.3, n=1e4)
D_mean <- apply(D_sim,2,mean)
D_PI <- apply(D_sim,2,PI)

plot(mu_mean ~ d$D, col=rangi2, ylim=range(mu_PI), xlab="Observed divorce",
     ylab="Predicted divorce")
abline(a=0, b=1, lty=2)
for(i in 1:nrow(d)) lines(rep(d$D[i],2), mu_PI[,i],  col=rangi2)

plot(D_mean ~ d$D, col=rangi2, ylim=range(D_PI), xlab="Observed divorce",
     ylab="Simulated divorce")
abline(a=0, b=1, lty=2)
for(i in 1:nrow(d)) lines(rep(d$D[i],2), D_PI[,i],  col=rangi2)

# Detect spurious association
N <- 100
x_real <- rnorm(N)
x_spur <- rnorm(N, x_real)
y <- rnorm(N,x_real)
d <- data.frame(y,x_real,x_spur)
pairs(d)

model5.spur <- quap(
  alist(
    y ~ dnorm(mu,sigma),
    mu <- a + br*x_real + bs*x_spur,
    a ~ dnorm(0,0.5),
    bs ~ dnorm(0,0.5),
    br ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model5.spur)

#5153 Counterfactual plots
data("WaffleDivorce")
d <- list()
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)

model5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M,sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0,0.2),
    bAM ~ dnorm(0,0.5),
    sigma_M ~ dexp(1)
  ), data=d
)
precis(model5.3_A)

A_seq <- seq(from=-2, to=2, length.out=40)

sim_data <- data.frame(A=A_seq)
s <- sim(model5.3_A, data=sim_data, vars=c("M","D"))

plot(sim_data$A, colMeans(s$D), ylim=c(-2,2), type="l", xlab="manipulated A", 
     ylab="counterfactual D")
shade(apply(s$D,2,PI), sim_data$A)

plot(sim_data$A, colMeans(s$M), ylim=c(-2,2), type="l", xlab="manipulated A", 
     ylab="counterfactual M")
shade(apply(s$M,2,PI), sim_data$A)

sim_data <- data.frame(M=seq(from=-2, to=2, length.out=40), A=0)
s <- sim(model5.3_A, data=sim_data, vars="D")
plot(sim_data$M, colMeans(s), ylim=c(-2,2), type="l", xlab="manipulated M", 
     ylab="counterfactual D")
shade(apply(s,2,PI), sim_data$M)

#5.2 Milk data - masked relationships
library(rethinking)
data("milk")
d <- milk

d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))

model_milk.1a <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bN*N,
    a ~ dnorm(0,1),
    bN ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)

dcc <- d[complete.cases(d$K,d$N,d$M),]

model_milk.1b <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bN*N,
    a ~ dnorm(0,1),
    bN ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=dcc
)
prior <- extract.prior(model_milk.1b)
xseq <- c(-2,2)
mu <- link(model_milk.1b,post=prior,data=list(N=xseq))
plot(NULL,xlim=xseq,ylim=xseq)
for(i in 1:50) lines(xseq,mu[i,],col=col.alpha("black",0.3))

model_milk.1c <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bN*N,
    a ~ dnorm(0,0.2),
    bN ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
prior <- extract.prior(model_milk.1c)
xseq <- c(-2,2)
mu <- link(model_milk.1c,post=prior,data=list(N=xseq))
plot(NULL,xlim=xseq,ylim=xseq)
for(i in 1:50) lines(xseq,mu[i,],col=col.alpha("black",0.3))

precis(model_milk.1c)

xseq <- seq(from=min(dcc$N)-0.15,to=max(dcc$N)+0.15,length.out=30)
mu <- link(model_milk.1c,data=list(N=xseq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.95)
plot(K~N, data=dcc)
lines(xseq,mu_mean,lwd=2)
shade(mu_PI,xseq)

model_milk.2 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
precis(model_milk.2)

xseq <- seq(from=min(dcc$M)-0.15,to=max(dcc$M)+0.15,length.out=30)
mu <- link(model_milk.2,data=list(M=xseq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.95)
plot(K~M, data=dcc)
lines(xseq,mu_mean,lwd=2)
shade(mu_PI,xseq)

model_milk.3 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0,0.2),
    bN ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dcc
)
precis(model_milk.3)

pairs(~K + M + N, dcc)

xseq <- seq(from=min(dcc$M)-0.15,to=max(dcc$M)+0.15,length.out=30)
mu <- link(model_milk.3,data=data.frame(M=xseq, N=0))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.95)
plot(NULL, xlim=range(dcc$M), ylim=range(dcc$K))
lines(xseq,mu_mean,lwd=2)
shade(mu_PI,xseq)

xseq <- seq(from=min(dcc$N)-0.15,to=max(dcc$N)+0.15,length.out=30)
mu <- link(model_milk.3,data=data.frame(N=xseq, M=0))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.95)
plot(NULL, xlim=range(dcc$N), ylim=range(dcc$K))
lines(xseq,mu_mean,lwd=2)
shade(mu_PI,xseq)

data(milk)
d <- milk
levels(d$clade)
d$clade_id <- as.integer(d$clade)
d$K <- standardize(d$kcal.per.g)

model_milk.4 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
labels <- paste("a[" , 1:4 , "]:", levels(d$clade), sep="")
plot(precis(model_milk.4,depth=2,pars="a"), labels=labels, xlab="Expected kcal")

#### PRACTICE (BOOK) ####
#5E1: only 2 and 4, as they have different beta & parameters
#5E2: D = u + sigma; mu = alpha + betaL*L + betaPD*PD
#5E3: PhD_time mu + sigma; u = apha + betaF*F + betaLS*LS
# with betaF & betaLS positive e.g. being from a log normal distribution
#5E4: 1, 3, 4 & 5 are equivalent, because with 3 indicators, you can always
# calculate the fourth

#5M1: Students following DST course ending up in dairy industry, relation
# will be gone when correcting for prior interest in dairy
#5M2: Students following DST and ending up in vegan industry would be masked
# by the negative correlation between being vegan and following the course
#5M3: People divorcing more will subsequently possible remarry, leading
# to a higher marriage rate. Use 1sr vs remarriage as a separate parameter?
#5M4: 

library(rethinking)
data("WaffleDivorce")
lsd <- read.csv("lsd_per_state.csv", header=TRUE, sep=";")
d <- WaffleDivorce
d$lsd <- lsd$LSD_perc

d$A <- standardize(d$MedianAgeMarriage)
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$L <- standardize(d$lsd)

model5M4 <- quap(
  alist(
    ## Combined effect of A, M & L on D
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A + bM*M + bL*L,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    bL ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model5M4)

#5M5: Obesity ~ dnorm(mu,sigma); mu <- alpa + bD*Driving + bE*Excrcise

#5H1: D dnd M cond. A

#5H2
data("WaffleDivorce")
d <- list()
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)

model5H2 <- quap(
  alist(
    ## M -> A
    A ~ dnorm(mu_A,sigma_A),
    mu_A <- a_A + bM*M,
    a_A ~ dnorm(0,1),
    bM ~ dnorm(0,1),
    sigma_A ~ dexp(1),
    ## A -> D
    D ~ dnorm(mu_D,sigma_D),
    mu_D <- a_D + bA*A,
    a_D ~ dnorm(0,1),
    bA ~ dnorm(0,1),
    sigma_D ~ dexp(1)
  ), data=d
)
precis(model5H2)

sim_data <- data.frame(A=seq(from=-2, to=2, length.out=40), M=0)
s <- sim(model5H2, data=sim_data, vars="D")
plot(sim_data$A, colMeans(s), ylim=c(-2,2), type="l", xlab="A at fixed M", 
     ylab="counterfactual D")
shade(apply(s,2,PI), sim_data$A)

#5H3
library(rethinking)
data("milk")
d <- milk

d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))

dcc <- d[complete.cases(d$K,d$N,d$M),]

model_milk.5H3 <- quap(
  alist(
    # M -> K <- N
    K ~ dnorm(mu,sigma),
    mu <- a + bN*N + bM*M,
    a ~ dnorm(0,0.2),
    bN ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1),
    # M -> N
    N ~ dnorm(muN,sigmaN),
    muN <- aN + bMN*M,
    aN ~ dnorm(0,0.2),
    bMN ~ dnorm(0,0.5),
    sigmaN ~ dexp(1)
  ), data=dcc
)

sim_data <- data.frame(M=c(mean(d$M),2*mean(d$M)))
sim <- sim(model_milk.5H3,data=sim_data,vars=c("N","K"))
mean(sim$K[,2]-sim$K[,1])

#5H4
data("WaffleDivorce")
d <- list()
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$S <- WaffleDivorce$South+1

library(dagitty)
dag5H4 <- dagitty("dag{S->A;A->D;S->D}")
drawdag(dag5H4)

model5H4NS <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA*A,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model5H4NS)

model5H4S <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a[S] + bA*A,
    a[S] ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model5H4S,depth=2)

post <- extract.samples(model5H4S)
dens(post$a[,1], col=2, xlim=c(-0.5,0.7))
dens(post$a[,2], col=3, add=TRUE)

D_NS <- rnorm(1e4,post$a[,1],post$sigma)
D_S <- rnorm(1e4,post$a[,2],post$sigma)
dens(D_NS, col=2)
dens(D_S, col=3, add=TRUE)

contrast_mu <- post$a[,2]-post$a[,1]
dens(contrast_mu, xlab="posterior mean weight contrast", lwd=3)