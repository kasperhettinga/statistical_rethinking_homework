#### Book Examples ####
library(rethinking)

# multicollinearity in leg and body length
N <- 100
height <- rnorm(N,10,2)
leg_prop <- runif(N,0.4,0.5)
leg_left <- leg_prop*height + rnorm(N,2,0.02)
leg_right <- leg_prop*height + rnorm(N,2,0.02)
d <- data.frame(height,leg_left,leg_right)

model6.1 <- quap(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)
  ), data=d
)
plot(precis(model6.1))

post <- extract.samples(model6.1)
plot(bl~br,post)

# MC in milk data
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat)
d$L <- standardize(d$perc.lactose)

model_milk1 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bF*F,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_milk1)

model_milk2 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bL*L,
    a ~ dnorm(0,0.2),
    bL ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_milk2)

model_milk3 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bF*F + bL*L,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    bL ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_milk3)

# Post-treatment bias
N <- 100
h0 <- rnorm(N,10,2)
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4)
h1 <- h0 + rnorm(N,5 - 3*fungus)
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d)

model_plant1 <- quap(
  alist(
    h1 ~ dnorm(mu,sigma),
    mu <- h0*p,
    p ~ dlnorm(1,0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_plant1)

model_plant2 <- quap(
  alist(
    h1 ~ dnorm(mu,sigma),
    mu <- h0*p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(1,0.25),
    bt ~ dnorm(0,0.5),
    bf ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_plant2)

library(dagitty)
plant_dag <- dagitty("dag{H0->H1 F->H1 T->F}")
drawdag(plant_dag)
impliedConditionalIndependencies(plant_dag)

# Post-treatment bias II
N <- 1000
h0 <- rnorm(N,10,2)
treatment <- rep(0:1, each=N/2)
moisture <- rbern(N)
fungus <- rbinom(N, size=1, prob=0.5-treatment*0.4+moisture*0.4)
h1 <- h0 + rnorm(N,5 - moisture*3)
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)

model_plant3 <- quap(
  alist(
    h1 ~ dnorm(mu,sigma),
    mu <- h0*p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(1,0.25),
    bt ~ dnorm(0,0.5),
    bf ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_plant3)

model_plant4 <- quap(
  alist(
    h1 ~ dnorm(mu,sigma),
    mu <- h0*p,
    p <- a + bt*treatment,
    a ~ dlnorm(1,0.25),
    bt ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_plant4)

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

# Haunted DAG
N <- 200
bGP <- 1
bGC <- 0
bPC <- 1
bU <- 2

U <- 2*rbern(N,0.5)-1
G <- rnorm(N)
P <- rnorm(N,bGP*G+bU*U)
C <- rnorm(N,bGC*G+bPC*P+bU*U)
d <- data.frame(C=C,P=P,G=G,U=U)

model_parentchild <- quap(
  alist(
    C ~ dnorm(mu,sigma),
    mu <- a+b_PC*P+b_GC*G,
    a ~ dnorm(0,1),
    c(b_PC,b_GC) ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_parentchild)

# Use a DAG for finding and selecting backdoors in confounded analyses
library(dagitty)
dag6.1 <- dagitty("dag{
                  U [unobserved]
                  X -> Y
                  X <- U <- A -> C -> Y
                  U -> B <- C
                  }")
adjustmentSets(dag6.1, exposure="X", outcome="Y")

#### Practices from book ####

#6E1: colinearity, post-treatment bias & colliding bias
#6E2a: colinearity: measuring clinical mastitis signs & SCC -> relate milk
#6E2b: post-treatment bias: Mastitis treatment, bacterial count -> milk
#6E2c: animal health and animal feed both influence milk composition, but
# are not related to each other
#6E3: if you have a fork, pipe, colliding or descendant relationship
# fork can be solved by conditioning on the joint point in the fork
# pipe makes start & end unrelated, ones mediator is included
# colliding makes a correlation arise due to an unseen parameter
# descendant can relate to any of the above, but less strong
#6E4: Example 6E2c is a collider bias. conditioning on milk can make one
# believe that a relation exists.

#6M1: 
library(dagitty)
dag6M1 <- dagitty("dag{
                  U [unobserved]
                  V [unobserved]
                  X -> Y
                  X <- U <- A -> C -> Y
                  U -> B <- C
                  Y <- V -> C
                  }")
coordinates(dag6M1) <- list(x=c(U=0, X=0, A=1, B=1, C=2, Y=2, V=2.5),
                            y=c(U=1, X=3, A=0, B=2, C=1, Y=3, V=2))
drawdag(dag6M1)

#X-Y; X-U-B-C-Y; X-U-B-C-V-Y; X-U-A-C-Y; X-U-A-C-V-Y (5 paths)
adjustmentSets(dag6M1, exposure="X", outcome="Y")
# You should condition on A (C and B are both collider points)

#6M2
N <- 1000
X <- rnorm(N)
Z <- rnorm(N,0.9*X)
Y <- rnorm(N,0.9*Z)
cor(X,Z)
cor(X,Y)
cor(Z,Y)
d <- data.frame(X=X, Y=Y, Z=Z)

model6M2 <- quap(
  alist(
    Y ~ dnorm(mu,sigma),
    mu <- a + bX*X + bZ*Z,
    a ~ dnorm(0,1),
    bX ~ dnorm(0,0.5),
    bZ ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
plot(precis(model6M2))

model6M2b <- quap(
  alist(
    Y ~ dnorm(mu,sigma),
    mu <- a + bX*X,
    a ~ dnorm(0,1),
    bX ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
plot(precis(model6M2b))
#This is different from the leg example in the book, as that was not a pipe
# but a fork

#6M3
#Upper left on Z (correct)
#Upper right on (A = wrong) nothing, as it is a collider situation
#Lower left on none (correct) as it is a collider situation
#Lower right on A (correct)

#6H1
library(dagitty)
library(rethinking)

dag6H1 <- dagitty("dag{
                  S -> W -> D
                  S -> A -> D
                  }")
coordinates(dag6H1) <- list(x=c(S=1, W=0, A=2, D=1),
                            y=c(S=0, W=1, A=1, D=2))
drawdag(dag6H1)

# Effect of W on D is influenced by confounder S

data(WaffleDivorce)
d <- WaffleDivorce
d$D <- standardize(d$Divorce)
d$S <- d$South + 1
d$W <- standardize(d$WaffleHouses)
d$A <- standardize(d$MedianAgeMarriage)

#Make a model for divorce rate vs waffle houses
model6H1a <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bW*W,
    a ~ dnorm(0,0.2),
    bW ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H1a)

model6H1b <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a[S] + bW*W,
    a[S] ~ dnorm(0,0.2),
    bW ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H1b,depth=2)

sim_W_change <- data.frame(W=seq(from=-3, to=3, length.out=30), S=1)
sim_W <- sim(model6H1b, data=sim_W_change)
plot(sim_W_change$W,colMeans(sim_W),type="l", ylim=c(-2,2))
shade(apply(sim_W,2,PI), sim_W_change$W)

#6H2
# S->A->D is a pipe, so conditioning on A & W would remove the effect of S
model6H2a <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a[S] + bA*A + bW*W,
    a[S] ~ dnorm(0,0.2),
    c(bA,bW) ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H2a,depth=2)

# Age is independent of south after correcting for waffle
model6H2b <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a[S] + bA*A,
    a[S] ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H2b,depth=2)

# Fox assignments

data(foxes)
d <- list()
d$W <- standardize(foxes$weight)
d$T <- standardize(foxes$area)
d$G <- standardize(foxes$groupsize)
d$F <- standardize(foxes$avgfood)

#6H3:
model6H3 <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bT*T,
    a ~ dnorm(0,0.2),
    bT ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

prior <- extract.prior(model6H3)
mu <- link(model6H3, post=prior, data=list(T=c(-2,2)))
plot(NULL, xlim=c(-2,2), ylim=c(-3,3), xlab="Territory area change (std)",
     ylab="Weight change (std)")
for(i in 1:50) lines(c(-2,2), mu[i,], col=col.alpha("black",0.4))

# Priors can go both ways, staying for the largest part within the expectation

T_seq <- seq(from=-3, to=3, length.out=30)
mu <- link(model6H3,data=list(T=T_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(W~T,data=d,col=rangi2)
lines(T_seq,mu.mean,lwd=2)
shade(mu.PI,T_seq)

precis(model6H3)

#Conclusions: no real effect of territory size on weight

#6H4:
# This is a pipe, so effect of area would be masked by including food
# By only including food, the total effect of food can be seen
# Also, there can be an indirect effect of group size, which we want to keep
model6H4 <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bF*F,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H4)

# Conclusion: no effect of average food available on weight through the 
# combination of both causal pathways

#6H5:

#First for a total effect of group size:
model6H5a <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bG*G,
    a ~ dnorm(0,0.2),
    bG ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H5a)

G_seq <- seq(from=-3, to=3, length.out=30)
mu <- link(model6H5a,data=list(G=G_seq))
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI,prob=0.91)
plot(G~W,data=d,col=rangi2)
lines(G_seq,mu.mean,lwd=2)
shade(mu.PI,G_seq)

# Conclusion small negative effect of group size

# Group size is confounded by food available, so use that as well

model6H5b <- quap(
  alist(
    W ~ dnorm(mu,sigma),
    mu <- a + bF*F + bG*G,
    a ~ dnorm(0,0.2),
    c(bF,bG) ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H5b)

sim_G_data <- data.frame(G=seq(from=-3, to=3, length.out=30), F=0)
sim_G <- sim(model6H5b,data=sim_G_data)
plot(sim_G_data$G,colMeans(sim_G),type="l",ylim=c(-2,2), xlab="simulated G",
     ylab="predicted weight")
shade(apply(sim_G,2,PI), sim_G_data$G)

# Counterfactual plot: at fixed F, large effect of G, with larger groups
# leading to smaller foxes

sim_F_data <- data.frame(F=seq(from=-3, to=3, length.out=30), G=0)
sim_F <- sim(model6H5b,data=sim_F_data)
plot(sim_F_data$F,colMeans(sim_F),type="l",ylim=c(-2,2), xlab="simulated F",
     ylab="predicted weight")
shade(apply(sim_F,2,PI), sim_F_data$F)

# Counterfactual plot: at fixed G, effect of F is strong as well with more
# food leading to heavier foxed

# This model does show an effect of both average food available & group size
cor(d$F,d$F)
# These two parameters are strongly correlated
# It thus seems that more food available leads to larger groups
# Weight is impacted negatively by group size and positively by food available

#6H6:
# The presence of mastitis (M) will impact proteolysis of caseins (P)
# Mastitis has three levels (health, subclinical, clinical; indexed 1, 2, 3)
# At any point in time on a farm, mastitis will be 0.8, 0.15, 0.05
# Mastitis will lead to higher somatic cell count (SCC)
# clinical mastitis will lead to more proteolysis & higher SSC
# Higher cell count will in general already lead to proteolysis
# Unobserved parameters: 
# days in lactation (Ud), with higher Ud leading to higher SSC & M
# b-casein ratio (Ub), with higher Ub leading to more proteolysis
library(rethinking)
library(dagitty)
dag6H6a <- dagitty("dag{
                  Ub [unobserved]
                  Ud [unobserved]
                  M -> P <- S
                  M <- Ud -> S <- M
                  Ub -> P
                  }")
coordinates(dag6H6) <- list(x=c(M=0, S=2, P=1, Ud=1, Ub=2),
                            y=c(M=1, S=1, P=2, Ud=0, Ub=2))
drawdag(dag6H6a)
adjustmentSets(dag6H6a, exposure="S", outcome="P")
impliedConditionalIndependencies(dag6H6a)

#What happens if we realize Ud is relevant and decide to measure it?
dag6H6b <- dagitty("dag{
                  Ub [unobserved]
                  M -> P <- S
                  M <- D -> S <- M
                  Ub -> P
                  }")
coordinates(dag6H6b) <- list(x=c(M=0, S=2, P=1, D=1, Ub=2),
                            y=c(M=1, S=1, P=2, D=0, Ub=2))
drawdag(dag6H6b)
adjustmentSets(dag6H6b, exposure="S", outcome="P")
impliedConditionalIndependencies(dag6H6b)

#What happens if we assume that the effect of days of lactation on cell count
# is only indirectly through mastitis risk?
dag6H6c <- dagitty("dag{
                  Ub [unobserved]
                  M -> P <- S
                  D -> S <- M
                  Ub -> P
                  }")
coordinates(dag6H6c) <- list(x=c(M=0, S=2, P=1, D=1, Ub=2),
                             y=c(M=1, S=1, P=2, D=0, Ub=2))
drawdag(dag6H6c)
adjustmentSets(dag6H6c, exposure="S", outcome="P")
impliedConditionalIndependencies(dag6H6c)

#6H7:
#Continuing making a simulated dataset based on dag6H6c
N <- 1e4
d <- as.data.frame(matrix(ncol=0,nrow=N))
d$D <- runif(N,min=0,max=1) # lactation stage as a fraction of total lactation
prob <- c(0.8,0.15,0.05) # mastitis probabilities
d$M <- sample(1:3,N,replace = T,prob=prob) # Mastitis probability turned to index
for(i in 1:nrow(d)){
  if(d$M[i] == 1){
    d$S[i] <- rnorm(1,50+80*d$D[i],20)
  } else if(d$M[i] == 2){
    d$S[i] <- rnorm(1,100+100*d$D[i],25)
  } else {
    d$S[i] <- rnorm(1,250+120*d$D[i],100)
  }
} # model for cell count based on mastitis status
dens(d$S)
for(i in 1:nrow(d)){
  if(d$M[i] == 1){
    d$P[i] <- rnorm(1,1.5+d$S[i]/100,0.5)
  } else if(d$M[i] == 2){
    d$P[i] <- rnorm(1,2+d$S[i]/200,0.6)
  } else {
    d$P[i] <- rnorm(1,2+d$S[i]/250,0.7)
  }
} # model for cell count based on mastitis status
dens(d$P)
pairs(d)
cor(d)
plot(d$S~d$P)
plot(d$S~d$P,col=d$M)

#Build a model to study the relation between SCC and proteolysis degree
# both for total effect as after adding mastitis as a confounder
# using the simulated data
model6H7a <- quap(
  alist(
    P ~ dnorm(mu,sigma),
    mu <- a + bs*S,
    a ~ dnorm(1,1),
    bs ~ dnorm(0.01,0.005),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H7a, digit=4)

model6H7b <- quap(
  alist(
    P ~ dnorm(mu,sigma),
    mu <- a[M] + bs*S,
    a[M] ~ dnorm(1,1),
    bs ~ dnorm(0.05,0.02),
    sigma ~ dexp(1)
  ), data=d
)
precis(model6H7b, depth=2, digit=4)
