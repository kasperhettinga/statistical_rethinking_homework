library(rethinking)
data(rugged)
d <- rugged

d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)

model_rugged_gdp1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a+b*(rugged_std-0.215),
    a ~ dnorm(1,0.1),
    b ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd
)

#Check prior predictions
set.seed(1234)
prior <- extract.prior(model_rugged_gdp1)
plot(NULL, xlim=c(0,1), ylim=c(0.5,1.5), xlab="ruggedness", ylab="log GDP")
abline(h=min(dd$log_gdp_std), lty=2)
abline(h=max(dd$log_gdp_std), lty=2)
rugged_seq <- seq(from=-0.1, to=1.1, length.out=30)
mu <- link(model_rugged_gdp1, post=prior, data=data.frame(rugged_std=rugged_seq))
for(i in 1:50) lines(rugged_seq,mu[i,],col=col.alpha("black",0.3))

dd$cid <- ifelse(dd$cont_africa==1,1,2)

model_rugged_gdp2 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd
)
precis(model_rugged_gdp2,depth=2)
compare(model_rugged_gdp1,model_rugged_gdp2)

post <- extract.samples(model_rugged_gdp2)
diff_1_2 <- post$a[,1]-post$a[,2]
PI(diff_1_2)

mu_Afr <- link(model_rugged_gdp2, data=data.frame(cid=1,rugged_std=rugged_seq))
mu_NotAfr <- link(model_rugged_gdp2, data=data.frame(cid=2,rugged_std=rugged_seq))
mu.Afr_mu <- apply(mu_Afr,2,mean)
mu.Afr_ci <- apply(mu_Afr,2,PI,prob=0.97)
mu.NotAfr_mu <- apply(mu_NotAfr,2,mean)
mu.NotAfr_ci <- apply(mu_NotAfr,2,PI,prob=0.97)

plot(dd$rugged_std,dd$log_gdp_std,xlab="ruggedness",ylab="log GDP")
lines(rugged_seq,mu.Afr_mu,col="red")
lines(rugged_seq,mu.NotAfr_mu,col="blue")
shade(mu.Afr_ci, rugged_seq, col=col.alpha("red",0.2))
shade(mu.NotAfr_ci, rugged_seq, col=col.alpha("blue",0.2))


model_rugged_gdp3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data=dd
)
precis(model_rugged_gdp3,depth=2)

compare(model_rugged_gdp1,model_rugged_gdp2,model_rugged_gdp3,func=PSIS)
plot(PSIS(model_rugged_gdp3,pointwise = TRUE)$k)

d.A1 <- dd[dd$cid==1,]
plot(d.A1$rugged_std, d.A1$log_gdp_std, pch=16, col="blue", 
     xlab="ruggedness std", ylab="log GDP std", xlim=c(0,1))
mu <- link(model_rugged_gdp3, data=data.frame(cid=1,rugged_std=rugged_seq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.97)
lines(rugged_seq,mu_mean)
shade(mu_PI,rugged_seq,col=col.alpha("blue",0.2))

d.A2 <- dd[dd$cid==2,]
plot(d.A2$rugged_std, d.A2$log_gdp_std, pch=16, col="red", 
     xlab="ruggedness std", ylab="log GDP std", xlim=c(0,1))
mu <- link(model_rugged_gdp3, data=data.frame(cid=2,rugged_std=rugged_seq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.97)
lines(rugged_seq,mu_mean)
shade(mu_PI,rugged_seq,col=col.alpha("red",0.2))

#8.2 Symmetry of interactions
rugged_seq <- seq(from=-0.2, to=1.2, length.out=30)
muA <- link(model_rugged_gdp3, data=data.frame(cid=1,rugged_std=rugged_seq))
muN <- link(model_rugged_gdp3, data=data.frame(cid=2,rugged_std=rugged_seq))
delta <- muA-muN
delta_mean <- apply(delta,2,mean)
delta_PI <- apply(delta,2,PI,prob=0.97)

plot(NULL, xlab="ruggedness std", ylab="log GDP std", xlim=c(0,1), ylim=c(-0.3,0.3))
lines(rugged_seq,delta_mean)
shade(delta_PI,rugged_seq,col=col.alpha("red",0.2))
abline(h=0, lty=2)

#8.3 Continuous interactions
data("tulips")
d <- tulips
d$blooms_std <- d$blooms/max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

model_tulips1 <- quap(
  alist(
    blooms_std ~ dnorm(mu,sigma),
    mu <- a + bW*water_cent + bS*shade_cent,
    a ~ dnorm(0.5,0.25),
    bW ~ dnorm(0,0.25),
    bS ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data=d
)

model_tulips2 <- quap(
  alist(
    blooms_std ~ dnorm(mu,sigma),
    mu <- a + bW*water_cent + bS*shade_cent + bWS*water_cent*shade_cent,
    a ~ dnorm(0.5,0.25),
    bW ~ dnorm(0,0.25),
    bS ~ dnorm(0,0.25),
    bWS ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data=d
)

precis(model_tulips2,depth=2)
compare(model_tulips1,model_tulips2,func=WAIC)

#Plot posterior predictions for tulip models
par(mfrow=c(1,3))
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx],d$blooms_std[idx],xlim=c(-1,1),ylim=c(0,1),
       xlab="water",ylab="blooms")
  mu <- link(model_tulips1,data=data.frame(shade_cent=s, water_cent=-1:1))
  for(i in 1:20) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx],d$blooms_std[idx],xlim=c(-1,1),ylim=c(0,1),
       xlab="water",ylab="blooms")
  mu <- link(model_tulips2,data=data.frame(shade_cent=s, water_cent=-1:1))
  for(i in 1:20) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}
for(w in -1:1){
  idx <- which(d$water_cent==s)
  plot(d$shade_cent[idx],d$blooms_std[idx],xlim=c(-1,1),ylim=c(0,1),
       xlab="shade",ylab="blooms")
  mu <- link(model_tulips2,data=data.frame(shade_cent=-1:1, water_cent=w))
  for(i in 1:20) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}

#Plotting priors
set.seed(1234)
priors <- extract.prior(model_tulips2)
for(w in -1:1){
  plot(NULL,xlim=c(-1,1),ylim=c(-1,2),xlab="shade",ylab="blooms")
  mu <- link(model_tulips2,post=priors,data=data.frame(shade_cent=-1:1, water_cent=w))
  for(i in 1:50) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}

#Practice
#8E1
# 1) temperature, 2) socio-economic status, 3) oxygen, gear

#8E2
# 1) heat * drying interaction 2) cyl & injector_type are independent, 
# 3) Two separate routes, no interaction, 4) no interaction, two are independent?

#8E3
# 1: mu <- a + bH*H + bD*D = bHD*H*D
# 3: mu <- a + bFam*Fam + bFr*Fr
# etc

#8M1
# The model with all the interactions will result in a mu of 0 if temperature = hot

#8M2
# mu <- a[temp] + ...
# In this way, when temp is hot, all parameters are zero
# Or just make an index variable that is 0 for hot and then multiply the formula
# by this index variable

#8M3
# Num_raven ~ dnorm(mu, sigma)
# mu <- a + bP*P + bW*W + bPW*P*W // with P being amount of prey and W being wolves

# Example based on: https://www.erikkusch.com/courses/rethinking/chapter-08/
N <- 1e5 # simulation size
rPW <- 0.2 # correlation between prey and wolf
bP <- 0.2 # regression coefficient for prey
bW <- -0.4 # regression coefficient for wolf
bPW <- 0.3 # regression coefficient for prey-by-wolf interaction
# Simulate data
prey <- as.integer(rnorm(N, mean = 100, sd = 15)) # as.integer, so we have "whole" animals
wolves <- as.integer(rnorm(N, mean = 10 + rPW * prey, sd = 7))
ravens <- as.integer(rnorm(N, mean = 5 + bP * prey + bW * wolves + bPW * wolves * prey, sd = 9))
d <- data.frame(prey = prey, wolves = wolves, ravens = ravens)
# plot the data
par(mfrow = c(1, 2))
plot(ravens ~ prey, data = d, main = "Ravens like prey!")
plot(ravens ~ wolves, data = d, main = "Ravens like wolves?")
#Visually correlated, at different correlation-values

#8M4
data("tulips")
d <- tulips
d$blooms_std <- d$blooms/max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

model_tulips8M4 <- quap(
  alist(
    blooms_std ~ dnorm(mu,sigma),
    mu <- a + bW*water_cent + bS*shade_cent + bWS*water_cent*shade_cent,
    a ~ dnorm(0.5,0.25),
    bW ~ dlnorm(0,0.25),
    bS ~ dlnorm(0,0.25),
    bWS ~ dlnorm(0,0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_tulips8M4)

#Plotting priors
set.seed(1234)
priors <- extract.prior(model_tulips8M4)
for(w in -1:1){
  plot(NULL,xlim=c(-1,1),ylim=c(-1,2),xlab="shade",ylab="blooms")
  mu <- link(model_tulips8M4,post=priors,data=data.frame(shade_cent=-1:1, water_cent=w))
  for(i in 1:50) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}

#8H1
data("tulips")
d <- tulips
d$blooms_std <- d$blooms/max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)
d$bed_id <- coerce_index(d$bed)

model_tulips8H1 <- quap(
  alist(
    blooms_std ~ dnorm(mu,sigma),
    mu <- a[bed_id] + bW*water_cent + bS*shade_cent + bWS*water_cent*shade_cent,
    a[bed_id] ~ dnorm(0.5,0.25),
    bW ~ dnorm(0,0.25),
    bS ~ dnorm(0,0.25),
    bWS ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_tulips8H1, depth=2)

post <- extract.samples(model_tulips8H1)
diff_a_b <- post$a[,1] - post$a[,2]
HPDI(diff_a_b)
diff_b_c <- post$a[,2] - post$a[,3]
HPDI(diff_b_c)

#8H2
compare(model_tulips2,model_tulips8H1,func=WAIC)

# Both models perform similarly. Bed thus has an effect of bloom (bed a being 
# worse), but knowing bed does not improve prediction of blooms

#Plot posterior predictions for tulip models
par(mfrow=c(1,3))
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx],d$blooms_std[idx],xlim=c(-1,1),ylim=c(0,1),
       xlab="water",ylab="blooms")
  mu <- link(model_tulips8H1,data=data.frame(shade_cent=s, water_cent=-1:1, bed_id=1))
  for(i in 1:20) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}
for(w in -1:1){
  idx <- which(d$water_cent==s)
  plot(d$shade_cent[idx],d$blooms_std[idx],xlim=c(-1,1),ylim=c(0,1),
       xlab="shade",ylab="blooms")
  mu <- link(model_tulips8H1,data=data.frame(shade_cent=-1:1, water_cent=w, bed_id=1))
  for(i in 1:20) lines (-1:1,mu[i,],col=col.alpha("blue",0.3))
}

# Plot predicted bloom size per bed
post <- extract.samples(model_tulips8H1)
post.a <- post$a[, 1]
post.b <- post$a[, 2]
post.c <- post$a[, 3]
dens(post.a, col = "red", xlim = c(0, 1), ylim = c(0, 15))
dens(post.b, col = "blue", add = TRUE)
dens(post.c, col = "black", add = TRUE)

#Bed a is clearly very different, whereas beds b&c are very similar. All have 
# lots of variation

#8H3
data(rugged)
d <- rugged

d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa==1,1,2)

model_rugged_8H3a <- quap(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.5),
    b[cid] ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dd
)
precis(model_rugged_8H3a, depth=2)
set.seed(1234)
PSIS_8H3a <- PSIS(model_rugged_8H3a,pointwise=TRUE)
set.seed(1234)
WAIC_8H3a <- WAIC(model_rugged_8H3a,pointwise=TRUE)
plot(PSIS_8H3a$k,WAIC_8H3a$penalty,xlab="PSIS k",ylab="WAIC penalty")
as.character(dd[WAIC_8H3a$penalty > .2, ]$country)
as.character(dd[PSIS_8H3a$k > .35, ]$country)
# All countries with either a high WAIC penalty or a high Pareto k-score
# are rich countries that are very rugged

model_rugged_8H3b <- quap(
  alist(
    log_gdp_std ~ dstudent(2,mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.5),
    b[cid] ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=dd
)
precis(model_rugged_8H3b, depth=2)
# Predictions for the robust model for alpha and beta values are very similar

set.seed(1234)
PSIS_8H3b <- PSIS(model_rugged_8H3b,pointwise=TRUE)
set.seed(1234)
WAIC_8H3b <- WAIC(model_rugged_8H3b,pointwise=TRUE)
plot(PSIS_8H3b$k,WAIC_8H3b$penalty,xlab="PSIS k",ylab="WAIC penalty")
as.character(dd[WAIC_8H3b$penalty > .2, ]$country)
as.character(dd[PSIS_8H3b$k > .2, ]$country)
# Seychelles is no longer such an extreme according to PSIS. In general
# penalty values have become smaller, showing the more robust prediction
# with PSIS no longer giving warnings.

#Now without Seychelles included
d2 <- dd[ dd$country!="Seychelles" , ]

model_rugged_8H3c <- quap(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.5),
    b[cid] ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d2
)
model_rugged_8H3d <- quap(
  alist(
    log_gdp_std ~ dstudent(2,mu,sigma),
    mu <- a[cid]+b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.5),
    b[cid] ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d2
)

compare(model_rugged_8H3a,model_rugged_8H3b,func=PSIS)

compare(model_rugged_8H3c,model_rugged_8H3d,func=PSIS)
#The more robust model based on the student distribution, gets a worse PSIS 
# score. Removing the Seychelles actually hardly changed the model score

precis(model_rugged_8H3a, depth=2)
precis(model_rugged_8H3c, depth=2)
# Slope in Africa has halved with the model without the Seychelles

#8H4
library(rethinking)
data("nettle")
d <- nettle
d$lang.per.cap <- d$num.lang/d$k.pop
d$log.lang.per.cap <- log(d$lang.per.cap)

plot(d$log.lang.per.cap~d$mean.growing.season)
#log.lang.per.cap averages around -6 with a variation of around +/- 2, so the
# prior SD will be set to 1 so that any values between -4 and -8 is within 2 SD
# There are 12 units of growing season over 6 units of log.lang.per.cap, so the 
# max slope to expect is around 0.5.
#I will build a naive model, with an expected slope of 0. With an SD of 0.25, 
# slopes up til 0.5 fall within 2 standard deviation units

model_lang_cap8H4a <- quap(
  alist(
    log.lang.per.cap ~ dnorm(mu,sigma),
    mu <- a + bM*mean.growing.season,
    a ~ dnorm(-6,1),
    bM ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_lang_cap8H4a)
#The model shows a weak but positive relation between mean growing season and 
# languages. 

plot(d$log.lang.per.cap~d$area)
#Are has a huge variation and will be log-scaled first
d$log.area <- log(d$area)
plot(d$log.lang.per.cap~d$log.area)
#There are 6 units of log.area with 6 units of log.lang.per.cap, so max 
# expected slope would be one. Chosen SD will therefore be 0.5.

model_lang_cap8H4b <- quap(
  alist(
    log.lang.per.cap ~ dnorm(mu,sigma),
    mu <- a + bM*mean.growing.season + bA*log.area,
    a ~ dnorm(-6,1),
    bM ~ dnorm(0,0.25),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_lang_cap8H4b)
#The model shows again a weak but positive relation between mean growing season
# and languages. Effect of area seems to be a small negative one.

plot(d$log.lang.per.cap~d$sd.growing.season)
#The SD of the growing season seems to be usable as is, with 6 units of 
# variation over 6 units of log.lang.per.cap.

model_lang_cap8H4c <- quap(
  alist(
    log.lang.per.cap ~ dnorm(mu,sigma),
    mu <- a + bS*sd.growing.season,
    a ~ dnorm(-6,1),
    bS ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_lang_cap8H4c)
#The model shows a negative relation between SD and log.lang, which is larger
# than the effects of either growing season or area.

model_lang_cap8H4d <- quap(
  alist(
    log.lang.per.cap ~ dnorm(mu,sigma),
    mu <- a + bS*sd.growing.season + bA*log.area,
    a ~ dnorm(-6,2),
    c(bS,bA) ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_lang_cap8H4d)
# Adding log area gives again a negative relation for both area and SD.

# As it is hard to predict the interaction factor, a wider prior is selected
model_lang_cap8H4e <- quap(
  alist(
    log.lang.per.cap ~ dnorm(mu,sigma),
    mu <- a + bM*mean.growing.season + bS*sd.growing.season + 
      bMS*mean.growing.season*sd.growing.season,
    a ~ dnorm(-6,1),
    bM ~ dnorm(0,0.25),
    bS ~ dnorm(0,0.5),
    bMS ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_lang_cap8H4e)
#Once an interaction term is added, SD becomes a factor with a positive relation
#Let's compare models and plot outcomes.

compare(model_lang_cap8H4a,model_lang_cap8H4b,model_lang_cap8H4c,
        model_lang_cap8H4d,model_lang_cap8H4e,func=WAIC)
#The models with SD and without interaction score worst, whereas the model
# with interactions scores best (although the dWAIC is in the same order as
# dSE)

#Plot posterior predictions for model e
par(mfrow=c(1,3))
for(s in c(0,2,4)){
  idx <- which(d$sd.growing.season==s)
  plot(d$mean.growing.season[idx],d$log.lang.per.cap[idx],xlim=c(0,12),ylim=c(-10,0),
       xlab="season length",ylab="log languages")
  mu <- link(model_lang_cap8H4e,data=data.frame(sd.growing.season=s, 
                                                mean.growing.season=0:12))
  for(i in 1:20) lines (0:12,mu[i,],col=col.alpha("blue",0.3))
}
# Only at low SD, there is a clear relation between season length and languages
# At high SD, the relation tends towards being negative

for(m in c(0,6,12)){
  idx <- which(d$mean.growing.season==s)
  plot(d$sd.growing.season[idx],d$log.lang.per.cap[idx],xlim=c(0,6),ylim=c(-10,0),
       xlab="variation",ylab="log languages")
  mu <- link(model_lang_cap8H4e,data=data.frame(mean.growing.season=m, 
                                                sd.growing.season=0:6))
  for(i in 1:20) lines (0:6,mu[i,],col=col.alpha("blue",0.3))
}
#At a long growing season, the model predicts a negative relation between SD
# and log languages, whereas at shorter seasons this effects becomes smaller
# or even disappears.

#8H5
library(rethinking)
data("Wines2012")
d <- Wines2012
d$score.std <- d$score/max(d$score) # to keep the 0 at 0, I divided by max
d$wine.index <- coerce_index(d$wine)
d$judge.index <- coerce_index(d$judge)
str(d)
plot(d$judge.index~d$score.std)
plot(d$wine.index~d$score.std)

model_wine8H5 <- quap(
  alist(
    score.std ~ dnorm(mu,sigma),
    mu <- w[wine.index] + j[judge.index],
    w[wine.index] ~ dnorm(0.7,0.25), # average std score is around 0.7
    j[judge.index] ~ dnorm(0.7,0.25), # and spread is unknown
    sigma ~ dexp(1)
  ), data=d
)
precis(model_wine8H5,depth=2)

# Wines 6 & 18 score rather low compared to the other wines
# judges 5&6 score high and judge 4 and 9 low on average

#8H6
library(rethinking)
data("Wines2012")
d <- Wines2012
d$score.std <- d$score/max(d$score) # to keep the 0 at 0, I divided by max
d$flight.index <- coerce_index(d$flight)
d$wine_am.index <- coerce_index(d$wine.amer) + 1
d$judge_am.index <- coerce_index(d$judge.amer) + 1
str(d)

model_wine8H6 <- quap(
  alist(
    score.std ~ dnorm(mu,sigma),
    mu <- f[flight.index] + w[wine_am.index] + j[judge_am.index],
    f[flight.index] ~ dnorm(0.5,0.25), # average std score is around 0.7 but
    w[wine_am.index] ~ dnorm(0.5,0.25), # effect of individual index will be
    j[judge_am.index] ~ dnorm(0.5,0.25), # lower & spread is unknown
    sigma ~ dexp(1)
  ), data=d
)
precis(model_wine8H6,depth=2)
# The color of the wine does not seem to have any influence on the score
# American wines score on average slightly lower
# American judges score on average slightly higher

#8H7
model_wine8H7 <- quap(
  alist(
    score.std ~ dnorm(mu,sigma),
    mu <- f[flight.index] + w[wine_am.index] + j[judge_am.index] +
      fw*flight.index*wine_am.index + fj*flight.index*judge_am.index + 
      wj*wine_am.index*judge_am.index,
    f[flight.index] ~ dnorm(0.5,0.25), # average std score is around 0.7 but
    w[wine_am.index] ~ dnorm(0.5,0.25), # effect of individual index will be
    j[judge_am.index] ~ dnorm(0.5,0.25), # lower & spread is unknown
    c(fw,fj,wj) ~ dnorm(0.5,0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(model_wine8H7,depth=2)
# Hard to interpret. What is the predicted posterior for different combis?

wine_predict <- data.frame(wine_am.index = rep(c(1, 2), 4), 
                      judge_am.index = rep(c(1, 2), each = 4), 
                      flight.index = c(1, 1, 2, 2, 1, 1, 2, 2))
mu <- link(model_wine8H7,wine_predict)
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI,prob=0.97)
plot(mu_mean)
# Highest score is for group 5, American judge drinking French red wine
# Lowest score is for group 2, French judge drinking American red wine