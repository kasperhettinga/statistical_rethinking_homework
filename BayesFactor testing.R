#https://richarddmorey.github.io/BayesFactor/
library(BayesFactor)

#One-sample t-test BF alternative
#Based on https://bayesfactor.blogspot.com/2014/02/bayes-factor-t-tests-part-1.html
improvements = with(sleep, extra[group == 2])
hist(improvements)
N = length(improvements)
t = mean(improvements)/sd(improvements) * sqrt(N)
t
deltaHat = t/sqrt(N)
deltaHat

#Shorter is to just calculate the normalized improvement directly
deltaHat <- mean(improvements)/sd(improvements)

dt(t, df = 9, ncp = 1 * sqrt(10))/dt(t, df = 9)

ttestBF(improvements, nullInterval = c(0, Inf))

#Two-sample t-test BF alternative
#Based on: https://bayesfactor.blogspot.com/2014/02/bayes-factor-t-tests-part-2-two-sample.html
randDotStereo <- read.csv("randomFusion.csv", header = TRUE, sep = ";")
boxplot(randDotStereo$logFuseTime,randDotStereo$condition)

#regular t-test
classical.test = t.test(logFuseTime ~ condition, data = randDotStereo, var.eq = TRUE)
classical.test

#two-sided BF test
bf = ttestBF(formula = logFuseTime ~ condition, data = randDotStereo)
bf

#one-sided BF test
bf.signed = ttestBF(formula = logFuseTime ~ condition, data = randDotStereo, 
                    nullInterval = c(-Inf, 0))
bf.signed

bf.pos.vs.neg = bf.signed[1]/bf.signed[2]
bf.pos.vs.neg

bf.negl = ttestBF(formula = logFuseTime ~ condition, data = randDotStereo, nullInterval = c(0, 
                                                                                            0.2))
bf.negl
bf.nonnegl = ttestBF(formula = logFuseTime ~ condition, data = randDotStereo, 
                     nullInterval = c(0.2, Inf))
bf.nonnegl
bf.nonnegl[2]/bf.negl[2]
