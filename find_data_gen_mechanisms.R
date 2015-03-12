# Try to maintain:
# 40% treated
# Unadjusted estimate is substantially biased and is larger than true estimate
# ATT > ATE

rm(list=ls())
n <- 1000

data <- list(NULL)

# Make the 6 datasets
for (covset in 1:2) {  
  set.seed(123)
  if (covset == 1) {
    W1 <- runif(n, min = 0.02, max = 0.7)  
    W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1)
  } else if (covset == 2) {  
    W1 <- rgamma(n, 1, 2)
    W2 <- rbinom(n, 1, 0.45+0.04*W1)
  }
  for (txmode in 1:3) {
    if (txmode == 1) A <- A1 <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
    if (txmode == 2) A <- A2 <- rbinom(n, 1, plogis(-.4 + 1.8*W1 -.8*I(W1^2) + 1.7*I(W1*W2) - 1.1*W2))
    if (txmode == 3) A <- A3 <- rbinom(n, 1, plogis(-.5 + 1.4*W1 - .7*I(W1^2) + I(W1*W2) - .7*W2))    
    Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
    Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
    Y <-  ifelse(A==1,Y.1,Y.0)  
    data <- data.frame(W1,W2,A,Y.0,Y.1,Y, covset = covset, txmode = txmode)
    if (covset == 1 & txmode == 1) { datasets <- data
    } else { datasets <- rbind(datasets, data) }
  }
}

# Create summary
summ <- expand.grid(covset=1:2, txmode=1:3, unadj=NA, ATE=NA, ATT=NA, txunits =NA)

for (covset in 1:2) {
  for (txmode in 1:3) {
    temp <- datasets[datasets$covset == covset & datasets$txmode == txmode, ]
    
    summ$unadj[summ$covset == covset & summ$txmode == txmode] <- round(summary(lm(Y ~ A, data = temp))$coefficients[2,1],2)
    summ$ATE[summ$covset == covset & summ$txmode == txmode] <- round(mean(temp$Y.1-temp$Y.0),2)
    summ$ATT[summ$covset == covset & summ$txmode == txmode] <-  round(mean(temp$Y.1[temp$A==1] - temp$Y.0[temp$A==1]),2)
    summ$txunits[summ$covset == covset & summ$txmode == txmode] <- table(temp$A)[2]
  }
}

summ

temp1 <- datasets[datasets$covset == 1 & datasets$txmode == 1, ]
m1 <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = temp1)
temp2 <- datasets[datasets$covset == 2 & datasets$txmode == 1, ]
m2 <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = temp2)

p1 <- hist(predict(m1, type="response")[temp1$A==1], breaks=30)
p2 <- hist(predict(m1, type="response")[temp1$A==0], breaks=30)
plot( p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 

summ
summary(predict(m1, type="response")[temp1$A==1])
summary(predict(m1, type="response")[temp1$A==0])
summary(predict(m2, type="response")[temp2$A==1])
summary(predict(m2, type="response")[temp2$A==0])



