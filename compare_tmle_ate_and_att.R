## Identify differences between tmle ATE and tmle ATT

rm(list=ls())
library(tmle)
library(tmlecte)

SL.library <- c("SL.glm","SL.step","SL.glm.interaction")

#data <- read.csv("C:/Users/kecolson/Google Drive/Intro to causal inference/Rlabs/RLab5.Fa2013.csv")

## Create data
set.seed(252)
expit <- function(x) { return(1/(1+exp(-x))) }
n <- 100000
UW1 <- runif(n,0,1); UW2 <- rnorm(n,1,2); UW3 <- runif(n,0,1); UW4 <- rnorm(n,2,1); UA <- runif(n,0,1); UY <- runif(n,0,1)
W1 <- as.numeric(UW1 <0.45)
W2 <- 0.75*UW2
W3 <- as.numeric(UW3 < expit(-1+1.3*W1+2.9*W2))
W4 <- as.numeric(UW4 <2.1)
A <- as.numeric(UA < expit(-1+2.6*W1+0.9*W2-2*W3) )
Y <- as.numeric(UY < expit(-2+A+0.5*A*W2+0.7*W1+W4))
Y0 <- as.numeric(UY < expit(-2+0+0.5*0*W2+0.7*W1+W4))
Y1 <- as.numeric(UY < expit(-2+1+0.5*1*W2+0.7*W1+W4))
cosW2<- cos(pi*W2)
sinW3<- sin(pi*W3)
W4sq<- W4*W4
data <- data.frame(W1,W2,W3,W4,A,Y,cosW2,sinW3,W4sq)
W <- subset(data, select=-c(A,Y))
X <- subset(data, select=c(W1,W2,W3,W4,A))
X1 <- X0 <- X
X1$A <- 1
X0$A <- 0
newdata <- rbind(X,X1,X0)

# Calculate true ATE/ATT
ATE <- mean(Y0-Y1)
ATT <- mean((Y0-Y1)[A==1])

## TMLE ATE by hand
  # 1. Create initial estimate of Q
Qinit <- SuperLearner(Y=data$Y, X=X, SL.library=SL.library,
                      family="binomial", newX=newdata)
  # 2: Estimate the treatment mechanism g_o(A|W) = Po(A|W)
gHatSL <- SuperLearner(Y=data$A, X=W, SL.library=SL.library,
                       family="binomial")
  # 3: create clever covariate
H.AW <- as.numeric(data$A==1)/gHatSL$SL.predict - as.numeric(data$A==0)/(1-gHatSL$SL.predict)
  # 4: update initial estimates
logitUpdate <- glm(data$Y ~ -1 + offset(qlogis(Qinit$SL.predict[1:n])) + H.AW, family="binomial")

tmle.ate.1 <- mean(plogis(qlogis(Qinit$SL.predict[(n+1):(2*n)]) + logitUpdate$coef*(1/gHatSL$SL.predict))) - mean(plogis(qlogis(Qinit$SL.predict[(2*n+1):(3*n)]) + logitUpdate$coef*(-1/(1-gHatSL$SL.predict))))

## TMLE ATE function
tmle.out<- tmle(Y=data$Y, A=data$A, W=W, Q.SL.library = SL.library, g.SL.library = SL.library, family="binomial")
tmle.ate.2 <- tmle.out$estimates$ATE$psi

## Confirm that TMLE ATE by hand and by function are the same
tmle.ate.1
tmle.ate.2
ATE

  # Close but not exactly the same-- same to ~3 decimal places. Estimate by hand is always slightly larger.
    # > ATE
    # [1] -0.29185
    # # Run 1
    # > tmle.ate.1
    # [1] 0.2921537
    # > tmle.ate.2
    # [1] 0.2920563
    # # Run 2
    # > tmle.ate.1
    # [1] 0.292137
    # > tmle.ate.2
    # [1] 0.2920485


## TMLE ATT by hand
# 4: update initial estimates
logitUpdate2 <- glm(data$Y[data$A==1] ~ -1 + offset(qlogis(Qinit$SL.predict[1:n][data$A==1])) + H.AW[data$A==1], family="binomial")

tmle.att.1 <- (mean(plogis(qlogis(Qinit$SL.predict[(n+1):(2*n)]) + logitUpdate$coef*(1/gHatSL$SL.predict))[data$A==1]) - mean(plogis(qlogis(Qinit$SL.predict[(2*n+1):(3*n)]) + logitUpdate$coef*(-1/(1-gHatSL$SL.predict)))[data$A==1]))

tmle.att.1b <- (mean(plogis(qlogis(Qinit$SL.predict[(n+1):(2*n)]) + logitUpdate2$coef*(1/gHatSL$SL.predict))[data$A==1]) - mean(plogis(qlogis(Qinit$SL.predict[(2*n+1):(3*n)]) + logitUpdate2$coef*(-1/(1-gHatSL$SL.predict)))[data$A==1]))

## TMLE ATT function
Y.scaled <- (data$Y-min(data$Y))/ (max(data$Y)-min(data$Y)) # Scale Y to be between 0 and 1 so it will play nice with tmle.att
tmle.out<- tmle.att(Y=Y.scaled, A=data$A, W=W, Q.SL.library = SL.library, g.SL.library = SL.library, family="binomial")
tmle.att.2<- tmle.out$psi * (max(data$Y)-min(data$Y)) # rescale at end

## Confirm that TMLE ATT by hand and by function are the same
tmle.att.1
tmle.att.1b
tmle.att.2
ATT
  # Close-- same to 2 decimal places
# > ATT
# [1] -0.3330542
# # Run 1
# > tmle.att.1
# [1] 0.3277985
# > tmle.att.1b
# [1] 0.3285735
# > tmle.att.2
# [1] 0.3263886

