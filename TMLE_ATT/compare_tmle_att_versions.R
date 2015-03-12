# Compare two different versions of weighted TMLE ATT 

rm(list=ls())
set.seed(1.31)
options(scipen=5)

setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/code/TMLE_ATT")

SL.library <- c("SL.glm", "SL.glm.interaction", "SL.step", "SL.step.interaction")

suppressPackageStartupMessages(library("tmle"))
suppressPackageStartupMessages(library("SuperLearner"))
suppressPackageStartupMessages(library("tmlecte"))
suppressPackageStartupMessages(library("MatchIt"))

# Simulate data
generateData<- function(n) {
  
  W1 <- runif(n, min=0.02, max=0.7)  
  W2 <- rnorm(n, mean=(0.2+0.0125*W1), sd=1) 
  A <- rbinom(n, 1, plogis(-.5 + W1 + I(W1^2) + I(W1*W2) - W2))
  W1sq <- W1*W1
  W2sq <- W2*W2
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
  Y <-  ifelse(A==1,Y.1,Y.0)
  Y.star<- Y
  Y.star.1<- Y.1
  Y.star.0<- Y.0
  data.frame(W1,W2,W1sq, W2sq, A,Y, Y.star, Y.1, Y.0, Y.star.1, Y.star.0)
}

O.true<- generateData(n=1000000)
ATE <- mean(O.true$Y.star.1 - O.true$Y.star.0)
ATT <- mean(O.true$Y.star.1[O.true$A==1] - O.true$Y.star.0[O.true$A==1])

# Create a sample dataset to analyze
sample <- generateData(n=500)

# Conduct nearest neighbor matching
match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="nearest", replace = T))
O <- match.data(match)

# Analyze for ATT using package + hacky adjustment
O$Y.star.scaled <- (O$Y.star-min(O$Y.star))/ (max(O$Y.star)-min(O$Y.star)) # Scale Y.star to be between 0 and 1 so it will play nice with tmle.att
tmle.out<- tmle.att(Y=c(O$Y.star.scaled,NA), A=c(O$A,mean(O$A)), 
                    W=cbind(c(O$W1,mean(O$W1)), c(O$W2,mean(O$W2)), c(O$W1sq,mean(O$W1sq)), c(O$W2sq,mean(O$W2sq))), 
                    family = "gaussian", 
                    Delta = c(rep(1,nrow(O)),0), 
                    gDelta.method = "user",
                    gDelta.1 = c(1/O$weights,NA),
                    g.method = "SL",
                    g.SL.library = SL.library, 
                    Q.method = "SL",
                    Q.SL.library = SL.library)
att1 <- (sum(O$weights * (tmle.out$Qstar[1:nrow(O),2] - tmle.out$Qstar[1:nrow(O),1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))

# Run it again just to see how much it varies stochastically
tmle.out<- tmle.att(Y=c(O$Y.star.scaled,NA), A=c(O$A,mean(O$A)), 
                    W=cbind(c(O$W1,mean(O$W1)), c(O$W2,mean(O$W2)), c(O$W1sq,mean(O$W1sq)), c(O$W2sq,mean(O$W2sq))), 
                    family = "gaussian", 
                    Delta = c(rep(1,nrow(O)),0), 
                    gDelta.method = "user",
                    gDelta.1 = c(1/O$weights,NA),
                    g.method = "SL",
                    g.SL.library = SL.library, 
                    Q.method = "SL",
                    Q.SL.library = SL.library)
att1b <- (sum(O$weights * (tmle.out$Qstar[1:nrow(O),2] - tmle.out$Qstar[1:nrow(O),1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))

# Analyze for ATT using adjusted source code from ATT package to prevent gDelta.1 from resetting to 1 for all if Y doesn't have any missingness
source("TMLE_ATT_edited.R")
tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A, 
                    W=cbind(O$W1, O$W2, O$W1sq, O$W2sq), 
                    family = "gaussian", 
                    Delta = rep(1,nrow(O)), 
                    gDelta.method = "user",
                    gDelta.1 = 1/O$weights,
                    g.method = "SL",
                    g.SL.library = SL.library, 
                    Q.method = "SL",
                    Q.SL.library = SL.library)
att2 <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))

# Run it again just to see how much it varies stochastically
tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A, 
                     W=cbind(O$W1, O$W2, O$W1sq, O$W2sq), 
                     family = "gaussian", 
                     Delta = rep(1,nrow(O)), 
                     gDelta.method = "user",
                     gDelta.1 = 1/O$weights,
                     g.method = "SL",
                     g.SL.library = SL.library, 
                     Q.method = "SL",
                     Q.SL.library = SL.library)
att2b <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))

# Analyze using package and ignoring weights
tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, 
                    W=cbind(O$W1, O$W2, O$W1sq, O$W2sq), 
                    family = "gaussian", 
                    g.method = "SL",
                    g.SL.library = SL.library, 
                    Q.method = "SL",
                    Q.SL.library = SL.library)
att3 <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star))

# Analyze using package with weights, but without hacky method-- this should be same as #3
tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, 
                    W=cbind(O$W1, O$W2, O$W1sq, O$W2sq), 
                    family = "gaussian", 
                    Delta = rep(1,nrow(O)), 
                    gDelta.method = "user",
                    gDelta.1 = 1/O$weights,
                    g.method = "SL",
                    g.SL.library = SL.library, 
                    Q.method = "SL",
                    Q.SL.library = SL.library)
att4 <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star))

# diff
ATT
att1
att1b
att2
att2b
att3
att4

# Run 1:
# [1] 1.701813
# [2] 1.692809
# [3] 1.700019
# [4] 1.697379

# Run 2: new sample, new match
# [1] 1.460922
# [2] 1.462239
# [3] 1.498367
# [4] 1.498256

# Run 3: another new sample, new match
# [1] 1.778959
# [1b] 1.776913
# [2] 1.754447
# [2b] 1.752553 varies stochastically a wee bit but not quite as much as the difference between 1 and 2
# [3] 1.811494
# [4] 1.811529

# perc diff -- each about ~ 0.5% different...
(att1-att2)/att2*100
(att3-att2)/att2*100
(att4-att2)/att2*100

# Conclusion:
# All methods are biased vs real ATT. This is probably a function of the sampling and matching and inability to correct for all bias
# Each of the methods is very similar. 
# The weighted and unweighted are different enough that it merits properly weighting
# The hacky and properly implemented are very similar in this scenario but they might not be similar in every scenario
# Safer to use the properly implemented one, even if it's a bit more annoying to code.










