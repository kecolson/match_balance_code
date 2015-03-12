# Compare regions of propensity score support for each matching method
# Ellie Colson
# 12/16/14


rm(list=ls())
library("MatchIt")
library("SuperLearner")
library("Matching")
library("tmle")
library("Hmisc")
library("sm")
options(scipen=5)

# Settings:
# Normal or non-normal covariates?
covars.norm <- T
# Propensity score overlap (1=good, 2= medium, 3=poor)
pscore.type <- 3

#--------------
# generateData:
#   input: n
# 	output: full data
#------------------
generateData<- function(n, covars.norm=T, pscore.type = 1) {
  
  if (covars.norm == T) {
    W1 <- runif(n, min = 0.02, max = 0.7)  
    #W2 <- rnorm(n, mean = (0.2 + 0.0125*W1), sd = 1) 
    W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1) 
  } else {
    W1 <- rgamma(n, 2, 10)
    W2 <- rbinom(n, 1, 0.55+0.05*W1)
  }
  
  if (pscore.type==1) {
    # A <- A1 <- rbinom(n, 1, plogis(-1.3 + W1 + I(W1^2) + I(W1*W2) - W2))
    # A <- rbinom(n, 1, plogis(-.5 + W1 + .1*I(W1^2) + .5*I(W1*W2) - W2))
    A <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
  } else if (pscore.type==2) {
    A <- A2 <- rbinom(n, 1, plogis(-.7 + 1.8*W1 -.1*I(W1^2) + 1.7*I(W1*W2) - 1.4*W2))
  } else if (pscore.type==3) {
    A <- A3 <- rbinom(n, 1, plogis(-.3 + 2*W1 -2*I(W1^2) + 2*I(W1*W2) - 2.5*W2))
  }
  
  ## Generate outcomes
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

# get the truth
set.seed(123)
O.true<- generateData(n=1000000, covars.norm = covars.norm, pscore.type = pscore.type)

ATE <- mean(O.true$Y.star.1 - O.true$Y.star.0)
ATT <- mean(O.true$Y.star.1[O.true$A==1] - O.true$Y.star.0[O.true$A==1])
ATE*1000
ATT*1000

m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = O.true)
summary(predict(m, type="response"))


#######
# plot.support:
# Plot propensity score overlap
# Input: data, title
# Output: graph
#######
plot.support <- function(data,title,g.pred,weights=NA) {

  if (!is.na(weights[1])) {
    plot(density(g.pred[data$A==1], adjust = 0.4, weights=weights[data$A==1]/sum(weights[data$A==1])), col=2, main=title, xlab="Propensity score", xlim=c(0,1))
    lines(density(g.pred[data$A==0], adjust = 0.4, weights=weights[data$A==0]/sum(weights[data$A==0])), col=3, lty=2)
  } else {
    plot(density(g.pred[data$A==1], adjust = 0.4), col=2, main=title, xlab="Propensity score", xlim=c(0,1))
    lines(density(g.pred[data$A==0], adjust = 0.4), col=3, lty=2)
  }
  legend("topright",legend=c("Treated units", "Controls"), cex=0.7, pch=c(15,3))
}



##########################
### Match and plot support
##########################

n=100000
ss <- 1000

# generate the population
set.seed(123)
pop <- generateData(n = n, covars.norm = covars.norm, pscore.type = pscore.type)

# Take a random (representative) sample from the population
set.seed(1235)
sample <- pop[sample(row.names(pop), size = ss, replace = F),]

# Estimate propensity score
g.model<- glm(A ~ poly(W1,2) + W2 + W1:W2, family='binomial', data = sample)
g.pred<- predict(g.model, type = 'resp')

pdf("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/pscore_support.pdf", height=7, width=5)

##############
# 1. SAMPLE AS IS
############

plot.support(sample,"Unmatched sample", g.pred)

#############
# 2. Nearest neighbor matching
#############
match <- Obs <- NULL
match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="nearest"))
Obs <- match.data(match)
plot.support(Obs,"Nearest neighbor", g.pred[row.names(sample) %in% row.names(Obs)])

# Obs dropped from sample: 
text(0.8, 3, paste("Treated units dropped:",(sum(sample$A)-sum(Obs$A))), cex=0.8)
text(0.8, 2.8, paste("Control units dropped:",((nrow(sample)-sum(sample$A))-(nrow(Obs)-sum(Obs$A)))), cex=0.8)

#############
# 3. Optimal matching 
#############
match <- Obs <- NULL
match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="optimal"))
Obs <- match.data(match)
plot.support(Obs,"Optimal matching", g.pred[row.names(sample) %in% row.names(Obs)])

# Obs dropped from sample: 
text(0.8, 3, paste("Treated units dropped:",(sum(sample$A)-sum(Obs$A))), cex=0.8)
text(0.8, 2.8, paste("Control units dropped:",((nrow(sample)-sum(sample$A))-(nrow(Obs)-sum(Obs$A)))), cex=0.8)

#############
# 4. Estimate pscore with superlearner, then nearest neighbor match
#############
SL.out <- match <- m <- Obs <- NULL
SL.out <- SuperLearner(Y=sample$A, X=sample[,c("W1","W2")], SL.library=c("SL.glm", "SL.glm.interaction","SL.step", "SL.step.interaction"), family='binomial', cvControl=list(V=10))
match <- suppressWarnings(Match(Y=sample$Y.star, Tr=sample$A, X=SL.out$SL.predict, estimand ="ATE", replace = F, ties = F, M=1, version="fast"))
m <- MatchBalance(A~W1+W2, data=sample, match.out = match, nboots=10)
Obs <- rbind(sample[unique(match$index.treated),] , sample[unique(match$index.control),])
plot.support(Obs,"SuperLearner + nearest neighbor", g.pred[row.names(sample) %in% row.names(Obs)])

# Obs dropped from sample: 
text(0.8, 3, paste("Treated units dropped:",(sum(sample$A)-sum(Obs$A))), cex=0.8)
text(0.8, 2.8, paste("Control units dropped:",((nrow(sample)-sum(sample$A))-(nrow(Obs)-sum(Obs$A)))), cex=0.8)

#############
# 5. genetic                           
#############
match <- Obs <- NULL
match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="genetic"))
Obs <- match.data(match)
plot.support(Obs,"Genetic matching", g.pred[row.names(sample) %in% row.names(Obs)], weights=Obs$weights) 
#plot.support(Obs,"Genetic matching")

# Obs dropped from sample: 
text(0.8, 3, paste("Treated units dropped:",(sum(sample$A)-sum(Obs$A))), cex=0.8)
text(0.8, 2.8, paste("Control units dropped:",((nrow(sample)-sum(sample$A))-(nrow(Obs)-sum(Obs$A)))), cex=0.8)

#############
# 6. subclassification
#############
match <- Obs <- NULL
match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="subclass", subclass=12))
Obs <- match.data(match)
plot.support(Obs,"Subclassification", g.pred[row.names(sample) %in% row.names(Obs)], weights=Obs$weights) 
#plot.support(Obs,"Subclassification")

# Obs dropped from sample: 
text(0.8, 3, paste("Treated units dropped:",(sum(sample$A)-sum(Obs$A))), cex=0.8)
text(0.8, 2.8, paste("Control units dropped:",((nrow(sample)-sum(sample$A))-(nrow(Obs)-sum(Obs$A)))), cex=0.8)

#############
# 7. full
#############
match <- Obs <- NULL
match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="full"))
Obs <- match.data(match)
plot.support(Obs,"Full matching", g.pred[row.names(sample) %in% row.names(Obs)], weights=Obs$weights) 
#plot.support(Obs,"Full matching")

# Obs dropped from sample: 
text(0.8, 3, paste("Treated units dropped:",(sum(sample$A)-sum(Obs$A))), cex=0.8)
text(0.8, 2.8, paste("Control units dropped:",((nrow(sample)-sum(sample$A))-(nrow(Obs)-sum(Obs$A)))), cex=0.8)

dev.off()


