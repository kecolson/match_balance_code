#--------------
# Simulations to examine matching, balance measures, and performance of different matching and analysis options
# Based on Laura Balzer's simulations created for Jen Ahern's high risk, high reward grant
# Version 2: adjusted treatment and outcome mechanisms; added IPTW and TMLE analysis methods where propensity score is not re-estimated; 
# incorporated distinctions and estimates of ATT and ATE, etc
# Ellie Colson
# 11/27/14
#----------------------
rm(list=ls())
options(scipen=5)

SL.library <- c("SL.glm", "SL.glm.interaction", "SL.step", "SL.step.interaction")

# Settings:
# Normal or non-normal covariates?
covars.norm <- T
# Propensity score overlap (1=good, 2= medium, 3=poor)
pscore.type <- 1
# Use cluster?
cluster <- T

n <- 100000
R <- 1000
index <- 1:R
ss <- 1000
n.analysis <- 12
n.match <- 7

# Cluster setup 
if (cluster==T) {
  library("Rmpi")
  
  # Spawn as many slaves as possible
  mpi.spawn.Rslaves()
  
  # In case R exits unexpectedly, have it automatically clean up
  # resources taken up by Rmpi (slaves, memory, etc...)
  .Last <- function(){
    if (is.loaded("mpi_initialize")){
      if (mpi.comm.size(1) > 0){
        print("Please use mpi.close.Rslaves() to close slaves.")
        mpi.close.Rslaves()
      }
      print("Please use mpi.quit() to quit R")
      .Call("mpi_finalize")
    }
  }
  
  # Tell all slaves to return a message identifying themselves
  mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
}

# Source edited version of tmle att so we can use weights properly
if (cluster == T) {
  source("TMLE_ATT_edited.R")
} else {
  setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse")
  source("code/TMLE_ATT/TMLE_ATT_edited.R") 
}

#--------------
# generateData:
# 	input: n
# 	output: full data
#------------------
generateData<- function(n, covars.norm=T, pscore.type = 1) {
  
  if (covars.norm == T) {
    W1 <- runif(n, min = 0.02, max = 0.7)  
    #W2 <- rnorm(n, mean = (0.2 + 0.0125*W1), sd = 1) 
    W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1) 
  } else {
    W1 <- sample(c(0,0.5,1), n, replace = T, prob=c(0.2,0.3,0.5))
    W2 <- rgamma(n, 2+2*W1, 20) 
  }
  
  if (pscore.type==1) {
    # A <- A1 <- rbinom(n, 1, plogis(-1.3 + W1 + I(W1^2) + I(W1*W2) - W2))
    # A <- A1 <- rbinom(n, 1, plogis(-.5 + W1 + 0.1*I(W1^2) + 0.5*I(W1*W2) - W2))
    A <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
  } else if (pscore.type==2) {
    A <- A2 <- rbinom(n, 1, plogis(1.3 + I(-5*W1) + I(2*W1^2) + I(3*W1*W2) - W2))
  } else if (pscore.type==3) {
    A <- A3 <- rbinom(n, 1, plogis(1.2 + I(-10*W1) + I(4*W1^2) + I(6*W1*W2) - I(-2*W2)))
  }
  
  ## Generate outcomes
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
  Y <-  ifelse(A==1,Y.1,Y.0)
  
  Y.star<- Y
  Y.star.1<- Y.1
  Y.star.0<- Y.0
  
  data.frame(W1,W2, A,Y, Y.star, Y.1, Y.0, Y.star.1, Y.star.0)
}

# get the truth
set.seed(123)
O.true<- generateData(n=10000000, covars.norm=T, pscore.type=1)

ATE <- mean(O.true$Y.star.1 - O.true$Y.star.0)
ATT <- mean(O.true$Y.star.1[O.true$A==1] - O.true$Y.star.0[O.true$A==1])

write.csv(data.frame(ATE), "ATE.csv")
write.csv(data.frame(ATT), "ATT.csv")

#---------
# ESTIMATORS: function to implement the estimators
#   input: O (obs data)
# 	output: point estimates from unadj, Gcomp for ATE, Gcomp for ATT, IPTW for ATE, IPTW for ATT, TMLE for ATE, and TMLE for ATT
#--------

ESTIMATORS<- function(O, pscore=NA) {
  
  unadj<- iptw.att.misp<- iptw.att.cor<- iptw.att.sl<- gcomp.att.misp<- gcomp.att.cor<- gcomp.att.sl<- tmle.att.misp.misp <-
    tmle.att.misp.cor<- tmle.att.cor.misp<- tmle.att.cor.cor<- tmle.att.sl<- NA
  
  txt<- control <- O
  txt$A <-1;  control$A <- 0
  n <- nrow(O)
  tx.indices <- O$A==1
  control.indices <- O$A==0
  
  ## Unadjusted
  unadj<- mean(O$Y.star[tx.indices]) - mean(O$Y.star[control.indices])
  
  ### IPTW for ATT, misspecified
  g.model<- glm(A~ W1 + W2, family='binomial', data=O)
  g.pred<- predict(g.model, type='resp')
  O$weight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y ~ A, weights = weight, data = O)
  iptw.att.misp <- mean(predict(mod, newdata = txt, type='response')) - 
    mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "weight"]
  
  ### IPTW for ATT, correctly specified
  g.model<- glm(A~ poly(W1,2) + W2 + W1:W2, family='binomial', data = O)
  g.pred<- predict(g.model, type = 'resp')
  O$weight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y ~ A, weights = weight, data = O)
  iptw.att.cor <- mean(predict(mod, newdata = txt, type='response')) - 
    mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "weight"]
  
  # IPTW for ATT, with superlearner
  library("SuperLearner")
  SL.out <- SuperLearner(Y = O$A, X = O[,c("W1","W2")], SL.library = SL.library, 
                         family = 'binomial', cvControl=list(V=10))
  g.pred <- SL.out$SL.predict
  O$weight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y ~ A, weights = weight, data = O)
  iptw.att.sl <- mean(predict(mod, newdata = txt, type='response')) - 
    mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "weight"]
  
  ### gcomp with misp linear- for ATT
  reg.model <- glm(Y.star ~ A + W2 + W1, data=O)
  gcomp.att.misp<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) -
   mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
  
  ### gcomp with correctly specified linear- for ATT
  reg.model <- glm(Y.star ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2, data=O)
  gcomp.att.cor<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) -
   mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
  
  ### gcomp with superlearner- for ATT
  library("SuperLearner")
  newX <- rbind(O, txt, control)[,c("A","W1","W2")]
  SL.out <- SuperLearner(Y = O$Y.star, X = O[,c("A","W1","W2")], SL.library = SL.library, newX = newX,
                        family = 'gaussian', cvControl=list(V=10))
  QbarAW <- SL.out$SL.predict[1:n]
  Qbar1W <- SL.out$SL.predict[(n+1):(2*n)]
  Qbar0W <- SL.out$SL.predict[(2*n+1):(3*n)]
  gcomp.att.sl<- mean(Qbar1W[tx.indices]) - mean(Qbar0W[tx.indices])

  # TMLE for ATT
  library("tmle")
  library("tmlecte")
  O$Y.star.scaled <- (O$Y.star-min(O$Y.star))/ (max(O$Y.star)-min(O$Y.star)) # Scale Y.star to be between 0 and 1 so it will play nice with tmle.att

  # ### TMLE for ATT, mispecified tx mech, mispecified outcome
  tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, W = O[,c("W1","W2")], 
                     family = "gaussian", 
                     g.method = "glm",
                     g.formula = A ~ W1 + W2,
                     Q.method = "glm", 
                     Q.formula = Y ~ A + W1 + W2)
  tmle.att.misp.misp <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end

  # ### TMLE for ATT, mispecified tx mech, correct outcome
  tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, W = O[,c("W1","W2")], 
                      family = "gaussian", 
                      g.method = "glm",
                      g.formula = A ~ W1 + W2,
                      Q.method = "glm", 
                      Q.formula = Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2)
  tmle.att.misp.cor <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end
  
  # ### TMLE for ATT, correct tx mech, mispecified outcome
  tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, W = O[,c("W1","W2")], 
                      family = "gaussian", 
                      g.method = "glm",
                      g.formula = A ~ poly(W1,2) + W2 + W1:W2,
                      Q.method = "glm", 
                      Q.formula = Y ~ A + W1 + W2)
  tmle.att.cor.misp <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end
  
  # ### TMLE for ATT, correct tx mech, correct outcome
  tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, W = O[,c("W1","W2")], 
                      family = "gaussian", 
                      g.method = "glm",
                      g.formula = A ~ poly(W1,2) + W2 + W1:W2,
                      Q.method = "glm", 
                      Q.formula = Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2)
  tmle.att.cor.cor <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end
  
  # ### TMLE for ATT, with super learner
  tmle.out<- tmle.att(Y=O$Y.star.scaled, A=O$A, W = O[,c("W1","W2")], 
                      family = "gaussian", 
                      g.method = "SL",
                      g.SL.library = SL.library,
                      Q.method = "SL", 
                      Q.SL.library = SL.library)
  tmle.att.sl <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end
  
  c(unadj, iptw.att.misp, iptw.att.cor, iptw.att.sl, gcomp.att.misp, gcomp.att.cor, gcomp.att.sl, tmle.att.misp.misp, tmle.att.misp.cor, tmle.att.cor.misp, tmle.att.cor.cor, tmle.att.sl)
  
}		



#---------
# ESTIMATORS.weighted: function to implement the 4 estimators with weights
#   input: O (obs data)
# 	output: point estimates from unadj, Gcomp, IPTW and TMLE
#--------

ESTIMATORS.weighted<- function(O, pscore=NA) {
  
  unadj<- iptw.att.misp<- iptw.att.cor<- iptw.att.sl<- gcomp.att.misp<- gcomp.att.cor<- gcomp.att.sl<- tmle.att.misp.misp <-
    tmle.att.misp.cor<- tmle.att.cor.misp<- tmle.att.cor.cor<- tmle.att.sl<- NA
  
  txt<- control <- O
  txt$A <-1;  control$A <- 0
  n <- nrow(O)
  tx.indices <- O$A==1
  control.indices <- O$A==0
  
  ## Unadjusted
  unadj<- weighted.mean(O$Y.star[tx.indices], w = O$weights[tx.indices]) - weighted.mean(O$Y.star[control.indices], w = O$weights[control.indices])
  
  ### gcomp with misp linear- for ATT
  reg.model <- glm(Y.star~ A +W2+W1, data=O, weights=O$weights)
  gcomp.att.misp<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) -
   mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
  
  ### gcomp with correctly specified linear- for ATT
  reg.model <- glm(Y.star~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2, data=O, weights=O$weights)
  gcomp.att.cor<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) -
   mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
  
  ### gcomp with superlearner- for ATT
  library("SuperLearner")
  newX <- rbind(O, txt, control)[,c("A","W1","W2")]
  SL.out <- SuperLearner(Y = O$Y.star, X = O[,c("A","W1","W2")], SL.library = SL.library, newX = newX,
                           family = 'gaussian', cvControl=list(V=10), obsWeights = O$weights)
  QbarAW <- SL.out$SL.predict[1:n]
  Qbar1W <- SL.out$SL.predict[(n+1):(2*n)]
  Qbar0W <- SL.out$SL.predict[(2*n+1):(3*n)]
  gcomp.att.sl<- mean(Qbar1W[tx.indices]) - mean(Qbar0W[tx.indices])
  
  # ### TMLE for ATT, 
  library("tmle")
  library("tmlecte") 
  O$Y.star.scaled <- (O$Y.star-min(O$Y.star))/ (max(O$Y.star)-min(O$Y.star)) # Scale Y.star to be between 0 and 1 so it will play nice with tmle.att
  
  # TMLE for ATT, mispecified treatment mechanism, mispecified outcome
  tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A,  # Using tmle.att2 here, which is my adjusted function to make weights work properly
                     W = O[,c("W1","W2")], 
                     family = "gaussian", 
                     Delta = rep(1,nrow(O)), 
                     gDelta.method = "user",
                     gDelta.1 = 1/O$weights,
                     g.method = "glm",
                     g.formula = A ~ W1 + W2, 
                     Q.method = "glm",
                     Q.formula = Y ~ A + W1 + W2)
  tmle.att.misp.misp <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))

  # TMLE for ATT, mispecified treatment mechanism, correct outcome
  tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A,  # Using tmle.att2 here, which is my adjusted function to make weights work properly
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "glm",
                       g.formula = A ~ W1 + W2, 
                       Q.method = "glm",
                       Q.formula = Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2)
  tmle.att.misp.cor <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))
  
  # TMLE for ATT, correct treatment mechanism, mispecified outcome
  tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A,  # Using tmle.att2 here, which is my adjusted function to make weights work properly
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "glm",
                       g.formula = A ~ poly(W1,2) + W2 + W1:W2, 
                       Q.method = "glm",
                       Q.formula = Y ~ A + W1 + W2)
  tmle.att.cor.misp <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))
  
  # TMLE for ATT, correct treatment mechanism, correct outcome
  tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A,  # Using tmle.att2 here, which is my adjusted function to make weights work properly
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "glm",
                       g.formula = A ~ poly(W1,2) + W2 + W1:W2, 
                       Q.method = "glm",
                       Q.formula = Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2)
  tmle.att.cor.cor <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))
  
  # TMLE for ATT, with superlearner
  tmle.out<- tmle.att2(Y=O$Y.star.scaled, A=O$A,  # Using tmle.att2 here, which is my adjusted function to make weights work properly
                       W = O[,c("W1","W2")], 
                       family = "gaussian", 
                       Delta = rep(1,nrow(O)), 
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "SL",
                       g.SL.library = SL.library, 
                       Q.method = "SL",
                       Q.SL.library = SL.library)
  tmle.att.sl <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))

  c(unadj, iptw.att.misp, iptw.att.cor, iptw.att.sl, gcomp.att.misp, gcomp.att.cor, gcomp.att.sl, tmle.att.misp.misp, tmle.att.misp.cor, tmle.att.cor.misp, tmle.att.cor.cor, tmle.att.sl)
  
}		


#---------
# outer.loop: function to run each matching method and set of estimators
#   input: pop (the population data from which to sample)
#   output: parameter estimates and balance metrics: estimands, bal
#--------

outer.loop <- function(iteration, pop, ss) {

  set.seed(as.numeric(iteration+1234))
  
  library("MatchIt")
  library("optmatch")
  library("Matching")
  library("Hmisc")
  
  # Make matrices to put the results in
  all<-  data.frame( matrix(NA, nrow=1, ncol = n.analysis))
  colnames(all)<- c('unadj','iptw.att.misp','iptw.att.cor','iptw.att.sl','gcomp.att.misp','gcomp.att.cor','gcomp.att.sl','tmle.att.misp.misp','tmle.att.misp.cor','tmle.att.cor.misp','tmle.att.cor.cor','tmle.att.sl')
  match.full <- match.sub <- match.gen <- match.sl <- match.nn <- match.opt<- all
  
  bal <- data.frame(matrix(NA,nrow=1, ncol = 2*n.match))
  colnames(bal)<-paste0( rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."),each=2), rep(c("W1","W2"), n.match) )
  
  # Take a random (representative) sample from the population
  sample <- pop[sample(row.names(pop), size=ss, replace=F),]
  
  # Balance denominators
  dW1 <- sqrt((var(sample$W1[sample$A==1]) + var(sample$W1[sample$A==0])) / 2)
  dW2 <- sqrt((var(sample$W2[sample$A==1]) + var(sample$W2[sample$A==0])) / 2)
  
  ##############
  # 1. ANALYZE THE SAMPLE AS IS
  ############
  print("all")
  all <- ESTIMATORS(sample, pscore=NA)
  bal$all.W1 <- (mean(sample$W1[sample$A==1]) - mean(sample$W1[sample$A==0])) / dW1
  bal$all.W2 <- (mean(sample$W2[sample$A==1]) - mean(sample$W2[sample$A==0])) / dW2
  
  ############# 
  # Different matching designs
  # Match on W1 (percent black) and W2 (baseline violence)
  # Keep default support
  ##############
  
  #############
  # 2. Nearest neighbor matching
  # Analysis does not need to be weighted
  #############
  print("nearest neighbor")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="nearest")))
  try(Obs <- match.data(match))
  try(match.nn<- ESTIMATORS(Obs, pscore=Obs$distance)) ## Estimation
  
  try(bal$match.nn.W1 <- (mean(Obs$W1[Obs$A==1]) - mean(Obs$W1[Obs$A==0])) / dW1 )
  try(bal$match.nn.W2 <- (mean(Obs$W2[Obs$A==1]) - mean(Obs$W2[Obs$A==0])) / dW2 )
  
  #############
  # 3. Optimal matching 
  # Analysis does not need to be weighted
  #############
  print("optimal")
  
  # Only works if more control units than treated units 
  if (sum(sample$A)< (nrow(sample)-sum(sample$A))) {
    match <- Obs <- NULL
    try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="optimal")))
    try(Obs <- match.data(match))
    try(match.opt<- ESTIMATORS(Obs, pscore=Obs$distance))
    
    try(bal$match.opt.W1 <- (mean(Obs$W1[Obs$A==1]) - mean(Obs$W1[Obs$A==0])) / dW1 )
    try(bal$match.opt.W2 <- (mean(Obs$W2[Obs$A==1]) - mean(Obs$W2[Obs$A==0])) / dW2 )
  }
  
  #############
  # 4. Estimate pscore with superlearner, then nearest neighbor match
  # Analysis does not need to be weighted
  #############
  print("sl")
  SL.out <- match <- m <- Obs <- NULL
  SL.out <- SuperLearner(Y=sample$A, X=sample[,c("W1","W2")], SL.library = SL.library, family = 'binomial', cvControl = list(V=10))
  try(match <- suppressWarnings(Match(Y=sample$Y.star, Tr=sample$A, X=SL.out$SL.predict, estimand ="ATT", replace = F, ties = F, M=1, version="fast")))
  try(m <- MatchBalance(A ~ W1 + W2, data = sample, match.out = match, nboots = 10))
  try(Obs <- rbind(sample[unique(match$index.treated),], sample[unique(match$index.control),]))
  
  try(match.sl<- ESTIMATORS(Obs, pscore=SL.out$SL.predict[row.names(sample) %in% row.names(Obs)]))
  
  try(bal$match.sl.W1 <- (mean(Obs$W1[Obs$A==1]) - mean(Obs$W1[Obs$A==0])) / dW1 )
  try(bal$match.sl.W2 <- (mean(Obs$W2[Obs$A==1]) - mean(Obs$W2[Obs$A==0])) / dW2 )
  
  #############
  # 5. genetic                           
  # Analysis much be weighted
  #############
  print("genetic")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="genetic")))
  try(Obs <- match.data(match))
  
  try(match.gen<- ESTIMATORS.weighted(Obs, pscore=Obs$distance))
  
  # Use weights in balance metric too
  try(bal$match.gen.W1 <- (weighted.mean(Obs$W1[Obs$A==1], w = Obs$weights[Obs$A==1]) - weighted.mean(Obs$W1[Obs$A==0], w = Obs$weights[Obs$A==0])) / dW1 )
  try(bal$match.gen.W2 <- (weighted.mean(Obs$W2[Obs$A==1], w = Obs$weights[Obs$A==1]) - weighted.mean(Obs$W2[Obs$A==0], w = Obs$weights[Obs$A==0])) / dW2 )
  
  #############
  # 6. subclassification
  # Analysis must be weighted
  # Analyses within each subclass, then average
  #############
  print("subclassification")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="subclass", subclass=12)))
  try(Obs <- match.data(match))
  
  # Estimate within each subclass and then average across the subclasses
  # Since we are estimating the ATT here, there is no need to weight this average across the subclasses, because each subclass has the same number of treated units.
  ests <- matrix(NA, nrow=12, ncol=ncol(match.sub))
  for (s in 1:12) {
    try(ests[s,]<- ESTIMATORS.weighted(Obs[Obs$subclass==s,], pscore= Obs$distance[Obs$subclass==s]))  
  }
  
  try(match.sub<- apply(ests,2, function(x) mean(x, na.rm = T)))
  
  # Use weights 
  try(bal$match.sub.W1 <- (weighted.mean(Obs$W1[Obs$A==1], w = Obs$weights[Obs$A==1]) - weighted.mean(Obs$W1[Obs$A==0], w=Obs$weights[Obs$A==0])) / dW1 )
  try(bal$match.sub.W2 <- (weighted.mean(Obs$W2[Obs$A==1], w = Obs$weights[Obs$A==1]) - weighted.mean(Obs$W2[Obs$A==0], w=Obs$weights[Obs$A==0])) / dW2  )
  
  #############
  # 7. full
  # Analysis needs to be weighted
  # Sample sizes within each subclass are very small
  # Estimate using weights
  #############
  print("full matching")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(A ~ W1 + W2, data=sample, method="full")))
  try(Obs <- match.data(match))
  try(match.full<- ESTIMATORS.weighted(Obs, pscore=Obs$distance))   
  
  # Use weights
  try(bal$match.full.W1 <- (weighted.mean(Obs$W1[Obs$A==1], w = Obs$weights[Obs$A==1]) - weighted.mean(Obs$W1[Obs$A==0], w=Obs$weights[Obs$A==0])) / dW1 )
  try(bal$match.full.W2 <- (weighted.mean(Obs$W2[Obs$A==1], w = Obs$weights[Obs$A==1]) - weighted.mean(Obs$W2[Obs$A==0], w=Obs$weights[Obs$A==0])) / dW2 )

  estimands<- as.data.frame(t(c(all, match.nn, match.opt, match.sl, match.gen, match.sub, match.full)))
  names(estimands) <- paste0(rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."), each = n.analysis),
                             c("unadj","iptw.att.misp","iptw.att.cor","iptw.att.sl","gcomp.att.misp","gcomp.att.cor","gcomp.att.sl","tmle.att.misp.misp","tmle.att.misp.cor","tmle.att.cor.misp","tmle.att.cor.cor","tmle.att.sl"))
  
  list(estimands = estimands, bal = bal)
}



########################
# compare different designs and different estimators

# Balance:
# Standardized mean difference = abs(mean treated-mean control) / SD pooled
# % of vars below 20%, 10%, 5%  SMD

# generate the population
set.seed(123)
pop <- generateData(n=n, covars.norm=T, pscore.type = 1)

# Run all the analyses and matching combinations
if (cluster==T) { 
  # Send all necessary info to the slaves
  mpi.bcast.Robj2slave(ESTIMATORS)  
  mpi.bcast.Robj2slave(ESTIMATORS.weighted)
  mpi.bcast.Robj2slave(outer.loop)
  mpi.bcast.Robj2slave(index)
  mpi.bcast.Robj2slave(pop)
  mpi.bcast.Robj2slave(ss) 
  mpi.bcast.Robj2slave(ATE) 
  mpi.bcast.Robj2slave(ATT)  
  mpi.bcast.Robj2slave(n.analysis)  
  mpi.bcast.Robj2slave(n.match)  
  mpi.bcast.Robj2slave(SL.library) 
  mpi.bcast.Robj2slave(.setColnames) # These last ones are functions needed to run my edited version of tmle.att2
  mpi.bcast.Robj2slave(.bound)
  mpi.bcast.Robj2slave(regress)
  mpi.bcast.Robj2slave(predict.regress)
  mpi.bcast.Robj2slave(tmle.att2)
  mpi.bcast.Robj2slave(print.cte)
  mpi.bcast.Robj2slave(tmle.cte)
  mpi.bcast.Robj2slave(tmle.nde)
  
  # Run our big function in parallel
  results <- mpi.parLapply(index, outer.loop, pop=pop, ss=ss) 
} else {
  results <- lapply(index, outer.loop, pop=pop, ss=ss) 
}

# Turn off slaves
if (cluster==T) { mpi.close.Rslaves() }

print(head(results))
print(class(results))
print(head(results[[1]]))
print(class(results[[1]]))
print(head(results[[1]]$estimands))
print(class(results[[1]]$estimands))

# Collapse results
est <- do.call(rbind, lapply(index, function(x) results[[x]]$estimands ))
balance <- do.call(rbind, lapply(index, function(x) results[[x]]$bal ))

# Save
write.csv(data.frame(lapply(est, as.character), stringsAsFactors=FALSE), "estimands1.csv")
write.csv(data.frame(lapply(balance, as.character), stringsAsFactors=FALSE), "balance1.csv")

# Close
if (cluster==T) { mpi.quit() }




