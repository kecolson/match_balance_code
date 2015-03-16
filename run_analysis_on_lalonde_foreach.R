# Run matching and analysis procedure on lalonde data
# Ellie Colson
# 3/13/2015


# ````````````````````````````````

# THIS CODE IS NOT FINISHED

# ``````````````````````````````````


rm(list=ls())
set.seed(1.31)
library("parallel")
library("doParallel")
library("foreach")
options(scipen=5)
#SL.library <- c("SL.glm", "SL.glm.interaction", "SL.step", "SL.step.interaction")
#SL.library <- c("SL.glm","SL.glm.interaction", "SL.step", "SL.step.interaction", "SL.earth", "SL.nnet", "SL.gam")
SL.library <- c("SL.glm","SL.glm.interaction","SL.mean","SL.step","SL.step.interaction")

n.analysis <- 10
n.match <- 7

cluster <- T
nCores <- 2

# Choose number of runs: 
index <- 1:4

# Choose data
# Tx choices: Experimental: full, dw, dw_big, dw_big_noisy
# Control choices: Experimental: full, dw, dw_big, dw_big_noisy 
# Observational: cps1, cps2, cps3, psid1, psid2, psid3
dn_tx <- "dw"  # full not working here for now
dn_ct <- "psid3"

if (cluster == F) {
  source("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/match_balance_code/TMLE_ATT/TMLE_ATT_edited.R") 
  setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/data/lalonde/")
} else {
  source("TMLE_ATT_edited.R") 
}

# Function for binomial standard error
bin_sd <- function(p,n) { sqrt(p*(1-p)/n) }

# Cluster setup 
if (cluster==T) {  
  registerDoParallel(nCores)
}

# Setup data
if (dn_tx == "full" & dn_ct == "full") {
  data <- data.frame(read.csv("nsw_full.csv"), re74=NA)
  data <- data[,c("data_id","treat","age","education","black","hispanic","married","nodegree","re74","re75","re78")]
  
} else if (dn_tx == "dw" & dn_ct == "dw") {
  data <- read.csv("nsw_dw.csv")
 
} else if (dn_tx == "dw_big" & dn_ct == "dw_big") {
  data <- read.csv("nsw_dw.csv")
  data <- data[sample(row.names(data), 2000, replace = T),]
 
} else if (dn_tx == "dw_big_noisy" & dn_ct == "dw_big_noisy") {
  data <- read.csv("nsw_dw.csv")
  n <- nrow(data)
  param <- data.frame(treat = bin_sd(mean(data$treat),n), age = sd(data$age), education = sd(data$education),
                      black = bin_sd(mean(data$black),n), hispanic = bin_sd(mean(data$hispanic),n), 
                      married = bin_sd(mean(data$married),n), nodegree = bin_sd(mean(data$nodegree),n), 
                      re74 = sd(data$re74), re75 = sd(data$re75), re78 = sd(data$re78))
  data <- data[sample(row.names(data), 2000, replace = T),]
  n <- nrow(data)
  # Binary vars: treat, black, hispanic, married, nodegree
  for (v in c("treat","black","hispanic","married","nodegree")) {
    data[[v]] <- ifelse(sample(c(T,F), n, replace = T, prob = c(param[[v]],1-param[[v]])), abs(data[[v]]-1), data[[v]])
  }
  # Semi-continuous whole numbers: age, education; Positive continuous: re74, re75, re78
  for (v in c("age","education","re74","re75","re78")) {
    data[[v]] <- data[[v]] + rnorm(n, mean=0, sd = param[[v]]/4)
    data[[v]] <- ifelse(data[[v]]<0, 0, data[[v]])
  }
  
} else if (substr(dn_ct,1,3)=="cps") {
  if (dn_tx == "full") {
    lalonde_tx <- data.frame(read.csv("nsw_full.csv"), re74=NA)
    lalonde_tx <- lalonde_tx[,c("data_id","treat","age","education","black","hispanic","married","nodegree","re74","re75","re78")]
   
  } else if (dn_tx == "dw") {
    lalonde_tx <- read.csv("nsw_dw.csv")
  }
  
  lalonde_ct <- read.csv(paste0("cps_controls",substr(dn_ct,4,4),".csv"))
  data <- rbind(lalonde_tx[lalonde_tx$treat==1,], lalonde_ct)
  
} else if (substr(dn_ct,1,4)=="psid") {
  if (dn_tx == "full") {
    lalonde_tx <- data.frame(read.csv("nsw_full.csv"), re74=NA)
    lalonde_tx <- lalonde_tx[,c("data_id","treat","age","education","black","hispanic","married","nodegree","re74","re75","re78")]
   
  } else if (dn_tx == "dw") {
    lalonde_tx <- read.csv("nsw_dw.csv")
    
  }
  lalonde_ct <- read.csv(paste0("psid_controls",substr(dn_ct,5,5),".csv"))
  data <- rbind(lalonde_tx[lalonde_tx$treat==1,], lalonde_ct)
}

if (dn_tx == "full") {
  qmodel <- "Y.star~A+W1+W2+W3+W4+W5+W6+W7"
  gmodel <- "A~W1+W2+W3+W4+W5+W6+W7"
  covs <- c("W1","W2","W3","W4","W5","W6","W7")
  ncovs <- length(covs)
} else if (dn_tx %in% c("dw","dw_big","dw_big_noisy")) {
  qmodel <- "Y.star~A+W1+W2+W3+W4+W5+W6+W7+W8"
  gmodel <- "A~W1+W2+W3+W4+W5+W6+W7+W8"
  covs <- c("W1","W2","W3","W4","W5","W6","W7","W8")
  ncovs <- length(covs)
}

names(data) <- c("data_id","A","W1","W2","W3","W4","W5","W6","W7","W8","Y.star")

#---------
# ESTIMATORS: function to implement the estimators
#   input: O (obs data)
# 	output: point estimates from unadj, Gcomp for ATE, Gcomp for ATT, IPTW for ATE, IPTW for ATT, TMLE for ATE, and TMLE for ATT
#--------

ESTIMATORS<- function(O) {
  
  unadj<- iptw.att.misp<- iptw.att.sl<- gcomp.att.misp<- gcomp.att.sl<- lalonde1<- lalonde2<- lalonde3<- tmle.att.misp.misp<- tmle.att.sl<- NA
  
  txt<- control <- O
  txt$A <- 1;  control$A <- 0
  n <- nrow(O)
  tx.indices <- O$A==1
  control.indices <- O$A==0

  library("SuperLearner")
  library("tmle")
  library("tmlecte")
  
  ## Unadjusted
  unadj<- mean(O$Y.star[tx.indices]) - mean(O$Y.star[control.indices])
  
  ### IPTW for ATT, misspecified parametric
  g.model<- glm(as.formula(gmodel), family='binomial', data=O)
  g.pred<- predict(g.model, type='resp')
  O$weight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y.star ~ A, weights = weight, data = O)
  iptw.att.misp <- mean(predict(mod, newdata = txt, type='response')) - mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "weight"]
  
  # IPTW for ATT, with superlearner
  SL.out <- SuperLearner(Y = O$A, X = O[,covs], SL.library = SL.library, 
                         family = 'binomial', cvControl = list(V=10))
  g.pred <- SL.out$SL.predict
  O$weight <- ifelse(O$A==1, 1, g.pred/(1-g.pred))
  mod <- glm(Y.star ~ A, weights = weight, data = O)
  iptw.att.sl <- mean(predict(mod, newdata = txt, type='response')) - mean(predict(mod, newdata = control, type='response'))
  O <- O[,names(O) != "weight"]
  
  ### gcomp for ATT with misp linear
  reg.model <- glm(as.formula(qmodel), data=O)
  gcomp.att.misp<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) - mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
 
  ### gcomp with superlearner- for ATT
  newX <- rbind(txt[,c("A",covs)], control[,c("A",covs)])
  SL.out <- SuperLearner(Y = O$Y.star, X = O[,c("A",covs)], SL.library = SL.library, newX = newX, family = 'gaussian', cvControl=list(V=10))
  tx.preds <- SL.out$SL.predict[1:n]
  ct.preds <- SL.out$SL.predict[(n+1):(2*n)]
  gcomp.att.sl <- mean(tx.preds[tx.indices]) - mean(ct.preds[tx.indices])  

  ## Lalonde model 1 (d): RE78 ~ treat + RE75
  lalonde.m1 <- glm(Y.star~A+W8, data=O)
  lalonde1 <- mean(predict(lalonde.m1, newdata = txt[tx.indices,], type='response')) - mean(predict(lalonde.m1, newdata = control[tx.indices,], type='response'))
  
  ## Lalonde model 2 (c): RE78 ~ treat + age + age^2 + yrs edu + no degree + black + hispanic
  lalonde.m2 <- glm(Y.star~A+W1+I(W1^2)+W2+W6+W3+W4, data=O)
  lalonde2 <- mean(predict(lalonde.m2, newdata = txt[tx.indices,], type='response')) - mean(predict(lalonde.m2, newdata = control[tx.indices,], type='response'))
  
  ## Lalonde model 3 (e): RE78 ~ treat + re75 + age + age^2 + yrs edu + no degree + black + hispanic
  lalonde.m3 <- glm(Y.star~A+W8+W1+I(W1^2)+W2+W6+W3+W4, data=O)
  lalonde3 <- mean(predict(lalonde.m3, newdata = txt[tx.indices,], type='response')) - mean(predict(lalonde.m3, newdata = control[tx.indices,], type='response'))
    
  # ### TMLE for ATT
  O$Y.star.scaled <- (O$Y.star-min(O$Y.star))/ (max(O$Y.star)-min(O$Y.star)) # Scale Y.star to be between 0 and 1 so it will play nice with tmle.att
  
  # TMLE for ATT, mispecified parametric
  tmle.out <- tmle.att(Y = O$Y.star.scaled, A = O$A, W = O[,covs], 
                      family="gaussian",
                      g.method = "glm",
                      g.formula = as.formula(gmodel),
                      Q.method = "glm",
                      Q.formula = as.formula(qmodel))
  tmle.att.misp.misp <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end
  
  # TMLE for ATT, with super learner
  tmle.out <- tmle.att(Y=O$Y.star.scaled, A=O$A, W=O[,covs], 
                      family="gaussian",
                      g.method = "SL",
                      g.SL.library = SL.library,
                      Q.method = "SL",
                      Q.SL.library = SL.library)
  tmle.att.sl <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star)) # rescale at end
  
  c(unadj, iptw.att.misp, iptw.att.sl, gcomp.att.misp, gcomp.att.sl, lalonde1, lalonde2, lalonde3, tmle.att.misp.misp, tmle.att.sl)
  
}		



#---------
# ESTIMATORS.weighted: function to implement the 4 estimators with weights
#   input: O (obs data)
# 	output: point estimates from unadj, Gcomp, IPTW and TMLE
#--------

ESTIMATORS.weighted<- function(O) {
  
  unadj<- iptw.att.misp<- iptw.att.sl<- gcomp.att.misp<- gcomp.att.sl<- lalonde1<- lalonde2<- lalonde3<- tmle.att.misp.misp<- tmle.att.sl<- NA
    # Have empty slots for IPTW here, but will never actually run IPTW for weighted analysis
  
  txt<- control <- O
  txt$A <-1;  control$A <- 0
  n <- nrow(O)
  tx.indices <- O$A==1
  control.indices <- O$A==0

  library("SuperLearner")
  library("tmle")
  library("tmlecte")
  
  ## Unadjusted
  unadj<- weighted.mean(O$Y.star[tx.indices], w = O$weights[tx.indices]) - weighted.mean(O$Y.star[control.indices], w = O$weights[control.indices])
  
  ### gcomp with misp linear- for ATT
  reg.model <- glm(as.formula(qmodel), data=O, weights=O$weights)
  gcomp.att.misp<- mean(predict(reg.model, newdata = txt[tx.indices,], type='response')) - mean(predict(reg.model, newdata = control[tx.indices,], type='response'))
  
  ### gcomp with superlearner- for ATT
  newX <- rbind(txt[,c("A",covs)], control[,c("A",covs)])
  SL.out <- SuperLearner(Y = O$Y.star, X = O[,c("A",covs)], SL.library = SL.library, newX = newX, family = 'gaussian', cvControl=list(V=10), obsWeights = O$weights)
  Qbar1W <- SL.out$SL.predict[1:n]
  Qbar0W <- SL.out$SL.predict[(n+1):(2*n)]
  gcomp.att.sl<- mean(Qbar1W[tx.indices]) - mean(Qbar0W[tx.indices])

  ## Lalonde model 1
  lalonde.m1 <- glm(Y.star~A+W8, data=O, weights=O$weights)
  lalonde1 <- mean(predict(lalonde.m1, newdata = txt[tx.indices,], type='response')) - mean(predict(lalonde.m1, newdata = control[tx.indices,], type='response'))
  
  ## Lalonde model 2
  lalonde.m2 <- glm(Y.star~A+W1+I(W1^2)+W2+W6+W3+W4, data=O, weights=O$weights)
  lalonde2 <- mean(predict(lalonde.m2, newdata = txt[tx.indices,], type='response')) - mean(predict(lalonde.m2, newdata = control[tx.indices,], type='response'))
  
  ## Lalonde model 3
  lalonde.m3 <- glm(Y.star~A+W8+W1+I(W1^2)+W2+W6+W3+W4, data=O, weights=O$weights)
  lalonde3 <- mean(predict(lalonde.m3, newdata = txt[tx.indices,], type='response')) - mean(predict(lalonde.m3, newdata = control[tx.indices,], type='response'))
  
  # ### TMLE for ATT
  O$Y.star.scaled <- (O$Y.star-min(O$Y.star))/ (max(O$Y.star)-min(O$Y.star)) # Scale Y.star to be between 0 and 1 so it will play nice with tmle.att
  
  # TMLE for ATT, mispecified parametric
  tmle.out<- tmle.att2(Y = O$Y.star.scaled, A = O$A, 
                      W = O[,covs], 
                      family = "gaussian",
                      Delta=rep(1,nrow(O)),
                      gDelta.method = "user",
                      gDelta.1 = 1/O$weights,
                      g.method = "glm",
                      g.formula = as.formula(gmodel), 
                      Q.method = "glm",
                      Q.formula = as.formula(qmodel))
  tmle.att.misp.misp <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star))
  #tmle.att.misp.misp <- mean(tmle.out$Qstar[,2] - tmle.out$Qstar[,1]) * (max(O$Y.star)-min(O$Y.star))
  #tmle.att.misp.misp <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))
  
  # TMLE for ATT, with superlearner
  tmle.out<- tmle.att2(Y = O$Y.star.scaled, A = O$A, 
                       W = O[,covs], 
                       family = "gaussian",
                       Delta=rep(1,nrow(O)),
                       gDelta.method = "user",
                       gDelta.1 = 1/O$weights,
                       g.method = "SL",
                       g.SL.library = SL.library, 
                       Q.method = "SL",
                       Q.SL.library = SL.library)
  tmle.att.sl <- tmle.out$psi * (max(O$Y.star)-min(O$Y.star))
  #tmle.att.sl <- mean(tmle.out$Qstar[,2] - tmle.out$Qstar[,1]) * (max(O$Y.star)-min(O$Y.star))
  #tmle.att.sl <- (sum(O$weights * (tmle.out$Qstar[,2] - tmle.out$Qstar[,1]))/sum(O$weights) ) * (max(O$Y.star)-min(O$Y.star))
  
  c(unadj, iptw.att.misp, iptw.att.sl, gcomp.att.misp, gcomp.att.sl, lalonde1, lalonde2, lalonde3, tmle.att.misp.misp, tmle.att.sl)
  
}		


#---------
# outer.loop: function to run each matching method and set of estimators
#   input: pop (the population data from which to sample)
#   output: parameter estimates and balance metrics: estimands, bal
#--------

outer.loop <- function(iteration, sample) {
  
  set.seed(as.numeric(iteration))
  
  library("MatchIt")
  library("optmatch")
  library("Matching")
  library("Hmisc")
  
  # Make matrices to put the results in
  all<- data.frame( matrix(NA, nrow = 1, ncol = n.analysis))
  colnames(all)<- c('unadj','iptw.att.misp','iptw.att.sl','gcomp.att.misp','gcomp.att.sl','lalonde1','lalonde2','lalonde3','tmle.att.misp.misp','tmle.att.sl')
  match.full <- match.sub <- match.gen <- match.sl <- match.nn <- match.opt<- all
  
  bal <- data.frame(matrix(NA,nrow=1, ncol=n.match*ncovs))
  colnames(bal)<-paste0( rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."),each = ncovs), rep(covs, n.match) )
  
  # Balance denominators
  denom <- NULL
  for (i in paste0("W",1:ncovs)) {
    denom <- c(denom, i= sqrt((var(sample[sample$A==1,i]) + var(sample[sample$A==0,i])) / 2) )
  }
  names(denom) <- paste0("W",1:ncovs)
  
  ##############
  # 1. ANALYZE THE SAMPLE AS IS
  ############
  print("all")
  all <- ESTIMATORS(sample)
  for (i in paste0("W",1:ncovs)) {
    bal[[paste0("all.",i)]] <- (mean(sample[sample$A==1, i]) - mean(sample[sample$A==0,i])) / denom[[i]]
  }
  
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
  try(match <- suppressWarnings(matchit(as.formula(gmodel), data= sample, method = "nearest", m.order = "random")))
  try(Obs <- match.data(match))
  try(match.nn<- ESTIMATORS(Obs)) ## Estimation
  
  for (i in paste0("W",1:ncovs)) {
    bal[[paste0("match.nn.",i)]] <- (mean(Obs[Obs$A==1, i]) - mean(Obs[Obs$A==0,i])) / denom[[i]]
  }
  
  #############
  # 3. Optimal matching 
  # Analysis does not need to be weighted
  #############
  print("optimal")
  
  # Only works if more control units than treated units 
  if (sum(sample$A)< (nrow(sample)-sum(sample$A))) {
    match <- Obs <- NULL
    try(match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="optimal")))
    try(Obs <- match.data(match))
    try(match.opt<- ESTIMATORS(Obs))
    
    for (i in paste0("W",1:ncovs)) {
      bal[[paste0("match.opt.",i)]] <- (mean(Obs[Obs$A==1, i]) - mean(Obs[Obs$A==0,i])) / denom[[i]]
    }
  }
  
  #############
  # 4. Estimate pscore with superlearner, then nearest neighbor match
  # Analysis does not need to be weighted
  #############
  print("sl")
  SL.out <- match <- m <- Obs <- NULL
  SL.out <- SuperLearner(Y=sample$A, X=sample[,covs], 
                         SL.library=SL.library, family='binomial', cvControl=list(V=10))
  try(match <- suppressWarnings(Match(Y=sample$Y.star, Tr=sample$A, X=SL.out$SL.predict, estimand ="ATT", replace = F, ties = F, M=1, version="fast")))
  try(m <- MatchBalance(as.formula(gmodel), data=sample, match.out = match, nboots=10))
  try(Obs <- rbind(sample[unique(match$index.treated),] , sample[unique(match$index.control),]))
  
  try(match.sl<- ESTIMATORS(Obs))
  
  for (i in paste0("W",1:ncovs)) {
    bal[[paste0("match.sl.",i)]] <- (mean(Obs[Obs$A==1, i]) - mean(Obs[Obs$A==0,i])) / denom[[i]]
  }
  
  #############
  # 5. genetic                           
  # Analysis much be weighted
  #############
  print("genetic")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="genetic", print.level = 0)))
  try(Obs <- match.data(match))
  
  try(match.gen<- ESTIMATORS.weighted(Obs))
  
  for (i in paste0("W",1:ncovs)) {
    bal[[paste0("match.gen.",i)]] <- (weighted.mean(Obs[Obs$A==1, i], w=Obs$weights[Obs$A==1]) - weighted.mean(Obs[Obs$A==0,i], w=Obs$weights[Obs$A==0])) / denom[[i]]
  }
  
  #############
  # 6. subclassification
  # Analysis must be weighted
  # Analyses within each subclass, then average
  #############
  print("subclassification")
  match <- Obs <- NULL
  nsubclass <- floor(sum(sample$A)/33)
  try(match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="subclass", subclass = nsubclass)))
  try(Obs <- match.data(match))
  
  # Estimate within each subclass and then average across the subclasses
  ests <- matrix(NA, nrow = nsubclass, ncol = ncol(match.sub))
  for (s in 1:nsubclass) {
    try(ests[s,]<- ESTIMATORS(Obs[Obs$subclass==s,]))  
  }

  try(match.sub<- apply(ests, 2, weighted.mean, na.rm = T, w = as.numeric(table(Obs$subclass[Obs$A == 1])) ))
  
  for (i in paste0("W",1:ncovs)) {
    bal[[paste0("match.sub.",i)]] <- (weighted.mean(Obs[Obs$A==1, i], w=Obs$weights[Obs$A==1]) - weighted.mean(Obs[Obs$A==0,i], w=Obs$weights[Obs$A==0])) / denom[[i]]
  }
  
  #############
  # 7. full
  # Analysis needs to be weighted
  # Sample sizes within each subclass are very small
  # Estimate using weights
  #############
  print("full matching")
  match <- Obs <- NULL
  try(match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="full")))
  try(Obs <- match.data(match))
  try(match.full<- ESTIMATORS.weighted(Obs))   
  
  for (i in paste0("W",1:ncovs)) {
    bal[[paste0("match.full.",i)]] <- (weighted.mean(Obs[Obs$A==1, i], w=Obs$weights[Obs$A==1]) - weighted.mean(Obs[Obs$A==0,i], w=Obs$weights[Obs$A==0])) / denom[[i]]
  }
  
  estimands<- as.data.frame(t(c(all, match.nn, match.opt, match.sl, match.gen, match.sub, match.full)))
  names(estimands) <- paste0(rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."), each=n.analysis),c("unadj","iptw.att.misp","iptw.att.sl","gcomp.att.misp","gcomp.att.sl","lalonde1","lalonde2","lalonde3","tmle.att.misp.misp","tmle.att.sl"))
  
  list(estimands = estimands, bal = bal)
}


#----------------------------------------------------
# compare different designs and different estimators
#----------------------------------------------------

# Run all the analyses and matching combinations
if (cluster==T) { 
  out <- foreach(i = index) %dopar% {
    cat('Starting ', i, 'th job.\n', sep = '')
    outSub <- outer.loop(i, pop = data)
    cat('Finishing ', i, 'th job.\n', sep = '')
    outSub # this will become part of the out object
  } 
  
} else {
  results <- lapply(index, outer.loop, pop = data) 
}

results <- out
print(head(results))
  
# Collapse results
est <- do.call(rbind, lapply(index, function(x) results[[x]]$estimands ))
balance <- do.call(rbind, lapply(index, function(x) results[[x]]$bal ))

# Coerce everything to make sure it is a data frame and not a list, so it will write.csv properly
est <- as.data.frame(est)
balance <- as.data.frame(balance)
est <- data.frame(lapply(est, as.character), stringsAsFactors=FALSE)
balance <- data.frame(lapply(balance, as.character), stringsAsFactors=FALSE)

# Save
if (cluster == F) {
  setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_02_20/")
  if (length(index) == 1 && index == 1) {
    bal2 <- data.frame(type=c("all","match.nn","match.opt","match.sl","match.gen","match.sub","match.full"), bal10=NA, bal5=NA)
    for (i in 1:n.match) {
      bal2$bal10[i] <- mean( as.numeric(results$bal[(((i-1)*ncovs)+1):(i*ncovs)]<0.10) )
      bal2$bal5[i]  <- mean( as.numeric(results$bal[(((i-1)*ncovs)+1):(i*ncovs)]<0.05) )
    }
    write.csv(bal2, paste0("lalonde_bal2_",dn_tx,"_",dn_ct,".csv"))
  }
  
  write.csv(est, paste0("lalonde_est_",dn_tx,"_",dn_ct,".csv"))
  write.csv(balance, paste0("lalonde_bal_",dn_tx,"_",dn_ct,".csv"))
  
} else {
  write.csv(est, paste0("many_runs_",dn_tx,"_",dn_ct,".csv"))
  write.csv(balance, paste0("many_runs_bal_",dn_tx,"_",dn_ct,".csv"))
}

# END