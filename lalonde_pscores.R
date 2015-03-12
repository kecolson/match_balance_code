# Estimate propensity scores for different datasets and matching methods
# Ellie Colson
# 3/5/2015

rm(list=ls())
set.seed(1.31)
library("Matching")
library("MatchIt")
library("optmatch")
library("SuperLearner")
library("Hmisc")
options(scipen=5)
SL.library <- c("SL.glm","SL.glm.interaction","SL.mean","SL.step","SL.step.interaction")
setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/data/lalonde/")

# Choose data
# Tx choices: Experimental: full, dw; 
# Control choices: Experimental: full, dw; Observational: cps1, cps2, cps3, psid1, psid2, psid3
dn_tx <- "dw"  # full not working here for now
dn_ct <- "psid3"

if (dn_tx == "full" & dn_ct == "full") {
  data <- data.frame(read.csv("nsw_full.csv"), re74=NA)
  data <- data[,c("data_id","treat","age","education","black","hispanic","married","nodegree","re74","re75","re78")]
  
} else if (dn_tx == "dw" & dn_ct == "dw") {
  data <- read.csv("nsw_dw.csv")
 
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
} else if (dn_tx == "dw") {
  qmodel <- "Y.star~A+W1+W2+W3+W4+W5+W6+W7+W8"
  gmodel <- "A~W1+W2+W3+W4+W5+W6+W7+W8"
  covs <- c("W1","W2","W3","W4","W5","W6","W7","W8")
  ncovs <- length(covs)
}

# Order of data names for dw, psid1, psid2, psid3, cps1, cps2, cps3
# "data_id"   "treat"     "age"       "education" "black"     "hispanic"  "married"   "nodegree"  "re74"      "re75"      "re78"  
names(data) <- c("data_id","A","W1","W2","W3","W4","W5","W6","W7","W8","Y.star")

est.pscore.para <- function(O) {
  set.seed(1)
  # Estimate parametric propensity scores
  g.model <- glm(as.formula(gmodel), family = 'binomial', data = O)
  g.pred <- predict(g.model, type='resp')
  return(g.pred)
}

est.pscore.sl <- function(O) {
  set.seed(1)
  # Estimate superlearner propensity scores
  SL.out <- SuperLearner(Y = O$A, X = O[,covs], SL.library = SL.library, 
                         family = 'binomial', cvControl = list(V=10))
  g.pred <- SL.out$SL.predict
  return(g.pred)
}

sample <- data

results <- data.frame(method = c("all","nn","opt","sl","gen","sub","full"), para.min = NA, para.max = NA, sl.min = NA, sl.max = NA)

# Sample as is
print("ALL DATA")
g.pred1 <- est.pscore.para(sample)
summary(g.pred1)
results[results$method=="all",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(sample)
summary(g.pred2)
results[results$method=="all",4:5] <- range(g.pred2)
warnings()

# NN
print("NEAREST NEIGHBOR")
match <- Obs <- NULL
match <- suppressWarnings(matchit(as.formula(gmodel), data= sample, method = "nearest"))
Obs <- match.data(match)
g.pred1 <- est.pscore.para(Obs)
summary(g.pred1)
results[results$method=="nn",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(Obs)
summary(g.pred2)
results[results$method=="nn",4:5] <- range(g.pred2)
warnings()
 
# Opt
print("OPTIMAL")
match <- Obs <- NULL
match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="optimal"))
Obs <- match.data(match)
g.pred1 <- est.pscore.para(Obs)
summary(g.pred1)
results[results$method=="opt",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(Obs)
summary(g.pred2)
results[results$method=="opt",4:5] <- range(g.pred2)
warnings()
    
# NN + SL
print("SL + NN")
SL.out <- match <- m <- Obs <- NULL
SL.out <- SuperLearner(Y=sample$A, X=sample[,covs], 
                       SL.library=SL.library, family='binomial', cvControl=list(V=10))
match <- suppressWarnings(Match(Y=sample$Y.star, Tr=sample$A, X=SL.out$SL.predict, estimand ="ATT", replace = F, ties = F, M=1, version="fast"))
m <- MatchBalance(as.formula(gmodel), data=sample, match.out = match, nboots=10)
Obs <- rbind(sample[unique(match$index.treated),] , sample[unique(match$index.control),])
g.pred1 <- est.pscore.para(Obs)
summary(g.pred1)
results[results$method=="sl",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(Obs)
summary(g.pred2)
results[results$method=="sl",4:5] <- range(g.pred2)
warnings()

# Gen
print("GENETIC")
match <- Obs <- NULL
match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="genetic", print.level = 0))
Obs <- match.data(match)
g.pred1 <- est.pscore.para(Obs)
summary(g.pred1)
results[results$method=="gen",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(Obs)
summary(g.pred2)
results[results$method=="gen",4:5] <- range(g.pred2)
warnings()

# Sub
print("SUBCLASS")
match <- Obs <- NULL
nsubclass <- floor(sum(sample$A)/33)
match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="subclass", subclass = nsubclass))
Obs <- match.data(match)
g.pred1 <- est.pscore.para(Obs)
summary(g.pred1)
results[results$method=="sub",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(Obs)
summary(g.pred2)
results[results$method=="sub",4:5] <- range(g.pred2)
warnings()
  
# Full 
print("FULL")
match <- Obs <- NULL
match <- suppressWarnings(matchit(as.formula(gmodel), data=sample, method="full"))
Obs <- match.data(match)
g.pred1 <- est.pscore.para(Obs)
summary(g.pred1)
results[results$method=="full",2:3] <- range(g.pred1)
warnings()
g.pred2 <- est.pscore.sl(Obs)
summary(g.pred2)
results[results$method=="full",4:5] <- range(g.pred2)
warnings()

results

write.csv(results, paste0("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/pscores/",dn_ct,".csv"))

results[,2:5] <- round(results[,2:5],3)

write.csv(results, paste0("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/pscores/",dn_ct,"_rounded.csv"))



