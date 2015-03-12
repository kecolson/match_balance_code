# Compile latest match/balance results

rm(list=ls())
setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/normalCovs_goodOverlap/2015_03_10")
n.batches <- 1
R <- 1000      # number of sim runs
n.match <- 7 # Number of matching methods
n.analysis <- 12

# Choose balance metric
# asmd_perc: percent of variables with absolute standardized mean difference less than 1%, 5%, 10%, 20%
# asmd_stats: mean, median, and max absolute standardized mean difference
# asmd_prog
# asmd_pscore
# tstat: mean t-statistic
# ks: mean K-S test statistic
bal_met <- "ks"

# Compile
estimands <- read.csv("estimands1.csv")[,2:(n.match*n.analysis+1)]
balance <- read.csv("balance1.csv")
if (n.batches > 1) {
for (k in 2:n.batches) {
  estimands <- rbind(estimands, read.csv(paste0("estimands",k,".csv")))
  balance <- rbind(balance, read.csv(paste0("balance",k,".csv")))
}
}

# Calculate and collapse
if (bal_met == "asmd_perc") {
  bal.collapsed <- data.frame(matrix(NA,nrow=R,ncol=n.match*4))
  names(bal.collapsed) <- paste0(rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."), each = 4),c("percb20","percb10","percb5","percb1"))
  for (i in c("all","match.nn","match.opt","match.sl","match.gen","match.sub","match.full")) {
    bal.collapsed[,names(bal.collapsed)==paste0(i,".percb20")] <- (as.numeric(abs(balance[,names(balance)==paste0(i,".W1")])<0.2) + as.numeric(abs(balance[,names(balance)==paste0(i,".W2")])<0.2))/2
    bal.collapsed[,names(bal.collapsed)==paste0(i,".percb10")] <- (as.numeric(abs(balance[,names(balance)==paste0(i,".W1")])<0.1) + as.numeric(abs(balance[,names(balance)==paste0(i,".W2")])<0.1))/2
    bal.collapsed[,names(bal.collapsed)==paste0(i,".percb5")] <-  (as.numeric(abs(balance[,names(balance)==paste0(i,".W1")])<0.05) + as.numeric(abs(balance[,names(balance)==paste0(i,".W2")])<0.05))/2
    bal.collapsed[,names(bal.collapsed)==paste0(i,".percb1")] <-  (as.numeric(abs(balance[,names(balance)==paste0(i,".W1")])<0.01) + as.numeric(abs(balance[,names(balance)==paste0(i,".W2")])<0.01))/2
  }
  bal.collapsed <- apply(bal.collapsed, 2, function(x) mean(x,na.rm=T) )

} else if (bal_met == "asmd_stats") {
  bal.collapsed <- data.frame("all" = NA, "match.nn" = NA, "match.opt" = NA,
                              "match.sl" = NA, "match.gen" = NA, 
                              "match.sub" = NA, "match.full" = NA)
  for (i in names(bal.collapsed)) {
    bal.collapsed[1,i] <- mean(apply(abs(balance[,grep(i,names(balance))]),1,mean))
    bal.collapsed[2,i] <- mean(apply(abs(balance[,grep(i,names(balance))]),1,median))
    bal.collapsed[3,i] <- mean(apply(abs(balance[,grep(i,names(balance))]),1,max))
  }
  row.names(bal.collapsed) <- c("mean","median","max")

} else if (bal_met %in% c("tstat","ks")) {
  bal.collapsed <- data.frame("all" = NA, "match.nn" = NA, "match.opt" = NA,
                              "match.sl" = NA, "match.gen" = NA, 
                              "match.sub" = NA, "match.full" = NA)
  for (i in names(bal.collapsed)) {
    bal.collapsed[1,i] <- mean(apply(abs(balance[,grep(i,names(balance))]),1,mean))
  }

} else if (bal_met %in% c("asmd_prog","asmd_pscore")) {
  bal.collapsed <- apply(balance, 2, function(x) mean(abs(x), na.rm=T))
}

###
MSE.funct<- function( X, Xtrue) {
  ave <- mean(X, na.rm = T)  	
  bias <- mean(X- Xtrue, na.rm = T)
  var <- var(X, na.rm = T)
  MSE <- mean( (X-Xtrue)^2, na.rm = T)
  c(ave,bias, var, MSE)
}
###

ATE <- read.csv("ATE.csv")[1,2]
ATT <- read.csv("ATT.csv")[1,2]

# ATE
done.ate <- apply(estimands, 2, MSE.funct, Xtrue=ATE)
row.names(done.ate)<- c('Ave', 'Bias', 'Var', 'MSE')
done.ate <- done.ate*1000

# ATT
done.att <- apply(estimands, 2, MSE.funct, Xtrue=ATT)
row.names(done.att)<- c('Ave', 'Bias', 'Var', 'MSE')
done.att <- done.att*1000

# Save
write.csv(done.ate, "done_ATE.csv", row.names = F)
write.csv(done.att, "done_ATT.csv", row.names = F)
write.csv(bal.collapsed, paste0("bal.collapsed_",bal_met,".csv"))

