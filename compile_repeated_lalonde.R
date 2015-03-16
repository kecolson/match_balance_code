# Ellie Colson
# 3/2/2015
# Compile results from repeatedly running analysis of lalonde experimental data

rm(list = ls())
options(scipen=10)
library(tidyr)
library(ggplot2)
library(ggvis)
n.analysis <- 10
n.match <- 7

match.methods <- c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full.")
an.methods <- c("unadj","iptw.att.misp","iptw.att.sl","gcomp.att.misp","gcomp.att.sl","lalonde1","lalonde2","lalonde3","tmle.att.misp.misp","tmle.att.sl")

# DW- DW
r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_02/many_runs_dw_dw_est.csv")
row.names(r) <- paste0(rep(match.methods, each = n.analysis),
                       an.methods)
r <- as.data.frame(r)
r <- r[!row.names(r) %in% paste0(rep(c("match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."), each = 2),c("iptw.att.misp","iptw.att.sl")),]
r <- as.data.frame(r)

# Examine distributions
pdf("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_02/distributions_of_estimates.pdf", height = 5, width = 8)
for (i in 1:nrow(r)) {
  temp <- round(as.numeric(r[i,]),6)
  dim <- min(c(30,length(unique(temp))))
  if (dim > 1) hist(temp, main = paste(row.names(r)[i]), breaks = dim)
}
dev.off()

# by match method
pdf("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_02/distributions_of_estimates_by_match.pdf", height = 5, width = 8)
for (m in match.methods) {
  if (m=="all.") par(mfrow = c(3,4), mar=c(3,2,2,1))
  if (m!="all.") par(mfrow = c(3,3), mar=c(3,2,2,1))
  for (a in an.methods) {
    if (m!="all." & (a=="iptw.att.misp" | a=="iptw.att.sl")) next
    temp <- round(as.numeric(r[row.names(r)==paste0(m,a),]),6)
    dim <- min(c(30,length(unique(temp))))
    if (dim > 1) hist(temp, main = paste(row.names(r)[row.names(r)==paste0(m,a)]), breaks = dim)
    if (dim == 1) plot (1:10, 1:10, type = "n", xlab = "", ylab = "", main = paste(row.names(r)[row.names(r)==paste0(m,a)]))
    if (dim == 1) text(5,7,paste("All",round(temp[1])))
  }
}
dev.off()

compile <- function(X) {
  mean <- mean(X, na.rm = T)
  med <- median(X, na.rm = T)
  min <- min(X, na.rm = T)
  max <- max(X, na.rm = T)
  Q1 <- quantile(X, probs = c(0.25), na.rm = T)
  Q3 <- quantile(X, probs = c(0.75), na.rm = T)
  sd <- sd(X, na.rm = T)
  n.miss <- sum(as.numeric(is.na(X)))
  c(min, Q1, med, mean, Q3, max, sd, n.miss)
}

compile2 <- function(X) {
  min <- min(X, na.rm = T)
  max <- max(X, na.rm = T)
  mean <- mean(X, na.rm = T)
  sd <- sd(X, na.rm = T)
  c(min, max, mean, sd)
}

summ <- apply(r, 1, compile)
row.names(summ) <- c("min", "Q1", "med", "mean", "Q3", "max", "sd", "n.miss")

summ2 <- apply(r, 1, compile2)
row.names(summ2) <- c("min","max","mean","sd")
summ2

summ3 <- data.frame(t(summ2))
summ3$sum <- paste0("range: ",round(summ3$min),"-",round(summ3$max),", mean: ",round(summ3$mean),", sd: ",summ3$sd)
write.csv(summ3, "C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_02/compiled.csv")


# PSID1, 2, 3

r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_06/many_runs_dw_psid1.csv")
r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_06/many_runs_dw_psid2.csv")
r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_11/many_runs_dw_psid3.csv")

compile2 <- function(X) {
  min <- min(X, na.rm = T)
  max <- max(X, na.rm = T)
  mean <- mean(X, na.rm = T)
  sd <- sd(X, na.rm = T)
  sum <- ifelse(round(sd(X, na.rm = T)) == 0, as.character(round(mean(X, na.rm = T))), paste0(round(min(X, na.rm = T))," - ",round(max(X, na.rm = T)),", mean: ",round(mean(X, na.rm = T)),", sd: ",round(sd(X, na.rm = T))))
  c(min, max, mean, sd, sum)
}

summ2 <- apply(r, 2, compile2)
row.names(summ2) <- c("min","max","mean","sd","summ")


summ2 <- as.data.frame(t(summ2))
summ2$summ[round(as.numeric(as.character(summ2$sd))) == 0] <- as.character(round(as.numeric(as.character(summ2$mean))))[round(as.numeric(as.character(summ2$sd))) == 0]

write.csv(summ2, "C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_11/compiled_psid3.csv")

# Plot distributions
r <- r[,2:ncol(r)] # Drop row names variable
r2 <- gather(r, method, est, all.unadj:match.full.tmle.att.sl) # Reshape long
r2$match <- r2$an <- "" # Split match method and analysis method
for (i in match.methods) {
  r2$match[grep(i, r2$method)] <- i
  r2$an   [grep(i, r2$method)] <- rep(an.methods, each = 500)
}
r2 <- r2[!(r2$method %in% paste0(rep(c("match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."), each = 2),c("iptw.att.misp","iptw.att.sl"))),] # Drop rows corresponding to matching + IPTW, because we would never do this
r2$an <- ordered(r2$an, levels = an.methods) # Make analysis method into an ordered factor
r2 <- r2[!is.na(r2$est),]
r2$est <- round(r2$est,2)

for (i in match.methods[match.methods!="match.opt."]) {
  temp <- r2[grep(i, r2$method),]
  for (j in an.methods) {
    if (dim(table(temp$est[temp$an == j]))==1) temp <- temp[temp$an!=j,]
  }
  p <- ggplot(temp, aes(x = est, colour = an)) + geom_density(size = 1) + ggtitle(i)
  if (i == "all.") { plots <- list(p) } else { plots <- list(plots, p) }
  ggsave(paste0("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_11/dist_psid3_",i,"jpg"), height = 5, width = 8)
} 
  

temp <- r2[grep("all.", r2$method),]
p1 <- ggplot(temp, aes(x = est, colour = an)) + geom_density(size = 1) + ggtitle("All data")
  ggsave(paste0("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_11/dist_psid3_",i,"jpg"), height = 5, width = 8)




