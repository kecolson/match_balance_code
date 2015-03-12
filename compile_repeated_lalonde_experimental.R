# Ellie Colson
# 3/2/2015
# Compile results from repeatedly running analysis of lalonde experimental data

rm(list = ls())
options(scipen=10)
n.analysis <- 10
n.match <- 7

# DW- DW
r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_02/many_runs_dw_dw_est.csv")
row.names(r) <- paste0(rep(c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full."), each = n.analysis),
                             c("unadj","iptw.att.misp","iptw.att.sl","gcomp.att.misp","gcomp.att.sl",
                               "lalonde1","lalonde2","lalonde3","tmle.att.misp.misp","tmle.att.sl"))
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
for (m in c("all.","match.nn.","match.opt.","match.sl.","match.gen.","match.sub.","match.full.")) {
  if (m=="all.") par(mfrow = c(3,4), mar=c(3,2,2,1))
  if (m!="all.") par(mfrow = c(3,3), mar=c(3,2,2,1))
  for (a in c("unadj","iptw.att.misp","iptw.att.sl","gcomp.att.misp","gcomp.att.sl",
                               "lalonde1","lalonde2","lalonde3","tmle.att.misp.misp","tmle.att.sl")) {
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


# PSID1, 2

r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_06/many_runs_dw_psid1.csv")
r <- read.csv("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_06/many_runs_dw_psid2.csv")

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

write.csv(summ2, "C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_03_06/compiled_psid2.csv")




  