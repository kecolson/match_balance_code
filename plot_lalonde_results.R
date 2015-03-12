# Plot lalonde results
# 1/27/15
# Ellie Colson

rm(list =ls())
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library("lattice")
n.analysis <- 10
n.match <- 7
colors <- brewer.pal(10,"Paired")

setwd("C:/Users/kecolson/Google Drive/simulation/questions and tests/match_balance_and_mse/results/lalonde/2015_02_20/")

# Bring in data
dw_dw <- read.csv("lalonde_est_dw_dw.csv")
dw_dw$dataset <- "exp"
dw_psid1 <- read.csv("lalonde_est_dw_psid1.csv")
dw_psid1$dataset <- "psid1"
dw_psid2 <- read.csv("lalonde_est_dw_psid2.csv")
dw_psid2$dataset <- "psid2"
dw_psid3 <- read.csv("lalonde_est_dw_psid3.csv")
dw_psid3$dataset <- "psid3"
dw_cps1 <- read.csv("lalonde_est_dw_cps1.csv")
dw_cps1$dataset <- "cps1"
dw_cps2 <- read.csv("lalonde_est_dw_cps2.csv")
dw_cps2$dataset <- "cps2"
dw_cps3 <- read.csv("lalonde_est_dw_cps3.csv")
dw_cps3$dataset <- "cps3"

# Plot estimates -----------------------

data <- rbind(dw_dw, dw_psid1, dw_psid2, dw_psid3, dw_cps1, dw_cps2, dw_cps3)
names(data) <- c("param","est","dataset")
data$match <- as.factor(rep(rep(c("all","nn","opt","sl","gen","subclass","full"), each = n.match), n.analysis))
data$an <- as.factor(rep(c("unadj","iptw.misp","iptw.sl","gcomp.misp","gcomp.sl","lalonde1","lalonde2","lalonde3","tmle.misp","tmle.sl"), n.match*n.analysis))

data$x <- as.factor(data$dataset)
data$y <- data$est

# drop certain analyses because they screws up the axis
data <- data[data$an!="unadj" & data$an!="lalonde2",]

# Reorder factor levels of an and match
data$an <- factor(data$an, levels(data$an)[c(7,1:6)])
data$match <- factor(data$match, levels(data$match)[rev(c(1,4,5,6,3,7,2))])

# Make tiled forest plot-- 
   # Grid where x boxes are analysis methods, y boxes are matching methods and points are estimates according to differnet datasets

xyplot(y ~ x | an + match, data=data, xlab = "Analysis method; dataset", ylab = "Matching method; estimate", scales = list(rot=90))



## Plot difference between experimental estimate and observational estimates ---------------

data2 <- rbind(dw_psid1, dw_psid2, dw_psid3, dw_cps1, dw_cps2, dw_cps3)
names(data2) <- c("param","est","dataset")
data2$match <- as.factor(rep(rep(c("all","nn","opt","sl","gen","subclass","full"), each = n.analysis), 6))
data2$an <- as.factor(rep(c("unadj","iptw.misp","iptw.sl","gcomp.misp","gcomp.sl","lalonde1","lalonde2","lalonde3","tmle.misp","tmle.sl"), n.match*6))

# Use as gold standard the estimate for the corresponding analysis in the experimental data
#data2$gold1 <- dw_dw$V1
#data2$est_diff <- data2$est - dw_dw$V1
#data2$est <- data2$est - dw_dw$V1

# Use as gold standard the unadjusted experimental estimate
#data2$gold2 <- dw_dw$V1[1]
#data2$est_diff2 <- data2$est - dw_dw$V1[1]
data2$est <- data2$est - dw_dw$V1[1]

data2$x <- as.factor(data2$dataset)
data2$y <- data2$est

# drop certain analyses because they screws up the axis
data2 <- data2[!is.na(data2$y) & abs(data2$y)<2000,]

# Drop IPTW for matched analyses
data2 <- data2[!((data2$an == "iptw.misp" | data2$an == "iptw.sl") & data2$match!="all"),]

# Reorder factor levels of an and match, then sort data
data2$an <- factor(data2$an, levels(data2$an)[c(10,3,4,1,2,5:9)]) 
data2$match <- factor(data2$match, levels(data2$match)[rev(c(1,4,5,6,3,7,2))])
data2 <- data2[order(data2$an, data2$match, data2$x),]


xyplot(y ~ x | match, groups = an, data = data2, xlab = "Analysis method; dataset", ylab = "Matching method; est difference", 
       key=list(text=list(unique(as.character(data2$an))), points=list(pch = rep(17,10), col=colors), columns=2),
       col = colors, pch = 17, type="b", lwd=3,
       scales = list(rot=90),
       layout = c(4,2),
       panel = function(x, y, ...) {
         panel.abline( h=0, lty = "dotted", col = "black")
         panel.abline( h=400, lty = "dotted", col = "black")
         panel.abline( h=-400, lty = "dotted", col = "black")
         panel.xyplot(x, y, ...)
       })

# 1 plot for each analysis and each matching method
# xyplot(y ~ x | an + match, data=data2, xlab = "Analysis method; dataset", ylab = "Matching method; est difference", 
#        scales = list(rot=90),
#        panel = function( x,y,...) {
#          panel.abline( h=0, lty = "dotted", col = "black")
#          panel.xyplot( x,y,...)
#        })

# Example code to change colors
# xyplot( result ~ time | location, data=dataset,groups=genotype,
#         fill.color = as.character(dataset$color),
#         panel = function(x, y,fill.color,...,subscripts) {
#           fill = fill.color [subscripts]
#           panel.xyplot(x, y,pch=19, col=fill, type ="b")}
# )





# Examine experimental estimates only
exp <- dw_dw
names(exp) <- c("param","est","dataset")
exp$match <- as.factor(rep(c("all","nn","opt","sl","gen","subclass","full"), each = n.analysis))
exp$an <- as.factor(rep(c("unadj","iptw.misp","iptw.sl","gcomp.misp","gcomp.sl","lalonde1","lalonde2","lalonde3","tmle.misp","tmle.sl"), n.match))
exp <- exp[,c(2,4,5)]
exp <- exp[!((exp$an == "iptw.misp" | exp$an == "iptw.sl") & exp$match!="all"),] # Drop IPTW for matched analyses
exp$an <- factor(exp$an, levels(exp$an)[c(10,3,4,1,2,5:9)]) # reorder factors
exp$match <- factor(exp$match, levels(exp$match)[c(1,4,5,6,3,7,2)])

exp2 <- dcast(exp, an ~ match, value.var = "est")
exp2

write.csv(exp2, "dw_dw_results_summary.csv", row.names=F)

hist(exp$est, breaks = 30, xlab = "Estimate", main = "Distribution of experimental estimates \n across analyses")


