# Check data generating mechansims
# Ellie Colson
# 1/7/14

library("tmle")

# Check bias, ATE, ATT
n<-1000
W1 <- runif(n, min=0.02, max=0.7)  
W2 <- rnorm(n, mean=(0.2+0.125*W1), sd=1) 

#W1 <- sample(c(-1,0,1), n, replace = T, prob=c(0.2,0.3,0.5))
#W2 <- rgamma(n, 2+2*W1, 20) 

A <- A1 <- rbinom(n, 1, plogis(-1.3 + W1 + I(W1^2) + I(W1*W2) - W2))
#A <- A2 <- rbinom(n, 1, plogis(0.5 + I(-5*W1) + I(2*W1^2) + I(3*W1*W2) - W2))
#A <- A3 <- rbinom(n, 1, plogis(0.5 + I(-10*W1) + I(4*W1^2) + I(6*W1*W2) - I(-2*W2)))
  
Y0<-rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
Y1<-rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
Y <-  ifelse(A==1,Y1,Y0)
mean(Y1-Y0)

summary(lm(Y ~ A))

mean(Y1[A==1] - Y0[A==1])

# Check positivity
m <- glm(A ~ W1 + W2, family = 'binomial')
summary(m$fitted.values)


# Examine correctly specified TMLE result

generateData<- function(n, covars.norm=T, pscore.type = 1) {
  
  if (covars.norm == T) {
    W1 <- runif(n, min = 0.02, max = 0.7)  
    W2 <- rnorm(n, mean = (0.2 + 0.0125*W1), sd = 1) 
  } else {
    W1 <- sample(c(-1,0,1), n, replace = T, prob=c(0.2,0.3,0.5))
    W2 <- rgamma(n, 2+2*W1, 20) 
  }
  
  if (pscore.type==1) {
    A <- A1 <- rbinom(n, 1, plogis(-1.2 + W1 + I(W1^2) + I(W1*W2) - W2))
  } else if (pscore.type==2) {
    A <- A2 <- rbinom(n, 1, plogis(0.8 + I(-5*W1) + I(2*W1^2) + I(3*W1*W2) - W2))
  } else if (pscore.type==3) {
    A <- A3 <- rbinom(n, 1, plogis(1 + I(-10*W1) + I(4*W1^2) + I(6*W1*W2) - I(-2*W2)))
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
O.true<- generateData(n=1000, covars.norm=T, pscore.type=1)
ATE <- mean(O.true$Y.star.1 - O.true$Y.star.0)
ATT <- mean(O.true$Y.star.1[O.true$A==1] - O.true$Y.star.0[O.true$A==1])
ATE 
ATT
table(O.true$A)

# compare to TMLE
O <- O.true
tmle.out<- tmle(Y=O$Y.star, A=O$A, W=cbind(O$W1, O$W2, O$W1sq, O$W2sq), 
                Qform = "Y ~ W1 + I(W1^2) + W2 + A + A:W1 + A:I(W1^2) + A:W2",
                gform = "A ~ W1 + I(W1^2) + W1:W2 + W2",
                gbound=0.001)
tmle.out$estimates$ATE$psi
ATE

# Unbiased-- good.

