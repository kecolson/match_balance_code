# Compare

n <- 1000


##############
# DGM 1: normal covs, good overlap
##############

# v1
set.seed(123)
W1 <- runif(n, min = 0.02, max = 0.7)  
W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1)
A <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
Y <-  ifelse(A==1,Y.1,Y.0)  
m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial")
gpred1 <- predict(m, type = 'response')

# v2
generateData<- function(n) {
  W1 <- runif(n, min = 0.02, max = 0.7)  
  W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1) 
  A <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
  Y <-  ifelse(A==1,Y.1,Y.0)
  Y.star<- Y
  Y.star.1<- Y.1
  Y.star.0<- Y.0
  data.frame(W1,W2, A,Y, Y.star, Y.1, Y.0, Y.star.1, Y.star.0)
}
set.seed(123)
data <- generateData(n)
m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = data)
gpred2 <- predict(m, type = 'response')

summary(gpred1)
summary(gpred2)
# Exactly the same. Good. Range of pscores = 0.19 to 0.63

# v3
set.seed(123)
pop <- generateData(100000)
set.seed(1235)
sample <- pop[sample(row.names(pop), size = n, replace = F),]

m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = pop)
gpred3 <- predict(m, type = 'response')

m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = sample)
gpred4 <- predict(m, type = 'response')

summary(gpred1) # Same
summary(gpred2) # Same
summary(gpred3) # Much wider range-- 0.09 to 0.78
summary(gpred4) # Close -- 0.16 to 0.67



##############
# DGM 2: normal covs, medium overlap
##############

generateData<- function(n) {
  W1 <- runif(n, min = 0.02, max = 0.7)  
  W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1) 
  A <- rbinom(n, 1, plogis(-.7 + 1.8*W1 -.1*I(W1^2) + 1.7*I(W1*W2) - 1.4*W2))
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
  Y <-  ifelse(A==1,Y.1,Y.0)
  Y.star<- Y
  Y.star.1<- Y.1
  Y.star.0<- Y.0
  data.frame(W1,W2, A,Y, Y.star, Y.1, Y.0, Y.star.1, Y.star.0)
}
set.seed(123)
pop <- generateData(100000)
set.seed(1235)
sample <- pop[sample(row.names(pop), size = 1000, replace = F),]

m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = pop)
gpred1 <- predict(m, type = 'response')
m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = sample)
gpred2 <- predict(m, type = 'response')

p1 <- hist(predict(m, type="response")[sample$A==1], breaks=30)
p2 <- hist(predict(m, type="response")[sample$A==0], breaks=30)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 

summary(gpred1) # 
summary(gpred2)

# True ATT:
round(mean(pop$Y.1[pop$A==1] - pop$Y.0[pop$A==1]),3)
# Unadjusted value:
round(summary(lm(Y ~ A, data = pop))$coefficients[2,1],3)
# Proportion of treated units
round(sum(pop$A==1)/nrow(pop)*100)



##############
# DGM 3: non-normal covs, good overlap
##############

generateData<- function(n) {
  W1 <- rgamma(n, 2, 10)
  W2 <- rbinom(n, 1, 0.55+0.05*W1)
  A <- rbinom(n, 1, plogis(-.5 + .5*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
  Y <-  ifelse(A==1,Y.1,Y.0)
  Y.star<- Y
  Y.star.1<- Y.1
  Y.star.0<- Y.0
  data.frame(W1,W2, A,Y, Y.star, Y.1, Y.0, Y.star.1, Y.star.0)
}
set.seed(123)
pop <- generateData(100000)
set.seed(1235)
sample <- pop[sample(row.names(pop), size = n, replace = F),]

m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = pop)
gpred1 <- predict(m, type = 'response')
m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = sample)
gpred2 <- predict(m, type = 'response')

summary(gpred1) # 
summary(gpred2)

# True ATT:
round(mean(pop$Y.1[pop$A==1] - pop$Y.0[pop$A==1]),3)
# Unadjusted value:
round(summary(lm(Y ~ A, data = pop))$coefficients[2,1],3)
# Proportion of treated units
round(sum(pop$A==1)/nrow(pop)*100)


##############
# DGM 4: normal covs, poor overlap
##############

generateData<- function(n) {
  W1 <- runif(n, min = 0.02, max = 0.7)  
  W2 <- rnorm(n, mean = (0.2 + 0.125*W1), sd = 1) 
  A <- rbinom(n, 1, plogis(-.3 + 2*W1 -2*I(W1^2) + 2*I(W1*W2) - 2.5*W2))
  Y.0 <- rnorm(n, mean=(-.5 + 2*poly(W1,2)  - W2), sd=1)
  Y.1 <- rnorm(n, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1)
  Y <-  ifelse(A==1,Y.1,Y.0)
  Y.star<- Y
  Y.star.1<- Y.1
  Y.star.0<- Y.0
  data.frame(W1,W2, A,Y, Y.star, Y.1, Y.0, Y.star.1, Y.star.0)
}
set.seed(123)
pop <- generateData(100000)
n <- 1000
set.seed(1235)
sample <- pop[sample(row.names(pop), size = n, replace = F),]

m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = pop)
gpred1 <- predict(m, type = 'response')
m <- glm(A ~ poly(W1,2) + W2 + W1:W2, family="binomial", data = sample)
gpred2 <- predict(m, type = 'response')

p1 <- hist(predict(m, type="response")[sample$A==1], breaks=30)
p2 <- hist(predict(m, type="response")[sample$A==0], breaks=30)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 

summary(gpred1) # 
summary(gpred2)

# True ATT:
round(mean(pop$Y.1[pop$A==1] - pop$Y.0[pop$A==1]),3)
# Unadjusted value:
round(summary(lm(Y ~ A, data = pop))$coefficients[2,1],3)

# Proportion of treated units
round(sum(pop$A==1)/nrow(pop)*100)




