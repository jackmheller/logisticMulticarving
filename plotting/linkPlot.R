
require(MASS)
require(glmnet)
require(Matrix)
require(tictoc)
require(hdi)
require(selectiveInference)
require(doSNOW)
require(parallel)
require(doRNG)
require(tmg)
require(git2r)

# toeplitz
n <- 100
p <- 200
rho <- 0
level<-0.05 #17/02/23 VK, setting significance level only once
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 2
sparsity <- length(sel.index) # 17/02/23 VK, changed so that value automatically updates
set.seed(42) # to make different methods comparable, fix the x-matrix
x <- mvrnorm(n, rep(0, p), Cov)
print (x[1,1])
# should create the right x on D-MATH server, x[1 ,1] = 0.958

xb <- x %*% beta
p.true.logit <- exp(xb) / (1 + exp(xb))

p.true.cloglog <- 1 - exp(-exp(xb))

plot(xb, p.true.logit, xlim=range(xb), ylim=range(p.true.logit), xlab=expression(paste(X^T, beta)), ylab="Probability of 1", 
     main = "Link Function Probability Differences",pch=16, col = 'red')

lines(xb[order(xb)], p.true.logit[order(xb)], xlim=range(xb), ylim=range(p.true.logit), pch=16, col = 'red')
# do cloglog
lines(xb[order(xb)], p.true.cloglog[order(xb)], xlim=range(xb), ylim=range(p.true.cloglog), pch=16, col = 'blue')
points(xb[order(xb)], p.true.cloglog[order(xb)], xlim=range(xb), ylim=range(p.true.cloglog), pch=16, col = 'blue')

# do beta binomial
require('VGAM')
p.true.bbinom <-rbetabinom(n, 1, p.true)

lines(xb[order(xb)], p.true.bbinom[order(xb)], xlim=range(xb), ylim=range(p.true.bbinom), pch=16, col = 'black')
points(xb[order(xb)], p.true.bbinom[order(xb)], xlim=range(xb), ylim=range(p.true.bbinom), pch=16, col = 'black')

legend('topleft', legend=c('Logit', 'Cloglog'), col = c('red', 'blue'), lwd = 3, ncol = 1)


