library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

# getwd()
# sourceCpp("src/tester.cpp")
compileAttributes()
devtools::build()
devtools::install()
library(psbpHMM)

# test1 functions
timesTwo(5)
addTwo(5)
minusTwo(5)
timesTwo(pi)

# test2 functions 
a1(diag(2))
a2(c(1,2))
a8(3,4,.5)

# test3 functions 
pilist2 = updatePi2(beta = beta.k, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = 288)
pilist2

# test4 functions 
returnPi()

# test5 functions
i = 1
t = 1
k = 1
y[[i]][1,]
mu[[k]]
Sigma[[k]]
p

mvndensity(y=y[[i]][1,], mu=mu[[k]], Sigma=Sigma[[k]], d=p)
dmvn(X = y[[i]][1,], mu = mu[[k]], sigma = Sigma[[k]], log = F, 1, isChol = FALSE)

state.list[[i]][[t]]
y
mu
Sigma
log.stuff
nu.df
detR.star
t.max
K
n
p
det(Sigma[[k]])

# test 6 functions 
upZ(stateList = state.list, y = y, mu = mu, Sigma = Sigma, logStuff = log.stuff, 
    nudf = nu.df, detRstar = detR.star, piz = pi.z, u = u, tmax = t.max, K = K, n = n, d = p)

logsum <- function(x) {
  log(sum(exp(x - max(x)))) + max(x)
}

x <- c(1,2,3,4,5)
prob <- c(.2, .1, .3, .1, .3); sum(prob)
csample_num(x, size = 1, replace = FALSE, prob = prob)
csample_int(x, 1, 0, prob)

(2*pi)^(-p/2)*det(Sigma[[k]])^(-0.5)*exp(-0.5 * (y[[i]][t,] - mu[[k]]) %*% solve(Sigma[[k]]) %*% t(y[[i]][t,] - mu[[k]]))
#pow( (2*M_PI), (-d/2))*temp3*exp(-0.5 * (yit - muK).t() * inv(sigK) * (yit - muK))


