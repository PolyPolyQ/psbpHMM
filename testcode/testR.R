rm(list=ls())
gc()

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

# getwd()
# sourceCpp("src/tester.cpp")
compileAttributes()
devtools::build()
devtools::install()
library(psbpHMM)

pisave1 = pi.z # pi from markovPSBP
pisave2 = pi.z[[1]] # pi from psbpHMM 
pisave3 = pi.z

range(pisave3[[5]][[214]]-pisave1[[214]])

