options(rgl.useNULL = TRUE)

# load package
library('this.path')
library('mgcv')
library('MGLM') 
library('BPST')
library('cmna')
library('LaplacesDemon')
library('Rcpp')
library('RcppArmadillo')
library('fastmatrix')
library('quantreg')
library('MultiRNG')

# working directory
setwd(this.path::here())

## load functions  
source('G.u.R')
source('ell.uh.R')
source('x1.R')
source('beta1.R')
source('tune.lambda.R') 
source('tune.lambda.qsvcm.R')
source('SQBiT.R')  
source('tpqr.plm.R')   
source('AInv.R')
source('rG.R')
source('smqsvcm_admm.wb.R') 
source('QBiT.R') 
source('stoper.R')
source('rhotau.R')
source('SQBiT_tune.R')
source('cv.pred.SQBiT.R')


# c++ matrix inversion
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }",
                  depends = "RcppArmadillo")
# c++ matrix multiplication
sourceCpp('test.cpp')


# population
pop.all <- as.matrix(read.csv('horse_tr1_n5K.csv', header = TRUE))
pop.all <- pop.all[, -1]
S <- pop.all[, c(1, 2)]
S <- as.matrix(as.data.frame(S))
# triangulation matrix
# nrow is the number of triangles
# each row is the indices of vertices 


#### Triangulation

Tr <- as.matrix(read.csv('Tr_1.csv', header = FALSE))     ########### change Tr_2.csv for Triangulation 2
# vertices of triangulation
V <- as.matrix(read.csv('V_1.csv', header = FALSE))       ########### change V_2.csv for Triangulation 2

# population
N.all <- nrow(pop.all)
# V0: vertices of a triangulation
# determine a point with (xx, yy) is inside a given triangulation
ind1 <- inVT(V0 = V, Tr0 = Tr, 
             xx = S[, 1], yy = S[, 2])
# indices of all the points which are inside the triangulation
ind1 <- ind1$ind.inside
ind2 <- (1:nrow(pop.all))[!is.na(pop.all[, 1])]
ind <- sort(intersect(ind1, ind2))
pop.r <- pop.all[ind, ]# sample size
Npop <- nrow(pop.r)


# coordinates
S.pop <- round(pop.r[, c(1, 2)], 2)
# unique covariates
u <- unique(round(pop.all[, 1], 2))
v <- unique(round(pop.all[, 2], 2))
pop.r <- cbind(pop.r, beta1(S.pop[, 1], S.pop[, 2]))


# Population basis
d <- 3
r <- 1
# bivariate spline basis function
# d = degree of piecewise polynomials
# r = smoothness parameters
# 0 <= r < d
# z = coordinates
B0.pop <- basis(V = V, Tr = Tr, 
                d = d, r = r, 
                Z = S.pop)
# matrix after QR decomposition of the smoothness matrix
Q2 <- B0.pop$Q2
B.pop <- B0.pop$B # new basis
BQ2.pop <- as.matrix(B.pop %*% Q2)
# thin-plate energy function
K <- B0.pop$K
P <- t(Q2) %*% K %*% Q2

# simulation
n <- 1000                               ########### change n = 1000 for n = 1000
nsim <- 200
set.seed(123)
seeds <- sample(1:1000, nsim)
tau <- 0.5                              ########### change tau = 0.75 for tau = 0.75

# set a seed
set.seed(123)

# simulating Y
eta <- matrix(c(1, 1, 1))
# covariates
Sigma <- matrix(0, 3, 3)
for (i in 1:3) {
  for (j in 1:3) {
    Sigma[i, j] <- 0.5 ^ (abs(i - j))
  }
}


# quantile
Finv <- qt(tau, df = 2)
colnames(pop.r) <- c('u', 'v', 'b0', 'b1') 

#### sensitivity analysis
cs <- seq(from = 1, to = 10, length.out = 10)
mses <- c()

for (c in cs) {
  
  # bandwidth for each c
  h.c <- c * tau * (1 - tau) * ((3 + dim(Q2)[2] * 2 + log(n)) / n) ^ (2/5)
  
  # Monte Carlo simulations
  res.c <- c()
  
  for (i in 1:nsim) {
    
    set.seed(seeds[i])
    ind.s <- sample(1:Npop, n)
    data <- as.matrix(pop.r[ind.s, ])
    beta0 <- data[, c('b0', 'b1')]
    
    # locations
    S <- data[, c(1, 2)]
    
    # constant covariates
    C <- draw.d.variate.uniform(no.row = n,
                                d = 3,
                                cov.mat = Sigma) * sqrt(3)
    C <- C - 0.5 * sqrt(3)
    X <- runif(n, -1, 1) 
    Y <- beta0[, 1] + X * beta0[, 2] + C %*% eta + (1 + .2 * C[, 3] + X * x1(S[, 'u'], S[, 'v'])) * rt(n = n, 2)
    X <- cbind(1, X)  
    # generate Bivariate spline basis 
    B <- B.pop[ind.s, ]
    # alpha
    true.alpha.samp <- beta0
    true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
    true.alpha.samp[, 2] <- true.alpha.samp[, 2] + x1(S[, 'u'], S[, 'v']) * Finv
    true.eta <- eta
    true.eta[3] <- eta[3] + .2 * Finv
    
    # tuning for SQBiT
    mod.SQBiT <- tune.lambda(Y = Y, X = X, C = C, max.iter = 100,
                             P = P, B = B, Q2 = Q2, var.j.2 = FALSE, kernel = 'unif',
                             h = h.c, tau = tau, lambda_start = 1e-4,  
                             eta.j1 = .32, eta.j2 = .32, nlambda = 10, new.nlambda = 10, 
                             lambda.scale = 100) 
    
    res.c <- rbind(res.c,  
                   c(norm(true.eta - mod.SQBiT$eta, '2') ^ 2 / 3,
                     norm(true.alpha.samp[, 1] - mod.SQBiT$beta[, 1], '2') ^ 2 / n,
                     norm(true.alpha.samp[, 2] - mod.SQBiT$beta[, 2], '2') ^ 2 / n))
    
    cat("iteration = ", i, "\n")
  }
  
  # mse for each c
  mses <- rbind(mses, apply(res.c, 2, mean))
  
}

###################### plot
par(mfrow = c(1, 2))
plot(cs, mses[, 1],  
     ylab = 'MSE', type = 'b', col = 'magenta',
     xlab = '$c$', main = '$\\eta_{0.50}$')
abline(v = 7, col = 'deepskyblue', lty = 2)
plot(cs, mses[, 3], 
     main = '$\\beta_{1,0.50}(\\mathbf{s})$', 
     type = 'b', col = 'indianred1',
     xlab = '$c$', ylab = 'MISE')
abline(v = 7, col = 'deepskyblue', lty = 2)


################# Estimation

# bandwidth
h <- 7 * tau * (1 - tau) * ((3 + dim(Q2)[2] * 2 + log(n)) / n) ^ (2/5)

res <- c()
res.QBiT <- c()
res.tp <- c() 
nsim <- 200
for (i in 1:nsim) {
  set.seed(seeds[i])
  ind.s <- sample(1:Npop, n)
  data <- as.matrix(pop.r[ind.s, ])
  beta0 <- data[, c('b0', 'b1')]
  
  # locations
  S <- data[, c(1, 2)]
  
  # constant covariates
  C <- draw.d.variate.uniform(no.row = n,
                              d = 3,
                              cov.mat = Sigma) * sqrt(3)
  C <- C - 0.5 * sqrt(3)
  X <- runif(n, -1, 1) 
  Y <- beta0[, 1] + X * beta0[, 2] + C %*% eta + (1 + .2 * C[, 3] + X * x1(S[, 'u'], S[, 'v'])) * rt(n = n, 2)
  X <- cbind(1, X)  
  # generate Bivariate spline basis 
  B <- B.pop[ind.s, ]
  # alpha
  true.alpha.samp <- beta0
  true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
  true.alpha.samp[, 2] <- true.alpha.samp[, 2] + x1(S[, 'u'], S[, 'v']) * Finv
  true.eta <- eta
  true.eta[3] <- eta[3] + .2 * Finv
  
  
  # tuning for QBiT
  tune.QBiT <- tune.lambda.qsvcm(y = Y, X = X, C = C, max.iter = 100,
                                 P = P, B = B, Q2 = Q2, var.j.2 = FALSE,
                                 tau = tau, 
                                 lambda_start = 1e-4, 
                                 lambda_end = 1e1, 
                                 nlambda = 10, new.nlambda = 10, 
                                 lambda.scale = 100) 
  
  # tuning for smoothed QSVCM
  tune.SQBiT <- tune.lambda(Y = Y, X = X, C = C, max.iter = 100,
                            P = P, B = B, Q2 = Q2, var.j.2 = FALSE,
                            h = h, tau = tau, 
                            lambda_start = 1e-4, 
                            lambda_end = 1e1, 
                            zeta = 14,
                            eta.j1 = .32, eta.j2 = .32, 
                            nlambda = 10, new.nlambda = 10, 
                            lambda.scale = 100, interval = FALSE, kernel = 'unif') 
  
  
  # tensor product QR
  tpqr.mod <- tpqr.plm(Y = Y, C = C, X = X, S = S, d = 3, tau = tau) 
  
  
  # output 
  res <- rbind(res, 
               c(tune.SQBiT$result_time,
                 mean(tune.SQBiT$iter),
                 norm(true.eta - tune.SQBiT$eta, '2') ^ 2 / 3,
                 norm(true.alpha.samp[, 2] - tune.SQBiT$beta[, 2], '2') ^ 2 / n))
  
  res.QBiT <- rbind(res.QBiT, 
                    c(tune.QBiT$time,
                      mean(tune.QBiT$iter),
                      norm(true.eta - tune.QBiT$eta, '2') ^ 2 / 3,
                      norm(true.alpha.samp[, 2] - tune.QBiT$beta[, 2], '2') ^ 2 / n ))
  
  res.tp <- rbind(res.tp, 
                  c(tpqr.mod$time,
                    0,
                    norm(true.eta - tpqr.mod$eta, '2') ^ 2 / 3,
                    norm(true.alpha.samp[, 2] - tpqr.mod$beta[, 2], '2') ^ 2 / n ))
  
  cat("iteration = ", i, "\n")
}

tab <- c()
for (j in 1:4) {
  if (j == 2) { 
    tab.1 <- paste0(round(mean(res[, j]), 2), ' (', round(sd(res[, j]) / sqrt(nsim), 2), ')')
    tab.2 <- paste0(round(mean(res.QBiT[, j]), 2), ' (', round(sd(res.QBiT[, j]) / sqrt(nsim), 2), ')')
    tab.3 <- paste0(round(mean(res.tp[, j]), 2), ' (', round(sd(res.tp[, j]) / sqrt(nsim), 2), ')') 
  } else {
    tab.1 <- paste0(round(mean(res[, j]) * 10 ^ 2, 2), ' (', round(sd(res[, j] * 10 ^ 2) / sqrt(nsim), 2), ')')
    tab.2 <- paste0(round(mean(res.QBiT[, j]) * 10 ^ 2, 2), ' (', round(sd(res.QBiT[, j] * 10 ^ 2) / sqrt(nsim), 2), ')')
    tab.3 <- paste0(round(mean(res.tp[, j]) * 10 ^ 2, 2), ' (', round(sd(res.tp[, j] * 10 ^ 2) / sqrt(nsim), 2), ')') 
  }
  tab <- cbind(tab, c(tab.1, tab.2, tab.3))
}
require('xtable')
tab <- cbind(c('SQBiT', 'QBiT', 'TPS'), tab)
tab <- as.data.frame(tab)
xtable(tab)



######################################### Coverage (5-fold CV)

## population varying coefficient
BQ2.pop <- as.matrix(B.pop %*% Q2)
beta.pop <- pop.r[, c('b0', 'b1')]
beta.pop[, 1] <- beta.pop[, 1] + Finv
beta.pop[, 2] <- beta.pop[, 2]

## coverage for VC 
all_all_Beta2 <- c() 

## coverage for constant
covs1 <- c()
lengths1 <- c() 
covs2 <- c()
lengths2 <- c()  

for (i in 1:nsim) {
  
  set.seed(seeds[i]) 
  ind.s <- sample(1:Npop, n)
  data <- as.matrix(pop.r[ind.s, ])
  beta0 <- data[, c('b0', 'b1')]
  
  # locations
  S <- data[, c(1, 2)]
  
  # constant covariates
  C <- draw.d.variate.uniform(no.row = n,
                              d = 3,
                              cov.mat = Sigma) * sqrt(3)
  C <- C - 0.5 * sqrt(3)
  X <- runif(n, -1, 1) 
  Y <- beta0[, 1] + X * beta0[, 2] + C %*% eta + (1 + .2 * C[, 3] + X * x1(S[, 'u'], S[, 'v'])) * rt(n = n, 2)
  X <- cbind(1, X)  
  # generate Bivariate spline basis 
  B <- B.pop[ind.s, ]
  # alpha
  true.alpha.samp <- beta0
  true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
  true.alpha.samp[, 2] <- true.alpha.samp[, 2] + x1(S[, 'u'], S[, 'v']) * Finv
  true.eta <- eta
  true.eta[3] <- eta[3] + .2 * Finv
  
  
  ####### tuning for SQBiT (5-fold CV)
  lambda_SQBiT <- SQBiT_tune(Y = Y, X = X, C = C,  
                             P = P, B = B, Q2 = Q2,  
                             h = h, tau = tau, 
                             lambda_start = 1e-4, lambda_end = 1e1,
                             nlambda = 10, new.nlambda = 10, 
                             lambda.scale = 100)  
  
  ###### Point Estimation
  mod_SQBiT <- SQBiT(Y = Y, X = X, C = C, 
                     P = P, B = B, Q2 = Q2, tau = tau,  
                     lambda = lambda_SQBiT$lambda, h = h, adaptive.h = FALSE,
                     gacv.compute = FALSE, interval = TRUE) 
  
  
  ##### asymptotic interval
  cis1 <- mod_SQBiT$cis
  
  ##### wild bootstrap
  mod_SQBiT.wb <- smqsvcm_admm.wb(h = h, tau = tau, Y = Y, C = C, X = X, P = P, B = B, Q2 = Q2, 
                                  eta.j = .32, var.j = FALSE, biascorr = TRUE,
                                  lambda = lambda_SQBiT$lambda, Br = 500, 
                                  compute.vc = TRUE, BQ2.eva = BQ2.pop, 
                                  gamma.hat = c(mod_SQBiT$eta, mod_SQBiT$gamma))
  cis2 <- mod_SQBiT.wb$cis
  
  # varying coefficient
  vc.lb <- mod_SQBiT.wb$betas.lb
  vc.ub <- mod_SQBiT.wb$betas.ub 
  
  all_Beta2 <- (vc.lb[, 2] < beta.pop[, 2]) & (vc.ub[,2] > beta.pop[, 2])
  all_all_Beta2 <- cbind(all_all_Beta2, all_Beta2)
  
  # constant
  cov1 <- c()  
  length1 <- c()
  cov2 <- c()  
  length2 <- c() 
  for (j in 1:3) { 
    cov1 <- cbind(cov1, ifelse(true.eta[j] >= cis1[j, 1] & true.eta[j] <= cis1[j, 2], 1, 0)) 
    cov2 <- cbind(cov2, ifelse(true.eta[j] >= cis2[j, 1] & true.eta[j] <= cis2[j, 2], 1, 0)) 
    length1 <- c(length1, cis1[j, 2] - cis1[j, 1]) 
    length2 <- c(length2, cis2[j, 2] - cis2[j, 1])  
  } 
  covs1 <- rbind(covs1, cov1) 
  lengths1 <- rbind(lengths1, length1) 
  covs2 <- rbind(covs2, cov2) 
  lengths2 <- rbind(lengths2, length2)  
  cat("iteration = ", i, "\n")
} 



## constant
c(apply(covs1, 2, mean) * 100, apply(lengths1, 2, mean) * 100) 
c(apply(covs1, 2, sd) * 100, apply(lengths1, 2, sd) * 100) / sqrt(nsim) 

c(apply(covs2, 2, mean) * 100, apply(lengths2, 2, mean) * 100) 
c(apply(covs2, 2, sd) * 100, apply(lengths2, 2, sd) * 100) / sqrt(nsim) 


## varying coefficient
mean(apply(all_all_Beta2, 1, mean)) 

######################################### Coverage (GACV)

## population varying coefficient
BQ2.pop <- as.matrix(B.pop %*% Q2)
beta.pop <- pop.r[, c('b0', 'b1')]
beta.pop[, 1] <- beta.pop[, 1] + Finv
beta.pop[, 2] <- beta.pop[, 2]

## coverage for VC
all_all_Beta2 <- c() 

## coverage for constant
covs1 <- c()
lengths1 <- c() 
covs2 <- c()
lengths2 <- c()  

for (i in 1:nsim) {
  
  
  set.seed(seeds[i]) 
  ind.s <- sample(1:Npop, n)
  data <- as.matrix(pop.r[ind.s, ])
  beta0 <- data[, c('b0', 'b1')]
  
  # locations
  S <- data[, c(1, 2)]
  
  # constant covariates
  C <- draw.d.variate.uniform(no.row = n,
                              d = 3,
                              cov.mat = Sigma) * sqrt(3)
  C <- C - 0.5 * sqrt(3)
  X <- runif(n, -1, 1) 
  Y <- beta0[, 1] + X * beta0[, 2] + C %*% eta + (1 + .2 * C[, 3] + X * x1(S[, 'u'], S[, 'v'])) * rt(n = n, 2)
  X <- cbind(1, X)  
  # generate Bivariate spline basis 
  B <- B.pop[ind.s, ]
  # alpha
  true.alpha.samp <- beta0
  true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
  true.alpha.samp[, 2] <- true.alpha.samp[, 2] + x1(S[, 'u'], S[, 'v']) * Finv
  true.eta <- eta
  true.eta[3] <- eta[3] + .2 * Finv
  
  ##### GACV
  mod_SQBiT <- tune.lambda(Y = Y, X = X, C = C, max.iter = 100,
                           P = P, B = B, Q2 = Q2, var.j.2 = FALSE,
                           h = h, tau = tau, lambda_start = 1e-4, zeta = 14,
                           eta.j1 = .32, eta.j2 = .32, nlambda = 10, new.nlambda = 10, 
                           lambda.scale = 100, interval = TRUE, kernel = 'unif')  
  cis1 <- mod_SQBiT$cis2
  
  ##### wild bootstrap
  mod_SQBiT.wb <- smqsvcm_admm.wb(h = h, tau = tau, Y = Y, C = C, X = X, P = P, B = B, Q2 = Q2, 
                                  eta.j = .32, var.j = FALSE, biascorr = TRUE,
                                  lambda = mod_SQBiT$lambdac, Br = 500, 
                                  compute.vc = TRUE, BQ2.eva = BQ2.pop, 
                                  gamma.hat = mod_SQBiT$gamma)
  cis2 <- mod_SQBiT.wb$cis 
  
  # varying coefficient
  vc.lb <- bs$betas.lb
  vc.ub <- bs$betas.ub 
  
  all_Beta2 <- (vc.lb[, 2] < beta.pop[, 2]) & (vc.ub[,2] > beta.pop[, 2])
  all_all_Beta2 <- cbind(all_all_Beta2, all_Beta2) 
  
  # constant 
  cov1 <- c()  
  length1 <- c()
  cov2 <- c()  
  length2 <- c() 
  for (j in 1:3) { 
    cov1 <- cbind(cov1, ifelse(true.eta[j] >= cis1[j, 1] & true.eta[j] <= cis1[j, 2], 1, 0)) 
    cov2 <- cbind(cov2, ifelse(true.eta[j] >= cis2[j, 1] & true.eta[j] <= cis2[j, 2], 1, 0))  
    length1 <- c(length1, cis1[j, 2] - cis1[j, 1]) 
    length2 <- c(length2, cis2[j, 2] - cis2[j, 1])  
  } 
  covs1 <- rbind(covs1, cov1) 
  lengths1 <- rbind(lengths1, length1) 
  covs2 <- rbind(covs2, cov2) 
  lengths2 <- rbind(lengths2, length2)  
  
  
  cat("iteration = ", i, "\n")
} 


## constant
c(apply(covs1, 2, mean) * 100, apply(lengths1, 2, mean) * 100) 
c(apply(covs1, 2, sd) * 100, apply(lengths1, 2, sd) * 100) / sqrt(nsim) 

c(apply(covs2, 2, mean) * 100, apply(lengths2, 2, mean) * 100) 
c(apply(covs2, 2, sd) * 100, apply(lengths2, 2, sd) * 100) / sqrt(nsim) 


## varying coefficient
mean(apply(all_all_Beta2, 1, mean)) 


############################################
############################################ heatmap
require("akima")
require("fields")
require('tikzDevice')
require('randomcoloR')

## plotting
options(tikzLatexPackages 
        = c(getOption( "tikzLatexPackages" ),
            "\\usepackage{amsmath,amsfonts,amsthm, palatino, mathpazo}"))

# boundary of horseshoe 
boundary <- read.csv('bb_horse.csv')

# prediction
p <- ncol(X)
q <- ncol(C)
Jn <- dim(B.pop %*% Q2)[2]
# population prediction

# tps basis
knots1 <- tpqr.mod$knots[[1]]
knots2 <- tpqr.mod$knots[[2]]
B1 <- bs(S.pop[, 1], df = length(knots1), degree = d)
B2 <- bs(S.pop[, 2], df = length(knots2), degree = d)
basis <- do.call(rbind, lapply(1:nrow(S.pop), function(i) kronecker(B1[i,], B2[i,])))
# X * basis
XB <- c()
for (j in 1:p) {
  XB <- cbind(XB, t(sapply(1:n, function(i) X[i, j] * basis[i, ])))
}
jn <- ncol(basis)


# VC estimates
mhat <- c()
mhat.sm <- c()
mhat.tp <- c() 
for(i in 1:p){
  mhat <- cbind(mhat, BQ2.pop %*% tune.qsvcm$gamma[(q + 1 + (i - 1) * Jn):(q + Jn * i)])
  mhat.sm <- cbind(mhat.sm, BQ2.pop %*% tune.sqsvcm$gamma[(q + 1 + (i - 1) * Jn):(q + Jn * i)])
  mhat.tp <- cbind(mhat.tp, basis %*% tpqr.mod$gamma[(jn * (i - 1) + 1):(jn * i)]) 
}

# horseshoe plotting
heat.2d.hs <- function(figData,  zlim, col, outside.color='white'){
  # figure data
  newz.outside <- zlim[2] + 2 * (zlim[2] - zlim[1]) / length(col) 
  # color 
  figData$z[which(figData$z < zlim[1] | figData$z > zlim[2])] <- newz.outside 
  zlim[2] <- zlim[2] + 2 * (zlim[2] - zlim[1]) / length(col) 
  col <- c(col, outside.color)
  # plot
  par(mar = c(2.1, 2.1, 2.1, 2.1))
  imagePlot(figData, frame = FALSE, axes = FALSE, 
            zlim = zlim, col = col)
}

# true coefficient beta0
tikz(paste0('plm_t2_heter_n', n, '_truth_beta0_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = pop.r[, 1], y = pop.r[, 2], duplicate = "mean",
              z = pop.r[, 3], nx = 140, ny = 140)
# heatmap
heat.2d.hs(int, zlim = c(-5, 5.5), 
           col = tim.colors(64),
           outside.color = 'black')
# horseshoe polygon
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# true coefficient beta1
tikz(paste0('plm_t2_heter_n', n, '_truth_beta1_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = pop.r[, 1], y = pop.r[, 2], duplicate = "mean",
              z = pop.r[, 4], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-1, 5.3), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# qplvcm beta0
tikz(paste0('plm_t2_heter_n', n, '_QBiT_beta0_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
              z = mhat[, 1], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-5, 5.5), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# qplvcm beta1
tikz(paste0('plm_t2_heter_n', n, '_QBiT_beta1_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
              z = mhat[, 2], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-1, 5.3), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# sqplvcm beta0
tikz(paste0('plm_t2_heter_n', n, '_SQBiT_beta0_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
              z = mhat.sm[, 1], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-5, 5.5), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# sqplvcm beta1
tikz(paste0('plm_t2_heter_n', n, '_SQBiT_beta1_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
              z = mhat.sm[, 2], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-1, 5.3), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# tps beta0
tikz(paste0('plm_t2_heter_n', n, '_tps_beta0_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
cols <- c(tim.colors(n = 64), 'black')

int <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
              z = mhat.tp[, 1], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-5, 5.5), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()

# tps beta1
tikz(paste0('plm_t2_heter_n', n, '_tps_beta1_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
              z = mhat.tp[, 2], nx = 140, ny = 140)
heat.2d.hs(int, zlim = c(-1, 5.3), 
           col = tim.colors(64),
           outside.color = 'black')
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()