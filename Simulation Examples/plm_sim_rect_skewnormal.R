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
library('sn')

### directory
setwd(this.path::here())

############################# source functions
file.sources <- paste0('functions/', list.files('functions', pattern = '*.R'))
sapply(file.sources, source, .GlobalEnv)

############################# source some Rcpp functions
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }",
                  depends = "RcppArmadillo")
sourceCpp('functions/test.cpp')


# population locations
set.seed(999)
s1 <- runif(10000)
s2 <- runif(10000)
S_pop <- data.frame(s1 = s1, s2 = s2)

# triangulation matrix
# nrow is the number of triangles
# each row is the indices of vertices 


#### Triangulation

Tr <- as.matrix(read.csv('Tr_rectangle.csv', header = TRUE))     ########### change Tr_rectangle_2.csv for Triangulation 2
# vertices of triangulation
V <- as.matrix(read.csv('V_rectangle.csv', header = TRUE))       ########### change V_rectangle_2.csv for Triangulation 2

# population
N.all <- nrow(S_pop)
# V0: vertices of a triangulation
# determine a point with (xx, yy) is inside a given triangulation
ind1 <- inVT(V0 = V, Tr0 = Tr, 
             xx = S_pop[, 1], yy = S_pop[, 2])
# indices of all the points which are inside the triangulation
ind1 <- ind1$ind.inside
ind2 <- (1:nrow(S_pop))[!is.na(S_pop[, 1])]
ind <- sort(intersect(ind1, ind2))
pop.r <- S_pop[ind, ]# sample size
Npop <- nrow(pop.r)


# coordinates
S.pop <- pop.r[, c(1, 2)]
# unique covariates
u <- unique(round(S.pop[, 1], 2))
v <- unique(round(S.pop[, 2], 2))

# coefficient functions
beta1 <- function(s1, s2) sin(pi * s1 * s2)
beta2 <- function(s1, s2) (1 - (1 - 2 * s1) ^ 2) * (1 - (1 - 2 * s2) ^ 2)
beta3 <- function(s1, s2) sin(pi * s1) * cos(pi * s2)

# zeta 1
xi2 <- function(s1, s2) (s1 ^ 2 + s2 ^ 2) / 6

pop.r <- cbind(S.pop, beta1(S.pop[, 1], S.pop[, 2]),
               beta2(S.pop[, 1], S.pop[, 2]),
               beta3(S.pop[, 1], S.pop[, 2]))
colnames(pop.r)[3:5] <- c('beta0', 'beta1', 'beta2')

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
                Z = as.matrix(S.pop))
# matrix after QR decomposition of the smoothness matrix
Q2 <- B0.pop$Q2
B.pop <- B0.pop$B # new basis
BQ2.pop <- as.matrix(B.pop %*% Q2)
# thin-plate energy function
K <- B0.pop$K
P <- t(Q2) %*% K %*% Q2

# simulation
n <- 2000                               ########### change n = 3000 for n = 3000
nsim <- 200
set.seed(123)
seeds <- sample(1:1000, nsim)
tau <- 0.5                              ########### change tau = 0.75, 0.90 for tau = 0.75, 0.90

# set a seed
set.seed(666)

########################### setup for constant covariates
eta <- matrix(c(1, 1, 1))
Sigma <- matrix(0, 3, 3)
for (i in 1:3) {
  for (j in 1:3) {
    Sigma[i, j] <- 0.5 ^ (abs(i - j))
  }
}

# quantile
Finv <- qsn(tau, alpha = 5)

#### sensitivity analysis
cs <- seq(from = 1, to = 10, length.out = 10)
mses <- c()

for (c in cs) {
  
  # bandwidth for each c
  h.c <- c * tau * (1 - tau) * ((3 + dim(Q2)[2] * 3 + log(n)) / n) ^ (2/5)
  
  # Monte Carlo simulations
  res.c <- c()
  
  for (i in 1:nsim) {
    
    set.seed(seeds[i])
    ind.s <- sample(1:Npop, n)
    data <- as.matrix(pop.r[ind.s, ])
    beta0 <- data[, c('beta0', 'beta1', 'beta2')]
    
    # locations
    S <- data[, c(1, 2)]
    
    # covariates that have constant effects
    C <- draw.d.variate.uniform(no.row = n,
                                d = 3,
                                cov.mat = Sigma) * sqrt(3)
    C <- C - 0.5 * sqrt(3)
    
    # covariates that have varying effects
    X <- matrix(runif(2 * n, -1, 1), n, 2)
    
    # simulate the response
    Y <- beta0[, 1] + X[, 1] * beta0[, 2] + X[, 2] * beta0[, 3] + 
      C %*% eta + (1 + X[, 2] * xi2(S[, 's1'], S[, 's2'])) * rsn(n = n, alpha = 5)
    
    
    # generate Bivariate spline basis 
    X <- cbind(1, X)  
    B <- B.pop[ind.s, ]
    
    # varying coefficients
    true.alpha.samp <- beta0
    true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
    true.alpha.samp[, 2] <- true.alpha.samp[, 2]
    true.alpha.samp[, 3] <- true.alpha.samp[, 3] + xi2(S[, 's1'], S[, 's2']) * Finv
    
    # constant coefficients
    true.eta <- eta
    
    ########## change lambda_start = 1e-4 to lambda_start = 1e-1 when tau = 0.90
    ########## change lambda.scale = 100 to lambda.scale = 10 when tau = 0.90
    mod.SQBiT <- SQBiT_gacv(Y = Y, X = X, C = C, max.iter = 100,
                            P = P, B = B, Q2 = Q2, var.j = FALSE,
                            h = h.c, tau = tau, lambda_start = 1e-4,  
                            eta.j1 = .32, eta.j2 = .32, nlambda = 10, new.nlambda = 10, 
                            lambda.scale = 100) 
    
    res.c <- rbind(res.c,  
                   c(norm(true.eta - mod.SQBiT$eta, '2') ^ 2 / 3,
                     norm(true.alpha.samp[, 1] - mod.SQBiT$beta[, 1], '2') ^ 2 / n,
                     norm(true.alpha.samp[, 2] - mod.SQBiT$beta[, 2], '2') ^ 2 / n,
                     norm(true.alpha.samp[, 3] - mod.SQBiT$beta[, 3], '2') ^ 2 / n))
    
    cat("iteration = ", i, "\n")
  }
  
  # mse for each c
  mses <- rbind(mses, apply(res.c, 2, mean))
  
}

################# Estimation (Table 1 in Supplementary Materials)

# bandwidth
h <- 7 * tau * (1 - tau) * ((3 + dim(Q2)[2] * 3 + log(n)) / n) ^ (2/5)

res <- c()
res.QBiT <- c()
res.tp <- c() 
nsim <- 200
for (i in 1:nsim) {
  set.seed(seeds[i])
  ind.s <- sample(1:Npop, n)
  data <- as.matrix(pop.r[ind.s, ])
  beta0 <- data[, c('beta0', 'beta1', 'beta2')]
  
  # locations
  S <- data[, c(1, 2)]
  
  # covariates that have constant effects
  C <- draw.d.variate.uniform(no.row = n,
                              d = 3,
                              cov.mat = Sigma) * sqrt(3)
  C <- C - 0.5 * sqrt(3)
  
  # covariates that have varying effects
  X <- matrix(runif(2 * n, -1, 1), n, 2)
  
  # simulate the response
  Y <- beta0[, 1] + X[, 1] * beta0[, 2] + X[, 2] * beta0[, 3] + 
    C %*% eta + (1 + X[, 2] * xi2(S[, 's1'], S[, 's2'])) * rsn(n = n, alpha = 5)
  
  
  # generate Bivariate spline basis 
  X <- cbind(1, X)  
  B <- B.pop[ind.s, ]
  
  # varying coefficients
  true.alpha.samp <- beta0
  true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
  true.alpha.samp[, 2] <- true.alpha.samp[, 2]
  true.alpha.samp[, 3] <- true.alpha.samp[, 3] + xi2(S[, 's1'], S[, 's2']) * Finv
  
  # constant coefficients
  true.eta <- eta
  
  
  ########## tuning for QBiT
  ########## change lambda_start = 1e-4 to lambda_start = 1e-1 when tau = 0.90
  ########## change lambda.scale = 100 to lambda.scale = 10 when tau = 0.90
  tune.QBiT <- tune.lambda.qsvcm(y = Y, X = X, C = C, max.iter = 100,
                                 P = P, B = B, Q2 = Q2, var.j.2 = FALSE,
                                 tau = tau, 
                                 lambda_start = 1e-4, 
                                 lambda_end = 1e1, 
                                 nlambda = 10, new.nlambda = 10, 
                                 lambda.scale = 100) 
  
  ########## tuning for SQBiT
  ########## change lambda_start = 1e-4 to lambda_start = 1e-1 when tau = 0.90
  ########## change lambda.scale = 100 to lambda.scale = 10 when tau = 0.90
  tune.SQBiT <- SQBiT_gacv(Y = Y, X = X, C = C, max.iter = 100,
                           P = P, B = B, Q2 = Q2, var.j = FALSE,
                           h = h, tau = tau, 
                           lambda_start = 1e-4, 
                           lambda_end = 1e1, 
                           zeta = 10, 
                           eta.j1 = .32, eta.j2 = .32, 
                           nlambda = 10, new.nlambda = 10, 
                           lambda.scale = 100, interval = FALSE)  
  
  
  # tensor product QR
  tpqr.mod <- tpqr.plm(Y = Y, C = C, X = X, S = S, d = 3, tau = tau) 
  
  
  # output 
  res <- rbind(res, 
               c(tune.SQBiT$result_time,
                 mean(tune.SQBiT$iter),
                 norm(true.eta - tune.SQBiT$eta, '2') ^ 2 / 3,
                 norm(true.alpha.samp[, 2] - tune.SQBiT$beta[, 2], '2') ^ 2 / n,
                 norm(true.alpha.samp[, 3] - tune.SQBiT$beta[, 3], '2') ^ 2 / n))
  
  res.QBiT <- rbind(res.QBiT, 
                    c(tune.QBiT$time,
                      mean(tune.QBiT$iter),
                      norm(true.eta - tune.QBiT$eta, '2') ^ 2 / 3,
                      norm(true.alpha.samp[, 2] - tune.QBiT$beta[, 2], '2') ^ 2 / n,
                      norm(true.alpha.samp[, 3] - tune.QBiT$beta[, 3], '2') ^ 2 / n ))
  
  res.tp <- rbind(res.tp, 
                  c(tpqr.mod$time,
                    0,
                    norm(true.eta - tpqr.mod$eta, '2') ^ 2 / 3,
                    norm(true.alpha.samp[, 2] - tpqr.mod$beta[, 2], '2') ^ 2 / n,
                    norm(true.alpha.samp[, 3] - tpqr.mod$beta[, 3], '2') ^ 2 / n))
  
  cat("iteration = ", i, "\n")
}

tab <- c()
for (j in 1:5) {
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



######################################### Coverage (Table 2 in Supplementary Materials)

# bandwidth
h <- 7 * tau * (1 - tau) * ((3 + dim(Q2)[2] * 3 + log(n)) / n) ^ (2/5)

## population varying coefficient
BQ2.pop <- as.matrix(B.pop %*% Q2)
beta.pop <- pop.r[, c('beta0', 'beta1', 'beta2')]
beta.pop[, 1] <- beta.pop[, 1] + Finv
beta.pop[, 2] <- beta.pop[, 2]
beta.pop[, 3] <- beta.pop[, 3] + xi2(pop.r[, 's1'], pop.r[, 's2']) * Finv

## coverage for VC 
all_all_Beta2 <- c() 
all_all_Beta3 <- c() 

## coverage for constant
covs1 <- c()
lengths1 <- c() 
covs2 <- c()
lengths2 <- c()  

for (i in 1:nsim) {
  
  set.seed(seeds[i])
  ind.s <- sample(1:Npop, n)
  data <- as.matrix(pop.r[ind.s, ])
  beta0 <- data[, c('beta0', 'beta1', 'beta2')]
  
  # locations
  S <- data[, c(1, 2)]
  
  # covariates that have constant effects
  C <- draw.d.variate.uniform(no.row = n,
                              d = 3,
                              cov.mat = Sigma) * sqrt(3)
  C <- C - 0.5 * sqrt(3)
  
  # covariates that have varying effects
  X <- matrix(runif(2 * n, -1, 1), n, 2)
  
  # simulate the response
  Y <- beta0[, 1] + X[, 1] * beta0[, 2] + X[, 2] * beta0[, 3] + 
    C %*% eta + (1 + X[, 2] * xi2(S[, 's1'], S[, 's2'])) * rsn(n = n, alpha = 5)
  
  
  # generate Bivariate spline basis 
  X <- cbind(1, X)  
  B <- B.pop[ind.s, ]
  
  # varying coefficients
  true.alpha.samp <- beta0
  true.alpha.samp[, 1] <- true.alpha.samp[, 1] + Finv
  true.alpha.samp[, 2] <- true.alpha.samp[, 2]
  true.alpha.samp[, 3] <- true.alpha.samp[, 3] + xi2(S[, 's1'], S[, 's2']) * Finv
  
  # constant coefficients
  true.eta <- eta
  
  
  ####### tuning for SQBiT (5-fold CV)
  lambda_SQBiT <- SQBiT_cv(Y = Y, X = X, C = C,  
                           P = P, B = B, Q2 = Q2,  
                           h = h, tau = tau, eta.j1 = .32, eta.j2 = .32,
                           lambda_start = 1e-4, lambda_end = 1e1,
                           nlambda = 10, new.nlambda = 10, 
                           lambda.scale = 100)  
  
  ###### Point Estimation
  mod_SQBiT <- tryCatch({
    SQBiT(Y = Y, X = X, C = C, 
          P = P, B = B, Q2 = Q2, tau = tau,  
          lambda = lambda_SQBiT$lambda, h = h, adaptive.h = FALSE,
          gacv.compute = FALSE, interval = TRUE) 
  }, error = function(e) {
    SQBiT(Y = Y, X = X, C = C, 
          P = P, B = B, Q2 = Q2, tau = tau,  
          lambda = lambda_SQBiT$lambda, h = h, adaptive.h = FALSE,
          gacv.compute = FALSE, interval = FALSE) 
  })
  
  ##### asymptotic interval
  cis1 <- mod_SQBiT$cis
  
  ##### wild bootstrap
  mod_SQBiT.wb <- SQBiT_wb(h = h, tau = tau, Y = Y, C = C, X = X, P = P, B = B, Q2 = Q2, 
                           eta.j = .32, var.j = FALSE, biascorr = TRUE,
                           lambda = lambda_SQBiT$lambda, Br = 500, 
                           compute.vc = TRUE, BQ2.eva = BQ2.pop, 
                           eta.hat = mod_SQBiT$eta, gamma.hat = mod_SQBiT$gamma)
  cis2 <- mod_SQBiT.wb$cis
  
  # varying coefficient
  vc.lb <- mod_SQBiT.wb$betas.lb
  vc.ub <- mod_SQBiT.wb$betas.ub 
  
  all_Beta2 <- (vc.lb[, 2] < beta.pop[, 2]) & (vc.ub[,2] > beta.pop[, 2])
  all_all_Beta2 <- cbind(all_all_Beta2, all_Beta2)
  
  all_Beta3 <- (vc.lb[, 3] < beta.pop[, 3]) & (vc.ub[,3] > beta.pop[, 3])
  all_all_Beta3 <- cbind(all_all_Beta3, all_Beta3)
  
  # constant
  cov1 <- c()  
  length1 <- c()
  cov2 <- c()  
  length2 <- c() 
  for (j in 1:3) { 
    cov1 <- cbind(cov1, ifelse(eta[j] >= cis1[j, 1] & eta[j] <= cis1[j, 2], 1, 0)) 
    cov2 <- cbind(cov2, ifelse(eta[j] >= cis2[j, 1] & eta[j] <= cis2[j, 2], 1, 0)) 
    length1 <- c(length1, cis1[j, 2] - cis1[j, 1]) 
    length2 <- c(length2, cis2[j, 2] - cis2[j, 1])  
  } 
  if (length(cov1) > 0) {
    covs1 <- rbind(covs1, cov1) 
  }
  lengths1 <- rbind(lengths1, length1) 
  covs2 <- rbind(covs2, cov2) 
  lengths2 <- rbind(lengths2, length2)  
  cat("iteration = ", i, "\n")
} 



## constant
c(apply(covs1, 2, function(x) mean(x, na.rm = TRUE)) * 100, 
  apply(lengths1, 2, function(x) mean(x, na.rm = TRUE)) * 100) 
c(apply(covs1, 2, function(x) sd(x, na.rm = TRUE)) * 100, 
  apply(lengths1, 2, function(x) sd(x, na.rm = TRUE)) * 100) / sqrt(nsim) 

c(apply(covs2, 2, mean) * 100, apply(lengths2, 2, mean) * 100) 
c(apply(covs2, 2, sd) * 100, apply(lengths2, 2, sd) * 100) / sqrt(nsim) 


## varying coefficient
mean(apply(all_all_Beta2, 1, mean)) 
mean(apply(all_all_Beta3, 1, mean)) 

############################################
############################################ heatmap
library("akima")
library("fields")
library('tikzDevice')
library('randomcoloR')
library('ggplot2')

## plotting
options(tikzLatexPackages 
        = c(getOption( "tikzLatexPackages" ),
            "\\usepackage{amsmath,amsfonts,amsthm, bm, palatino, mathpazo}"))


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
  mhat <- cbind(mhat, BQ2.pop %*% tune.QBiT$gamma[(q + 1 + (i - 1) * Jn):(q + Jn * i)])
  mhat.sm <- cbind(mhat.sm, BQ2.pop %*% tune.SQBiT$gamma[(q + 1 + (i - 1) * Jn):(q + Jn * i)])
  mhat.tp <- cbind(mhat.tp, basis %*% tpqr.mod$gamma[(jn * (i - 1) + 1):(jn * i)]) 
}

# horseshoe plotting
heat.2d.hs <- function(figData,  zlim, col, outside.color='white', legend = TRUE){
  # figure data
  newz.outside <- zlim[2] + 2 * (zlim[2] - zlim[1]) / length(col) 
  # color 
  figData$z[which(figData$z < zlim[1] | figData$z > zlim[2])] <- newz.outside 
  zlim[2] <- zlim[2] + 2 * (zlim[2] - zlim[1]) / length(col) 
  col <- c(col, outside.color)
  # plot
  par(mar = c(2.1, 2.1, 2.1, 2.1))
  if (legend == TRUE) {
    imagePlot(figData, frame = FALSE, axes = FALSE,
              zlim = zlim, col = col)
  } else {
    image(figData, frame = FALSE, axes = FALSE,
          zlim = zlim, col = col)
  }
}

# true coefficient beta0
beta1 <- function(s1, s2) sin(pi * s1 * s2)
beta2 <- function(s1, s2) (1 - (1 - 2 * s1) ^ 2) * (1 - (1 - 2 * s2) ^ 2)
beta3 <- function(s1, s2) sin(2 * pi * s1) * sin(2 * pi * s2)




################################################## plots for estimates of beta1

################################# True coefficient 
tikz('beta1.tex', width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = beta2(S.pop[, 1], S.pop[, 2]), nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(0, 1), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_1(\\boldsymbol{s})$",
    fill = "$\\beta_1(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()

################################# SQBiT

tikz(paste0('beta1_SQBiT_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = mhat.sm[, 2], nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(0, 1), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_1(\\boldsymbol{s})$",
    fill = "$\\beta_1(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()


################################# QBiT

tikz(paste0('beta1_QBiT_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = mhat[, 2], nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(0, 1), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_1(\\boldsymbol{s})$",
    fill = "$\\beta_1(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()

################################# TPS

tikz(paste0('beta1_TPS_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = mhat.tp[, 2], nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(0, 1), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_1(\\boldsymbol{s})$",
    fill = "$\\beta_1(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()


################################################## plots for estimates of beta2

tikz(paste0('beta2_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = beta3(S.pop[, 1], S.pop[, 2]) + 
                          xi2(S.pop[, 1], S.pop[, 2]) * Finv, nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(-1, 1.1), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_{2,\\tau}(\\boldsymbol{s})$",
    fill = "$\\beta_{2,\\tau}(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()

################################# SQBiT

tikz(paste0('beta2_SQBiT_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = mhat.sm[, 3], nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(-1, 1.25), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_2(\\boldsymbol{s})$",
    fill = "$\\beta_2(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()


################################# QBiT

tikz(paste0('beta2_QBiT_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = mhat[, 3], nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(-1, 1.25), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_2(\\boldsymbol{s})$",
    fill = "$\\beta_2(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()


################################# TPS

tikz(paste0('beta2_TPS_tau', tau * 100, '.tex'), width = 4, height = 3, standAlone = TRUE)

interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
                        z = mhat.tp[, 3], nx = 140, ny = 140)

interp_df <- with(interp_result, expand.grid(x = x, y = y))
interp_df$z <- as.vector(interp_result$z)

# Plot smooth heatmap
ggplot(interp_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient(limits = c(-1, 1.25), 
                      low = "deepskyblue", 
                      high = "magenta", 
                      oob = scales::squish) +
  labs(
    x = "$s_1$",
    y = "$s_2$",
    z = "$\\beta_2(\\boldsymbol{s})$",
    fill = "$\\beta_2(\\boldsymbol{s})$"
  ) + 
  coord_fixed() +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

dev.off()

############ histogram
library('sn')
x <- rsn(5000, xi = 0, omega = 1, alpha = 5)

tikz('hist_skewnormal.tex', width = 4, height = 3, standAlone = TRUE)
hist(x, border = 'magenta', col = 'white',
     xlab = '$\\varepsilon$',
     main = "Right-skewed Skew-Normal Distribution")
dev.off()



