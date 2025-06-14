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
source('cv.pred.tps.R') 

# c++ matrix inversion
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }",
                  depends = "RcppArmadillo")
# c++ matrix multiplication
sourceCpp('test.cpp')

# Quantile level
tau <- 0.75

# Load matrices for triangles (Tr) and vertices (V)
Tr <- as.matrix(read.csv('Tr_mort2.csv', header = TRUE))
V <- as.matrix(read.csv('V_mort2.csv', header = TRUE))
bbound <- as.matrix(read.csv('brusa.csv', header = FALSE))
head(bbound) 

# Tr = VT$Tr
# V = VT$V

# Load data
# set working directory
setwd('~/Desktop/Research/Quantile Regression/mortality_app/')
Mort_data <- read.csv('Data_mortality.csv', header = T)
S_pop <- as.matrix(read.csv('Mortal_pop_location_d025.csv', header = T))

# degree and smoothness
d <- 3
r <- 1

# population basis and index
B0_pop <- basis(V, Tr, d, r, S_pop)
Q2 <- B0_pop$Q2
nq <- dim(Q2)[2]
BQ2_pop <- as.matrix(B0_pop$B %*% Q2)
Ind_all <- B0_pop$Ind.inside
K <- B0_pop$K
P <- t(Q2) %*% K %*% Q2

# Load response and covariates
Y <- as.matrix(Mort_data[,7])
X <- Mort_data[,1:4]
X <- as.matrix(scale(X, scale = TRUE)) 
X <- as.matrix(cbind(1, X))
S <- cbind(Mort_data$Longtitude, Mort_data$Latitude)


# Generate sample Bivariate spline basis 
B0 <- basis(V, Tr, d, r, S)
Ind <- B0$Ind.inside
Y <- Y[Ind]
X <- X[Ind, ]
C <- as.matrix(X[, 3])
X <- X[, -3]
X <- as.matrix(X)
S <- S[Ind,]
B <- B0$B
BQ2 <- as.matrix(B %*% Q2) 

# sample size
n <- length(Y)

### sensitivity analysis
cs <- seq(from = 1, to = 10, length.out = 10)
lengths <- c()
for (c in cs) {
  
  ### cross-validation for h
  h.de <- c * tau * (1 - tau) * ((1 + dim(Q2)[2] * 4 + log(n)) / n) ^ (2/5)
  
  
  ### tuning parameter selection for smoothed QSVCM
  lambda_SQBiT <- SQBiT_tune(Y = Y, X = X, C = C,  
                             P = P, B = B, Q2 = Q2,  
                             h = h.de, tau = tau, 
                             lambda_start = 1e-4, lambda_end = 1e1,
                             nlambda = 10, new.nlambda = 10, 
                             lambda.scale = 100) 
  
  mod.sm.ie <- SQBiT(h = h.de, C = C, Y = Y, X = X, B = B, P = P, Q2 = Q2,
                     tau = tau, lambda = matrix(lambda_SQBiT$lambda, ncol = 4), 
                     var.j = FALSE, adaptive.h = FALSE,
                     interval = TRUE) 
  
  # wild bootstrap for interval estimation of varying coefficient
  set.seed(666)
  
  mod.sm.wb <- smqsvcm_admm.wb(h = h.de, tau = tau, Y = Y, C = C, X = X, P = P, B = B, Q2 = Q2, eta.j = .32,
                               lambda = matrix(lambda_SQBiT$lambda, ncol = 4), 
                               biascorr = TRUE,
                               Br = 500, BQ2.eva = BQ2_pop, max.iter = 200, 
                               compute.vc = TRUE,
                               gamma.hat = c(mod.sm.ie$eta, mod.sm.ie$gamma))
  
  lengths <- c(lengths, mod.sm.wb$cis2[2] - mod.sm.wb$cis2[1]) # constant coefficient
  
}

plot(cs, lengths, xlab = '$c$', ylab = 'Length',
     type = 'b', col = 'deepskyblue',
     main = 'Sensitivity Analysis ($\\tau = 0.75$)')

### bandwidth
h.de <- 2 * tau * (1 - tau) * ((1 + dim(Q2)[2] * 4 + log(n)) / n) ^ (2/5)


### tuning parameter selection for smoothed QSVCM
lambda_SQBiT <- SQBiT_tune(Y = Y, X = X, C = C,  
                           P = P, B = B, Q2 = Q2,  
                           h = h.de, tau = tau, 
                           lambda_start = 1e-4, lambda_end = 1e1,
                           nlambda = 10, new.nlambda = 10, 
                           lambda.scale = 100) 

# tensor product QR
tpqr.mod <- tpqr.plm(Y = Y, C = C, X = X, S = S, d = 3, tau = tau,
                     c0 = 1.5, c1 = 3.5, nc = 6) 



# boundary of US
SS_pop <- S_pop[Ind_all, ]
boundary <- read.csv("main_boundusa.csv",header=TRUE)
S1b <- cbind(min(boundary[, 1]), max(boundary[, 1])) 
S2b <- cbind(min(boundary[, 2]), max(boundary[, 2])) 

# tps basis
knots1 <- tpqr.mod$knots[[1]]
knots2 <- tpqr.mod$knots[[2]]
B1 <- bs(SS_pop[, 1], df = length(knots1), degree = d)
B2 <- bs(SS_pop[, 2], df = length(knots2), degree = d)
basis <- do.call(rbind, lapply(1:nrow(SS_pop), function(i) kronecker(B1[i,], B2[i,])))
jn <- ncol(basis)

# coefficient plot (population)
p <- ncol(X)
Jn <- dim(BQ2)[2] 


# interval estimation
mod.sm.ie <- SQBiT(h = h.de, C = C, Y = Y, X = X, B = B, P = P, Q2 = Q2,
                   tau = tau, lambda = matrix(lambda_SQBiT$lambda, ncol = 4), 
                   var.j = FALSE, adaptive.h = FALSE,
                   interval = TRUE) 


# wild bootstrap for interval estimation of varying coefficient
set.seed(666)

mod.sm.wb <- smqsvcm_admm.wb(h = h.de, tau = tau, Y = Y, C = C, X = X, P = P, B = B, Q2 = Q2, eta.j = .32,
                             lambda = matrix(lambda_SQBiT$lambda, ncol = 4), 
                             biascorr = FALSE,
                             Br = 500, BQ2.eva = BQ2_pop, max.iter = 200, 
                             compute.vc = TRUE,
                             gamma.hat = c(mod.sm.ie$eta, mod.sm.ie$gamma))

# interval estimates
mod.sm.wb$cis2 # constant coefficient
mhat_all.sm_lb <- mod.sm.wb$betas.lb # lower bound of varying coefficient
mhat_all.sm_ub <- mod.sm.wb$betas.ub # upper bound of varying coefficient



################### heat-map

mhat_all.sm <- c()
mhat.tp <- c()
for(i in 1:p){ 
  mhat_all.sm <- cbind(mhat_all.sm, BQ2_pop %*% mod.sm.ie$gamma[(1 + (i - 1) * Jn):(Jn * i)])
  mhat.tp <- cbind(mhat.tp, basis %*% tpqr.mod$gamma[(jn * (i - 1) + 1):(jn * i)])
}


# longitude and latitude
dist <- 0.25
uu <- seq(S1b[1], S1b[2], dist)
vv <- seq(S2b[1], S2b[2], dist)
n1 <- length(uu)
n2 <- length(vv)
u <- rep(uu, n2)
v <- rep(vv, rep(n1, n2))
uvpop <- cbind(u, v)

# match
Index <- match(data.frame(t(SS_pop)), data.frame(t(uvpop)))

# limits
zlims <- rbind(c(.5, 1.5),
               c(-.15, .15),
               c(-.15, .15),
               c(-.15, .15))

# coefficient plots
options(tikzLatexPackages 
        = c(getOption( "tikzLatexPackages" ),
            "\\usepackage{amsmath,amsfonts,amsthm, palatino, mathpazo}"))

# SQPLSVCM
for (i in 1:4) {
  tikz(paste0('heat75beta', i, 'sm.tex'), width = 6, height = 3, standAlone = TRUE)
  mhat_all.sm <- as.matrix(mhat_all.sm)
  Zpop <- matrix(NaN, dim(uvpop)[1], dim(mhat_all.sm)[2])
  Zpop[Index, ] <- mhat_all.sm
  Y0 <- matrix(rep(NA, n1 * n2), ncol = 1)
  Y0[Index] <- Zpop[Index, i] # for the first one
  index <- point.in.polygon(u, v, boundary[, 1], boundary[, 2])
  Y0[index == 0] <- NA
  
  
  Y0[which(Y0 > zlims[i, 2])] <- zlims[i, 2]
  Y0[which(Y0 < zlims[i, 1])] <- zlims[i, 1]
  
  # coefficient plot
  alpha_fitted <- data.frame(u, v, Y0)
  coordinates(alpha_fitted) <- ~u + v
  gridded(alpha_fitted) <- TRUE
  rasterDF <- raster(alpha_fitted)
  colours <- colorRamps::matlab.like(100)  
  
  par(mar = c(1, 2, 1, 2) + 0.1)
  plot(rasterDF, col = colours, axes = FALSE, box = FALSE, 
       legend.width = 0.75, legend.shrink = 0.85, zlim = zlims[i, ])
  US(xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = TRUE, lty = 2)
  dev.off()
}

for (i in 1:4) {
  tikz(paste0('heat75beta', i, 'tps.tex'), width = 6, height = 3, standAlone = TRUE)
  mhat.tp <- as.matrix(mhat.tp)
  Zpop <- matrix(NaN, dim(uvpop)[1], dim(mhat.tp)[2])
  Zpop[Index, ] <- mhat.tp
  Y0 <- matrix(rep(NA, n1 * n2), ncol = 1)
  Y0[Index] <- Zpop[Index, i] # for the first one
  index <- point.in.polygon(u, v, boundary[, 1], boundary[, 2])
  Y0[index == 0] <- NA
  
  Y0[which(Y0 > zlims[i, 2])] <- zlims[i, 2]
  Y0[which(Y0 < zlims[i, 1])] <- zlims[i, 1]
  
  # coefficient plot
  alpha_fitted <- data.frame(u, v, Y0)
  coordinates(alpha_fitted) <- ~u + v
  gridded(alpha_fitted) <- TRUE
  rasterDF <- raster(alpha_fitted)
  colours <- colorRamps::matlab.like(100) 
  
  par(mar = c(1, 2, 1, 2) + 0.1)
  plot(rasterDF, col = colours, axes = FALSE, box = FALSE, 
       legend.width = 0.75, legend.shrink = 0.85, zlim = zlims[i, ])
  US(xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = TRUE, lty = 2)
  dev.off()
}

# lower bound
for (i in 1:4) {
  tikz(paste0('heat75beta', i, 'smlb.tex'), width = 6, height = 3, standAlone = TRUE)
  mhat_all.sm_lb <- as.matrix(mhat_all.sm_lb)
  Zpop <- matrix(NaN, dim(uvpop)[1], dim(mhat_all.sm_lb)[2])
  Zpop[Index, ] <- mhat_all.sm_lb
  Y0 <- matrix(rep(NA, n1 * n2), ncol = 1)
  Y0[Index] <- Zpop[Index, i] # for the first one
  index <- point.in.polygon(u, v, boundary[, 1], boundary[, 2])
  Y0[index == 0] <- NA
  
  
  Y0[which(Y0 > zlims[i, 2])] <- zlims[i, 2]
  Y0[which(Y0 < zlims[i, 1])] <- zlims[i, 1]
  
  
  # coefficient plot
  alpha_fitted <- data.frame(u, v, Y0)
  coordinates(alpha_fitted) <- ~u + v
  gridded(alpha_fitted) <- TRUE
  rasterDF <- raster(alpha_fitted)
  colours <- colorRamps::matlab.like(100) 
  
  par(mar = c(1, 2, 1, 2) + 0.1)
  plot(rasterDF, col = colours, axes = FALSE, box = FALSE, 
       legend.width = 0.75, legend.shrink = 0.85, zlim = zlims[i, ])
  US(xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = TRUE, lty = 2)
  dev.off()
}

# upper bound
for (i in 1:4) {
  tikz(paste0('heat75beta', i, 'smub.tex'), width = 6, height = 3, standAlone = TRUE)
  mhat_all.sm_ub <- as.matrix(mhat_all.sm_ub)
  Zpop <- matrix(NaN, dim(uvpop)[1], dim(mhat_all.sm_ub)[2])
  Zpop[Index, ] <- mhat_all.sm_ub
  Y0 <- matrix(rep(NA, n1 * n2), ncol = 1)
  Y0[Index] <- Zpop[Index, i] # for the first one
  index <- point.in.polygon(u, v, boundary[, 1], boundary[, 2])
  Y0[index == 0] <- NA
  
  
  Y0[which(Y0 > zlims[i, 2])] <- zlims[i, 2]
  Y0[which(Y0 < zlims[i, 1])] <- zlims[i, 1]
  
  # coefficient plot
  alpha_fitted <- data.frame(u, v, Y0)
  coordinates(alpha_fitted) <- ~u + v
  gridded(alpha_fitted) <- TRUE
  rasterDF <- raster(alpha_fitted)
  colours <- colorRamps::matlab.like(100) 
  
  par(mar = c(1, 2, 1, 2) + 0.1)
  plot(rasterDF, col = colours, axes = FALSE, box = FALSE, 
       legend.width = 0.75, legend.shrink = 0.85, zlim = zlims[i, ])
  US(xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = TRUE, lty = 2)
  dev.off()
}

