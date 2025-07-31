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


# population
pop.all <- as.matrix(read.csv('horse_tr1_n5K.csv', header = TRUE))
pop.all <- pop.all[, -1]
S <- pop.all[, c(1, 2)]
S <- as.matrix(as.data.frame(S))
# triangulation matrix
# nrow is the number of triangles
# each row is the indices of vertices 


#### Triangulation

Tr <- as.matrix(read.csv('Tr_1.csv', header = TRUE))     ########### change Tr_2.csv for Triangulation 2
# vertices of triangulation
V <- as.matrix(read.csv('V_1.csv', header = TRUE))       ########### change V_2.csv for Triangulation 2


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
colnames(pop.r) <- c('u', 'v', 'b0', 'b1')

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
n <- 2000                               
nsim <- 200
set.seed(123)
seeds <- sample(1:1000, nsim)
tau <- 0.50                             ########### change tau = 0.75, 0.90 for tau = 0.75, 0.90

########################### setup for constant covariates
eta <- matrix(c(1, 1, 1))
Sigma <- matrix(0, 3, 3)
for (i in 1:3) {
  for (j in 1:3) {
    Sigma[i, j] <- 0.5 ^ (abs(i - j))
  }
}
 
################# varying coefficient function
a1 <- 0.10               ########### signal strength for beta2
a2 <- 0.35               ########### signal strength for beta3
beta2 <- function(s1, s2) 2 * sin(a1 * pi * (s1 ^ 2 + s2 ^ 2))
beta3 <- function(s1, s2) -sin(a2 * pi * (s1 + 0.15)) + sin(a2 * pi * s2)



################# Accuracy for model selection

#### oracle percentage
op <- c()
nv <- c()
acc.c <- c()
acc.x <- c()

#### simulation studies

for (i in 1:nsim) {
  set.seed(seeds[i])
  
  ### sample locations
  ind.s <- sample(1:Npop, n)
  data <- as.matrix(pop.r[ind.s, ])
  beta0 <- data[, c('b0', 'b1')]
  
  #### sample locations 
  S <- data[, c(1, 2)]
  
  #### covariates with constant effects
  C <- draw.d.variate.uniform(no.row = n,
                              d = 3,
                              cov.mat = Sigma) * sqrt(3)
  C <- C - 0.5 * sqrt(3)
  
  ### covariates with varying effects
  X <- matrix(runif(2 * n, -1, 1) , n, 2)
  Y <- beta0[, 1] + X[, 1] * beta2(S[, 1], S[, 2]) + X[, 2] * beta3(S[, 1], S[, 2])  + C %*% eta + 
    rsn(n = n, alpha = 5)
  X <- cbind(1, X)  
  
  #### Bivariate spline basis 
  B <- B.pop[ind.s, ]
  
  #### forward selection based on GACV 
  cv.sets <- SQBiT_forward(tau = tau, Y = Y, C = cbind(C, X[, c(2, 3)]), Q2 = Q2,
                           P = P, B = B,
                           lambda_start = 1e-4, 
                           lambda_end = 1e1,
                           lambda.scale = 100,
                           nlambda = 10, new.nlambda = 10)
  print(list(constant = cv.sets$constant, varying = cv.sets$varying))
  
  #### indices of covariates identified for constant coefficients
  constant.sets <- cv.sets$constant
  constant.sets <- constant.sets[order(constant.sets)]
  constant.sets <- as.numeric(constant.sets)
  
  #### indices of covariates identified for constant coefficients
  if (is.null(cv.sets$varying) == TRUE) {
    varying.sets <- 0
  } else {
    varying.sets <- cv.sets$varying
    varying.sets <- varying.sets[order(varying.sets)]
    varying.sets <- as.numeric(varying.sets)
  }
  
  #### accuracy for identifying them correctly
  acc.c <- c(acc.c, mean(c(1, 2, 3) %in% constant.sets))
  acc.x <- c(acc.x, mean(c(4, 5) %in% varying.sets))
  
  #### oracle percentage
  op[i] <- ifelse(identical(constant.sets, c(1, 2, 3)) & identical(varying.sets, c(4, 5)), 1, 0) 
  
  cat("iteration = ", i, "\n")
}


#### accuracy for the constant
mean(acc.c) * 100
sd(acc.c) / sqrt(100) * 100

#### accuracy for the varying
mean(acc.x)
sd(acc.x) / sqrt(100) * 100

#### oracle percentage
mean(op) * 100
sd(op) / sqrt(100) * 100 
 
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
tikz(paste0('sim_ms_beta2_a1_', a1, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = pop.r[, 1], y = pop.r[, 2], duplicate = "mean",
              z = beta2(pop.r[, 1], pop.r[, 2]), nx = 140, ny = 140)
# heatmap
heat.2d.hs(int, zlim = c(-2, 2), 
           col = tim.colors(64),
           outside.color = 'black')
# horseshoe polygon
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()


tikz(paste0('sim_ms_beta3_a2_', a2, '.tex'), width = 4, height = 3, standAlone = TRUE)
int <- interp(x = pop.r[, 1], y = pop.r[, 2], duplicate = "mean",
              z = beta3(pop.r[, 1], pop.r[, 2]), nx = 140, ny = 140)
# heatmap
heat.2d.hs(int, zlim = c(-2, 2), 
           col = tim.colors(64),
           outside.color = 'black')
# horseshoe polygon
polygon(rbind(boundary, c(3.5, 0.9),
              c(-0.9, 0.9),
              c(-0.9, -0.9),
              c(3.4, -0.9),
              c(3.5, 0.9)), border = 'transparent', lty = 2, col = 'white')
dev.off()