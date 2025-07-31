#' Smoothed Quantile Bivariate Triangulation using penalization parameters selected by k-fold Cross Validation (SQBiT-CV)
#'
#' This function performs the point estimation of quantile spatial model (QSM) with penalization parameters selected by k-fold CV using SQBiT. It returns the selected penalization parameters.
#'
#' \strong{What is Quantile Spatial Model (QSM)?}
#'
#' The proposed Quantile Spatial Model (QSM) is given by:
#' \deqn{
#' Q_{\tau}(Y_i \mid \mathbf{C}_i,\mathbf{X}_i, \mathbf{S}_i) = \mathbf{C}_i^{\top} \boldsymbol{\eta}_{\tau}
#' + \sum_{j=0}^p X_{ij} \beta_{j,\tau}(\mathbf{S}_i), \quad i = 1, \ldots, n.
#' }
#' where \eqn{Q_Y(\tau \mid \cdot) = \inf\{y: P(Y \leq y \mid \cdot) \geq \tau\}} denotes the
#' \eqn{\tau}th (\eqn{\tau \in (0, 1)}) conditional quantile of \eqn{Y},
#' \eqn{\boldsymbol{\eta}_{\tau} = (\eta_{1,\tau}, \ldots, \eta_{q,\tau})^{\top}} captures the constant effects,
#' and \eqn{\beta_{j,\tau} \colon \Omega \mapsto \mathbb{R}} is the unknown spatially varying coefficient function.
#' \cr
#'
#' This semi-parametric model captures spatial heterogeneity through the spatially varying
#' coefficients \eqn{\beta_{j,\tau}(\cdot)}, and non-spatial covariate effects through the
#' constant coefficients \eqn{\eta_{\ell,\tau}}. The partially linear nature of the model offers both flexibility and interpretability.
#'
#' \strong{How to use QSM?}
#'
#' In general, there are four steps to analyze spatial data using QSM.
#'
#' 1. Create a triangulation mesh for your spatial domain. This can be done using "Triangulation" package. It can be downloaded in github using install_github("funstatpackages/Triangulation").
#' 2. Generate the bivariate triangulation spline given a degree of polynomials \eqn{d} and a smoothness parameter \eqn{r}. This can be done using "BPST". It can be downloaded in github using install_github("FIRST-Data-Lab/BPST").
#' 3. Fit QSM using SQBiT.
#' 4. Plot the estimates for the varying coefficient.
#'
#' \strong{Optimization:}
#'
#' The optimization is based on alternating direction method of multipliers (ADMM). For the algorithmic details, we refer users to Boyd et al (2011). Note that it has been shown that the algorithm will converge under mild conditions even if a fixed step size is used for ADMM.
#'
#' \strong{How the tuning parameters are selected?}
#'
#' To reduce computational cost, we use a two-step tuning approach to select \eqn{\boldsymbol{\lambda}} using k-fold CV.
#' \cr
#' First, we examine 10 equally spaced points between -4 and 1 for \eqn{\log_{10}(\lambda_j)},
#' where \eqn{j = 0, \ldots, p}, to obtain an initial penalization parameter \eqn{\lambda_j^*}.
#' \cr
#'
#' Then, we refine \eqn{\lambda_j^*} by selecting from 10 equally spaced points between
#' \eqn{\log_{10}(10^{-2} \lambda_j^*)} and \eqn{\log_{10}(10^2 \lambda_j^*)},
#' allowing a broader initial search and finer adjustment around \eqn{\lambda_j^*}.
#'
#' The current implementation assumes the same penalization parameters for different varying coefficient functions for simplicity. The above procedure uses lambda_start = 1e-4, lambda_end = 1e1, and lambda.scale = 100. But they can be modified according to the user's needs.
#'
#' @param h smoothing bandwidth
#' @param tau quantile level
#' @param Y response variable
#' @param C covariates with constant effects
#' @param X covariates with varying coefficients
#' @param B Bernstein basis polynoimials. An object from basis() in package BPST.
#' @param Q2 Q2 matrix from the QR decomposition of the smoothness matrix of triangles
#' @param P Penalty matrix. See examples. It is based on output from basis() in package BPST.
#' @param nfold the number of equal-sized subsets (called folds) into which the dataset is split for the purpose of model validation.
#' @param eta.j1 step size of the first round. Set it to be 1 if var.j = TRUE
#' @param eta.j2 step size of the second round. Set it to be 1 if var.j = TRUE
#' @param var.j whether to vary the step size in the ADMM
#' @param eps.abs first convergence criterion of ADMM
#' @param eps.rel second convergence criterion of ADMM
#' @param lambda_start the lower end of the range for the search of penalization parameters in the first round
#' @param lambda_end the upper end of the range for the search of penalization parameters in the first round
#' @param nlambda the number of penalization parameters to be searched in the first round
#' @param new.nlambda the number of penalization parameters to be searched in the second round
#' @param lambda.scale the scale of the search of penalization parameters around the penalization parameters selected in the first round. This is for the second round.
#' @return A list with components:
#' \describe{
#'   \item{\code{lambdac}}{The selected penalization parameter.}
#'   \item{\code{cverr.min}}{The k-fold cross validation quantile check loss for the selected penalization parameter.}
#'   \item{\code{lambdas}}{The set of penalization parameters in the second round.}
#'   \item{\code{cverr}}{The k-fold cross validation quantile check loss}
#'   \item{\code{time}}{Computational time used.}
#' }
#' @examples
#' \donttest{
#'###### Step 1. Create a Triangulation Mesh
#' boundaries <- matrix(c(0, 0,
#'                        0, 1,
#'                        1, 1,
#'                        1, 0), nrow = 4, byrow = TRUE)
#' tri <- TriMesh(boundaries, n = 8)
#' Tr <- tri$Tr
#' V <- round(tri$V, 3)
#'
#' # Population locations
#' s1 <- runif(10000)
#' s2 <- runif(10000)
#' S_pop <- data.frame(s1 = s1, s2 = s2)
#'
#' # Determine which points are inside triangulation
#' ind1 <- inVT(V0 = V, Tr0 = Tr, xx = S_pop[,1], yy = S_pop[,2])
#' ind1 <- ind1$ind.inside
#' ind2 <- (1:nrow(S_pop))[!is.na(S_pop[,1])]
#' ind <- sort(intersect(ind1, ind2))
#' pop.r <- S_pop[ind, ]
#' Npop <- nrow(pop.r)
#'
#' # Coordinates and coefficient functions
#' S.pop <- pop.r[, c(1, 2)]
#' beta1 <- function(s1, s2) sin(pi * s1 * s2)
#' beta2 <- function(s1, s2) (1 - (1 - 2 * s1)^2) * (1 - (1 - 2 * s2)^2)
#' pop.r <- cbind(S.pop, beta1(S.pop[,1], S.pop[,2]), beta2(S.pop[,1], S.pop[,2]))
#' colnames(pop.r)[3:4] <- c('beta0', 'beta1')
#'
#' ####### Step 2. Generate bivariate triangulation basis
#' d <- 3
#' r <- 1
#' B0.pop <- basis(V = V, Tr = Tr, d = d, r = r, Z = as.matrix(S.pop))
#' Q2 <- B0.pop$Q2
#' B.pop <- B0.pop$B
#' K <- B0.pop$K
#' P <- t(Q2) %*% K %*% Q2
#' BQ2.pop <- as.matrix(B.pop %*% Q2) ## basis function for population locations
#'
#' # Simulation parameters
#' n <- 2000
#' tau <- 0.5
#' eta <- matrix(c(1, 1, 1))
#' Sigma <- diag(3)
#'
#' # Sample population and simulate covariates
#' ind.s <- sample(1:Npop, n)
#' data <- as.matrix(pop.r[ind.s, ])
#' S <- data[, c(1, 2)]
#' C <- matrix(runif(3 * n, -1, 1), ncol = 3)
#' X <- matrix(runif(n, -1, 1), n, 1)
#' beta0 <- data[, c('beta0', 'beta1')]
#' Y <- beta0[,1] + X[,1] * beta0[,2] + C %*% eta + rnorm(n = n)
#'
#' ####### Step 3. Estimate QSM using SQBiT
#' X <- cbind(1, X)
#' B <- B.pop[ind.s, ]
#' lambda.cv <- SQBiT_cv(Y = Y, X = X, C = C, P = P, B = B, Q2 = Q2, tau = tau)
#' mod <- SQBiT(Y = Y, X = X, C = C, P = P, B = B, Q2 = Q2, lambda = lambda.cv$lambda, tau = tau)
#'
#' ####### Step 4. Plot varying coefficient
#' p <- ncol(X)
#' q <- ncol(C)
#' Jn <- dim(B.pop %*% Q2)[2]
#' mhat.sm <- c()
#' for(i in 1:p){
#'  mhat.sm <- cbind(mhat.sm, BQ2.pop %*% mod$gamma[(1 + (i - 1) * Jn):(Jn * i)])
#'}
#' interp_result <- interp(x = S.pop[, 1], y = S.pop[, 2], duplicate = "mean",
#'                        z = mhat.sm[, 2], nx = 140, ny = 140)
#' interp_df <- with(interp_result, expand.grid(x = x, y = y))
#' interp_df$z <- as.vector(interp_result$z)
#'
#' # Plot smooth heatmap
#' ggplot(interp_df, aes(x = x, y = y, fill = z)) +
#'  geom_raster(interpolate = TRUE) +
#'  scale_fill_gradient(limits = c(0, 1),
#'                      low = "deepskyblue",
#'                      high = "magenta",
#'                      oob = scales::squish) +
#'  labs(
#'    x = "$s_1$",
#'    y = "$s_2$",
#'    z = "$\\beta_1(\\boldsymbol{s})$",
#'    fill = "$\\beta_1(\\boldsymbol{s})$"
#'  ) +
#'  coord_fixed() +
#'  theme_minimal() +
#'  theme(panel.grid = element_blank(),
#'        plot.title = element_text(hjust = 0.5))
#' }
#'
#' @export
#' @references
#' Kim, M., Wang, L., and Wang, H. J. (2025). Estimation and inference of quantile spatially varying coefficient models over complicated domains. *Journal of the American Statistical Association*, 1--15. Taylor & Francis. Forthcoming.
#' @references
#' He, X., Pan, X., Tan, K.M., & Zhou, W.-X. (2023). Smoothed quantile regression with large-scale inference. *Journal of Econometrics*, 232(2), 367–388.
#' @references
#' Boyd, S., Parikh, N., Chu, E., Peleato, B., Eckstein, J., and others (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers. *Foundations and Trends in Machine Learning*, 3(1), 1--122. Now Publishers.
SQBiT_cv <- function(h = 7 * tau * (1 - tau) * ((ncol(C) + dim(Q2)[2] * ncol(X) + log(length(Y))) / length(Y)) ^ (2/5),
                     tau = .5, Y, C, X, P, B, Q2,
                     nlambda = 10, new.nlambda = 10, lambda.scale = 100,
                     eps.abs = 1e-4, eps.rel = 1e-2, nfold = 5,
                     eta.j1 = .32, eta.j2 = .32, var.j = FALSE,
                     lambda_start = 10 ^ (-4), lambda_end = 10 ^ 1) {
  start_time <- Sys.time()
  lambda_start <- log10(lambda_start)
  lambda_end <- log10(lambda_end)
  lambda <- 10 ^ (seq(lambda_start, lambda_end, length.out = nlambda))
  lambda <- as.matrix(lambda)

  cverr <- c()
  for (i in 1:nlambda) {
    cverr[i] <- mean(cv.pred.SQBiT(y = Y, C, X, B, Q2, P, lambda = lambda[i], h = h,
                                   eta.j = eta.j1, nfold = nfold, tau = tau,
                                   var.j = var.j))
  }
  lambdac <- lambda[which.min(cverr)]

  new.lambda <- 10 ^ (seq(log10(lambdac / lambda.scale),
                          log10(lambdac * lambda.scale),
                          length.out = new.nlambda))

  cverr <- c()
  for (i in 1:new.nlambda) {
    cverr[i] <- mean(cv.pred.SQBiT(y = Y, C, X, B, Q2, P, lambda = new.lambda[i],
                                   eta.j = eta.j2, nfold = nfold, tau = tau,
                                   h = h, var.j = var.j))
  }
  lambdac <- new.lambda[which.min(cverr)]
  end_time <- Sys.time()

  return(list(lambda = lambdac, lambdas = new.lambda, cverr = cverr, cverr.min = min(cverr), time = end_time - start_time))
}
