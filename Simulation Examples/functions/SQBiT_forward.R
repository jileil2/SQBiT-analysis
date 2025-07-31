#' Forward Selection of Quantile Spatial Model (QSM) based on Generalized Approximate Cross Validation using Smoothed Quantile Bivariate Triangulation
#'
#' This function performs forward selection of quantile spatial model (QSM) using GACV. It returns the indices (column numbers) for the constant coefficients, varying coefficients.
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
#' To reduce computational cost, we use a two-step tuning approach to select \eqn{\boldsymbol{\lambda}} using GACV.
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
#' @param C covariates to be decided for constant or varying coefficients
#' @param rlt the rule of thumb for \eqn{\tau (1-\tau)  \{[q+J_n p +\log n]/n\}^{2/5}}, where \eqn{J_n} is the number of basis functions
#' @param nlambda the number of penalization parameters to be searched in the first round
#' @param new.nlambda the number of penalization parameters to be searched in the second round
#' @param lambda_start the lower end of the range for the search of penalization parameters in the first round
#' @param lambda_end the upper end of the range for the search of penalization parameters in the first round
#' @param lambda.scale the scale of the search of penalization parameters around the penalization parameters selected in the first round. This is for the second round.
#' @return A list with components:
#' \describe{
#'   \item{\code{constant}}{The indices for covariates that are selected for constant coefficients.}
#'   \item{\code{varying}}{The indices for covariates that are selected for varying coefficients.}
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
#' QSM <- SQBiT_forward(tau = tau, Y = Y, C = cbind(C, X[, -1]), P = P, B = B, Q2 = Q2)
#' XC <- cbind(C, X[, -1])
#' new.C <- XC[, QSM$constant]
#' new.X <- cbind(1, XC[, QSM$varying])
#' mod <- SQBiT_gacv(Y = Y, X = new.X, C = new.C, P = P, B = B, Q2 = Q2, tau = tau)
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
#' @references
#' Yuan, M. (2006). GACV for quantile smoothing splines. *Computational Statistics & Data Analysis*, 50(3), 813-829.
SQBiT_forward <- function(tau, Y, C, B, Q2, P, rlt = 7,
                          lambda_start = 1e-4, lambda_end = 1e1,
                          lambda.scale = 100,
                          nlambda = 10, new.nlambda = 10) {

  # dimension
  p <- ncol(C)
  n <- length(Y)

  # base model with all constant effects
  h.base <- rlt * tau * (1 - tau) * ((p + dim(Q2)[2] * 1 + log(n)) / n) ^ (2/5)

  mod.base <- SQBiT_gacv(Y = Y, X = matrix(1, nrow = n), C = C,
                         P = P, B = B, Q2 = Q2,
                         nlambda = nlambda, new.nlambda = new.nlambda,
                         lambda_start = lambda_start, lambda_end = lambda_end,
                         lambda.scale = lambda.scale,
                         h = h.base, tau = tau)
  cri.base <- mod.base$gacv.min

  # forward selection
  index.C <- 1:p
  index.X <- c()
  ms.cri <- vector('list', length = p + 1)
  ms.cri[[1]] <- cri.base
  t <- 1

  while(p >= 1) {
    cris <- c()
    for (j in index.C) {

      C.j <- C[, setdiff(index.C, j), drop = FALSE]
      X.j <- cbind(1, C[, c(index.X, j)])

      # make bandwidth adaptive to the dimensions
      h.j <- rlt * tau * (1 - tau) * ((ncol(C.j) + dim(Q2)[2] * ncol(X.j) + log(n)) / n) ^ (2/5)


      mod.j <- SQBiT_gacv(Y = Y, X = X.j, C = C.j,
                          P = P, B = B, Q2 = Q2,
                          nlambda = nlambda, new.nlambda = nlambda,
                          lambda_start = lambda_start, lambda_end = lambda_end,
                          lambda.scale = lambda.scale,
                          h = h.j, tau = tau)
      cri.j <- mod.j$gacv.min


      cris <- c(cris, cri.j)
    }

    t <- t + 1
    ms.cri[[t]] <- cris

    if (min(cris) > cri.base) {
      break
    }  else {
      index.X <- c(index.X, index.C[which.min(cris)])
      index.C <- index.C[-which.min(cris)]
    }

    cri.base <- min(cris)
    p <- length(index.C)
  }

  return(list(constant = index.C, varying = index.X))
}
