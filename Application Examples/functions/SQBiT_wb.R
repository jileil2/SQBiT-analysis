#' Wild Bootstrap of Smoothed Quantile Bivariate Triangulation (SQBiT-WB)
#'
#' This function performs wild bootstrap inference for quantile spatial model (QSM) with initial estimates from output from SQBiT. It returns the confidence intervals for constant coefficients and varying coefficients. If the basis function for the locations to be predicted is provided, it also computes the confidence intervals for the given locations.
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
#' \strong{Inference:}
#'
#' The inference for SQBiT can be conducted using wild bootstrap. We refer to Feng et al (2011) for the details of wild bootstrap for quantile regression. Bias correction can be applied by specifying bias.corr = TRUE, which employs the Bahadur representation of SQBiT. Note that in our investigation, we found that the coverage using penalization parameters selected k-fold CV tends to be more stable than those by GACV.
#'
#' @param h smoothing bandwidth
#' @param tau quantile level
#' @param Y response variable
#' @param C covariates with constant effects
#' @param X covariates with varying coefficients
#' @param B Bernstein basis polynoimials. An object from basis() in package BPST.
#' @param Q2 Q2 matrix from the QR decomposition of the smoothness matrix of triangles
#' @param BQ2.eva the basis to be evaluated. It can be the basis function for the locations to be predicted.
#' @param lambda the penalization parameter
#' @param compute.vc whether to compute the pointwise confidence intervals for the varying coefficient. If TRUE, BQ2.eva must be provided.
#' @param Br the number of bootstrap replications
#' @param level significance level
#' @param biascorr whether to apply bias correction using Bahadur representation of the SQBiT estimator
#' @param eta.j step size. Set it to be 1 if var.j = TRUE
#' @param var.j whether to vary the step size in the ADMM
#' @param eps.abs first convergence criterion of ADMM
#' @param eps.rel second convergence criterion of ADMM
#' @param incr factor of increment for step size if stepsize is set to vary
#' @param zeta criterion for the step size to increment if stepsize is set to vary
#' @return A list with components:
#' \describe{
#'   \item{\code{etas}}{The estimated constant coefficients from bootstrap replications.}
#'   \item{\code{cis}}{The wild bootstrap confidence intervals for the constant coefficients.}
#'   \item{\code{betas.lb}}{The lower bound for the varying coefficient.}
#'   \item{\code{betas.ub}}{The upper bound for the varying coefficient.}
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
#' wb.est <- SQBiT_wb(tau = tau, Y = Y, C = C, X = X, P = P, B = B, Q2 = Q2, lambda = lambda.cv$lambda,
#'                    compute.vc = TRUE, BQ2.eva = BQ2.pop, eta.hat = mod$eta, gamma.hat = mod$gamma)
#' cis <- wb.est$cis        ### confidence intervals for constant coefficients
#' vc.lb <- wb.est$betas.lb ### lower bound for varying coefficients
#' vc.ub <- wb.est$betas.ub ### upper bound for varying coefficients
#'
#' ####### Step 4. Plot varying coefficient
#' p <- ncol(X)
#' q <- ncol(C)
#' Jn <- dim(B.pop %*% Q2)[2]
#' mhat.sm <- c()
#' for(i in 1:p){
#'  mhat.sm <- cbind(mhat.sm, BQ2.pop %*% mod$gamma[(1 + (i - 1) * Jn):(Jn * i)])
#'}
#'
#' ####### replace mhat.sm[, 2] with vc.lb[, 2] or vc.ub[, 2] for lower bound or upper bound of confidence intervals
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
#' @references
#' Feng, X., He, X., & Hu, J. (2011). Wild bootstrap for quantile regression. *Biometrika*, 98(4), 995-999.
SQBiT_wb <- function(h = 7 * tau * (1 - tau) * ((ncol(C) + dim(Q2)[2] * ncol(X) + log(length(Y))) / length(Y)) ^ (2/5),
                     tau, Y, C, X, P, B, Q2, max.iter = 50,
                     lambda, eps.abs = 1e-4, eps.rel = 1e-2,
                     var.j = FALSE, zeta = 10, incr = 2, eta.j = .32,
                     BQ2.eva, compute.vc = TRUE,
                     Br = 500, level = 0.05, eta.hat, gamma.hat, biascorr = TRUE) {
  start_time <- Sys.time()
  P <- as.matrix(P)
  Q2 <- as.matrix(Q2)
  B <- as.matrix(B)

  coef.hat <- c(eta.hat, gamma.hat)

  if (tau < 0 | tau > 1) {
    warning("Quantile level should be between 0 and 1. The default quantile level of 0.5 is used.")
    tau <- 0.5
  }

  n <- length(Y)
  np <- ncol(X)
  nq <- ncol(C)
  J <- dim(Q2)[2]

  lambda <- as.matrix(lambda)
  nl <- nrow(lambda)
  if (ncol(lambda) == 1) {
    lambda <- matrix(rep(lambda, times = np), nl, np)
  }

  # nl <- dim(lambda)[1]
  BQ2 <- eigenMapMatMult(B, Q2)
  W <- as.matrix(kr(X, BQ2, byrow = TRUE))
  W <- cbind(C, W)
  cp.W <- eigenMapMatMult(t(W), W)

  psi.eq <- function(r, eta.j) {

    psi1 <- (eta.j * r - (tau - 1 / 2)) / (eta.j + 1 / (2 * h))
    psi2 <- -(tau - 1) / eta.j + r
    psi3 <- -tau / eta.j + r
    psi1.r <- psi1 * ifelse(abs(psi1) <= h, yes = 1, no = 0)
    psi2.r <- psi2 * ifelse(psi2 < -h, yes = 1, no = 0)
    psi3.r <- psi3 * ifelse(psi3 > h, yes = 1, no = 0)
    psi <- psi1.r + psi2.r + psi3.r

    return(psi)
  }

  kp <- cppkp(diag(lambda[1,]), P)
  zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
  zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
  kp <- rbind(zm2, cbind(zm1, kp))
  cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
  mat <- eigenMapMatMult(cp.inv, t(W))

  # compute residual
  resid <- as.numeric(Y - eigenMapMatMult(W, coef.hat))
  Sigma.h <- diag(dunif(resid / h, -1, 1) / h)

  # Bahadur representation (bias correction)
  J.h <- eigenMapMatMult(t(W), W) / n
  f0 <- akj(resid, 0)$dens
  pgrad <- t(eigenMapMatMult( t(matrix(tau - punif(-resid / h, -1, 1) )), W) / n)
  if (biascorr == FALSE) {
    bias <- 0
  } else {
    # bias <- as.numeric(eigenMapMatMult(W, eigenMapMatMult(AInv(J.h), pgrad)))
    bias <- as.numeric(eigenMapMatMult(W, eigenMapMatMult(AInv(J.h) * 1 / f0, pgrad)))
  }

  # wild bootstrap
  eta.B <- c()
  gamma.B <- c()
  for (b in 1:Br) {

    # bootstrap response
    xis <- rG(n = n, tau = tau)
    Y.b <- eigenMapMatMult(W, coef.hat) + xis * abs(resid + bias)

    # old
    new.u <- (1 / eta.j) * Y.b
    new.gamma <- coef.hat
    new.fit <- eigenMapMatMult(W, new.gamma)

    ###############
    # ADMM iterations

    iter.admm <- 0
    repeat{
      iter.admm <- iter.admm + 1

      # old
      old.u <- new.u
      old.gamma <- new.gamma
      old.fit <- new.fit


      # start updating
      r <- Y.b - old.fit - old.u
      new.psi <- psi.eq(r, eta.j = eta.j)

      # equation (10)
      new.gamma <- eigenMapMatMult(mat, Y.b - new.psi - old.u)

      # equation (9)
      new.fit <- eigenMapMatMult(W, new.gamma)
      res.new <- Y.b - new.fit
      new.u <- old.u + (new.psi - res.new)

      # r = psi in the paper
      r.prime <- res.new - new.psi
      r.dual <- eta.j * (new.fit - old.fit)

      # sqrt norm
      eps.prime <- sqrt(length(r.prime)) * eps.abs +
        eps.rel * max(norm(new.psi, type = "2"),
                      norm(new.fit, type = "2"),
                      norm(Y, type = "2"))
      eps.dual <- sqrt(length(r.dual)) * eps.abs +
        eps.rel * norm(new.u, type = "2")

      # ################################
      n.prime <- norm(r.prime, type = "2")
      n.dual <- norm(r.dual, type = "2")


      if (var.j == TRUE) {
        old.eta.j <- eta.j
        if(n.prime > zeta * norm(r.dual, type = "2")){
          eta.j <- incr * eta.j
        } else if (zeta * n.prime < n.dual){
          eta.j <- eta.j / incr
        } else {
          eta.j <- eta.j
        }
      }


      if((n.prime < eps.prime) & (n.prime < eps.dual)){
        break
      }
      else if (iter.admm >= max.iter) {
        break
      }
      else {
        if (var.j == TRUE) {
          if (eta.j != old.eta.j) {
            cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
            mat <- eigenMapMatMult(cp.inv, t(W))
          }
        }
      }

    }
    gamma <- new.gamma
    eta.B <- rbind(eta.B, gamma[1:nq])
    gamma.B <- rbind(gamma.B, gamma[-(1:nq)])
  }

  # constant
  eta.hat <- coef.hat[1:nq]
  cis2 <- c()
  cis2.vc <- c()
  for (j in 1:nq) {
    ci <- quantile(eta.B[,j], c(level / 2, 1 - level / 2))
    cis2 <- rbind(cis2, c(2 * eta.hat[j] - ci[2], 2 * eta.hat[j] - ci[1]))
  }

  if (compute.vc == TRUE) {
    # VC
    gamma.hat <- coef.hat[-(1:nq)]
    betas.lb <- c()
    betas.ub <- c()
    for (j in 1:np) {
      beta.hat.j <- BQ2.eva %*% gamma.hat[(1 + (j - 1) * J):(J * j)]
      betas.B.j <- BQ2.eva %*% t(gamma.B[, (1 + (j - 1) * J):(J * j)])
      lb.j <- apply(betas.B.j, 1, function(x) quantile(x, level / 2))
      ub.j <- apply(betas.B.j, 1, function(x) quantile(x, 1 - level / 2))
      lb.j <- matrix(lb.j); ub.j <- matrix(ub.j)
      betas.lb <- cbind(betas.lb, 2 * beta.hat.j - ub.j)
      betas.ub <- cbind(betas.ub, 2 * beta.hat.j - lb.j)
    }
  } else {
    betas.lb <- NULL
    betas.ub <- NULL
  }


  ###################################
  end_time <- Sys.time()
  result_time <- end_time - start_time

  return(list(etas = eta.B, cis = cis2,
              betas.lb = betas.lb, betas.ub = betas.ub,
              time = result_time))
}
