#' Smoothed Quantile Bivariate Triangulation using penalization parameters selected by Generalized Approximate Cross Validation (SQBiT-GACV)
#'
#' This function performs the point estimation of quantile spatial model (QSM) with penalization parameters selected by GACV using SQBiT. It returns estimates for the constant coefficients, varying coefficients, and confidence intervals for the constant coefficients. SQBiT employs convolution smoothing based on uniform kernel and estimates the asymptotic covariance matrix using kernel density estimation based on uniform kernel.
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
#' The current implementation assumes the same penalization parameters for different varying coefficient functions for simplicity. The above procedure uses lambda_start = 1e-4, lambda_end = 1e1, and lambda.scale = 100. But they can be modified according to the user's needs. Note that at tail quantile levels, penalization parameter selection using GACV may be unstable for estimation. For this, we recommend using a narrower search by specifying lambda_start = 1e-1 and lambda.scale = 10.
#'
#' @param h smoothing bandwidth
#' @param tau quantile level
#' @param Y response variable
#' @param C covariates with constant effects
#' @param X covariates with varying coefficients
#' @param B Bernstein basis polynoimials. An object from basis() in package BPST.
#' @param Q2 Q2 matrix from the QR decomposition of the smoothness matrix of triangles
#' @param P Penalty matrix. See examples. It is based on output from basis() in package BPST.
#' @param interval whether to compute asymptotic confidence intervals for constant coefficients (TRUE/FALSE)
#' @param level significance level
#' @param max.iter maximum iterations for ADMM to converge
#' @param eta.j1 step size of the first round. Set it to be 1 if var.j = TRUE
#' @param eta.j2 step size of the second round. Set it to be 1 if var.j = TRUE
#' @param var.j whether to vary the step size in the ADMM
#' @param eps.abs first convergence criterion of ADMM
#' @param eps.rel second convergence criterion of ADMM
#' @param incr factor of increment for step size if stepsize is set to vary
#' @param zeta criterion for the step size to increment if stepsize is set to vary
#' @param lambda_start the lower end of the range for the search of penalization parameters in the first round
#' @param lambda_end the upper end of the range for the search of penalization parameters in the first round
#' @param nlambda the number of penalization parameters to be searched in the first round
#' @param new.nlambda the number of penalization parameters to be searched in the second round
#' @param lambda.scale the scale of the search of penalization parameters around the penalization parameters selected in the first round. This is for the second round.
#' @return A list with components:
#' \describe{
#'   \item{\code{lambdac}}{The selected penalization parameter.}
#'   \item{\code{gacv.min}}{The GACV of the estimation using the selected penalization parameter.}
#'   \item{\code{lambdas}}{The set of penalization parameters in the second round.}
#'   \item{\code{gacv}}{The GACVs of the estimation using the penalization parameters in the second round.}
#'   \item{\code{gamma}}{Spline coefficients.}
#'   \item{\code{theta}}{Basis coefficients (\code{Q2 \%*\% gamma}).}
#'   \item{\code{beta}}{Estimated varying coefficients.}
#'   \item{\code{eta}}{Estimated constant coefficients.}
#'   \item{\code{cis}}{Asymptotic intervals for the constant coefficients.}
#'   \item{\code{ses}}{Estimated asymptotic standard errors.}
#'   \item{\code{iter}}{Number of iterations used for convergence for the first and second rounds.}
#'   \item{\code{time}}{Computational time used.}
#'   \item{\code{gacv}}{The GACV of the estimated QSM; returns multiple values if multiple penalization parameters are used.}
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
#' mod <- SQBiT_gacv(Y = Y, X = X, C = C, P = P, B = B, Q2 = Q2, tau = tau)
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
SQBiT_gacv <- function(h = 7 * tau * (1 - tau) * ((ncol(C) + dim(Q2)[2] * ncol(X) + log(length(Y))) / length(Y)) ^ (2/5),
                       tau = .5, Y, C, X, P, B, Q2, var.j = FALSE,
                       nlambda = 10, new.nlambda = 10, lambda.scale = 100,
                       eps.abs = 1e-4, eps.rel = 1e-2, zeta = 10, incr = 2,
                       max.iter = 500, eta.j1 = .32, eta.j2 = .32, lambda_start = 10 ^ (-4),
                       lambda_end = 10 ^ 1, interval = FALSE, level = .05) {
  start_time <- Sys.time()

  P <- as.matrix(P)
  Q2 <- as.matrix(Q2)
  B <- as.matrix(B)
  if (is.null(ncol(C)) == FALSE) {
    C <- as.matrix(C)
  }

  if (tau < 0 | tau > 1) {
    warning("Quantile level should be between 0 and 1. The default quantile level of 0.5 is used.")
    tau <- 0.5
  }

  if (tau <= 0.10 | tau >= 0.90) {
    warning("Estimation at tail quantile levels using GACV is more challenging. Please consider narrowing the search of penalization parameters.")
  }

  n <- length(Y)
  np <- ncol(X)
  if (is.null(ncol(C)) == FALSE) {
    nq <- ncol(C)
  } else {
    nq <- 0
  }
  J <- dim(Q2)[2]

  lambda_start <- log10(lambda_start)
  lambda_end <- log10(lambda_end)
  lambda <- 10 ^ (seq(lambda_start, lambda_end, length.out = nlambda))
  lambda <- as.matrix(lambda)
  nl <- nrow(lambda)
  if (ncol(lambda) == 1) {
    lambda <- matrix(rep(lambda, times = np), nl, np)
  }

  # nl <- dim(lambda)[1]
  BQ2 <- eigenMapMatMult(B, Q2)
  W <- as.matrix(kr(X, BQ2, byrow = TRUE))
  if (is.null(ncol(C)) == FALSE) {
    W <- cbind(C, W)
  }
  cp.W <- eigenMapMatMult(t(W), W)

  gacv_all <- rep(0, nl)
  iter1 <- rep(0, nl)

  psi.eq <- function(r, eta.j, kernel = kernel) {
    psi1 <- (eta.j * r - (tau - 1 / 2)) / (eta.j + 1 / (2 * h))
    psi2 <- -(tau - 1) / eta.j + r
    psi3 <- -tau / eta.j + r
    psi1.r <- psi1 * ifelse(abs(psi1) <= h, yes = 1, no = 0)
    psi2.r <- psi2 * ifelse(psi2 < -h, yes = 1, no = 0)
    psi3.r <- psi3 * ifelse(psi3 > h, yes = 1, no = 0)
    psi <- psi1.r + psi2.r + psi3.r
    return(psi)
  }

  for(il in 1:nl){


    #####
    if (length(lambda[il,]) == 1) {
      kp <- lambda[il,] * P
    } else {
      kp <- cppkp(diag(lambda[il,]), P)
    }
    if (is.null(ncol(C)) == FALSE) {
      zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
      zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
      kp <- rbind(zm2, cbind(zm1, kp))
    }
    cp.inv <- Matrix::chol2inv(chol((2 / eta.j1) * kp + cp.W))
    mat <- eigenMapMatMult(cp.inv, t(W))

    # old
    new.u <- (1 / eta.j1) * Y
    new.gamma <- eigenMapMatMult(mat, Y)
    new.fit <- eigenMapMatMult(W, new.gamma)


    ###############

    ###############

    iter.admm <- 0
    repeat{
      iter.admm <- iter.admm + 1

      # old
      old.u <- new.u
      old.gamma <- new.gamma
      old.fit <- new.fit


      # start updating
      r <- Y - old.fit - old.u
      new.psi <- psi.eq(r, eta.j = eta.j1, kernel = kernel)

      # equation (10)
      new.gamma <- eigenMapMatMult(mat, Y - new.psi - old.u)
      # equation (9)
      new.fit <- eigenMapMatMult(W, new.gamma)
      res.new <- Y - new.fit
      new.u <- old.u + (new.psi - res.new)

      # r = psi in the paper
      r.prime <- res.new - new.psi
      r.dual <- eta.j1 * (new.fit - old.fit) #
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
        old.eta.j1 <- eta.j1
        if(n.prime > zeta * norm(r.dual, type = "2")){
          eta.j1 <- incr * eta.j1
        } else if (zeta * n.prime < n.dual){
          eta.j1 <- eta.j1 / incr
        } else {
          eta.j1 <- eta.j1
        }
      }


      if((n.prime < eps.prime) & (n.dual < eps.dual)){
        break
      } else if (iter.admm >= max.iter) {
        break
      } else {
        if (var.j == TRUE) {
          if (old.eta.j1 != eta.j1) {
            cp.inv <- Matrix::chol2inv(chol((2 / eta.j1) * kp + cp.W))
            mat <- eigenMapMatMult(cp.inv, t(W))
          }
        }
      }
    }
    # calculate gcv / new.b -> as.matrix
    #### ADMM df
    # kronecker(diag(lambda[il, ]), P) = D_{Lambda}
    df <- sum_cpp(diag(eigenMapMatMult(cp.inv, cp.W)))
    # sic <- log(mean(rhotau(y - yhat, tau))) + (1 / 2) * (1 / length(y)) * df * log(length(y))
    # gcv <- 2 * log(mean(rhotau(y - yhat, tau))) - 2 * log(1 - df / length(y))

    #############################################################
    # psi_hat <- rhotau(y - yhat, tau) / (y - yhat)
    # weight1 <- psi_hat / (2 * (y - yhat))
    # weight1.max <- diag(as.vector(weight1))
    # Hmat <- W %*% solve(t(W) %*% weight1.max %*% W + lambda[il] * D) %*% t(W) %*% weight1.max
    # df <- sum(diag(Hmat))

    # equation (11)
    gacv <- sum_rhotau(res.new, tau) / (n - df)
    #############################################################
    gacv_all[il] <- gacv
    iter1[il] <- iter.admm
  }

  candidates <- which(iter1 < max.iter)
  j <- which.min(gacv_all[candidates])
  if (ncol(lambda) == 1) {
    lambdac <- lambda[candidates][j]
  } else {
    lambdac <- lambda[candidates, ][j, ]
  }


    # nlambda #
    new.lambda <- 10 ^ (seq(log10(lambdac[1] / lambda.scale),
                            log10(lambdac[1] * lambda.scale),
                            length.out = new.nlambda))
    lambda <- as.matrix(new.lambda)
    nl <- nrow(lambda)
    if (ncol(lambda) == 1) {
      lambda <- matrix(rep(lambda, times = np), nl, np)
    }

    alpha_all <- matrix(rep(0, (np * J + nq) * nl), ncol = nl)
    iter2 <- rep(0, nl)
    gacv_all <- rep(0, nl)
    for(il in 1:nl){

      #####
      if (length(lambda[il,]) == 1) {
        kp <- lambda[il,] * P
      } else {
        kp <- cppkp(diag(lambda[il,]), P)
      }
      if (is.null(ncol(C)) == FALSE) {
        zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
        zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
        kp <- rbind(zm2, cbind(zm1, kp))
      }
      cp.inv <- Matrix::chol2inv(chol((2 / eta.j2) * kp + cp.W))
      mat <- eigenMapMatMult(cp.inv, t(W))

      # old
      new.u <- (1 / eta.j2) * Y
      new.gamma <- eigenMapMatMult(mat, Y)
      new.fit <- eigenMapMatMult(W, new.gamma)


      ###############
    iter.admm <- 0
    repeat{
      iter.admm <- iter.admm + 1

      # old
      old.u <- new.u
      old.gamma <- new.gamma
      old.fit <- new.fit


      # start updating
      r <- Y - old.fit - old.u
      new.psi <- psi.eq(r, eta.j = eta.j2, kernel = kernel)

      # equation (10)
      new.gamma <- eigenMapMatMult(mat, Y - new.psi - old.u)
      # equation (9)
      new.fit <- eigenMapMatMult(W, new.gamma)
      res.new <- Y - new.fit
      new.u <- old.u + (new.psi - res.new)

      # r = psi in the paper
      r.prime <- res.new - new.psi
      r.dual <- eta.j2 * (new.fit - old.fit) #
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
        old.eta.j2 <- eta.j2
        if(n.prime > zeta * norm(r.dual, type = "2")){
          eta.j2 <- incr * eta.j2
        } else if (zeta * n.prime < n.dual){
          eta.j2 <- eta.j2 / incr
        } else {
          eta.j2 <- eta.j2
        }
      }


      if((n.prime < eps.prime) & (n.dual < eps.dual)){
        break
      }  else if (iter.admm >= max.iter) {
        break
      } else {
        if (var.j == TRUE) {
          if (old.eta.j2 != eta.j2) {
            cp.inv <- Matrix::chol2inv(chol((2 / eta.j2) * kp + cp.W))
            mat <- eigenMapMatMult(cp.inv, t(W))
          }
        }
      }
    }
    alpha_all[, il] <- as.matrix(new.gamma)
    # calculate gcv / new.b -> as.matrix
    #### ADMM df
    # kronecker(diag(lambda[il, ]), P) = D_{Lambda}
    df <- sum_cpp(diag(eigenMapMatMult(cp.inv, cp.W)))
    # sic <- log(mean(rhotau(y - yhat, tau))) + (1 / 2) * (1 / length(y)) * df * log(length(y))
    # gcv <- 2 * log(mean(rhotau(y - yhat, tau))) - 2 * log(1 - df / length(y))

    #############################################################
    # psi_hat <- rhotau(y - yhat, tau) / (y - yhat)
    # weight1 <- psi_hat / (2 * (y - yhat))
    # weight1.max <- diag(as.vector(weight1))
    # Hmat <- W %*% solve(t(W) %*% weight1.max %*% W + lambda[il] * D) %*% t(W) %*% weight1.max
    # df <- sum(diag(Hmat))

    # equation (11)

    gacv <- sum_rhotau(res.new, tau) / (n - df)
    #############################################################
    gacv_all[il] <- gacv
    iter2[il] <- iter.admm
  }

    candidates <- which(iter2 < max.iter)
    j <- which.min(gacv_all[candidates])
    gcv <- gacv_all[candidates][j]
    if (ncol(lambda) == 1) {
      lambdac <- lambda[candidates][j]
    } else {
      lambdac <- lambda[candidates, ][j, ]
    }
    alpha_hat <- alpha_all[, candidates][, j]
    gamma <- alpha_hat


    ###################################
    if (is.null(ncol(C)) == FALSE) {
      eta <- gamma[1:nq]
      gamma.mtx <- matrix(gamma[-(1:nq)], J, np)
    } else {
      eta <- NULL
      gamma.mtx <- matrix(gamma, J, np)
    }
    theta <- eigenMapMatMult(Q2, gamma.mtx)
    beta <- eigenMapMatMult(B, theta)

    if (interval == TRUE) {
      resid <- as.numeric(Y - eigenMapMatMult(W, gamma))
      if (kernel == 'norm') {
        kd <- dnorm(resid / h) / h
      } else if (kernel == 'unif') {
        kd <- dunif(resid / h, -1, 1) / h
      }
      index <- which(kd > 0)
      sub.C <- C[index, ]
      sub.Z <- Z[index, ]
      resid <- resid[index]
      sn <- length(index)

      if (kernel == 'norm') {
        Sigma.h <- diag(dnorm(resid / h) / h)
        Sigma.h.sqrt <- diag(sqrt(dnorm(resid / h) / h))
        Sigma.h.sqrt.inv <- diag(1 / sqrt(dnorm(resid / h) / h))
      } else if (kernel == 'unif') {
        Sigma.h <- diag(dunif(resid / h, -1, 1) / h)
        Sigma.h.sqrt <- diag(sqrt(dunif(resid / h, -1, 1) / h))
        Sigma.h.sqrt.inv <- diag(1 / sqrt(dunif(resid / h, -1, 1) / h))
      }


      # weighted projection matrix Z
      Z.w <- eigenMapMatMult(Sigma.h.sqrt, sub.Z)
      Z.w.prod <- eigenMapMatMult(t(Z.w), Z.w)
      Z.w.prod.inv <- AInv(Z.w.prod)
      Z.w.Z.w.prod.inv <- eigenMapMatMult(Z.w, Z.w.prod.inv)
      P.Z.w <- eigenMapMatMult(Z.w.Z.w.prod.inv, t(Z.w))

      # sandwich matrices
      C.w <- eigenMapMatMult(Sigma.h.sqrt, sub.C)
      C.w.proj <- eigenMapMatMult(diag(sn) - P.Z.w, C.w)
      Sigma.c.h <- eigenMapMatMult(t(C.w), C.w.proj) / n
      Sigma.c.h.inv <- AInv(Sigma.c.h)

      # computing center matrix
      C.proj <- eigenMapMatMult(diag(sn) - P.Z.w, C.w)
      C.proj <- eigenMapMatMult(Sigma.h.sqrt.inv, C.proj)
      Sigma.c.e.h <- eigenMapMatMult(t(C.proj), C.proj) / n

      # final matrix
      Xi.n.h <- tau * (1 - tau) * eigenMapMatMult(Sigma.c.h.inv, Sigma.c.e.h)
      Xi.n.h <- eigenMapMatMult(Xi.n.h, Sigma.c.h.inv)

      # confidence intervals
      cv.n <- qnorm(1 - level / 2)
      cis <- c()
      ses <- c()
      for (j in 1:nq) {
        eta.hat.j <- eta[j]
        se.j <-  Xi.n.h[j, j] / n
        ci <- c(eta.hat.j - cv.n * sqrt(se.j), eta.hat.j + cv.n * sqrt(se.j))
        ses <- c(ses, sqrt(se.j))
        cis <- rbind(cis, ci)
      }
    } else {
      cis <- NULL
      ses <- NULL
    }

    ###################################
    gamma <- gamma[-(1:nq)]
    gamma.mtx <- matrix(gamma, J, np)
    theta <- eigenMapMatMult(Q2, gamma.mtx)
    beta <- eigenMapMatMult(B, theta)

    ###################################
    end_time <- Sys.time()
    result_time <- end_time - start_time

  return(list(gamma = gamma, theta = theta, beta = beta, iter = c(iter1, iter2), eta = eta,
              lambdac = lambdac, lambdas = lambda, time = result_time, gacv = gacv_all,
              cis = cis, ses = ses, gacv.min = gcv))
}
