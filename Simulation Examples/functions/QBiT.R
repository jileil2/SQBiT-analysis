######################################################################
# Myungjin Kim and Lily Wang  
# 01/11/2019
#
# Input Arguments:
#	(1) y: response variable
#	(2) X: explanatory variable
#	(3) BQ2: bivariate spline basis times Q2 matrix in the paper, that is B%*%Q2.
#	(4) lambda: smoothness penalty parameter candidates
# (5) tau:
# (6) eta:
# (7) eps.abs:
# (8) eps.rel:
#
# Output Arguments:
#	(1) theta_hat: estimator for bivariate spine coefficients 
#	(2) lambdac: roughness penalty parameter chosen
#	(3) gacv: generalized cross-validation criterion based on chosen estimator
# (4) df: effective degrees of freedom
# (5) theta: estimated parameter theta for the negative binomial (only GSVCM_est_nb function)
######################################################################
# updated:
# v<-old.u-y+W%*%old.b-(2*tau-1)/(2*eta)
# for simulation eta<-1.25
# new.b[1:J]<-new.b[1:J]-Finv

# updated 2/27/2020:
# theta, beta is updated and argument BQ2 -> B and Q2.


#updated 3/9/2020:
#   lambda<-as.matrix(lambda)
# nl<-nrow(lambda)
#   if(ncol(lambda)<-<-1){
#    lambda<-matrix(rep(lambda,times<-np),nl,np)
#} # ADDED
# lambda[il]*D --> kronecker(diag(lambda[il,]),P)
# Argument B and Q2 used instead of BQ2 
# Output (fucntion) is different
# np and nl should be before if...
QBiT <- function(Y, X, B, Q2, P, C, lambda, tau = 0.5, eta.j = 1.0,
                      var.j = FALSE, zeta = 10, 
                      eps.abs = 10 ^ (-4), eps.rel = 10 ^ (-2),
                      method = "ADMM", step = 1, # step <- 1: non-varying, 2 <- varying
                      incr = 2, alpha = 1.5,    ###### tuning for varying step/ those might be arguments later
                      s.method = "GACV", max.iter = 50,
                      b.initial = FALSE, screen = FALSE,
                      gacv.compute = TRUE){
  start_time <- Sys.time()
  
  st.1 <- Sys.time()
  
  P <- as.matrix(P)
  Q2 <- as.matrix(Q2)
  B <- as.matrix(B)
  
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
  
  alpha_all <- matrix(rep(0, (np * J + nq) * nl), ncol = nl)
  gacv_all <- rep(0, nl)
  df_all <- rep(0, nl) 
  iter <- rep(0, nl)
  
  et.1 <- Sys.time()
  
  time.admm <- 0
  time.gacv <- 0
  time.inv <- 0

  for(il in 1:nl){ 
    
    # time admm
    st.inverse <- Sys.time()
    
    #####
    kp <- kronecker.prod(diag(lambda[il,]), P)
    zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
    zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
    kp <- rbind(zm2, cbind(zm1, kp))
    cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
    mat <- eigenMapMatMult(cp.inv, t(W))
    
    et.inverse <- Sys.time()
    
    time.inv <- time.inv + et.inverse - st.inverse
    
    st.l <- Sys.time()
    
    # initials
    a <- 1 / (2 * eta.j)
    # initialization
    old.u <- new.u <- (1 / eta.j) * Y
    if (b.initial == FALSE){ # initial for b
      #####################
      old.b <- new.b <- eigenMapMatMult(mat, Y) 

    } else {
      old.b <- new.b <- b.initial
    }
    # v <- old.u - y + W %*% old.b - (2 * tau - 1) / (2 * eta)
    v <- -old.u + Y - eigenMapMatMult(W, old.b) - (2 * tau - 1) / (2 * eta.j)
    old.r <- new.r <- stoper(v, a)
    new.fit <- eigenMapMatMult(W, new.b)

    ###############

    ###############
    
    iter.admm <- 0
    all.var.eta <- c() ##
    repeat{
      iter.admm <- iter.admm + 1
      
      old.u <- new.u
      old.b <- new.b
      old.r <- new.r
      old.fit <- new.fit

      # v <- old.u - y + W %*% old.b - (2 * tau - 1) / (2 * var.eta) #
      v <- -old.u + Y - eigenMapMatMult(W, old.b) - (2 * tau - 1) / (2 * eta.j) #
      new.r <- stoper(v, a)
      
      if (method == "ADMM"){
        # equation (10)
        new.b <- eigenMapMatMult(mat, Y - new.r - old.u)
        new.fit <- eigenMapMatMult(W, new.b)
        res.new <- Y - new.fit
        # equation (9)
        new.u <- old.u + (new.r - res.new)
        
      } else if(method == "ORADMM"){
        new.b <- eigenMapMatMult(mat, y - (alpha * new.r - (1 - alpha) * (old.fit - Y)) - old.u)
        new.fit <- eigenMapMatMult(W, new.b)
        new.u <- old.u + (alpha * new.r - (1 - alpha) * (old.fit - y) + new.fit - Y)
      }

      # r = psi in the paper
      r.prime <- Y - new.fit - new.r
      r.dual <- eta.j * (new.fit - old.fit) #
      # sqrt norm
      eps.prime <- sqrt(length(r.prime)) * eps.abs + 
        eps.rel * max(norm(new.r, type = "2"),
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
          a <- 1 / (2 * eta.j)
        } else if (zeta * n.prime < n.dual){
          eta.j <- eta.j / incr  
          a <- 1 / (2 * eta.j)
        } else {
          eta.j <- eta.j 
          a <- 1 / (2 * eta.j)
        }
      }
      
      if ((n.prime < eps.prime) & (n.prime < eps.dual)){
        break
      } else if (iter.admm >= max.iter) {
        break
      } else {
        if (var.j == TRUE) {
          if (old.eta.j != eta.j) {
            cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
            mat <- eigenMapMatMult(cp.inv, t(W))
          }
        }
      }
    }
    
    et.l <- Sys.time()
    time.admm <- time.admm + et.l - st.l
    
    # gacv time
    st.l.gacv <- Sys.time()
    
    
    
    alpha_all[, il] <- as.matrix(new.b)
    
    if (gacv.compute == TRUE) {
      df <- sum_cpp(diag(eigenMapMatMult(cp.inv, cp.W)))
      gacv <- sum_rhotau(res.new, tau) / (n - df) 
      gacv_all[il] <- gacv
      df_all[il] <- df 
    } else {
      gacv_all <- gacv <- NULL
    }
    
    #############################################################   
    iter[il] <- iter.admm
    
    et.l.gacv <- Sys.time()
    time.gacv <- time.gacv + et.l.gacv - st.l.gacv
    
  }
  if (nrow(lambda) == 1) {
    if (gacv.compute == TRUE) {
      gacv <- gacv_all
    } 
    lambdac <- lambda
    alpha_hat <- alpha_all
    gamma <- alpha_hat
  } else {
    candidates <- which(iter < max.iter)
    
    if (gacv.compute == TRUE) {
      j <- which.min(gacv_all[candidates])
      gacv <- gacv_all[candidates][j]
    } 
    
    lambdac <- lambda[candidates, ][j, ]
    alpha_hat <- alpha_all[, candidates][, j]
    gamma <- alpha_hat
    iter <- iter[candidates][j]
  }
  
  # time for computing beta 
  st.beta <- Sys.time()
  ###################################
  eta <- gamma[1:nq]
  gamma <- gamma[-(1:nq)]
  gamma.mtx <- matrix(gamma, J, np)
  theta <- eigenMapMatMult(Q2, gamma.mtx)
  beta <- eigenMapMatMult(B, theta)
  # 
  et.beta <- Sys.time()
  
  ###################################
  end_time <- Sys.time()
  result_time <- end_time - start_time  


  list(gamma = gamma, theta = theta, beta = beta, eta = eta,
       sic = gacv, lambdac = lambdac, time = result_time, 
       iter = iter,
       time.decop = c(et.1 - st.1, time.inv, time.admm, time.gacv, et.beta - st.beta))
}

