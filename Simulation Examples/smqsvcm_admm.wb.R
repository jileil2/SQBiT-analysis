# smoothed QSVCM (wild bootstrap)
smqsvcm_admm.wb <- function(h, tau, Y, C, X, P, B, Q2, max.iter = 50,
                            lambda, eps.abs = 1e-4, eps.rel = 1e-2, kernel = 'unif',
                            var.j = FALSE, zeta = 10, incr = 2, eta.j = .32,
                            BQ2.eva, compute.vc = TRUE,
                            Br = 200, level = 0.05, gamma.hat, biascorr = FALSE) {
  start_time <- Sys.time() 
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
  
  psi.eq <- function(r, eta.j, kernel = kernel) {
    
    if (kernel == "unif") {
      psi1 <- (eta.j * r - (tau - 1 / 2)) / (eta.j + 1 / (2 * h))
      psi2 <- -(tau - 1) / eta.j + r
      psi3 <- -tau / eta.j + r
      psi1.r <- psi1 * ifelse(abs(psi1) <= h, yes = 1, no = 0)
      psi2.r <- psi2 * ifelse(psi2 < -h, yes = 1, no = 0)
      psi3.r <- psi3 * ifelse(psi3 > h, yes = 1, no = 0)
      psi <- psi1.r + psi2.r + psi3.r
    } else if (kernel == "norm") {
      psi <- sapply(r, function(ri) nleqslv(0, function(u) tau - pnorm(-u / h) + eta.j * (u - ri))$x)
    }
    return(psi)
  } 
  
  kp <- cppkp(diag(lambda[1,]), P)
  zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
  zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
  kp <- rbind(zm2, cbind(zm1, kp))
  cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
  mat <- eigenMapMatMult(cp.inv, t(W))
  
  # compute residual
  resid <- as.numeric(Y - eigenMapMatMult(W, gamma.hat))
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
    Y.b <- eigenMapMatMult(W, gamma.hat) + xis * abs(resid + bias)
    
    # old
    new.u <- (1 / eta.j) * Y.b
    new.gamma <- gamma.hat
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
      new.psi <- psi.eq(r, kernel = kernel, eta.j = eta.j)
      
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
  eta.hat <- gamma.hat[1:nq]
  cis2 <- c()
  cis2.vc <- c()
  for (j in 1:nq) {
    ci <- quantile(eta.B[,j], c(level / 2, 1 - level / 2))
    cis2 <- rbind(cis2, c(2 * eta.hat[j] - ci[2], 2 * eta.hat[j] - ci[1]))
  }
  
  if (compute.vc == TRUE) {
    # VC
    gamma.hat <- gamma.hat[-(1:nq)]
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
  
  return(list(etas = eta.B, cis2 = cis2, 
              betas.lb = betas.lb, betas.ub = betas.ub,
              result_time = result_time))
}
