tune.lambda <- function(tau = .5, Y, C, X, P, B, Q2, h, var.j.1 = FALSE, var.j.2 = TRUE,
                        nlambda = 10, new.nlambda = 10, lambda.scale = 100,
                        kernel = c('unif', 'norm'), 
                        eps.abs = 1e-4, eps.rel = 1e-2, zeta = 10, incr = 2,
                        max.iter = 500, eta.j1 = .45, eta.j2 = 1, lambda_start = 10 ^ (-4),
                        lambda_end = 10 ^ 1, interval = FALSE, level = .05) {
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
  Z <- as.matrix(kr(X, BQ2, byrow = TRUE))
  W <- cbind(C, Z)
  cp.W <- eigenMapMatMult(t(W), W)
   
  gacv_all <- rep(0, nl) 
  iter1 <- rep(0, nl) 
  
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
  
  for(il in 1:nl){
     
    
    #####
    kp <- cppkp(diag(lambda[il,]), P)
    zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
    zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
    kp <- rbind(zm2, cbind(zm1, kp))
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
      
      
      if (var.j.1 == TRUE) {
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
        if (var.j.1 == TRUE) { 
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
  lambdac <- lambda[candidates, ][j, ]
  
  
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
      kp <- cppkp(diag(lambda[il,]), P)
      kp <- rbind(zm2, cbind(zm1, kp))
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
      
      
      if (var.j.2 == TRUE) {
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
        if (var.j.2 == TRUE) { 
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
    lambdac <- lambda[candidates, ][j, ]
    alpha_hat <- alpha_all[, candidates][, j]
    gamma <- alpha_hat 
    ###################################
    eta <- gamma[1:nq]
    gamma.mtx <- matrix(gamma[-(1:nq)], J, np)
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
    end_time <- Sys.time()
    result_time <- end_time - start_time  
  
  return(list(gamma = gamma, theta = theta, beta = beta, iter = c(iter1, iter2), eta = eta, 
              lambdac = lambdac, lambda = lambda, result_time = result_time, gacv = gacv_all, 
              cis = cis, ses = ses, gacv.min = gcv))
}
