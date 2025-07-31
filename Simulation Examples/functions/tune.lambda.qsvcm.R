tune.lambda.qsvcm <- function(tau = .5, y, C, X, P, B, Q2, 
                        nlambda = 10, new.nlambda = 10, max.iter = 50, 
                        lambda_start = 10 ^ (-4), var.j.1 = FALSE, var.j.2 = TRUE,
                        lambda_end = 10 ^ 1, eta.j1 = 1, eta.j2 = 1, 
                        zeta = 10, incr = 2,
                        b.initial = FALSE, eps.abs = 10 ^ (-4), eps.rel = 10 ^ (-2),
                        lambda.scale = 25) {
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
  
  n <- length(y)
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
  iter1 <- c()
  for(il in 1:nl){ 
    
    #####
    kp <- kronecker.prod(diag(lambda[il,]), P)
    if (is.null(ncol(C)) == FALSE) {
      zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
      zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
      kp <- rbind(zm2, cbind(zm1, kp))
    } 
    cp.inv <- Matrix::chol2inv(chol((2 / eta.j1) * kp + cp.W))
    mat <- eigenMapMatMult(cp.inv, t(W))
     
    # initials
    a <- 1 / (2 * eta.j1)
    # initialization
    old.u <- new.u <- (1 / eta.j1) * y 
    old.b <- new.b <- eigenMapMatMult(mat, y) 
    
    # v <- old.u - y + W %*% old.b - (2 * tau - 1) / (2 * eta)
    v <- -old.u + y - eigenMapMatMult(W, old.b) - (2 * tau - 1) / (2 * eta.j1)
    old.r <- new.r <- stoper(v, a)
    new.fit <- eigenMapMatMult(W, new.b)
    
    ###############
    
    ###############
    
    iter.admm <- 0 
    repeat{
      iter.admm <- iter.admm + 1
      
      old.u <- new.u
      old.b <- new.b
      old.r <- new.r
      old.fit <- new.fit
      
      # v <- old.u - y + W %*% old.b - (2 * tau - 1) / (2 * var.eta) #
      v <- -old.u + y - eigenMapMatMult(W, old.b) - (2 * tau - 1) / (2 * eta.j1) #
      new.r <- stoper(v, a)
      
      # equation (10)
      new.b <- eigenMapMatMult(mat, y - new.r - old.u)
      new.fit <- eigenMapMatMult(W, new.b)
      res.new <- y - new.fit
      # equation (9)
      new.u <- old.u + (new.r - res.new)
      
      # r = psi in the paper
      r.prime <- y - new.fit - new.r
      r.dual <- eta.j1 * (new.fit - old.fit) #
      # sqrt norm
      eps.prime <- sqrt(length(r.prime)) * eps.abs + 
        eps.rel * max(norm(new.r, type = "2"),
                      norm(new.fit, type = "2"),
                      norm(y, type = "2"))
      eps.dual <- sqrt(length(r.dual)) * eps.abs + 
        eps.rel * norm(new.u, type = "2")
      
      # ################################
      n.prime <- norm(r.prime, type = "2")
      n.dual <- norm(r.dual, type = "2")
      
      
      if (var.j.1 == TRUE) {
        old.eta.j1 <- eta.j1
        if(n.prime > zeta * norm(r.dual, type = "2")){
          eta.j1 <- incr * eta.j1  
          a <- 1 / (2 * eta.j1)
        } else if (zeta * n.prime < n.dual){
          eta.j1 <- eta.j1 / incr  
          a <- 1 / (2 * eta.j1)
        } else {
          eta.j1 <- eta.j1 
          a <- 1 / (2 * eta.j1)
        }
      }
      
      if ((n.prime < eps.prime) & (n.prime < eps.dual)){
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
    iter1 <- c(iter1, iter.admm)
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
  gacv_all <- rep(0, nl) 
  iter2 <- rep(0, nl) 
  
  for(il in 1:nl){  
    
    #####
    kp <- kronecker.prod(diag(lambda[il,]), P)
    if (is.null(ncol(C)) == FALSE) {
      zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
      zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
      kp <- rbind(zm2, cbind(zm1, kp))
    } 
    cp.inv <- Matrix::chol2inv(chol((2 / eta.j2) * kp + cp.W))
    mat <- eigenMapMatMult(cp.inv, t(W))
    
    # initials
    a <- 1 / (2 * eta.j2)
    # initialization
    old.u <- new.u <- (1 / eta.j2) * y
    old.b <- new.b <- eigenMapMatMult(mat, y) 
    # v <- old.u - y + W %*% old.b - (2 * tau - 1) / (2 * eta)
    v <- -old.u + y - eigenMapMatMult(W, old.b) - (2 * tau - 1) / (2 * eta.j2)
    old.r <- new.r <- stoper(v, a)
    new.fit <- eigenMapMatMult(W, new.b)
    
    ###############
    
    ###############
    
    iter.admm <- 0
    repeat{
      iter.admm <- iter.admm + 1
      
      old.u <- new.u
      old.b <- new.b
      old.r <- new.r
      old.fit <- new.fit
      
      # v <- old.u - y + W %*% old.b - (2 * tau - 1) / (2 * var.eta) #
      v <- -old.u + y - eigenMapMatMult(W, old.b) - (2 * tau - 1) / (2 * eta.j2) #
      new.r <- stoper(v, a)
      
      # equation (10)
      new.b <- eigenMapMatMult(mat, y - new.r - old.u)
      new.fit <- eigenMapMatMult(W, new.b)
      res.new <- y - new.fit
      # equation (9)
      new.u <- old.u + (new.r - res.new)
      
      # r = psi in the paper
      r.prime <- y - new.fit - new.r
      r.dual <- eta.j2 * (new.fit - old.fit) #
      # sqrt norm
      eps.prime <- sqrt(length(r.prime)) * eps.abs + 
        eps.rel * max(norm(new.r, type = "2"),
                      norm(new.fit, type = "2"),
                      norm(y, type = "2"))
      eps.dual <- sqrt(length(r.dual)) * eps.abs + 
        eps.rel * norm(new.u, type = "2")
      
      # ################################
      n.prime <- norm(r.prime, type = "2")
      n.dual <- norm(r.dual, type = "2")
      
      
      if (var.j.2 == TRUE) {
        old.eta.j2 <- eta.j2
        if(n.prime > zeta * norm(r.dual, type = "2")){
          eta.j2 <- incr * eta.j2  
          a <- 1 / (2 * eta.j2)
        } else if (zeta * n.prime < n.dual){
          eta.j2 <- eta.j2 / incr  
          a <- 1 / (2 * eta.j2)
        } else {
          eta.j2 <- eta.j2 
          a <- 1 / (2 * eta.j2)
        }
      }
      
      if ((n.prime < eps.prime) & (n.prime < eps.dual)){
        break
      } else if (iter.admm >= max.iter) {
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
    alpha_all[, il] <- as.matrix(new.b)
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
  lambdac <- lambda[candidates, ][j, ]
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
  
  ###################################
  end_time <- Sys.time()
  result_time <- end_time - start_time  
   
  list(gamma = gamma, theta = theta, beta = beta, eta = eta, 
       lambdac = lambdac, time = result_time, 
       iter = c(iter1, iter2))
}
