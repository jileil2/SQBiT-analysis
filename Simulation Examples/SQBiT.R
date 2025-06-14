# smoothed QSVCM
SQBiT <- function(h, tau, Y, C, X, P, B, Q2, max.iter = 50,
                         lambda, eps.abs = 1e-4, eps.rel = 1e-2,
                         var.j = FALSE, zeta = 10, incr = 2, ie = FALSE, gamma.hat,
                         eta.j = .32, interval = FALSE, level = 0.05,
                         cond.est = FALSE, adaptive.h = FALSE, cc = .8,
                         delta = .01, gacv.compute = TRUE) {
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
  Z <- as.matrix(kr(X, BQ2, byrow = TRUE))
  W <- cbind(C, Z)
  cp.W <- eigenMapMatMult(t(W), W)
  
  alpha_all <- matrix(rep(0, (np * J + nq) * nl), ncol = nl)
  gacv_all <- rep(0, nl)
  df_all <- rep(0, nl)
  iter <- rep(0, nl)
  fit <- c()
  
  psi.eq <- function(r) {
    psi1 <- (eta.j * r - (tau - 1 / 2)) / (eta.j + 1 / (2 * h))
    psi2 <- -(tau - 1) / eta.j + r
    psi3 <- -tau / eta.j + r
    psi1.r <- psi1 * ifelse(abs(psi1) <= h, yes = 1, no = 0)
    psi2.r <- psi2 * ifelse(psi2 < -h, yes = 1, no = 0)
    psi3.r <- psi3 * ifelse(psi3 > h, yes = 1, no = 0)
    psi <- psi1.r + psi2.r + psi3.r
    return(psi)
  }
  
  et.1 <- Sys.time()
  
  time.admm <- 0
  time.gacv <- 0
  time.inv <- 0
  
  for(il in 1:nl){
    
    # time admm
    st.inverse <- Sys.time()
    
    #####
    kp <- cppkp(diag(lambda[il,]), P)
    zm1 <- matrix(0, nrow = nrow(kp), ncol = nq)
    zm2 <- matrix(0, nrow = nq, ncol = nrow(kp) + nq)
    kp <- rbind(zm2, cbind(zm1, kp))
    cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
    mat <- eigenMapMatMult(cp.inv, t(W))
    
    et.inverse <- Sys.time()
    
    time.inv <- time.inv + et.inverse - st.inverse
    
    st.l <- Sys.time()
    
    # old
    new.u <- (1 / eta.j) * Y
    if (ie == TRUE) {
      new.gamma <- gamma.hat
    } else {
      new.gamma <- eigenMapMatMult(mat, Y)
    }
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
      new.psi <- psi.eq(r)
      
      # equation (10)
      new.gamma <- eigenMapMatMult(mat, Y - new.psi - old.u)
      # equation (9)
      new.fit <- eigenMapMatMult(W, new.gamma)
      res.new <- Y - new.fit
      new.u <- old.u + (new.psi - res.new)
      
      # r = psi in the paper
      r.prime <- res.new - new.psi
      r.dual <- eta.j * (new.fit - old.fit) #
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
      } else if (iter.admm >= max.iter) {
        break
      } else {
        if (var.j == TRUE) {
          if (eta.j != old.eta.j) {
            st.inverse <- Sys.time()
            cp.inv <- Matrix::chol2inv(chol((2 / eta.j) * kp + cp.W))
            mat <- eigenMapMatMult(cp.inv, t(W))
            et.inverse <- Sys.time()
            time.inv <- time.inv + et.inverse - st.inverse
          }
        }
      }
    }
    
    et.l <- Sys.time()
    time.admm <- time.admm + et.l - st.l
    
    
    alpha_all[, il] <- as.matrix(new.gamma)
    
    # gacv time
    st.l.gacv <- Sys.time()
    
    if (gacv.compute == TRUE) {
      df <- sum_cpp(diag(eigenMapMatMult(cp.inv, cp.W)))
      gacv <- sum_rhotau(res.new, tau) / (n - df) 
      gacv_all[il] <- gacv
      df_all[il] <- df 
    } else {
      gacv_all <- NULL
      df <- NULL
    }
    
    
    
    #############################################################      
    iter[il] <- iter.admm
    fit[[il]] <- new.fit
    
    et.l.gacv <- Sys.time()
    time.gacv <- time.gacv + et.l.gacv - st.l.gacv
    
  }
  if (nrow(lambda) == 1) {
    if (gacv.compute == TRUE) {
      gcv <- gacv_all
    } 
    lambdac <- lambda
    alpha_hat <- alpha_all
    gamma <- alpha_hat
    fit <- new.fit
  } else {
    candidates <- which(iter < max.iter)
    
    if (gacv.compute == TRUE) {
      j <- which.min(gacv_all[candidates])
      gcv <- gacv_all[candidates][j]
    } 
    
    lambdac <- lambda[candidates, ][j, ]
    alpha_hat <- alpha_all[, candidates][, j]
    gamma <- alpha_hat
    iter <- iter[candidates][j]
    fit <- fit[candidates][j]
  }
  
  # time for computing beta 
  st.beta <- Sys.time()
  ###################################
  eta <- gamma[1:nq]
  
  et.beta <- Sys.time()
  
  if (interval == TRUE) {
    resid <- as.numeric(Y - eigenMapMatMult(W, gamma))
    
    if (cond.est[1] == FALSE) {
      
      if (adaptive.h == FALSE) {
        kd <- dunif(resid / h, -1, 1) / h
      } else {
        bounds <- quantile(resid, c(.1, .9))
        rr <- bounds[2] - bounds[1] 
        delta <- rr * (abs(tau - 0.5) + 0.5) * cc * n ^ (-1 / 3)
        kd <- dunif(resid / delta, -1, 1) / delta
      }
      
    } else {
      kd <- cond.est
    }
    
    
    sub.Z <- Z[which(kd > 0), ]
    sub.C <- C[which(kd > 0), ]
    kd <- kd[which(kd > 0)]
    
    Sigma.h <- diag(kd)
    Sigma.h.sqrt <- diag(sqrt(kd))
    Sigma.h.sqrt.inv <- diag(1 / sqrt(kd))
    
    
    # weighted projection matrix Z
    Z.w <- eigenMapMatMult(Sigma.h.sqrt, sub.Z)
    Z.w.prod <- eigenMapMatMult(t(Z.w), Z.w)
    Z.w.prod.inv <- AInv(Z.w.prod)
    Z.w.Z.w.prod.inv <- eigenMapMatMult(Z.w, Z.w.prod.inv)
    P.Z.w <- eigenMapMatMult(Z.w.Z.w.prod.inv, t(Z.w))
    
    # sandwich matrices
    C.w <- eigenMapMatMult(Sigma.h.sqrt, sub.C)
    C.w.proj <- eigenMapMatMult(diag(nrow(sub.Z)) - P.Z.w, C.w)
    Sigma.c.h <- eigenMapMatMult(t(C.w), C.w.proj) / n
    Sigma.c.h.inv <- AInv(Sigma.c.h)
    
    # computing center matrix
    C.proj <- eigenMapMatMult(diag(nrow(sub.Z)) - P.Z.w, C.w)
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
      ses <- c(ses, se.j)
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
  
  end_time <- Sys.time()
  result_time <- end_time - start_time  
  
  return(list(gamma = gamma, theta = theta, beta = beta, iter = iter, eta = eta, 
              time = result_time, gacv = gacv_all, c = c, delta = delta, df = df,
              time.decop = c(et.1 - st.1, time.inv, time.admm, time.gacv, et.beta - st.beta),
              cis = cis, ses = ses))
}
