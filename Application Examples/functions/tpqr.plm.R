# Tensor Product Spline Regression, He and Shi (1996)
tpqr.plm <- function(Y, C, X, S, d, tau, c0 = 2, c1 = 4.5, nc = 6) {
  start_time <- Sys.time()
  # load packages
  require("splines")
  require("quantreg")
  # setup
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(C)
  s1 <- S[, 1]
  s2 <- S[, 2]
  
  # knots
  ord <- n ^ (1 / (2 * d + 2))
  c <- seq(from = c0, to = c1, length.out = nc)
  nknots <- round(c * ord)
  
  # collect AIC
  AICs <- matrix(0, nc, nc)
  mods <- vector('list', nc * nc)
  Bss <- vector('list', nc * nc)
  knots <- vector('list', nc * nc)
  time <- matrix(0, nc, nc)
  # models
  for (a in 1:nc) {
    for (b in 1:nc) {
      start_time <- Sys.time()
      
      # basis construction
      k1 <- nknots[a]
      k2 <- nknots[b]
      knots.1 <- seq(from = min(s1), to = max(s1), length.out = k1)
      knots.2 <- seq(from = min(s2), to = max(s2), length.out = k2)
      B1 <- bs(s1, df = k1, degree = d)
      B2 <- bs(s2, df = k2, degree = d)
      
      # basis
      basis <- do.call(rbind, lapply(1:n, function(i) kronecker(B1[i,], B2[i,])))
      # X * basis
      XB <- c()
      for (j in 1:p) {
        XB <- cbind(XB, t(sapply(1:n, function(i) X[i, j] * basis[i, ])))
      }
      
      # design
      design <- cbind(C, XB) 
      # QR
      mod <- rq(Y ~ design - 1, tau = tau)
      end_time <- Sys.time()
      # AIC
      AICs[a, b] <- AIC(mod)
      mods[[nc * (a - 1) + b]] <- mod
      Bss[[nc * (a - 1) + b]] <- basis
      time[a, b] <- end_time - start_time
      knots[[nc * (a - 1) + b]] <- list(knots.1, knots.2)
    }
  }
  
  # minimizer
  sets <- as.numeric(which(AICs == min(AICs), arr.ind = TRUE))
  mod.best <- mods[[nc * (sets[1] - 1) + sets[2]]]
  Bs.best <- Bss[[nc * (sets[1] - 1) + sets[2]]]
  time.best <- time[sets[1], sets[2]]
  knots.best <- knots[[nc * (sets[1] - 1) + sets[2]]]
  
  # estimates
  gamma <- coef(mod.best)
  eta <- gamma[1:q]
  beta <- Bs.best %*% matrix(gamma[-(1:q)], ncol = p)
  
  end_time <- Sys.time()
  result_time <- end_time - start_time  
  return(list(eta = eta, beta = beta, knots = knots.best,
              time = result_time, 
              gamma = gamma[-(1:q)],
              mod = mod.best))
}




