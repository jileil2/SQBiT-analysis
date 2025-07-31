QBiT_tune <- function(tau = .5, Y, C, X, P, B, Q2, 
                      nlambda = 10, new.nlambda = 10, lambda.scale = 100,
                      eps.abs = 1e-4, eps.rel = 1e-2,
                      eta.j1 = 1, eta.j2 = 1, var.j = FALSE,
                      lambda_start = 10 ^ (-4), lambda_end = 10 ^ 1) {
  
  start_time <- Sys.time()
  lambda_start <- log10(lambda_start)
  lambda_end <- log10(lambda_end)
  lambda <- 10 ^ (seq(lambda_start, lambda_end, length.out = nlambda))
  lambda <- as.matrix(lambda)
  
  cverr <- c()
  for (i in 1:10) {
    cverr[i] <- mean(cv.pred.QBiT(y = Y, C, X, B, Q2, P, lambda = lambda[i], var.j = var.j))
  }
  lambdac <- lambda[which.min(cverr)]
  
  new.lambda <- 10 ^ (seq(log10(lambdac / lambda.scale),
                          log10(lambdac * lambda.scale),
                          length.out = new.nlambda))
  
  cverr <- c()
  for (i in 1:10) {
    cverr[i] <- mean(cv.pred.QBiT(y = Y, C, X, B, Q2, P, lambda = new.lambda[i], var.j = var.j))
  }
  lambdac <- lambda[which.min(cverr)]
  end_time <- Sys.time()
  
  return(list(lambda = lambdac, time = end_time - start_time))
}
