#' k-fold cross-validation MSPE for spatially varying coefficient regression
#'
#' This function implements k-fold cross-validation MSPE for spartially varying coefficient regression, and returns the mean squared prediction error.
#'
#' @importFrom BPST basis
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param S The cooridinates of dimension \code{n} by two. Each row is the coordinates of an observation.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 2.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param lambda The vector of the candidates of penalty parameter -- default is grid points of 10 to the power of a sequence from -6 to 6 by 0.5.
#' \cr
#' @param nfolds The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' \cr
#' @param initial The seed used for cross-validation -- default is 666
#' \cr
#' @return The mean squared prediction error (MSPE) based on k-fold cross-validation
#' @details
#'
#'
#' @export
#'
cv.pred.SQBiT = function(y, C, X, B, Q2, P, d = 3, r = 1, lambda, nfold = 5, initial = 2024,
                          tau = 0.50, eta.j = 0.32, h,
                          var.j = FALSE) {

  ####################################################
  n <- length(y)
  sfold <- round(n / nfold)
  set.seed(initial)
  Test <- sample(1:n)
  cv.error <- c()
  
  for(ii in 1:nfold){
    if(ii < nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):(ii * sfold)])
    }
    if(ii == nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):n])
    }
    Train.set = setdiff(1:n, Test.set)

    # consider univariate case / after using an index, it becomes a vector again above.
    if(is.vector(X) == 1){
      
      C.test <- as.matrix(C[Test.set, drop = FALSE])
      C.train <- as.matrix(C[Train.set, drop = FALSE])
      
      X.test <- as.matrix(X[Test.set])
      X.train <- as.matrix(X[Train.set])
      
    } else {
      
      C.test <- C[Test.set, , drop = FALSE]
      C.train <- C[Train.set, , drop = FALSE]
      
      X.test <- X[Test.set, ]
      X.train <- X[Train.set, ]
    }
    ###############################################
    B.test <- as.matrix(B[Test.set, ])
    B.train <- B[Train.set, ]
    # BQ2.test=BQ2[Test.set,]
    # BQ2.train=BQ2[Train.set,]
    ###############################################
    y.test <- y[Test.set]
    y.train <- y[Train.set]
    
    # model fit
    mfit.ii <- SQBiT(Y = y.train, C = C.train, X = X.train, B = B.train, Q2 = Q2, P = P, 
                           lambda = lambda, tau = tau, eta.j = eta.j, var.j = var.j,
                           h = h, gacv.compute = FALSE)
    
    W.test <- as.matrix(kr(X.test, eigenMapMatMult(B.test, Q2), byrow = TRUE)) #B.test
    ypred.ii <- eigenMapMatMult(C.test, as.vector(mfit.ii$eta)) + eigenMapMatMult(W.test, as.vector(mfit.ii$gamma))
    
    ## prediction check loss
    pred.error <- mean(rhotau(u = y.test - ypred.ii, tau = tau)) 
    cv.error <- c(cv.error, pred.error)
  }
  cv.error
 }
