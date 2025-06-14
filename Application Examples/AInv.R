AInv <- function(A) {
  svd.A <- svd(A)
  d <- svd.A$d
  u <- svd.A$u
  v <- svd.A$v
  d[which(d == 0)] <- 1e-20
  D.m1 <- diag(1 / d)
  eigenMapMatMult(eigenMapMatMult(v, D.m1), t(u))
}
