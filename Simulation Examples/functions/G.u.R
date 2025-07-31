# smoothed quantile loss
G.u <- function(u) {
  sqrt(2 / pi) * exp(-u ^ 2 / 2) + u * (1 - 2 * pnorm(-u))
}