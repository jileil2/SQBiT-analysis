rG <- function(n, tau) {
  Is <- rbinom(n, size = 1, prob = tau)
  Is * (- 2 * tau) + 2 * (1 - tau) * (1 - Is)
}
