ell.uh <- function(u, h, tau) {
  h / 2 * G.u(u / h) + (tau - 1 / 2) * u
}