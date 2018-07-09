#' Unbiased central moment estimates
#'
#' Calculate unbiased estimates of central moments and their powers and
#' products.
#'
#' @family unbiased estimates (one-sample)
#'
#' @param m2 naive biased variance estimate \eqn{m_2 = 1/n \sum_{i = 1}^n ((X_i
#'   - \bar{X})^2}{m[2] = mean((X - X-bar)^2)} for a vector \code{X}.
#' @param m3 naive biased third central moment estimate \eqn{m_3 = 1/n \sum_{i =
#'   1}^n ((X_i - \bar{X})^3}{m[3] = mean((X - X-bar)^3)} for a vector \code{X}.
#' @param m4 naive biased fourth central moment estimate \eqn{m_4 = 1/n \sum_{i
#'   = 1}^n ((X_i - \bar{X})^4}{m[4] = mean((X - X-bar)^4)} for a vector
#'   \code{X}.
#' @param m6 naive biased sixth central moment estimate \eqn{m_6 = 1/n \sum_{i =
#'   1}^n ((X_i - \bar{X})^6}{m[6] = mean((X - X-bar)^6)} for a vector \code{X}.
#' @param n sample size.
#' @return Unbiased estimate of a sixth central moment.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:6) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM6(m[2], m[3], m[4], m[6], n)
#' @export
uM6 <- function(m2, m3, m4, m6, n) {
  15*m2^3*(3*n - 10)*n^2/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) - 40*(n^2 - 6*n + 10)*m3^2*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) - 15*(n^3 - 8*n^2 + 29*n - 40)*m2*m4*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) + (n^4 - 9*n^3 + 31*n^2 - 39*n + 40)*m6*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased variance estimate.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' m <- c(m, mean((smp - m[1])^2))
#' uM2(m[2], n) - var(smp)
#' @export
uM2 <- function(m2, n) {
  m2*n/(n - 1)
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased estimate of a third central moment.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:3) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM3(m[3], n)
#' @export
uM3 <- function(m3, n) {
  m3*n^2/((n - 1)*(n - 2))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased estimate of squared variance \eqn{\mu_2^2}{\mu[2]^2}, where
#'   \eqn{\mu_2}{\mu[2]} is a variance.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:4) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM2pow2(m[2], m[4], n)
#' @export
uM2pow2 <- function(m2, m4, n) {
  (n^2 - 3*n + 3)*m2^2*n/((n - 1)*(n - 2)*(n - 3)) - m4*n/((n - 2)*(n - 3))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased estimate of a fourth central moment.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:4) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM4(m[2], m[4], n)
#' @export
uM4 <- function(m2, m4, n) {
  -3*m2^2*(2*n - 3)*n/((n - 1)*(n - 2)*(n - 3)) + (n^2 - 2*n + 3)*m4*n/((n - 1)*(n - 2)*(n - 3))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @param m5 naive biased fifth central moment estimate \eqn{m_5 = \sum_{i =
#'   1}^n ((X_i - \bar{X})^5}{m[5] = mean((X - X-bar)^5)} for a vector \code{X}.
#' @return Unbiased estimate of a product of second and third central moments
#'   \eqn{\mu_2 \mu_3}{\mu[2] \mu[3]}, where \eqn{\mu_2}{\mu[2]} and
#'   \eqn{\mu_3}{\mu[3]} are second and third central moments respectively.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:5) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM2M3(m[2], m[3], m[5], n)
#' @export
uM2M3 <- function(m2, m3, m5, n) {
  (n^2 - 2*n + 2)*m2*m3*n^2/((n - 1)*(n - 2)*(n - 3)*(n - 4)) - m5*n^2/((n - 2)*(n - 3)*(n - 4))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @inheritParams uM2M3
#' @return Unbiased estimate of a fifth central moment.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:5) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM5(m[2], m[3], m[5], n)
#' @export
uM5 <- function(m2, m3, m5, n) {
  -10*m2*m3*n^2/((n - 1)*(n - 3)*(n - 4)) + (n^2 - 5*n + 10)*m5*n^2/((n - 1)*(n - 2)*(n - 3)*(n - 4))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased estimate of cubed variance  central moment
#'   \eqn{\mu_2^3}{\mu[2]^3}, where \eqn{\mu_2}{\mu[2]} is a variance.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:6) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM2pow3(m[2], m[3], m[4], m[6], n)
#' @export
uM2pow3 <- function(m2, m3, m4, m6, n) {
  (n^2 - 7*n + 15)*m2^3*n^2/((n - 1)*(n - 3)*(n - 4)*(n - 5)) - 3*(n^2 - 5*n + 10)*m2*m4*n/((n - 1)*(n - 3)*(n - 4)*(n - 5)) + 2*m6*n/((n - 3)*(n - 4)*(n - 5)) - 2*(3*n^2 - 15*n + 20)*m3^2*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased estimate of squared third central moment
#'   \eqn{\mu_3^2}{\mu[3]^2}, where \eqn{\mu_3}{\mu[3]} is a third central
#'   moment.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:6) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM3pow2(m[2], m[3], m[4], m[6], n)
#' @export
uM3pow2 <- function(m2, m3, m4, m6, n) {
  -3*(3*n^2 - 15*n + 20)*m2^3*n^2/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) - (n^2 - n + 4)*m6*n/((n - 2)*(n - 3)*(n - 4)*(n - 5)) + (n^4 - 8*n^3 + 25*n^2 - 10*n - 40)*m3^2*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) + 3*(2*n^3 - 5*n^2 - 5*n + 20)*m2*m4*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5))
}

#' @family unbiased estimates (one-sample)
#' @inherit uM6 title description params
#' @return Unbiased estimate of a product of second and fourth central moments
#'   \eqn{\mu_2 \mu_4}{\mu[2] \mu[4]}, where \eqn{\mu_2}{\mu[2]} and
#'   \eqn{\mu_4}{\mu[4]} are second and fourth central moments respectively.
#' @examples
#' n <- 10
#' smp <- rgamma(n, shape = 3)
#' m <- mean(smp)
#' for (j in 2:6) {
#'   m <- c(m, mean((smp - m[1])^j))
#' }
#' uM2M4(m[2], m[3], m[4], m[6], n)
#' @export
uM2M4 <- function(m2, m3, m4, m6, n) {
  -3*m2^3*(2*n - 5)*n^2/((n - 1)*(n - 3)*(n - 4)*(n - 5)) + 4*(n^2 - 5*n + 10)*m3^2*n/((n - 1)*(n - 3)*(n - 4)*(n - 5)) - (n^2 - 3*n + 8)*m6*n/((n - 2)*(n - 3)*(n - 4)*(n - 5)) + (n^4 - 9*n^3 + 53*n^2 - 135*n + 120)*m2*m4*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5))
}
