#' @family pooled estimates (two-sample)
#' @inherit uM6pool title description params
#' @return Pooled variance estimate.
#' @examples
#' nx <- 10
#' ny <- 8
#' shp <- 3
#' smpx <- rgamma(nx, shape = shp) - shp
#' smpy <- rgamma(ny, shape = shp)
#' m2 <- mean(c((smpx - mean(smpx))^2, (smpy - mean(smpy))^2))
#' uM2pool(m2, nx, ny)
#' @export
uM2pool <- function(m2, n_x, n_y) {
  m2*(n_x + n_y)/(n_x + n_y - 2)
}

#' @family pooled estimates (two-sample)
#' @inherit uM6pool title description params
#' @return Pooled estimate of a third central moment.
#' @examples
#' nx <- 10
#' ny <- 8
#' shp <- 3
#' smpx <- rgamma(nx, shape = shp) - shp
#' smpy <- rgamma(ny, shape = shp)
#' mx <- mean(smpx)
#' my <- mean(smpy)
#' m  <- numeric(3)
#' for (j in 2:3) {
#'   m[j] <- mean(c((smpx - mx)^j, (smpy - my)^j))
#' }
#' uM3pool(m[3], nx, ny)
#' @export
uM3pool <- function(m3, n_x, n_y) {
  m3*n_x*n_y*(n_x + n_y)/(n_x^2*n_y + n_x*n_y^2 - 6*n_x*n_y + 2*n_x + 2*n_y)
}

#' @family pooled estimates (two-sample)
#' @inherit uM6pool title description params
#' @return Pooled estimate of squared variance \eqn{\mu_2^2}{\mu[2]^2}, where
#'   \eqn{\mu_2}{\mu[2]} is a variance.
#' @examples
#' nx <- 10
#' ny <- 8
#' shp <- 3
#' smpx <- rgamma(nx, shape = shp) - shp
#' smpy <- rgamma(ny, shape = shp)
#' mx <- mean(smpx)
#' my <- mean(smpy)
#' m  <- numeric(4)
#' for (j in 2:4) {
#'   m[j] <- mean(c((smpx - mx)^j, (smpy - my)^j))
#' }
#' uM2pow2pool(m[2], m[4], nx, ny)
#' @export
uM2pow2pool <- function(m2, m4, n_x, n_y) {
  n_x*n_y*(-m2^2*(n_x^2 + 2*n_x*n_y + n_y^2)*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 6*n_x^2*n_y - 3*n_x^2 + 6*n_x*n_y^2 - 3*n_y^2) + m4*n_x*n_y*(n_x + n_y)*(n_x^2*n_y + n_x*n_y^2 - 4*n_x*n_y + n_x + n_y))/(3*(n_x^2*n_y + n_x*n_y^2 - 4*n_x*n_y + n_x + n_y)*(4*n_x^2*n_y^2 - 5*n_x^2*n_y + 3*n_x^2 - 5*n_x*n_y^2 + 3*n_y^2) - (n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 6*n_x^2*n_y - 3*n_x^2 + 6*n_x*n_y^2 - 3*n_y^2)*(n_x^3*n_y + 2*n_x^2*n_y^2 - 5*n_x^2*n_y + n_x*n_y^3 - 5*n_x*n_y^2 + 12*n_x*n_y - 3*n_x - 3*n_y))
}

#' @family pooled estimates (two-sample)
#' @inherit uM6pool title description params
#' @return Pooled estimate of a fourth central moment.
#' @examples
#' nx <- 10
#' ny <- 8
#' shp <- 3
#' smpx <- rgamma(nx, shape = shp) - shp
#' smpy <- rgamma(ny, shape = shp)
#' mx <- mean(smpx)
#' my <- mean(smpy)
#' m  <- numeric(4)
#' for (j in 2:4) {
#'   m[j] <- mean(c((smpx - mx)^j, (smpy - my)^j))
#' }
#' uM4pool(m[2], m[4], nx, ny)
#' @export
uM4pool <- function(m2, m4, n_x, n_y) {
  n_x*n_y*(3*m2^2*(n_x^2 + 2*n_x*n_y + n_y^2)*(4*n_x^2*n_y^2 - 5*n_x^2*n_y + 3*n_x^2 - 5*n_x*n_y^2 + 3*n_y^2) - m4*n_x*n_y*(n_x + n_y)*(n_x^3*n_y + 2*n_x^2*n_y^2 - 5*n_x^2*n_y + n_x*n_y^3 - 5*n_x*n_y^2 + 12*n_x*n_y - 3*n_x - 3*n_y))/(3*(n_x^2*n_y + n_x*n_y^2 - 4*n_x*n_y + n_x + n_y)*(4*n_x^2*n_y^2 - 5*n_x^2*n_y + 3*n_x^2 - 5*n_x*n_y^2 + 3*n_y^2) - (n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 6*n_x^2*n_y - 3*n_x^2 + 6*n_x*n_y^2 - 3*n_y^2)*(n_x^3*n_y + 2*n_x^2*n_y^2 - 5*n_x^2*n_y + n_x*n_y^3 - 5*n_x*n_y^2 + 12*n_x*n_y - 3*n_x - 3*n_y))
}

#' @family pooled estimates (two-sample)
#' @inherit uM6pool title description params
#' @param m5 naive biased fifth central moment estimate \eqn{m_5 = 1/(n_x + n_y)
#'   \sum_{i = 1}^{n_x} ((X_i - \bar{X})^5 + \sum_{i = 1}^{n_y} ((Y_i -
#'   \bar{Y})^5}{m[5] = mean(c((X - X-bar)^5, (Y - Y-bar)^5))} for vectors
#'   \code{X} and \code{Y}.
#' @return Pooled estimate of a product of second and third central moments
#'   \eqn{\mu_2 \mu_3}{\mu[2] \mu[3]}, where \eqn{\mu_2}{\mu[2]} and
#'   \eqn{\mu_3}{\mu[3]} are second and third central moments respectively.
#' @examples
#' nx <- 10
#' ny <- 8
#' shp <- 3
#' smpx <- rgamma(nx, shape = shp) - shp
#' smpy <- rgamma(ny, shape = shp)
#' mx <- mean(smpx)
#' my <- mean(smpy)
#' m  <- numeric(5)
#' for (j in 2:5) {
#'   m[j] <- mean(c((smpx - mx)^j, (smpy - my)^j))
#' }
#' uM2M3pool(m[2], m[3], m[5], nx, ny)
#' @export
uM2M3pool <- function(m2, m3, m5, n_x, n_y) {
  n_x^2*n_y^2*(m2*m3*(n_x^2 + 2*n_x*n_y + n_y^2)*(n_x^4*n_y^3 + n_x^3*n_y^4 - 10*n_x^3*n_y^3 + 10*n_x^3*n_y^2 - 10*n_x^3*n_y + 4*n_x^3 + 10*n_x^2*n_y^3 - 10*n_x*n_y^3 + 4*n_y^3) - m5*n_x*n_y*(n_x + n_y)*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 5*n_x^2*n_y - 2*n_x^2 + 5*n_x*n_y^2 - 2*n_y^2))/(10*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 5*n_x^2*n_y - 2*n_x^2 + 5*n_x*n_y^2 - 2*n_y^2)*(-2*n_x^3*n_y^3 + 5*n_x^3*n_y^2 - 8*n_x^3*n_y + 4*n_x^3 + 5*n_x^2*n_y^3 - 8*n_x*n_y^3 + 4*n_y^3) + (n_x^4*n_y^3 + n_x^3*n_y^4 - 10*n_x^3*n_y^3 + 10*n_x^3*n_y^2 - 10*n_x^3*n_y + 4*n_x^3 + 10*n_x^2*n_y^3 - 10*n_x*n_y^3 + 4*n_y^3)*(n_x^4*n_y^2 + 2*n_x^3*n_y^3 - 12*n_x^3*n_y^2 + 2*n_x^3*n_y + n_x^2*n_y^4 - 12*n_x^2*n_y^3 + 60*n_x^2*n_y^2 - 42*n_x^2*n_y + 20*n_x^2 + 2*n_x*n_y^3 - 42*n_x*n_y^2 + 20*n_y^2))
}

#' @family pooled estimates (two-sample)
#' @inherit uM6pool title description params
#' @inheritParams uM2M3pool
#' @return Pooled estimate of a fifth central moment.
#' @examples
#' nx <- 10
#' ny <- 8
#' shp <- 3
#' smpx <- rgamma(nx, shape = shp) - shp
#' smpy <- rgamma(ny, shape = shp)
#' mx <- mean(smpx)
#' my <- mean(smpy)
#' m  <- numeric(5)
#' for (j in 2:5) {
#'   m[j] <- mean(c((smpx - mx)^j, (smpy - my)^j))
#' }
#' uM5pool(m[2], m[3], m[5], nx, ny)
#' @export
uM5pool <- function(m2, m3, m5, n_x, n_y) {
  n_x^2*n_y^2*(10*m2*m3*(n_x^2 + 2*n_x*n_y + n_y^2)*(-2*n_x^3*n_y^3 + 5*n_x^3*n_y^2 - 8*n_x^3*n_y + 4*n_x^3 + 5*n_x^2*n_y^3 - 8*n_x*n_y^3 + 4*n_y^3) + m5*n_x*n_y*(n_x + n_y)*(n_x^4*n_y^2 + 2*n_x^3*n_y^3 - 12*n_x^3*n_y^2 + 2*n_x^3*n_y + n_x^2*n_y^4 - 12*n_x^2*n_y^3 + 60*n_x^2*n_y^2 - 42*n_x^2*n_y + 20*n_x^2 + 2*n_x*n_y^3 - 42*n_x*n_y^2 + 20*n_y^2))/(10*(n_x^3*n_y^2 + n_x^2*n_y^3 - 8*n_x^2*n_y^2 + 5*n_x^2*n_y - 2*n_x^2 + 5*n_x*n_y^2 - 2*n_y^2)*(-2*n_x^3*n_y^3 + 5*n_x^3*n_y^2 - 8*n_x^3*n_y + 4*n_x^3 + 5*n_x^2*n_y^3 - 8*n_x*n_y^3 + 4*n_y^3) + (n_x^4*n_y^3 + n_x^3*n_y^4 - 10*n_x^3*n_y^3 + 10*n_x^3*n_y^2 - 10*n_x^3*n_y + 4*n_x^3 + 10*n_x^2*n_y^3 - 10*n_x*n_y^3 + 4*n_y^3)*(n_x^4*n_y^2 + 2*n_x^3*n_y^3 - 12*n_x^3*n_y^2 + 2*n_x^3*n_y + n_x^2*n_y^4 - 12*n_x^2*n_y^3 + 60*n_x^2*n_y^2 - 42*n_x^2*n_y + 20*n_x^2 + 2*n_x*n_y^3 - 42*n_x*n_y^2 + 20*n_y^2))
}


