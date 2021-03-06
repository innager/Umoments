% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expectStr.R
\name{one_combination}
\alias{one_combination}
\title{Generate symbolic expression for expectation}
\usage{
one_combination(powvect, smpsize = "n")
}
\arguments{
\item{powvect}{vector of non-negative integers representing exponents
\eqn{j_1, \dots, j_m}{j[1], ..., j[m]} of non-central moments in
expectation (see "Details"). The position (index) of an element of this
vector indicates a corresponding moment, e.g. for \eqn{E(\overline{X}^5
\overline{X^4})}, \code{powvect = c(5, 0, 0, 1)}. Thus the vector will have
\code{m} elements if \code{m}'th is the highest moment.}

\item{smpsize}{symbol to be used for sample size. Defaults to \code{"n"}.}
}
\value{
A string representing a symbolic expression for further processing
  using computer algebra (e.g. with \emph{Sage} or \emph{SymPy}), for
  calculating numeric values, or to be rendered with \emph{Latex}.
}
\description{
Generate a string with symbolic expression for expectation of powers and
products of non-central (raw) sample moments of an arbitrary order.
}
\details{
For a zero-mean random variable \code{X} and a sample \eqn{X_1, \dots,
X_n}{X[1], ..., X[n]}, find \eqn{E(\bar{X}^{j_1} \overline{X^2}^{j_2}
\overline{X^3}^{j_3} \cdots \overline{X^m}^{j_m})}{E(X-bar^{j_1}
{X^2}-bar^{j[2]} \cdots {X^m}-bar^{j[m]})}, where \eqn{overline{X^k} = 1/n
\sum_{i = 1}^n X_i^{k}}{{X^k}-bar = sum_{i = 1}^n X[i]^k} is a \eqn{k}'th
non-central sample moment. The expression is given in terms of sample size
and true moments \eqn{\mu_k}{\mu[k]} of \eqn{X}. These expectations can
subsequently be used for generating unbiased central moment estimators of an
arbitrary order, Edgeworth expansions, and possibly solving other
higher-order problems.
}
\examples{
one_combination(c(5, 0, 2, 1))

}
