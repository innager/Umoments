#=========================================================================#
#                     E(Xbar^k1*X2bar^k2*X3bar^k3*...)                    #
#                 powvect: c(k1, k2, k3, ...), k's can be 0,              #
#                  e.g. E(Xbar^k1*X3bar^k3) -> c(k1, 0, k3)               #
#=========================================================================#

# sorting the groups in a specific way;
# this function creates a next grouping from a given one
# from 4 3 3 2 2 1  get 4 3 3 3 1 1
one_forward <- function(vect) {
  m <- sum(vect)
  if (length(vect) == 1) return(vect)
  for (i in (length(vect) - 1):1) {
    if (i == 1) {
      vect[i] <- vect[i] + 1
      return(c(vect[1:i], rep(1, m - sum(vect[1:i]))))
    }
    else if (vect[i - 1] > vect[i]) {
      vect[i] <- vect[i] + 1
      return(c(vect[1:i], rep(1, m - sum(vect[1:i]))))
    }
  }
}

# creates a list of all possible groupings using our sorting rule
groups <- function(nelem) { # number of elements, sum(powvect)
  l <- list(vect <- rep(1, nelem))
  while (length(vect) > 1) {
    vect <- one_forward(vect)
    l <- c(l, list(vect))
  }
  return(l)
}

# all the combinations for letters where "a" is Xbar, "b" is X2bar and so on
# input is a vector of powXbar (# of a's), powX2bar (# of b's), powX3bar, ...
perm_categories <- function(powvect) {
  K <- length(powvect)
  nall <- sum(powvect)
  prev_vect <- paste(rep(0, nall), collapse = "")
  spaces_left <- c(nall - c(0, cumsum(powvect)[-K]))
  for (k in 1:K) {
    ind_next <- combn(1:spaces_left[k], powvect[k])
    prev_vect <- add_letter(ind_next, prev_vect, letters[k])
  }
  return(prev_vect)
}

# helper for perm_categories()
add_letter <- function(ind_next, prev_vect, letter) {
  next_vect <- matrix(nrow = length(prev_vect), ncol = ncol(ind_next))
  for (i in 1:length(prev_vect)) {
    for (j in 1:ncol(ind_next)) {
      old <- strsplit(prev_vect[i], split = "")[[1]]
      not_filled <- which(old == "0")  # 0's originally
      old[not_filled[ind_next[, j]]] <- letter
      next_vect[i, j] <- paste(old, collapse = "")
    }
  }
  return(as.character(next_vect))
}

# sort letters in a string - needed to compare strings
sort_str <- function(str) {
  return(paste(sort((strsplit(str, split = "")[[1]])), collapse = ""))
}

# calculate coefficient for one grouping
# elements sorted, K is length(powvect); the biggest letter is letters[K]
# excl is a vector (grep class) to remove letters and leave only one
calculate_coef <- function(str_vect, K, excl) {
  coef_numer(str_vect, K, excl)/coef_denom(str_vect)
}

# numerator for calculate_coef(), use choose() to get integer coefs
coef_numer <- function(str_vect, K, excl) {
  coefs <- matrix(nrow = K, ncol = length(str_vect))
  for (k in 1:K) {
    nn <- nchar(gsub(excl[k], "", str_vect))
    cn <- c(0, cumsum(nn))[- (length(nn) + 1)]  # 0 in front, remove last
    coefs[k, ] <- choose(sum(nn) - cn, nn)
  }
  return(prod(coefs))
}

# denominator for calculate_coef()
# permutations of repeated groups
# str_vect has to be sorted - strings and groups
coef_denom <- function(str_vect) {
  if (length(str_vect) == 1) return(1)
  repeats <- duplicated(str_vect)
  k <- 1
  multiples <- NULL
  for (i in 2:(length(str_vect))) {
    if (repeats[i]) k <- k + 1
    else {
      multiples <- c(multiples, k)
      k <- 1
    }
  }
  multiples <- c(multiples, k)
  return(prod(factorial(multiples)))
}

# table of all the groups and their coefficients (number of occurences)
all_groups <- function(powvect) {
  K <- length(powvect)
  m <- sum(powvect)
  groupings <- groups(m)                            # list
  perms     <- perm_categories(powvect)             # vector of strings
  excl <- paste("[^", letters[1:K], "]", sep = "")  # letters to exclude

  mat <- matrix("", ncol = m + 2)
  for (g in 1:length(groupings)) {
    grp <- groupings[[g]]
    for (p in 1:length(perms)) {
      permgroup <- substring(perms[p],          # break the string
                             c(1, cumsum(grp[-length(grp)]) + 1),
                             cumsum(grp))
      if (any(permgroup == "a")) next           # check order
      newrow <- sort(sapply(permgroup, sort_str))
      # sorted strings and sorted vector
      mat <- rbind(mat, c(newrow, rep("", m - length(newrow)),
                          length(newrow),
                          calculate_coef(newrow, K, excl)))
    }
  }
  all_empty <- apply(mat, 2, function(s) all(s == ""))
  mat <- unique(mat[, !all_empty])[-1, ]
  if (!inherits(mat, "matrix"))
    mat <- matrix(mat, nrow = 1)                  # if vector
  colnames(mat) <- c(paste("group", 1:(ncol(mat) - 2)), "d", "coef")
  return(mat)
}

# takes output of all_groups()
# subscript of mu is the order (moment)
# combines same order combinations and adds coefficients
combine_coef <- function(char_mat) {
  if (all(char_mat == "")) return(char_mat)
  J <- ncol(char_mat)
  if (all(dim(char_mat) == c(1, 3)))
    return(data.frame(mu = str_to_moment(char_mat[1, 1]), k = 1, coef = 1))
  num_mat <- t(apply(char_mat[, 1:(J-2)], 1, str_to_moment))
  mu_vect <- apply(num_mat, 1, vect_to_one)
  coef_vect <- tapply(char_mat[, "coef"], mu_vect,     # add coefficients
                      function(x) sum(as.numeric(x)))
  d_vect    <- tapply(char_mat[, "d"], mu_vect, unique)
  df_d <- data.frame(mu = names(d_vect), d = d_vect)
  df_coef <- data.frame(mu = names(coef_vect), coef = coef_vect)
  return(merge(df_d, df_coef, by = "mu"))
}

# helper for combine_coef()
# order of each group, eg "aacdd" has order 13 (mu13)
str_to_moment <- function(str_vect) {
  sapply(str_vect, function(str) {
    sum(match(strsplit(str, split = "")[[1]], letters))
  })
}
# helper for combine_coef()
# representation of moments (orders) in a grouping in single string
vect_to_one <- function(vect) {
  vect <- sort(vect[vect != 0])
  return(paste(vect, collapse = ":"))
}

#-------------------------------------------------------------------------#
#              Generate an expression string for Sage/SymPy:              #
#                    SR(expr_str) or sympify(expr_str)                    #
#-------------------------------------------------------------------------#

# inputs: vector of three elements (mu, d, and coef),
#         sample size notation (choose_from)
one_grouping <- function(vect3elem, smpsize) {
  if (all(vect3elem == "")) return("")
  mu <- strsplit(vect3elem[1], split = ":")
  mutab <- rev(table(mu))
  muexpr <- paste("mu", names(mutab), "^", mutab, "*",
                  sep = "", collapse = "")
  muexpr <- substr(muexpr, start = 1, stop = nchar(muexpr) - 1)

  k <- as.numeric(vect3elem[2])
  nexpr <- smpsize
  if (k > 1) {
    for (i in 1:(k-1)) {
      nexpr <- paste(nexpr, "*(", smpsize, "-", i, ")",
                     sep = "", collapse = "")
    }
  }
  coef_n_mu <- paste(vect3elem[3], "*", nexpr, "*", muexpr, " + ",
                     sep = "", collapse = "")
  return(coef_n_mu)
}

#' Generate symbolic expression for expectation
#'
#' Generate a string with symbolic expression for expectation of powers and
#' products of non-central (raw) sample moments of an arbitrary order.
#'
#' For a zero-mean random variable \code{X} and a sample \eqn{X_1, \dots,
#' X_n}{X[1], ..., X[n]}, find \eqn{E(\bar{X}^{j_1} \overline{X^2}^{j_2}
#' \overline{X^3}^{j_3} \cdots \overline{X^m}^{j_m})}{E(X-bar^{j_1}
#' {X^2}-bar^{j[2]} \cdots {X^m}-bar^{j[m]})}, where \eqn{overline{X^k} = 1/n
#' \sum_{i = 1}^n X_i^{k}}{{X^k}-bar = sum_{i = 1}^n X[i]^k} is a \eqn{k}'th
#' non-central sample moment. The expression is given in terms of sample size
#' and true moments \eqn{\mu_k}{\mu[k]} of \eqn{X}. These expectations can
#' subsequently be used for generating unbiased central moment estimators of an
#' arbitrary order, Edgeworth expansions, and possibly solving other
#' higher-order problems.
#'
#' @param powvect vector of non-negative integers representing exponents
#'   \eqn{j_1, \dots, j_m}{j[1], ..., j[m]} of non-central moments in
#'   expectation (see "Details"). The position (index) of an element of this
#'   vector indicates a corresponding moment, e.g. for \eqn{E(\overline{X}^5
#'   \overline{X^4})}, \code{powvect = c(5, 0, 0, 1)}. Thus the vector will have
#'   \code{m} elements if \code{m}'th is the highest moment.
#' @param smpsize symbol to be used for sample size. Defaults to \code{"n"}.
#'
#' @return A string representing a symbolic expression for further processing
#'   using computer algebra (e.g. with \emph{Sage} or \emph{SymPy}), for
#'   calculating numeric values, or to be rendered with \emph{Latex}.
#'
#' @examples
#' one_combination(c(5, 0, 2, 1))
#'
#' @export
#'
one_combination <- function(powvect, smpsize = "n") {
  if (!length(powvect))                        return("1")
  if (powvect[1] == 1 & sum(powvect[-1]) == 0) return("0")
  if (!sum(powvect))                           return("1")  # all 0's

  res <- all_groups(powvect)
  combined_res <- combine_coef(res)
  combined_res[, 1] <- as.character(combined_res[, 1])
  vect <- apply(combined_res, 1, one_grouping, smpsize = smpsize)
  str <- paste(vect, collapse = "")
  return(paste(" (", substr(str, start = 1, stop = nchar(str) - 2),
               ") / ", smpsize, "^", sum(powvect), sep = ""))
}











