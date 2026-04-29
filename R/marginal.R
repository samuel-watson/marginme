# =============================================================================
# Marginal effects for GLMMs (binomial-logit and Poisson-log only).
# R port of glmmr::Model::marginal.
#
# Computes one of:
#   - DyDx : dE[Y]/dx evaluated at xvals[1]
#   - Diff : E[Y | x = xvals[1]] - E[Y | x = xvals[2]]
#   - Ratio: log(E[Y | x = xvals[1]]) - log(E[Y | x = xvals[2]])
#
# Fixed effects can be conditioned on (`at`, `atmeans`) or averaged
# (`average`).  Random effects can be set to a value (`At`), set to their
# posterior-mean Zu (`AtEstimated`), set to zero (`AtZero`), or averaged over
# MCMC draws (`Average`).
# =============================================================================


# ---- inverse-link helpers ---------------------------------------------------
# For each family we need the value mu = g^{-1}(eta), the first derivative
# dmu/deta, and the second derivative d2mu/deta2.

.link_derivs <- function(eta, family) {
  if (family == "binomial") {
    mu <- stats::plogis(eta)
    dmu  <- mu * (1 - mu)
    d2mu <- dmu * (1 - 2 * mu)
  } else if (family == "poisson") {
    mu   <- exp(eta)
    dmu  <- mu
    d2mu <- mu
  } else {
    stop("Only 'binomial' (logit) and 'poisson' (log) supported.")
  }
  list(mu = mu, dmu = dmu, d2mu = d2mu)
}


# ---- core: one batch of (X-row, Zu) evaluations -----------------------------
# `Xrow` is either a 1xP matrix (single_row) or an nxP matrix.
# `zu` is either a vector (length matching nrow(Xrow)) or an n-by-iter matrix
# (only used in re_type = "Average").  When `Xrow` has 1 row but `zu` is a
# matrix, the single X row is reused for every Zu entry.
#
# Returns sums (not means) of the value and the P-vector gradient w.r.t. beta.
# This lets the caller divide by n*iter (or whatever denominator applies).

.sum_eval <- function(Xrow, beta, zu, family, mode, xcol = NULL) {

  if (is.matrix(zu)) {
    # Average branch: zu is n x iter.
    if (nrow(Xrow) == 1L) {
      xb  <- sum(Xrow[1, ] * beta)             # scalar
      eta <- xb + zu                           # n x iter
    } else {
      stopifnot(nrow(Xrow) == nrow(zu))
      xb  <- as.numeric(Xrow %*% beta)         # length n
      eta <- sweep(zu, 1, xb, "+")             # n x iter
    }
  } else {
    # zu is a length-N vector matching rows of Xrow.
    eta <- as.numeric(Xrow %*% beta) + zu      # length N
  }

  d <- .link_derivs(eta, family)

  if (mode == "BetaFirst") {
    # value = mu, grad_p = X[,p] * dmu
    value_sum <- sum(d$mu)
    if (is.matrix(zu)) {
      if (nrow(Xrow) == 1L) {
        grad_sum <- as.numeric(Xrow[1, ]) * sum(d$dmu)
      } else {
        grad_sum <- as.numeric(t(Xrow) %*% rowSums(d$dmu))
      }
    } else {
      grad_sum <- as.numeric(t(Xrow) %*% d$dmu)
    }

  } else if (mode == "XBeta") {
    # value = dmu * beta[xcol]
    # grad_p = d2mu * X[,p] * beta[xcol] + dmu * 1{p == xcol}
    bx <- beta[xcol]
    value_sum <- bx * sum(d$dmu)
    if (is.matrix(zu)) {
      if (nrow(Xrow) == 1L) {
        grad_sum <- bx * as.numeric(Xrow[1, ]) * sum(d$d2mu)
      } else {
        grad_sum <- bx * as.numeric(t(Xrow) %*% rowSums(d$d2mu))
      }
    } else {
      grad_sum <- bx * as.numeric(t(Xrow) %*% d$d2mu)
    }
    grad_sum[xcol] <- grad_sum[xcol] + sum(d$dmu)

  } else {
    stop("Unknown mode")
  }

  list(value_sum = value_sum, grad_sum = grad_sum)
}


#' Marginal effect from a binomial-logit or Poisson-log GLMM
#'
#' @param X            n x P design matrix.  Must have column names.
#' @param beta         length-P fitted coefficient vector.
#' @param M            P x P variance-covariance matrix of beta (caller chooses
#'                     model-based / sandwich / KR etc.).
#' @param family       "binomial" (logit) or "poisson" (log).
#' @param x            Name (character) or index (integer) of the column whose
#'                     marginal effect we want.
#' @param at           Character vector of column names to fix at user values.
#' @param atvals       Numeric vector matching `at`.
#' @param atmeans      Character vector of column names to fix at column mean.
#' @param average      Character vector of column names to average over (their
#'                     empirical distribution in X is preserved).
#' @param re_type      Random-effect treatment:
#'                       "At"          - fix Zu at user value(s) via `Zu`
#'                       "AtEstimated" - posterior-mean Zu = rowMeans(zu_samples)
#'                       "AtZero"      - Zu = 0
#'                       "Average"     - integrate over MCMC draws in zu_samples
#' @param type         "DyDx", "Diff", or "Ratio".
#' @param xvals        Length-2 numeric.  DyDx is evaluated at xvals[1];
#'                     Diff and Ratio compare xvals[1] vs xvals[2].
#' @param Zu           Used only when re_type = "At".  Either a single scalar
#'                     (when no covariates are averaged) or a length-n vector
#'                     of pre-computed Z u values (when `average` is non-empty
#'                     or only x is named with re_type Average/AtEstimated).
#' @param zu_samples   n x iter matrix of MCMC draws of Zu.  Required for
#'                     re_type = "Average" or "AtEstimated".
#' @param has_intercept TRUE if the first / an intercept column is in X and
#'                     should not be counted among "named" variables.  Used
#'                     only for the `at` + `atmeans` + `average` + 1 == P - int
#'                     consistency check.
#'
#' @return list with elements `estimate`, `se`, plus the input choices.
#' @export
glmm_marginal <- function(X, beta, M, family,
                          x,
                          at = character(0), atvals = numeric(0),
                          atmeans = character(0),
                          average = character(0),
                          re_type = c("AtZero", "AtEstimated", "At", "Average"),
                          type    = c("DyDx", "Diff", "Ratio"),
                          xvals   = c(0, 1),
                          Zu = NULL, zu_samples = NULL,
                          has_intercept = TRUE) {

  re_type <- match.arg(re_type)
  type    <- match.arg(type)
  family  <- match.arg(family, c("binomial", "poisson"))

  # --- validation --------------------------------------------------------
  if (!is.matrix(X)) stop("X must be a matrix.")
  if (is.null(colnames(X))) stop("X must have column names.")
  n <- nrow(X); P <- ncol(X); cn <- colnames(X)

  if (length(beta) != P)  stop("length(beta) != ncol(X)")
  if (any(dim(M) != c(P, P))) stop("M must be P x P")
  if (length(xvals) < 2L) xvals <- c(xvals, xvals)  # tolerate scalar for DyDx

  if (is.character(x)) {
    xcol <- match(x, cn)
    if (is.na(xcol)) stop(sprintf("Column '%s' not in X.", x))
  } else {
    xcol <- as.integer(x)
  }

  total_p    <- length(at) + length(atmeans) + length(average) + 1L
  expected_p <- P - as.integer(has_intercept)
  if (total_p != expected_p) {
    stop(sprintf("All non-intercept variables must be named: %d named, expected %d.",
                 total_p, expected_p))
  }
  if (length(at) != length(atvals)) stop("length(at) != length(atvals)")

  # --- decide single_row vs full n ---------------------------------------
  single_row <- !(length(average) > 0 ||
                    (total_p == 1L && re_type %in% c("Average", "AtEstimated")))
  N <- if (single_row) 1L else n

  # --- build newX --------------------------------------------------------
  if (single_row) {
    newX <- matrix(0, 1L, P, dimnames = list(NULL, cn))
    if (has_intercept) {
      icol <- match("(Intercept)", cn)
      if (is.na(icol)) icol <- 1L          # fall back to first column
      newX[1, icol] <- 1
    }
  } else {
    newX <- matrix(0, n, P, dimnames = list(NULL, cn))
    if (has_intercept) {
      icol <- match("(Intercept)", cn)
      if (is.na(icol)) icol <- 1L
      newX[, icol] <- 1
    }
    # preserve x column and any "average" columns from the original X
    newX[, xcol] <- X[, xcol]
    for (p in average) {
      idx <- match(p, cn)
      if (is.na(idx)) stop(sprintf("Variable '%s' not in column names.", p))
      newX[, idx] <- X[, idx]
    }
  }
  if (length(at) > 0) {
    for (k in seq_along(at)) {
      idx <- match(at[k], cn)
      if (is.na(idx)) stop(sprintf("Variable '%s' not in column names.", at[k]))
      newX[, idx] <- atvals[k]
    }
  }
  if (length(atmeans) > 0) {
    for (p in atmeans) {
      idx <- match(p, cn)
      if (is.na(idx)) stop(sprintf("Variable '%s' not in column names.", p))
      newX[, idx] <- mean(X[, idx])
    }
  }

  # --- assemble Zu (for At / AtEstimated / AtZero branch) ----------------
  if (re_type != "Average") {
    zu_vec <- numeric(N)
    if (re_type == "At") {
      if (is.null(Zu)) stop("`Zu` required when re_type = 'At'.")
      if (single_row) {
        if (length(Zu) != 1L) stop("Need a single Zu value when single_row.")
        zu_vec[1] <- Zu
      } else {
        if (length(Zu) != N) stop(sprintf("Need length-%d Zu vector.", N))
        zu_vec <- as.numeric(Zu)
      }
    } else if (re_type == "AtEstimated") {
      if (single_row) stop("AtEstimated requires non-single-row evaluation.")
      if (is.null(zu_samples)) stop("`zu_samples` required for AtEstimated.")
      if (nrow(zu_samples) != n) stop("zu_samples must have n rows.")
      zu_vec <- rowMeans(zu_samples)
    } # AtZero: leave zeros
  } else {
    if (is.null(zu_samples)) stop("`zu_samples` required for re_type = 'Average'.")
    if (nrow(zu_samples) != n) stop("zu_samples must have n rows.")
  }

  # --- core dispatch -----------------------------------------------------
  delta <- numeric(P)

  if (re_type != "Average") {
    # ---------- At / AtEstimated / AtZero ----------
    if (type == "DyDx") {
      newX[, xcol] <- xvals[1]
      s <- .sum_eval(newX, beta, zu_vec, family,
                     mode = "XBeta", xcol = xcol)
      est   <- s$value_sum / N
      delta <- s$grad_sum  / N

    } else if (type == "Diff") {
      newX[, xcol] <- xvals[1]
      s1 <- .sum_eval(newX, beta, zu_vec, family, mode = "BetaFirst")
      newX[, xcol] <- xvals[2]
      s2 <- .sum_eval(newX, beta, zu_vec, family, mode = "BetaFirst")
      est   <- (s1$value_sum - s2$value_sum) / N
      delta <- (s1$grad_sum  - s2$grad_sum)  / N

    } else { # Ratio
      newX[, xcol] <- xvals[1]
      s1 <- .sum_eval(newX, beta, zu_vec, family, mode = "BetaFirst")
      newX[, xcol] <- xvals[2]
      s2 <- .sum_eval(newX, beta, zu_vec, family, mode = "BetaFirst")
      # log E[mu1] - log E[mu2]; the 1/N factors cancel in each ratio.
      est   <- log(s1$value_sum / N) - log(s2$value_sum / N)
      delta <- s1$grad_sum / s1$value_sum - s2$grad_sum / s2$value_sum
    }

  } else {
    # ---------- Average ----------
    iter <- ncol(zu_samples)
    denom <- n * iter   # corrected normalisation (C++ used N*iter for delta)

    if (type == "DyDx") {
      newX[, xcol] <- xvals[1]
      s <- .sum_eval(newX, beta, zu_samples, family,
                     mode = "XBeta", xcol = xcol)
      est   <- s$value_sum / denom
      delta <- s$grad_sum  / denom

    } else if (type == "Diff") {
      newX[, xcol] <- xvals[1]
      s1 <- .sum_eval(newX, beta, zu_samples, family, mode = "BetaFirst")
      newX[, xcol] <- xvals[2]
      s2 <- .sum_eval(newX, beta, zu_samples, family, mode = "BetaFirst")
      est   <- (s1$value_sum - s2$value_sum) / denom
      delta <- (s1$grad_sum  - s2$grad_sum)  / denom

    } else { # Ratio
      newX[, xcol] <- xvals[1]
      s1 <- .sum_eval(newX, beta, zu_samples, family, mode = "BetaFirst")
      newX[, xcol] <- xvals[2]
      s2 <- .sum_eval(newX, beta, zu_samples, family, mode = "BetaFirst")
      est   <- log(s1$value_sum) - log(s2$value_sum)
      delta <- s1$grad_sum / s1$value_sum - s2$grad_sum / s2$value_sum
    }
  }

  se <- sqrt(as.numeric(crossprod(delta, M %*% delta)))

  list(estimate = est, se = se,
       type = type, re_type = re_type, family = family,
       xvals = xvals, x = if (is.character(x)) x else cn[xcol])
}
