#' Runs Bayesian instrumental variables regression using the same syntax as AER::ivreg
#'
#' @inheritParams AER::ivreg
#' @export


stan_ivreg <- function(formula, instruments, data, subset, na.action, weights, offset, ...) {

  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights",
               "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (!missing(instruments)) {
    formula <- Formula::as.Formula(formula, instruments)
    cl$instruments <- NULL
    cl$formula <- formula(formula)
  }
  else {
    formula <- Formula::as.Formula(formula)
  }
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in%
              1:2)
  has_dot <- function(formula) inherits(try(terms(formula),
                                            silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2))
      formula <- Formula::as.Formula(f1, update(formula(formula,
                                                        lhs = 0, rhs = 1), f2))
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf, contrasts)
  if (length(formula)[2] < 2L) {
    mtZ <- NULL
    Z <- NULL
  }
  else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf, contrasts)
  }
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if (is.null(offset))
    offset <- 0
  if (length(offset) == 1)
    offset <- rep(offset, NROW(Y))
  offset <- as.vector(offset)

  rval <- stan_iv(Y, X, Z, ...)

  class(rval) <- "stan_ivreg"
  return(rval)
}

#' Runs Bayesian IV
#'
#' @param Y a vector of outcomes
#' @param X_endog a matrix of endogenous regressors
#' @param Z a matrix of instruments
#' @param offset a vector offset
#' @param ... additional arguments to be passed to rstan::sampling
#' @export

stan_iv <- function(Y, X, Z, ...) {

  Y <- as.data.frame(Y)
  X <- as.data.frame(X)
  Z <- as.data.frame(Z)

  X_exog <- Z[,names(Z) %in% names(X)]
  X_exog[,"(Intercept)"] <- NULL
  Z <- Z[,!names(Z) %in% names(X_exog)]
  X_endog <- X[,!names(X) %in% names(X_exog)]
  starts_at <- "(Intercept)" %in% names(Z) + 1


  data_list <- list(N  = length(Y[,1]),
                    PN = ncol(as.matrix(X_endog)),
                    PZ = ncol(as.matrix(Z)),
                    PX = ncol(as.matrix(X_exog)),
                    X_endog = as.matrix(X_endog),
                    X_exog = as.matrix(X_exog),
                    Z = as.matrix(Z),
                    Y_outcome = Y[,1],
                    starts_at = starts_at)

  model_fit <- rstan::sampling(stanmodels$classic_iv, data = data_list, ...)
  out <- list(model_fit = model_fit,
              X_exog = X_exog,
              X_endog = X_endog,
              Z = Z,
              Y = Y)

  return(out)
}

#' print.stan_ivreg a print method for stan_ivreg
#' @param model_fit a model fit by
#' @import broom
#' @export

print.stan_ivreg <- function(model_fit) {

  par_table <- tidy(model_fit$model_fit, conf.int = T, pars = c("gamma1", "gamma2"))
  par_table[,1] <- c(names(cbind(model_fit$X_endog, model_fit$X_exog)), names(cbind(model_fit$Z, model_fit$X_exog)))
  par_table$stage <- c(rep("Second Stage", ncol(cbind(model_fit$X_endog, model_fit$X_exog))), rep("First stage",ncol(cbind(model_fit$Z, model_fit$X_exog) )))

  par_table
}
