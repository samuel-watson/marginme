#' Estimates Relative Risks, Risk Differences, and Marginal Effects from Mixed Models
#'
#' @description
#' Calculates the marginal effect of variable x. There are several options for
#' marginal effect and several types of conditioning or averaging. The type of marginal
#' effect can be the derivative of the mean with respect to x (`dydx`), the expected
#' difference E(y|x=a)-E(y|x=b) (`diff`), or the expected log ratio log(E(y|x=a)/E(y|x=b)) (`ratio`).
#' Other fixed effect variables can be set at specific values (`at`), set at their mean values
#' (`atmeans`), or averaged over (`average`). Averaging over a fixed effects variable here means
#' using all observed values of the variable in the relevant calculation.
#' The random effects can similarly be set at their
#' estimated value (`re="estimated"`), set to zero (`re="zero"`), set to a specific value
#' (`re="at"`), or averaged over (`re="average"`). The standard errors are calculated using the delta method with one
#' of several options for the variance matrix of the fixed effect parameters.
#' Several of the arguments require the names
#' of the variables as given to the model object. Most variables are as specified in the formula,
#' factor variables are specified as the name of the `variable_value`, e.g. `t_1`. To see the names
#' of the stored parameters and data variables see the member function `names()`.
#'
#' @param fit Either a lme4, glmmTMB, or glmmrBase model fit.
#' @param x String. Name of the variable to calculate the marginal effect for.
#' @param type String. Either `dydx` for derivative, `diff` for difference, or `ratio` for log ratio. See description.
#' @param re String. Either `estimated` to condition on estimated values, `zero` to set to zero, `at` to
#' provide specific values, or `average` to average over the random effects.
#' @param at Optional. A vector of strings naming the fixed effects for which a specified value is given.
#' @param atmeans Optional. A vector of strings naming the fixed effects that will be set at their mean value.
#' @param average Optional. A vector of strings naming the fixed effects which will be averaged over.
#' @param xvals Optional. A vector specifying the values of `a` and `b` for `diff` and `ratio`. The default is (1,0).
#' @param atvals Optional. A vector specifying the values of fixed effects specified in `at` (in the same order).
#' @param revals Optional. If `re="at"` then this argument provides a vector of values for the random effects.
#' @param oim Logical. If TRUE use the observed information matrix, otherwise use the expected information matrix for standard error and related calculations.
#' @param sampling Integer. Number of MCMC samples to use.
#' @return A named vector with elements `margin` specifying the point estimate and `se` giving the standard error.
#' @importFrom glmmrBase lme4_to_glmmr
#' @importFrom glmmrBase Model
#' @importFrom methods is
#' @importFrom stats binomial
#' @examples
#' ## fit a model using glmmTMB
#' fit <- glmmTMB::glmmTMB(y ~ Treatment + x1 + x2 + x3 + x4 + (1|Cluster),
#'   data = trial_data, family = binomial(link="logit"),REML = TRUE)
#' ## relative risk, average over random effects and fixed effects
#' m1 <- margin(fit,
#'        x = "Treatment",
#'        type = "ratio",
#'        average = c("x1","x2","x3","x4"),
#'        re = "average")
#' summary(m1)
#' ## stata default for margins command is to set random effects to zero
#' m2 <- margin(fit,
#'        x = "Treatment",
#'        type = "ratio",
#'        average = c("x1","x2","x3","x4"),
#'        re = "zero")
#' summary(m2)
#' ## finally estimate a risk difference, with random effects at zero and fixed effects
#' ## at mean values
#' m3 <- margin(fit,
#'        x = "Treatment",
#'        type = "diff",
#'        atmeans = c("x1","x2","x3","x4"),
#'        re = "zero")
#' summary(m3)
#' @export
margin <- function(fit, x,type,re,at = c(),atmeans = c(),average=c(),
                   xvals=c(1,0),atvals=c(),revals=c(),oim = FALSE, sampling = 250){

  if(is(fit,"glmmTMB")){
    df <- fit$frame
    f1 <- tryCatch(glmmrBase::lme4_to_glmmr(fit$call$formula,names(df)),error = function(e)return(list(NA)))
    if(!is(f1,"formula"))stop("Complex model formula in glmmTMB cannot be easily converted, please create a glmmrBase Model object")
    if(fit$modelInfo$REML){
      mean_pars <- fit$fit$parfull[grepl("beta",names(fit$fit$parfull))]
    } else {
      mean_pars <- fit$fit$par[which(names(fit$fit$par)=="beta")]
    }
    model <- glmmrBase::Model$new(
      f1,
      data = fit$frame,
      mean = mean_pars,
      covariance = exp(fit$fit$par[which(names(fit$fit$par)=="theta")])^2,
      family = binomial()
    )
  } else if(is(fit,"glmerMod")){
    df <- fit@frame
    f1 <- glmmrBase::lme4_to_glmmr(fit@call$formula,names(df))
    model <- glmmrBase::Model$new(
      f1,
      data = fit@frame,#df[!is.na(df[,as.character(f1[[2]])]),],
      mean = fit@beta,
      covariance = fit@theta^2,
      family = binomial()
    )
  } else if(is(fit,"Model")) {
    model <- Model$new(fit)
  } else {
    stop("Fit should be glmerMod, glmmTMB, or Model class")
  }

  suppressMessages(suppressWarnings(model$sample_u(sampling)))
  X <- model$mean$X
  colnames(X) <- gsub("b_","",names(model$mean$parameters))
  out <- glmm_marginal(X, model$mean$parameters,
                solve(model$information_matrix()), "binomial",
                "Treatment", atmeans = c("x1","x2","x3","x4"),
                type = "Diff", re_type = "Average",
                zu_samples = model$covariance$Z %*% model$u())
  out$family <- binomial()
  out$formula <- f1
  out$sampling <- sampling
  class(out) <- "margin"

  return(out)
}

#' Prints the marginal output
#'
#' Print method for class "`margin`"
#'
#' @param x An object of class "`margin`" resulting from a call to margin
#' @param ... Further arguments passed from other methods
#'
#' @return No return, called for effects
#' @examples
#' ## fit a model using glmmTMB
#' fit <- glmmTMB::glmmTMB(y ~ Treatment + x1 + x2 + x3 + x4 + (1|Cluster),
#'   data = trial_data, family = binomial(link="logit"),REML = TRUE)
#' ## relative risk, average over random effects and fixed effects
#' margin(fit,
#'        x = "Treatment",
#'        type = "ratio",
#'        average = c("x1","x2","x3","x4"),
#'        re = "average")
#' @export
print.margin <- function(x, ...){
  digits = 4
  cat("Marginal Effects from Mixed Model")
  f1 <- as.character(x$formula)
  cat("\nFormula: ",f1[[2]]," ~ ",f1[[3]])
  res <- round(x$estimate,digits)
  if(x$type == "Ratio"){
    cat("\nLog RR: ",res)
  } else if(x$type == "Diff"){
    cat("\nRisk diff.: ",res)
  } else {
    cat("\nMarginal effect (dydx): ",res)
  }
  cat(" SE: ",round(x$se,digits))
}


#' Summarises the marginal output
#'
#' Summary method for "`margin`" class
#'
#' @param object An object of class "`margin`" resulting from a call to margin
#' @param ... Further arguments passed from other methods
#'
#' @return No return, called for effects
#' @examples
#' ## fit a model using glmmTMB
#' fit <- glmmTMB::glmmTMB(y ~ Treatment + x1 + x2 + x3 + x4 + (1|Cluster),
#'   data = trial_data, family = binomial(link="logit"),REML = TRUE)
#' ## relative risk, average over random effects and fixed effects
#' m1 <- margin(fit,
#'        x = "Treatment",
#'        type = "ratio",
#'        average = c("x1","x2","x3","x4"),
#'        re = "average")
#' summary(m1)
#' @export
summary.margin <- function(object, ...){
  digits = 3
  cat("Marginal Effects from Mixed Model")
  f1 <- as.character(object$formula)
  cat("\nFormula: ",f1[[2]]," ~ ",f1[[3]])
  if(object$type == "Ratio"){
    name <- "Log RR "
  } else if(object$type == "Diff"){
    name <- "Risk diff. "
  } else {
    name <- "Marginal effect (dydx) "
  }
  name <- paste0(name, "(",object$x,")")
  cat("\n\n")
  out <- data.frame(est = round(object$estimate, digits),
                    se = round(object$se, digits),
                    z = round(object$estimate/object$se,digits))
  colnames(out) <- c("Estimate","Std. Err.","z value")
  rownames(out) <- name
  print(out)
  cat("\nRE type: ",object$re,", SE: ",object$se,", MCMC samples: ",object$sampling)
  return(invisible(out))
}

#' Confidence interval for marginal effect
#'
#' Confidence interval method for class "`margin`"
#'
#' @details
#' Computes confidence intervals using a standard Wald test for the marginal effect. If argument
#' `df` is used, then a t-statistic is use to construct the interval with `df` degrees of freedom,
#' otherwise a z-statistic is used.
#'
#'
#' @param object An object of class "`margin`" resulting from a call to margin
#' @param parm No effect.
#' @param level the confidence level required
#' @param ... additional argument(s) for methods
#'
#' @return A named vector giving lower and upper confidence limits for the marginal effect. They
#' will be labelled as (1-level)/2 and 1-(1-level).2 An argument df can be provided to use a t-test
#' to construct the confidence interval.
#' @examples
#' ## fit a model using glmmTMB
#' fit <- glmmTMB::glmmTMB(y ~ Treatment + x1 + x2 + x3 + x4 + (1|Cluster),
#'   data = trial_data, family = binomial(link="logit"),REML = TRUE)
#' ## relative risk, average over random effects and fixed effects
#' m1 <- margin(fit,
#'        x = "Treatment",
#'        type = "ratio",
#'        average = c("x1","x2","x3","x4"),
#'        re = "average")
#' confint(m1)
#' @importFrom stats qnorm qt
#' @export
 confint.margin <- function(object, parm, level = 0.95,...){
  x <- object
  args <- list(...)
  df <- if("df"%in%names(args)) NULL else args[["df"]]
  if(is.null(df)){
    cint <- c(x$estimate - qnorm(1 - (1-level)/2)*x$se, x$estimate + qnorm(1 - (1-level)/2)*x$se)
  } else {
    cint <- c(x$estimate - qt(1 - (1-level)/2, df = df)*x$se, x$estimate + qt(1 - (1-level)/2, df = df)*x$se)
  }
   names(cint) <- c(paste0((1-level)*100/2,"%"),paste0(100-(1-level)*100/2,"%"))
   return(cint)
}

#' Simulated trial data
#'
#' Simulated trial data used to demonstrate the estimation of relative risk from an adjusted mixed logistic regression model.
#' See \link[marginme]{margin}.
#' @name trial_data
#' @docType data
NULL
