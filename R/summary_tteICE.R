#' @title Summary method for 'tteICE' objects
#'
#' @description This function summarizes the results
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param ... Other arguments in function \code{\link{summary}}
#'
#' @importFrom stats quantile
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' bmt$A = A
#' X = as.matrix(bmt[,c('z1','z3','z5')])
#' 
#' ## Composite variable strategy,
#' ## nonparametric estimation without covariates
#' fit1 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#' summary(fit1)
#' 
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
#' predict(fit2)
#' 
#' library(survival)
#' fit3 = tteICE(Surv(t2, factor(d4))~A|z1+z3+z5, 
#'               data=bmt, strategy="composite", method='eff')
#' summary(fit3)
#'
#' @seealso
#' \code{\link[tteICE]{surv.tteICE}},
#' \code{\link[tteICE]{scr.tteICE}},
#' \code{\link[tteICE]{tteICE}},
#' \code{\link[tteICE]{print.tteICE}}
#'
#' @method summary tteICE
#' @return A list that consists of summaries of a tteICE object: data type, strategy, estimation method, maximum follow-up time,
#' sample size, treated sample size, controlled sample size, p-value, and predicted risks at quartiles
#' @export

summary.tteICE <- function(object, ...) {
  res = list(call=object$call,dtype=object$dtype, strategy=object$strategy, method=object$method, maxt=max(object$time),
             n=object$n, n1=object$n1, n0=object$n0, p.val=object$p.val, coef=coef(object), est=predict(object))
  class(res) <- "summary.tteICE"
  print(object)
  cat("-----------------------------------------------------------------------\n")
  if (!is.null(object$coef)) {
    coef1 = rbind(object$coef$coef11, object$coef$coef10)
    if (!is.null(coef1)){
      rownames(coef1) = c("A=1", "A=0")
      cat("Coefficients of covariates in the Cox model for event 1\n")
      cat(coef1,"\n")
    }
    coef2 = rbind(object$coef$coef21, object$coef$coef20)
    if (!is.null(coef2)){
      rownames(coef2) = c("A=1", "A=0")
      cat("Coefficients of covariates in the Cox model for event 2\n")
      cat(coef2,"\n")
    }
    coef = rbind(object$coef$coef1, object$coef$coef0)
    if (!is.null(coef)){
      rownames(coef) = c("A=1", "A=0")
      cat("Coefficients of covariates in the Cox model\n")
      cat(coef,"\n")
    }
  }
  cat("-----------------------------------------------------------------------\n")
  cat("The estimated cumulative incidences and treatment effects at quartiles:\n")
  print(round(object$est, digits))
  cat("\n")
  invisible(res)
}
