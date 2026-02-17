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
#' fit3 = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5, 
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
             n=object$n, n1=object$n1, n0=object$n0, p.val=object$p.val, est=predict(object))
  class(res) <- "summary.tteICE"
  res
}
