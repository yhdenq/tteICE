#' @title
#' Predict method for 'tteICE' objects at specific time points
#'
#' @description
#' This function predicts the potential cumulative incidence function and treatment effect at
#' specific time points.
#'
#' @param object
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param timeset
#' Time at which to predict the risk.
#' If \code{timeset=NULL}, risks will be predict at the quartiles of the maximum follow-up time.
#'
#' @param ... Other arguments in function \code{\link[stats]{predict}}
#'
#' @return
#' A matrix with each row being time points, potential cumulative incidences (under
#' treated and under control), treatment effects, standard errors, and P-values.
#'
#'
#' @seealso 
#' \code{\link[tteICE]{scr.tteICE}}, \code{\link[tteICE]{surv.tteICE}}, \code{\link[tteICE]{tteICE}}
#' \code{\link[tteICE]{surv.boot}}
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' bmt$A = A
#' X = as.matrix(bmt[,c('z1','z3','z5')])
#'
#' ## predict results at specified time points
#' ## model fitting using semicompeting risk data
#' fit1 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#' predict(fit1, timeset=c(670,2000))
#' 
#' ## predict results without specifying any time points
#' ## model fitting using competing risk data
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
#' predict(fit2)
#'
#' ## a simpler way
#' library(survival)
#' fit3 = tteICE(Surv(t2, factor(d4))~A|z1+z3+z5,
#'               data=bmt, strategy="composite", method='eff')
#' predict(fit3, timeset=c(670,2000))
#' predict(fit3)
#'
#'
#' @seealso 
#' \code{\link[tteICE]{surv.tteICE}}, \code{\link[tteICE]{scr.tteICE}}, 
#' \code{\link[tteICE]{tteICE}}
#' 
#' @method predict tteICE
#' @return predict a tteICE object.
#' The meanings of each row are: time points, potential cumulative incidences (under
#' treated and under control), treatment effects, standard errors, and P-values.
#' @export

predict.tteICE <- function(object, timeset=NULL, ...){

  # .riskpredict_validate(fit, timeset)
  if (is.null(timeset)) {
    maxt = max(object$time)
    timeset = c(0.25,0.5,0.75,1)*maxt
  }
  cif1 = .matchy(object$cif1,object$time,timeset)
  cif0 = .matchy(object$cif0,object$time,timeset)
  ate = cif1 - cif0
  se1 = .matchy(object$se1,object$time,timeset)
  se0 = .matchy(object$se0,object$time,timeset)
  se = .matchy(object$se,object$time,timeset)
  p = 2*pnorm(-abs(ate/se))
  res = rbind(cif1,se1,cif0,se0,ate,se,p)
  colnames(res) = timeset
  rownames(res) = c('CIF1','se1','CIF0','se0','ATE','se','p.val')
  return(res)
}
