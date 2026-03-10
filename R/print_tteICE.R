#' @title
#' Print method for 'tteICE' objects
#'
#' @description This function summarizes the results
#'
#' @param x
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param digits
#' The digits of the results
#'
#' @param ... Other arguments in function \code{\link{print.default}}
#'
#' @importFrom stats quantile
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' bmt$A = A
#'
#' ## print the results
#' fit1 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
#' print(fit1)
#'
#' fit2 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#' print(fit2)
#'
#' fit3 = tteICE(Surv(t2, fator(d4))~A|z1+z3+z5,
#'               data=bmt, strategy="composite", method='eff')
#' print(fit3, digits=4)
#'
#' @seealso
#' \code{\link[tteICE]{surv.tteICE}},
#' \code{\link[tteICE]{scr.tteICE}},
#' \code{\link[tteICE]{tteICE}}
#'
#' @method print tteICE
#' @return Print the summary of a tteICE object
#' @export

print.tteICE <- function(x, digits=3, ...){
  # if (!inherits(x, "tteICE"))  stop("Must be an object returned by `surv.tteICE` or `scr.tteICE`.", call. = FALSE)
  if(is.null(x$p.val)) p=NA else p=x$p.val
  dtype = c(cmprsk="competing risks", smcmprsk="semicompeting risks")
  strat = c(treatment="treatment policy strategy", composite="composite variable strategy",
           natural="hypothetical strategy (controlling the hazard of ICEs)",
           removed="hypothetical strategy (removing ICEs)",
           whileon="while on treatment strategy", principal="principal stratum strategy")
  meth = c(np="nonparametric estimation", eff="semiparametrically efficient estimation",
           ipw="inverse probability weighting")
  if(!is.null(x$call)){
   cat("Input:\n")
  print(x$call)
  }
  cat("-----------------------------------------------------------------------\n")
  cat("Data type:", dtype[x$dtype], "\n")
  cat("Strategy:", strat[x$strategy], "\n")
  cat("Estimation method:", meth[x$method], "\n")
  cat("Observations:", x$n, '(including', x$n1, 'treated and', x$n0, 'control)\n')
  cat("Maximum follow-up time:", max(x$time), '\n')
  cat("P-value of the average treatment effect:", round(p, digits), "\n")
  invisible(x)
}
