#' @title Print the summary of 'tteICE'
#' @description Print the summary of 'tteICE'
#' @param x
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#' @param digits
#' The digits of the results
#' @param ...
#' Other arguments in function \code{\link{print.default}}
#' @method print summary.tteICE
#' @return Print the summary of a tteICE object
#' @export
#'
print.summary.tteICE <- function(x, digits=3, ...) {
  if(!is.null(x$call)){
    cat("Input:\n")
    print(x$call)
  }
  if(is.null(x$p.val)) p=NA else p=x$p.val
  cat("-----------------------------------------------------------------------\n")
  cat("Data type:", x$dtype, "\n")
  cat("Strategy:", x$strategy, "\n")
  cat("Estimation method:", x$method, "\n")
  cat("Observations:", x$n, '(including', x$n1, 'treated and', x$n0, 'control)\n')
  cat("Maximum follow-up time:", round(x$maxt, digits), '\n')
  cat("P-value of the average treatment effect:", round(p, digits), "\n")
  if (!is.null(x$coef)) {
    cat("-----------------------------------------------------------------------\n")
    coef1 = rbind(x$coef$coef11, x$coef$coef10)
    if (!is.null(coef1)){
      rownames(coef1) = c("A=1", "A=0")
      cat("Coefficients of covariates in the Cox model for event 1\n")
      cat(coef1,"\n")
    }
    coef2 = rbind(x$coef$coef21, x$coef$coef20)
    if (!is.null(coef2)){
      rownames(coef2) = c("A=1", "A=0")
      cat("Coefficients of covariates in the Cox model for event 2\n")
      cat(coef2,"\n")
    }
    coef = rbind(x$coef$coef1, x$coef$coef0)
    if (!is.null(coef)){
      rownames(coef) = c("A=1", "A=0")
      cat("Coefficients of covariates in the Cox model\n")
      cat(coef,"\n")
    }
  }
  cat("-----------------------------------------------------------------------\n")
  cat("The estimated cumulative incidences and treatment effects at quartiles:\n")
  print(round(x$est, digits))
  cat("\n")
  invisible(x)
}
