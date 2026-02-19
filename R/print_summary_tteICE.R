#' @title Print the summary of 'tteICE'
#' @description Print the summary of 'tteICE'
#' @param object A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#' @export

print.summary.tteICE <- function(x) {
  if(!is.null(x$call)){
    cat("Input:\n")
    print(x$call)
  }
  cat("-----------------------------------------------------------------------\n")
  cat("Data type:", x$dtype, "\n")
  cat("Strategy:", x$strategy, "\n")
  cat("Estimation method:", x$method, "\n")
  cat("Observations:", x$n, '(including', x$n1, 'treated and', x$n0, 'control)\n')
  cat("Maximum follow-up time:", x$maxt, '\n')
  cat("P-value of the average treatment effect:", round(x$p.val, 3), "\n")
  cat("-----------------------------------------------------------------------\n")
  cat("The estimated cumulative incidences and treatment effects at quartiles:\n")
  print(round(x$est, 3))
  cat("\n")
  invisible(x)
}
