#' @title Checking proportional hazards of 'tteICE' objects
#'
#' @description This function checks the proportional hazards assumption in the Cox models using Schoenfeld residuals
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @keywords internal

coef.tteICE <- function(object) {
  object$ph
}
