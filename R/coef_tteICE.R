#' @title Coefficients of 'tteICE' objects
#'
#' @description This function extracts the coefficients in the Cox models
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @keywords internal

coef.tteICE <- function(object) {
  object$coef
}
