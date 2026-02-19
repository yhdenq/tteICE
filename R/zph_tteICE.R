#' @title Checking proportional hazards of 'tteICE' objects
#'
#' @description This function checks the proportional hazards assumption in the Cox models using Schoenfeld residuals
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @return
#' P-values of testing the proportional hazards (PH) assumption in the working Cox models, stratified by treatment groups.
#'
#' @export

coef.tteICE <- function(object) {
  object$ph
}
