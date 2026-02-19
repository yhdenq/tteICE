#' @title Coefficients of 'tteICE' objects
#'
#' @description This function extracts the coefficients in the Cox models
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @return
#' A list of coefficients of covariates in the working Cox models, stratified by treatment groups. 
#' For the treatment policy strategy and composite variable strategy, only one Cox model is fit (for the primary 
#' outcome event or the composite event). In these two strategies, \code{coef1} is the coefficients in the treated 
#' group, \code{coef0} is the coefficients in the control group. For other strategies, Cox models are fitted for 
#' each event (primary outcome event and intercurrent event). In these strategies, \code{coef11} is the coefficients 
#' for the primary outcome event in the treatment group, \code{coef10} is the coefficients for the primary outcome 
#' event in the control group, \code{coef21} is the coefficients for the inercurrent event in the treated group, 
#' \code{coef20} is the coefficients for the intercurrent in the control group.
#'
#' @export

coef.tteICE <- function(object) {
  object$coef
}
