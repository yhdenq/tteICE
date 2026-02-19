#' @title Baseline hazards of 'tteICE' objects
#'
#' @description This function extracts the baseline cumulative hazards in the survival models.
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @return
#' A data frame of baseline cumulative hazards in the working Kaplan-Meier or Cox models, stratified by treatment groups. 
#' The first column is time, the following columns are baseline cumulative hazards.
#'
#' @export

basehaz.tteICE <- function(object) {
  object$cumhaz
}
