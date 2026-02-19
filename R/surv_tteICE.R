#' @title 
#' Fit CIFs for competing risks time-to-event data with intercurrent events.
#'
#' @description This function estimates the potential cumulative incidence function
#' for time-to event data under ICH E9 (R1) to address intercurrent events. The input data 
#' should be of a competing risks structure.
#'
#' @param A Treatment indicator, 1 for treatment and 0 for control.
#'
#' @param Time Time to event.
#'
#' @param cstatus Indicator of event, 1 for the primary event, 2 for the intercurrent event, 0 for censoring.
#'
#' @param strategy Strategy to address intercurrent events, \code{"treatment"} indicating treatment policy strategy,
#' \code{"composite"} indicating composite variable strategy, \code{"natural"} indicating hypothetical strategy
#' (Scenario I, controlling the hazard of intercurrent events), \code{"removed"} indicating hypothetical strategy
#' (Scenario II, removing intercurrent events), \code{"whileon"} indicating while on treatment strategy, and
#' \code{"principal"} indicating principal stratum strategy.
#'
#' @param cov1 Baseline covariates.
#'
#' @param method Estimation method, \code{"np"} indicating nonparametric estimation, \code{"np"} indicating inverse 
#' treatment probability weighting, \code{"eff"} indicating semiparametrically efficient estimation based on efficient 
#' influence functions.
#'
#' @param weights Weight for each subject.
#'
#' @param subset Subset, either numerical or logical.
#'
#' @param na.rm Whether to remove missing values.
#'
#' @param nboot Number of resamplings in the boostrapping method. If \code{nboot} is 0 or 1, then
#' asymptotic standard error based on the explicit form is calculated instead of bootstrapping.
#'
#' @param seed Seed for bootstrapping.
#'
#'
#' @return A list including the fitted object and input variables.
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' X = as.matrix(bmt[,c('z1','z3','z5')])
#' 
#' ## Composite variable strategy, 
#' ## nonparametric estimation without covariates
#' fit1 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
#' 
#' ## Hypothetical strategy (natural effects), 
#' ## nonparametric estimation with inverse probability weighting
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "natural", X, method='ipw')
#' 
#' ## nonparametric estimation with weights as inverse propensity score
#' ps = predict(glm(A ~ X, family='binomial'), type='response')
#' w = A/ps + (1-A)/(1-ps)
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "natural", weights=w)
#' 
#' ## Hypothetical strategy (removing intercurrent events),
#' ## semiparametrically efficient estimation with covariates
#' fit3 = surv.tteICE(A, bmt$t2, bmt$d4, "removed", X, method='eff')
#'
#' @details
#' \describe{
#' \item{Background}{Intercurrent events refer to the events occurring after treatment initiation of clinical trials that
#' affect either the interpretation of or the existence of the measurements associated with the clinical
#' question of interest. The International Conference on Harmonization (ICH) E9 (R1) addendum proposed
#' five strategies to address intercurrent events, namely, treatment policy strategy,
#' composite variable strategy, while on treatment strategy, hypothetical strategy, and principal stratum
#' strategy. To answer a specific scientific question, a strategy with a particular estimand is chosen
#' before the study design.}
#' \item{Model}{We adopt the potential outcomes framework that defines a causal estimand as the contrast between
#' functionals of potential outcomes. Consider a randomized controlled trial with \eqn{n} individuals
#' randomly assigned to one of two treatment conditions, denoted by \eqn{w}, where \eqn{w = 1} represents
#' the active treatment (a test drug) and \eqn{w = 0} represents the control (placebo). Assume that all
#' patients adhere to their treatment assignments and do not discontinue treatment. Associated with individual
#' \eqn{i = 1, ..., n} are two potential time-to-event primary outcomes \eqn{T_i(1)} and \eqn{T_i(0)},
#' if any, which represent the time durations from treatment initiation to the primary outcome event under
#' two treatment assignments respectively. Let \eqn{R_i(1)} and \eqn{R_i(0)} denote the occurrence time of
#' potential intercurrent events, if any, under the two treatment assignments, respectively. Intercurrent
#' events are considered as absent if no post-treatment intercurrent events occur until the end of study.}
#' \item{Estimand}{We adopt the potential cumulative incidences under both treatment assignments as the target estimands.
#' Potential cumulative incidences describe the probability of time-to-event outcomes occurring at each
#' time point. We define the treatment effect as the contrast of two potential cumulative incidences.
#' Cumulative incidences are model-free and collapsible, enjoying causal interpretations.}
#' }
#'
#' @references
#' Deng, Y., Han, S., & Zhou, X. H. (2025).
#' Inference for Cumulative Incidences and Treatment Effects in Randomized Controlled Trials With Time-to-Event Outcomes Under ICH E9 (R1).
#' \emph{Statistics in Medicine}. \doi{10.1002/sim.70091}
#' 
#' @seealso \code{\link[tteICE]{surv.boot}}, \code{\link[tteICE]{scr.tteICE}}
#'
#' @keywords internal

surv.tteICE <- function(A,Time,cstatus,strategy='composite',cov1=NULL,method='np',
                        weights=NULL,subset=NULL,na.rm=FALSE,nboot=0,seed=0){
  
  # strategy <- match.arg(strategy, c('treatment','composite','natural','removed','whileon','principal'))
  # method <- match.arg(method, c('np','ipw','eff'))
  
  if (!strategy %in% c('treatment','composite','natural','removed','whileon','principal')){
    warning("Please choose a strategy from the following:\n treatment, composite, natural, removed, whileon, principal\n
            composite variable strategy is used by default", call. = FALSE)
    strategy = 'composite'
  }
  if (!method %in% c('np','ipw','eff')){
    warning("Please choose a method from the following:\n np, ipw, eff\n
            nonparametric estimation is used by default", call. = FALSE)
    method = 'np'
  }
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
  if (is.null(subset)) subset = rep(TRUE,N)
  if (inherits(subset,"logical")) subset = (1:N)[subset]
  if (na.rm){
    cc = complete.cases(data.frame(A,Time,cstatus,weights,subset,cov1))
    cc = (1:N)[cc]
    subset = subset[subset%in%cc]
  }
  if (length(unique(A))!=2) {
    stop('Treatment should be binary!', call. = FALSE)
  } else {
    A = as.numeric(A)
    if (min(A)!=0 | max(A)!=1) {
      A = as.numeric(A==max(A))
      warning(paste0('Treatment should be either 0 or 1! A=1 if A=',max(A)), call. = FALSE)
    }
  }
  
  A = A[subset]
  Time = Time[subset]
  cstatus = cstatus[subset]
  weights = weights[subset]
  if (!is.null(cov1)) cov1 = as.matrix(cov1)[subset,]
  
  n = length(A); n1 = sum(A==1); n0 = sum(A==0)
  if (method=='ipw') {
    ips = .ipscore(A,cov1,TRUE,weights)
  } else {
    ips = rep(1, n)
  }
  if (method=='np' | method=='ipw') {
    if (strategy=='treatment') fit = surv.treatment(A,Time,cstatus,weights*ips)
    if (strategy=='composite') fit = surv.composite(A,Time,cstatus,weights*ips)
    if (strategy=='natural') fit = surv.natural(A,Time,cstatus,weights*ips)
    if (strategy=='removed') fit = surv.removed(A,Time,cstatus,weights*ips)
    if (strategy=='whileon') fit = surv.whileon(A,Time,cstatus,weights*ips)
    if (strategy=='principal') fit = surv.principal(A,Time,cstatus,weights*ips)
  } else if (method=='eff') {
    if (strategy=='treatment') fit = surv.treatment.eff(A,Time,cstatus,cov1)
    if (strategy=='composite') fit = surv.composite.eff(A,Time,cstatus,cov1)
    if (strategy=='natural') fit = surv.natural.eff(A,Time,cstatus,cov1)
    if (strategy=='removed') fit = surv.removed.eff(A,Time,cstatus,cov1)
    if (strategy=='whileon') fit = surv.whileon.eff(A,Time,cstatus,cov1)
    if (strategy=='principal') fit = surv.principal.eff(A,Time,cstatus,cov1)
  }
  fit = c(fit,list(A=A,Time=Time,cstatus=cstatus,strategy=strategy,cov1=cov1,
                   method=method,weights=weights,na.rm=FALSE,dtype='cmprsk'))
  if (nboot>-1) fit = surv.boot(fit,nboot,seed)
  fit = c(fit, list(n=n, n1=n1, n0=n0, call= match.call()))
  
  class(fit) = "tteICE"
  return(fit)
}
