#' @title Using formula to fit CIFs for time-to-event data with intercurrent events
#'
#' @description
#' This function estimates the potential cumulative incidence function for time-to event data under ICH E9 (R1) to address intercurrent events.
#' The input data can be competing or semicompeting risks data structure.
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class).
#' A symbolic description of the model to be fitted.
#' For example, \code{formula=Surv(time, status)~treatment | baseline.covariate}.
#' The details of model specification are given under ‘Details’.
#'
#' @param add.scr Required for semicompeting data.
#' An object of class "Surv" (or one that can be coerced to that class).
#' For example, \code{add.scr=~Surv(time.intercurrent, status.intercurrent)}.
#' The details of model specification are given under ‘Details’.
#'
#' @param data Data or object coercible by as.data.frame to a data frame,
#' containing the variables in the model.
#'
#' @param strategy Strategy to address intercurrent events, \code{"treatment"} indicating treatment policy strategy,
#' \code{"composite"} indicating composite variable strategy, \code{"natural"} indicating hypothetical strategy
#' (Scenario I, controlling the hazard of intercurrent events), \code{"removed"} indicating hypothetical strategy
#' (Scenario II, removing intercurrent events), \code{"whileon"} indicating while on treatment strategy, and
#' \code{"principal"} indicating principal stratum strategy.
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
#' @return A list including the fitted object and input variables.
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' X = as.matrix(bmt[,c('z1','z3','z5')])
#' bmt$A = A
#'
#' library(survival)
#' ## Composite variable strategy,
#' ## nonparametric estimation without covariates
#' ## Composite variable strategy,
#' ## nonparametric estimation without covariates
#'
#' ## model fitting for competing risk data without covariates
#' fit1 = tteICE(Surv(t2, factor(d4)) ~ A,
#'  data=bmt, strategy="composite", method='np')
#' print(fit1)
#'
#' ## model fitting for competing risk data without covariates
#' ## with bootstrap confidence intervals
#' fit.bt1 = tteICE(Surv(t2, factor(d4)) ~ A,
#'  data=bmt, strategy="composite", method='eff', nboot=20, seed=2)
#' print(fit.bt1)
#'
#' ## model fitting for competing risk data with covariates
#' fit2 = tteICE(Surv(t2, factor(d4)) ~ A | z1 + z3 + z5,
#'  data=bmt, strategy="composite", method='eff')
#' print(fit2)
#'
#' ## model fitting for semicompeting risk data without covariates
#' fitscr1 = tteICE(Surv(t1, d1) ~ A, ~Surv(t2, d2),
#'  data=bmt, strategy="composite", method='np')
#' print(fitscr1)
#'
#' ## model fitting for semicompeting risk data without covariates
#' fitscr2 = tteICE(Surv(t1, d1) ~ A | z1 + z3 + z5, ~Surv(t2, d2),
#'  data=bmt, strategy="composite", method='eff')
#' print(fitscr2)
#'
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
#' \item{Formula specifications}{
#' The formula should be set in the following two ways.
#'
#' When data take format of competing risk data, set the first argument \code{formula = Surv(time, status) ~ treatment | covariate1+covariate2}
#' or \code{formula = Surv(time, status)~ A} without any baseline covariates, where \code{status} is a factor variable with levels 0,1,2 
#' (1 for the primary event, 2 for the intercurrent event, and 0 for censoring).
#'
#' When data take the format of semicompeting risk data, set the first argument \code{formula = Surv(time, status) ~ treatment | covariate1+covariate2}
#' or \code{formula = Surv(time, status) ~ A} without any baseline covariates, where \code{status}=0,1 (1 for the primary event and 0 for censoring).
#' In addition, the second argument \code{add.scr = ~ Surv(time.intercurrent, status.intercurrent)} is required.
#' }
#' }
#'
#' @references
#' Deng, Y., Han, S., & Zhou, X. H. (2025).
#' Inference for Cumulative Incidences and Treatment Effects in Randomized Controlled Trials With Time-to-Event Outcomes Under ICH E9 (R1).
#' \emph{Statistics in Medicine}. \doi{10.1002/sim.70091}
#'
#' @seealso \code{\link[tteICE]{surv.boot}}, \code{\link[tteICE]{scr.tteICE}}
#'
#' @importFrom stats model.frame model.response terms as.formula
# ' @importFrom survival Surv
#' @export

tteICE <- function(formula, add.scr=NULL, data, strategy='composite', method='np',
                     weights=NULL,subset=NULL,na.rm=FALSE,nboot=0,seed=0){

  # extract A, Time, cstatus
  if (missing(formula)) stop("`formula` is required,
    e.g. competing risk data type: `formula = Surv(time, status) ~ A | X`,
    or semicompeting risk data type: `formula = Surv(time.p, status.p) ~ A | X , add.scr = ~Surv(time.i, status.i)`")

  # if (missing(treatment)) stop("`treatment` is required.")
  if (missing(data)) stop("`data` is required.")

  fm <- model.frame(formula, data = data)

  ## extract outcomes
  Yp  <- model.response(fm)
  if (!inherits(Yp, "Surv")) stop("Use `Surv(time, status)` or `Surv(time, status)` as dependent variable")
  Time   <- Yp[, 1]
  cstatus <- Yp[, 2]

 if(!is.null(add.scr)) {
  Y2 = stats::model.frame(add.scr, data = data)[[1]]
  if (!inherits(Y2,"Surv")) stop("Use `~Surv(time, status)` ")
    Time_int   <- Y2[, 1]
    status_int <- Y2[, 2]
 }

# extract variables
rhs <- deparse(formula[[3]])
# treatment
parts <- strsplit(rhs, "\\|")[[1]]
A_name <- trimws(parts[1])
# A <- data$A_name
if ((length(A_name) < 1) || !A_name %in% names(data)) stop("Check the treatment variable.")
A <- unlist(subset(data, select=A_name),use.names = FALSE)

## covariates
if (length(parts) > 1) {
  rhs2 <- trimws(parts[2])
  cov <- all.vars(as.formula(paste("~", rhs2)))
  # cov1 <- data[,cov, drop = FALSE]
  cov1 <- subset(data, select=cov)
  } else cov1 <- NULL

  # tt <- terms(fm, data = data)
  # rhs_vars <- attr(tt, "term.labels")  # e.g. c("X1", "X2")
  # # if (length(rhs_vars) < 1) stop("Formula must include at least one treatment variable.")
  # X_name <- rhs_vars        # first term is A
  # X <- mf[,X_name]

  # if(missing(A)) stop("Must include at least one treatment variable.")
  # A <- treatment
  # cov1 <- cov1


  # strategy <- match.arg(strategy, c('treatment','composite','natural','removed','whileon','principal'))
  # method <- match.arg(method, c('np','ipw','eff'))
  #   if(!is.null(cov.formula)){
  #   if (!inherits(cov.formula, "formula") || length(cov.formula) != 2) {
  #   stop("`cov.formula` must be a one-sided formula like ~ X1 + X2 + X3")
  # }
  #   mf.cov <- model.frame(cov.formula, data = data)
  #   tt <- terms(mf.cov, data = data)
  #   rhs <- attr(tt, "term.labels")
  #   cov1 <-data[,rhs, drop = FALSE]
  #   cov1 <- as.matrix(cov1)[subset,]
  # } else cov1 <- NULL


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
  if (!is.null(add.scr)){
    Time_int = Time_int[subset]
    status_int = status_int[subset]
  }
  if (!is.null(cov1)) cov1 = as.matrix(cov1)[subset,]

  # extract covariates
  # if covariates exist
  # if(!if(0 %in% dim(X))){
  # #   if (!inherits(cov.formula, "formula") || length(cov.formula) != 2) {
  # #   stop("`cov.formula` must be a one-sided formula like ~ X1 + X2 + X3")
  # # }
  #   # mf.cov <- model.frame(cov.formula, data = data)
  #   # tt <- terms(mf.cov, data = data)
  #   # rhs <- attr(tt, "term.labels")
  #   # cov1 = data[,rhs, drop = FALSE]
  #   cov1 = as.matrix(cov1)[subset,]
  # }


  ## estimation
  if(is.null(add.scr)){
    if (method=='ipw') {
        weights = weights*.ipscore(A,cov1,TRUE,weights)
      }
      if (method=='np' | method=='ipw') {
        if (strategy=='treatment') fit = surv.treatment(A,Time,cstatus,weights)
        if (strategy=='composite') fit = surv.composite(A,Time,cstatus,weights)
        if (strategy=='natural') fit = surv.natural(A,Time,cstatus,weights)
        if (strategy=='removed') fit = surv.removed(A,Time,cstatus,weights)
        if (strategy=='whileon') fit = surv.whileon(A,Time,cstatus,weights)
        if (strategy=='principal') fit = surv.principal(A,Time,cstatus,weights)
      } else if (method=='eff') {
        if (strategy=='treatment') fit = surv.treatment.eff(A,Time,cstatus,cov1)
        if (strategy=='composite') fit = surv.composite.eff(A,Time,cstatus,cov1)
        if (strategy=='natural') fit = surv.natural.eff(A,Time,cstatus,cov1)
        if (strategy=='removed') fit = surv.removed.eff(A,Time,cstatus,cov1)
        if (strategy=='whileon') fit = surv.whileon.eff(A,Time,cstatus,cov1)
        if (strategy=='principal') fit = surv.principal.eff(A,Time,cstatus,cov1)
      }

      fit = c(fit,list(A=A,Time=Time,cstatus=cstatus,strategy=strategy,cov1=cov1,
                        method=method,weights=weights,na.rm=na.rm,dtype='cmprsk'))
  } else{
    if (method=='np' | method=='ipw') {
    if (strategy=='treatment') fit = scr.treatment(A,Time,cstatus,Time_int,status_int,weights)
    if (strategy=='composite') fit = scr.composite(A,Time,cstatus,Time_int,status_int,weights)
    if (strategy=='natural') fit = scr.natural(A,Time,cstatus,Time_int,status_int,weights)
    if (strategy=='removed') fit = scr.removed(A,Time,cstatus,Time_int,status_int,weights)
    if (strategy=='whileon') fit = scr.whileon(A,Time,cstatus,Time_int,status_int,weights)
    if (strategy=='principal') fit = scr.principal(A,Time,cstatus,Time_int,status_int,weights)
  } else if (method=='eff') {
    if (strategy=='treatment') fit = scr.treatment.eff(A,Time,cstatus,Time_int,status_int,cov1)
    if (strategy=='composite') fit = scr.composite.eff(A,Time,cstatus,Time_int,status_int,cov1)
    if (strategy=='natural') fit = scr.natural.eff(A,Time,cstatus,Time_int,status_int,cov1)
    if (strategy=='removed') fit = scr.removed.eff(A,Time,cstatus,Time_int,status_int,cov1)
    if (strategy=='whileon') fit = scr.whileon.eff(A,Time,cstatus,Time_int,status_int,cov1)
    if (strategy=='principal') fit = scr.principal.eff(A,Time,cstatus,Time_int,status_int,cov1)
  }
  fit = c(fit, list(A=A,Time=Time,status=cstatus,Time_int=Time_int,status_int=status_int,
                    strategy=strategy,cov1=cov1,method=method,
                    weights=weights,na.rm=FALSE,dtype='smcmprsk'))

  }

  if (nboot>-1) fit = surv.boot(fit,nboot,seed)
  n = length(A); n1 = sum(A==1); n0 = sum(A==0)
  fit = c(fit, list(n=n, n1=n1, n0=n0, call= match.call()))
  class(fit) = "tteICE"

  return(fit)
}

