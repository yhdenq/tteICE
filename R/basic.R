# library(survival)
## Hidden functions

.matchy <- function(yvec,xvec,newx,exact=FALSE){
  # options(warn = -1)
  if (exact) {
    ivec = sapply(newx, function(x) suppressWarnings(max(which(xvec==x))))
  } else {
    ivec = sapply(newx, function(x) suppressWarnings(max(which(xvec<=x))))
  }
  if (is.vector(yvec)) {
    newy = yvec[ivec]
  } else {
    newy = yvec[ivec,]
  }
  newy[is.infinite(ivec)] = 0
  return(newy)
}

.ipscore <- function(A, X=NULL, standardize=TRUE, weights=rep(1,length(A))){
  if (is.null(X)) {
    ps = mean(A)
  } else {
    fps = glm.fit(cbind(1,X), A, family='binomial'(link='logit'), weights=weights)
    ps = fitted(fps)
  }
  if (standardize){
    ips = A/ps*mean(A/ps) + (1-A)/(1-ps)*mean((1-A)/(1-ps))
  } else {
    ips = A/ps + (1-A)/(1-ps)
  }
  return(ips)
}

.generatedata <- function(N=500, a1=0.05, a0=0.03, c1=0.04, c0=0.05){
  X1 = rbinom(N,1,0.5)
  X2 = rbinom(N,1,0.5)
  X = cbind(X1,X2)
  A = rbinom(N,1,1/(1+exp(0.5-0.3*X1-0.6*X2)))
  a1x = a1 + X1*0.02 - X2*0.03
  a0x = a0 + X1*0.02 - X2*0.01
  c1x = c1 + X1*0.01 - X2*0.02
  c0x = c0 - X1*0.01 - X2*0.01
  T1 = rweibull(N, 2, sqrt(2/a1x))
  T0 = rweibull(N, 2, sqrt(2/a0x))
  R1 = rexp(N, c1x)
  R0 = rexp(N, c0x)
  C = runif(N,4,8)
  C[C>7] = 7
  T = T1*A + T0*(1-A)
  R = R1*A + R0*(1-A)
  R[R>=T] = 99
  dT = as.numeric(T <= C)
  dR = as.numeric(R <= C)
  T = T*dT + C*(1-dT)
  R = R*dR + C*(1-dR)
  Time = (T+R-abs(T-R))/2
  cstatus = dT + 2*dR
  cstatus[cstatus>2] = 2
  return(list(Z=A,T=T,R=R,dT=dT,dR=dR,Time=Time,cstatus=cstatus,X=X))
}

.plot_ate_validate <- function(fit, decrease, conf.int, xlab, xlim, ylim) {
  # ---- fit ----
  if (!inherits(fit, "tteICE"))
    stop("`fit` must be an object returned by `surv.tteICE` or `scr.tteICE`.", call. = FALSE)
  # needed <- c("time", "est_treat", "est_control")  # customize to your object
  # missing_fields <- setdiff(needed, names(fit))
  # if (length(missing_fields))
  #   stop("`fit` is missing required components: ", paste(missing_fields, collapse = ", "), call. = FALSE)
  # if (!is.numeric(fit$time) || !isTRUE(all(is.finite(fit$time))) || !is.unsorted(fit$time, strictly = TRUE))
  #   stop("`fit$time` must be a strictly increasing numeric vector without NA/Inf.", call. = FALSE)

  # ---- decrease ----
  if (!(is.logical(decrease) && length(decrease) == 1 && !is.na(decrease)))
    stop("`decrease` must be a single non-NA logical.", call. = FALSE)

  # ---- conf.int ----
  if (!is.null(conf.int)) {
    if (!is.numeric(conf.int) || length(conf.int) != 1 || !is.finite(conf.int))
      stop("`conf.int` must be NULL or a single finite numeric value.", call. = FALSE)

    # Accept users passing 95 instead of 0.95
    if (conf.int > 1 && conf.int <= 100) {
      warning("Interpreting `conf.int = ", conf.int, "` as ", conf.int/100, ".")
      conf.int <- conf.int / 100
    }
    if (conf.int <= 0 || conf.int >= 1)
      stop("`conf.int` must be in (0, 1) when not NULL.", call. = FALSE)
  }

  # ---- labels & limits ----
  if (!(is.character(xlab) && length(xlab) == 1 && !is.na(xlab)))
    stop("`xlab` must be a single non-NA character string.", call. = FALSE)

  if (!is.null(xlim)) {
    if (!is.numeric(xlim) || length(xlim) != 2 || any(!is.finite(xlim)))
      stop("`xlim` must be NULL or a numeric length-2 vector of finite values.", call. = FALSE)
    if (!(xlim[1] < xlim[2]))
      stop("`xlim[1]` must be strictly less than `xlim[2]`.", call. = FALSE)
  }

  if (!is.null(ylim)) {
    if (!is.numeric(ylim) || length(ylim) != 2 || any(!is.finite(ylim)))
      stop("`ylim` must be NULL or a numeric length-2 vector of finite values.", call. = FALSE)
    if (!(ylim[1] < ylim[2]))
      stop("`ylim[1]` must be strictly less than `ylim[2]`.", call. = FALSE)
  }

  # Optionally constrain CIF/S(t) difference range if that’s a design assumption:
  # if (!is.null(ylim) && (ylim[1] < -1 || ylim[2] > 1))
  #   warning("`ylim` extends beyond [-1, 1]; ensure this is intentional for difference curves.")

  # Return possibly adjusted conf.int (e.g., 95 -> 0.95)
  list(conf.int = conf.int)
}


.plot_inc_validate <- function(fit, decrease, conf.int, xlab, xlim, ylim) {
  # ---- fit ----
  if (!inherits(fit, "tteICE"))
    stop("`fit` must be an object returned by `surv.tteICE` or `scr.tteICE`.", call. = FALSE)
  # needed <- c("time", "est_treat", "est_control")  # customize to your object
  # missing_fields <- setdiff(needed, names(fit))
  # if (length(missing_fields))
  #   stop("`fit` is missing required components: ", paste(missing_fields, collapse = ", "), call. = FALSE)
  # if (!is.numeric(fit$time) || !isTRUE(all(is.finite(fit$time))) || !is.unsorted(fit$time, strictly = TRUE))
  #   stop("`fit$time` must be a strictly increasing numeric vector without NA/Inf.", call. = FALSE)

  # ---- decrease ----
  if (!(is.logical(decrease) && length(decrease) == 1 && !is.na(decrease)))
    stop("`decrease` must be a single non-NA logical.", call. = FALSE)

  # ---- conf.int ----
  if (!is.null(conf.int)) {
    if (!is.numeric(conf.int) || length(conf.int) != 1 || !is.finite(conf.int))
      stop("`conf.int` must be NULL or a single finite numeric value.", call. = FALSE)

    # Accept users passing 95 instead of 0.95
    if (conf.int > 1 && conf.int <= 100) {
      warning("Interpreting `conf.int = ", conf.int, "` as ", conf.int/100, ".")
      conf.int <- conf.int / 100
    }
    if (conf.int <= 0 || conf.int >= 1)
      stop("`conf.int` must be in (0, 1) when not NULL.", call. = FALSE)
  }

  # ---- labels & limits ----
  if (!(is.character(xlab) && length(xlab) == 1 && !is.na(xlab)))
    stop("`xlab` must be a single non-NA character string.", call. = FALSE)

  if (!is.null(xlim)) {
    if (!is.numeric(xlim) || length(xlim) != 2 || any(!is.finite(xlim)))
      stop("`xlim` must be NULL or a numeric length-2 vector of finite values.", call. = FALSE)
    if (!(xlim[1] < xlim[2]))
      stop("`xlim[1]` must be strictly less than `xlim[2]`.", call. = FALSE)
  }

  if (!is.null(ylim)) {
    if (!is.numeric(ylim) || length(ylim) != 2 || any(!is.finite(ylim)))
      stop("`ylim` must be NULL or a numeric length-2 vector of finite values.", call. = FALSE)
    if (!(ylim[1] < ylim[2]))
      stop("`ylim[1]` must be strictly less than `ylim[2]`.", call. = FALSE)
  }

    if (!is.null(ylim) && (ylim[1] < 0 || ylim[2] > 1))
    warning("`ylim` extends beyond [0, 1]; ensure this is intentional for the curves.")

  # if (!is.character(legend) || length(legend) != 2)
  #   stop("'legend' must be a character vector of length 2.")

  # # col: character length 2
  # if (!is.character(col) || length(col) != 2)
  #   stop("'col' must be a character vector of length 2 (colors).")

  list(conf.int = conf.int)
}


.riskpredict_validate <- function(fit, timeset) {
  # ---- fit ----
  if (!inherits(fit, "tteICE"))
    stop("`fit` must be an object returned by `surv.tteICE` or `scr.tteICE`.", call. = FALSE)

  # # ---- timeset ----
  # if (!is.null(timeset)) {
  #   if (!is.numeric(timeset) || length(timeset) != 1 || !is.finite(timeset))
  #     stop("`timeset` must be NULL or a single finite numeric value.", call. = FALSE)
  # }
  }


#' Internal wrapper of survival:::event
#' @noRd
.event <- function(time, value=NULL, censor=NULL) {
        x <- list(time=time, value=value, censor=censor);
        class(x) <-"event"; x}

#' Internal wrapper of survival:::tdc
#' @noRd
.tdc <- function(time, value=NULL, init=NULL) {
        x <- list(time=time, value=value, default= init);
        class(x) <- "tdc"; x}




