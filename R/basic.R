# library(survival)
## Hidden functions

.matchy <- function(yvec,xvec,newx,exact=FALSE){
  options(warn = -1)
  if (exact) {
    ivec = sapply(newx, function(x) max(which(xvec==x)))
  } else {
    ivec = sapply(newx, function(x) max(which(xvec<=x)))
  }
  if (is.vector(yvec)) {
    newy = yvec[ivec]
  } else {
    newy = yvec[ivec,]
  }
  newy[is.infinite(ivec)] = 0
  return(newy)
}

.ipscore <- function(A, X, standardize=TRUE, weights=rep(1,length(A)), subset=rep(TRUE,length(A))){
  if (is.null(X)) {
    return(rep(1, length(A)))
  }
  fps = glm(A~X, family='binomial', weights=weights, subset=subset)
  ps = predict(fps, type='response')
  ips = rep(1, length(A))
  if (standardize){
    ips0 = A/ps*mean(A/ps) + (1-A)/(1-ps)*mean((1-A)/(1-ps))
  } else {
    ips0 = A/ps + (1-A)/(1-ps)
  }
  ips[subset] = ips0
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

.phfit <- function(Tr,Dr,Td,Dd,A,X,a){
  X0 = X
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  Dr[Tr==max(Tr)] = 0
  Dd[Td==max(Td)] = 0
  tt = sort(unique(c(Td,Tr[Dr==1])))
  maxt = max(tt)*0.9
  L = length(tt)
  N = length(Td)
  if (!is.null(X)) {
    X = as.matrix(as.matrix(X)[A==a,])
    p = ncol(X)
    betad = rep(0, p)
    betar = rep(0, p)
  } else {
    p = 0
    betad = betar = 0
  }
  lamd = rep(1/L, L)
  lamr = rep(1/L, L)
  # d
  eXb = rep(1, N)
  delta_r = 0
  iter = 0; tol = 1
  while(iter<100){
    iter = iter+1
    est0 = c(betad,delta_r,lamd[tt<maxt])
    S0 = sapply(Td, function(l) sum((Td>=l)*eXb*exp(delta_r*Dr*(Tr<l))))
    if (!is.null(X)){
      S1 = t(sapply(Td, function(l) colSums((Td>=l)*eXb*exp(delta_r*Dr*(Tr<l))*X)))
      S2 = t(sapply(Td, function(l) t(X)%*%diag((Td>=l)*eXb*exp(delta_r*Dr*(Tr<l)))%*%X))
      if (p==1) {S1=t(S1);S2=t(S2)}
      dbeta = as.numeric(t(X-S1/S0)%*%Dd)
      ddbeta = t(S1/S0)%*%diag(Dd)%*%(S1/S0) - matrix(colSums(Dd*S2/S0),p,p)
      betad = betad - ginv(ddbeta) %*% dbeta
      eXb = exp(as.numeric(X%*%betad))
    }
    if (sum(Dd*(Tr<Td))>0) {
      S1 = sapply(Td, function(l) sum((Td>=l)*eXb*exp(delta_r*Dr*(Tr<l))*Dr*(Tr<l)))
      ddelta_r = sum(Dd*((Td>Tr)*Dr-S1/S0))
      dddelta_r = sum(Dd*(S1/S0)^2) - sum(Dd*S1/S0)
      delta_r = delta_r - ddelta_r/dddelta_r
    }
    delta_r = sign(delta_r)*min(abs(delta_r),5)
    lamd = sapply(tt, function(t) sum(Dd*(Td==t))/sum((Td>=t)*eXb*exp((Tr<t)*delta_r)))
    lamd[is.nan(lamd)] = 0
    est = c(betad,delta_r,lamd[tt<maxt])
    tol = max(abs(est-est0))
    if (tol<0.00001) break
  }
  delta_rd = delta_r
  # r
  eXb = rep(1, N)
  iter = 0; tol = 1
  while(iter<100){
    iter = iter+1
    est0 = c(betar,lamr[tt<maxt])
    S0 = sapply(Tr, function(l) sum((Tr>=l)*eXb))
    if (!is.null(X)){
      S1 = t(sapply(Tr, function(l) colSums((Tr>=l)*eXb*X)))
      S2 = t(sapply(Tr, function(l) t(X)%*%diag((Tr>=l)*eXb)%*%X))
      if (p==1) {S1=t(S1);S2=t(S2)}
      dbeta = as.numeric(t(X-S1/S0)%*%Dr)
      ddbeta = t(S1/S0)%*%diag(Dr)%*%(S1/S0) - matrix(colSums(Dr*S2/S0),p,p)
      betar = betar - ginv(ddbeta) %*% dbeta
      eXb = exp(as.numeric(X%*%betar))
    }
    lamr = sapply(tt, function(t) sum(Dr*(Tr==t))/sum((Tr>=t)*eXb))
    lamr[is.nan(lamr)] = 0
    est = c(betar,lamr[tt<maxt])
    tol = max(abs(est-est0))
    if (tol<0.00001) break
  }
  if (!is.null(X)){
    betad = as.numeric(betad)
    betar = as.numeric(betar)
    Xbd = as.numeric(X0%*%betad)
    Xbr = as.numeric(X0%*%betar)
  } else {
    betad = betar = 0
    Xbd = Xbr = rep(0,N)
  }
  return(list(betad=betad, betar=betar, delta=delta_rd,
              Xbd=Xbd, Xbr=Xbr, tt=tt, lamd=lamd, lamr=lamr))
}

.phfit_c <- function(Tr,Dr,Td,Dd,A,X,a){
  X0 = X
  Td = Td[A==a]
  Dd = Dd[A==a]
  Dd[Td==max(Td)] = 1
  Dc = 1 - Dd
  tt = sort(unique(Td[Dc==1]))
  maxt = max(tt)*0.9
  L = length(tt)
  N = length(Td)
  if (!is.null(X)) {
    X = as.matrix(as.matrix(X)[A==a,])
    p = ncol(X)
    beta = rep(0, p)
  } else {
    p = 0
    beta = 0
  }
  lam = rep(1/L, L)
  eXb = rep(1, N)
  iter = 0; tol = 1
  while(iter<100){
    iter = iter+1
    est0 = c(beta,lam[tt<maxt])
    S0 = sapply(Td, function(l) sum((Td>=l)*eXb))
    if (!is.null(X)){
      S1 = t(sapply(Td, function(l) colSums((Td>=l)*eXb*X)))
      S2 = t(sapply(Td, function(l) t(X)%*%diag((Td>=l)*eXb)%*%X))
      if (p==1) {S1=t(S1);S2=t(S2)}
      dbeta = as.numeric(t(X-S1/S0)%*%Dc)
      ddbeta = t(S1/S0)%*%diag(Dc)%*%(S1/S0) - matrix(colSums(Dc*S2/S0),p,p)
      beta = beta - ginv(ddbeta) %*% dbeta
      eXb = exp(as.numeric(X%*%beta))
    }
    lam = sapply(tt, function(t) sum(Dc*(Td==t))/sum((Td>=t)*eXb))
    lam[is.nan(lam)] = 0
    est = c(beta,lam[tt<maxt])
    tol = max(abs(est-est0))
    if (tol<0.00001) break
  }
  if (!is.null(X)){
    beta = as.numeric(beta)
    Xb = as.numeric(X0%*%beta)
  } else {
    beta = 0
    Xb = rep(0,N)
  }
  return(list(beta=beta, Xb=Xb, tt=tt, lam=lam))
}


.plot_ate_validate <- function(fit, decrease, conf.int, nboot, seed, xlab, xlim, ylim) {
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

  # ---- nboot ----
  if (!is.numeric(nboot) || length(nboot) != 1 || nboot < 0 || !is.finite(nboot) || nboot != as.integer(nboot))
    stop("`nboot` must be a single nonnegative integer.", call. = FALSE)

  if (nboot == 0 && !is.null(conf.int)) {
    # Analytical SE intended; fine.
    # No action.
  }
  if (nboot > 0 && is.null(conf.int)) {
    warning("`nboot > 0` but `conf.int = NULL`. Bootstrap SEs will be computed but no CI will be drawn.")
  }

  # ---- seed ----
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed != as.integer(seed))
      stop("`seed` must be NULL or a single integer.", call. = FALSE)
    if (nboot == 0)
      warning("`seed` is supplied but `nboot = 0`; seed will be ignored.")
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


.plot_inc_validate <- function(fit, decrease, conf.int, nboot, seed, xlab, xlim, ylim) {
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

  # ---- nboot ----
  if (!is.numeric(nboot) || length(nboot) != 1 || nboot < 0 || !is.finite(nboot) || nboot != as.integer(nboot))
    stop("`nboot` must be a single nonnegative integer.", call. = FALSE)

  if (nboot == 0 && !is.null(conf.int)) {
    # Analytical SE intended; fine.
    # No action.
  }
  if (nboot > 0 && is.null(conf.int)) {
    warning("`nboot > 0` but `conf.int = NULL`. Bootstrap SEs will be computed but no CI will be drawn.")
  }

  # ---- seed ----
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed != as.integer(seed))
      stop("`seed` must be NULL or a single integer.", call. = FALSE)
    if (nboot == 0)
      warning("`seed` is supplied but `nboot = 0`; seed will be ignored.")
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


.riskpredict_validate <- function(fit, timeset, nboot, seed) {
  # ---- fit ----
  if (!inherits(fit, "tteICE"))
    stop("`fit` must be an object returned by `surv.tteICE` or `scr.tteICE`.", call. = FALSE)

  # ---- timeset ----
  if (!is.null(timeset)) {
    if (!is.numeric(timeset) || !is.finite(timeset))
      stop("`timeset` must be NULL or finite numeric values.", call. = FALSE)
  }

  # ---- nboot ----
  if (!is.numeric(nboot) || length(nboot) != 1 || nboot < 0 || !is.finite(nboot) || nboot != as.integer(nboot))
    stop("`nboot` must be a single nonnegative integer.", call. = FALSE)

  # ---- seed ----
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed != as.integer(seed))
      stop("`seed` must be NULL or a single integer.", call. = FALSE)
    if (nboot == 0)
      warning("`seed` is supplied but `nboot = 0`; seed will be ignored.")
  }
  }
