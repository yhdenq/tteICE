#' @title
#' Plot estimated cumulative incidence functions (CIFs)
#'
#' @description
#' This function plots the estimated potential cumulative incidence function, along with pointwise confidence intervals.
#'
#' @param fit
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param decrease
#' A logical variable indicating the type of curve to display.
#' If \code{decrease = FALSE} (default), cumulative incidence functions (CIFs) are plotted.
#' If \code{decrease = TRUE}, survival functions are plotted instead.
#'
#' @param conf.int
#' Confidence level for the pointwise confidence intervals
#' If \code{conf.int = NULL}, no confidence intervals are provided.
#'
#' @param xlab
#' Label for the x-axis.
#'
#' @param xlim
#' A numeric vector of length 2 specifying the limits of the x-axis.
#' If \code{xlim = NULL} (default), the limits are determined automatically from the data.
#'
#' @param ylim
#' A numeric vector of length 2 specifying the limits of the y-axis.
#' Defaults to \code{ylim = c(0, 1)}.
#'
#' @param plot.configs
#' A named \code{list} of additional plot configurations. Common entries include:
#'
#' \itemize{
#'     \item \code{ylab}: character, label for the y-axis (default: \code{ylab=NULL}, use the default label).
#'     \item \code{main}: character, title for the plot (default: \code{main=NULL}, use the default label).
#'     \item \code{lty}: line type for the curve (default: \code{lty=1}).
#'     \item \code{lwd}: line width for the curve (default: \code{lwd=2}).
#'     \item \code{ci.lty}: line type for confidence interval curves (default: \code{ci.lty=5, ci.lwd=1.5}).
#'     \item \code{ci.lwd}: line width for confidence interval curves (default: \code{ci.lwd=1.5}).
#'     \item \code{legend}: legend for the two group (default: \code{legend=c('Treated','Control')}).
#'     \item \code{col}: color of the curve for the two group (default: \code{col=c('brown','darkcyan')}).
#'     \item \code{legend.cex}: font size for the legend (default: \code{legend.cex=0.9}).
#'     \item \code{show.p.value}: whether to show the p-value between two groups (default: \code{show.p.value=TRUE}, show the p-value)
#'   }
#'
#' @param ... Additional graphical arguments passed to function \code{\link{plot.default}} or function \code{\link{curve}}
#'
#' @importFrom graphics plot points
#' @importFrom utils modifyList
#'
#' @seealso
#' \code{\link[graphics]{plot.default}},
#' \code{\link[graphics]{points}},
#' \code{\link[graphics]{curve}},
#' \code{\link[tteICE]{plot.tteICE}}
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' bmt$A = A
#'
#' ## simple model fitting and plotting
#' library(survival)
#' fit = tteICE(Surv(t2,d4,type = "mstate")~A, data=bmt)
#' plot_inc(fit)
#'
#' ## model fitting using competing risk data
#' fit1 = surv.tteICE(A, bmt$t2, bmt$d4, 'treatment')
#'
#' ## plot asymptotic confidence intervals based on explicit formulas
#' plot_inc(fit1, ylim=c(0,1),
#'          plot.configs=list(legend=c('AML','ALL'), show.p.value=FALSE) )
#'
#' ## plot bootstrap confidence intervals
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, 'treatment', nboot=50)
#' plot_inc(fit2, ylim=c(0,1),
#'          plot.configs=list(legend=c('AML','ALL')))
#'
#' ## model fitting using semicompeting risk data
#' fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#'
#' ## plot asymptotic confidence intervals based on explicit formulas
#' plot_inc(fit3, ylim=c(0,1), plot.configs=list(add.null.line=FALSE))
#'
#' ## plot bootstrap confidence intervals
#' fit4 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2,
#'                   "composite", nboot=50)
#' plot_inc(fit4, ylim=c(0,1),
#'          plot.configs=list(lty=2, lwd=3, main="My title"))
#'
#' @return Plot the cumulative incidence function results from a tteICE object
#' @keywords internal

plot_inc <- function(fit,decrease=FALSE,conf.int=.95,
  xlab='Time',xlim=NULL,ylim=c(0,1),
  plot.configs=list(ylab=NULL, main=NULL,
                    lty=1, lwd=2,
                    ci.lty=5, ci.lwd=1.5,
                    legend=c('Treated','Control'),
                    col=c('brown','darkcyan'),
                    legend.cex=0.9,
                    show.p.value=TRUE
                    ),...){

  #---- validate input ----#
  adj <- .plot_inc_validate(fit, decrease, conf.int, xlab, xlim, ylim)
  conf.int <- adj$conf.int

   #---- set configurations ----#
  default.configs <- list(ylab=NULL, main=NULL,
                    lty=1, lwd=2,
                    ci.lty=5, ci.lwd=1.5,
                    legend=c('Treated','Control'),
                    col=c('brown','darkcyan'),
                    legend.cex=0.9,
                    show.p.value=TRUE)
  plot.configs <- modifyList(default.configs, plot.configs)

  #---- plot ----#
  if(is.null(plot.configs[["main"]])){
  if (fit$strategy=='treatment') stname = 'Treatment policy strategy'
  if (fit$strategy=='composite') stname = 'Composite variable strategy'
  if (fit$strategy=='natural') stname = 'Hypothetical I (natural) strategy'
  if (fit$strategy=='removed') stname = 'Hypothetical II (removed) strategy'
  if (fit$strategy=='whileon') stname = 'While on treatment strategy'
  if (fit$strategy=='principal') stname = 'Principal stratum strategy'
} else stname=plot.configs[["main"]]

  if(is.null(fit$p.val)) p=NA else p = fit$p.val

  tt = fit$time
  if (decrease==TRUE){
    cif1 = 1-fit$cif1
    cif0 = 1-fit$cif0
    x = 'bottomleft'
    if(is.null(plot.configs[["ylab"]])) ylab = 'Survival probability' else ylab=plot.configs[["ylab"]]
  } else {
    cif1 = fit$cif1
    cif0 = fit$cif0
    x = 'topleft'
    if(is.null(plot.configs[["ylab"]])) ylab = 'Cumulative incidence' else ylab=plot.configs[["ylab"]]
  }
  col1 = plot.configs[["col"]][1]
  col0 = plot.configs[["col"]][2]

  if (!is.null(xlim)){
    ti = (tt>=xlim[1])&(tt<=xlim[2])
    tt = tt[ti]
    cif1 = cif1[ti]
    cif0 = cif0[ti]
  } else {
    ti = rep(TRUE,length(tt))
  }

  #---- main plot ----#
  plot(tt,cif1,type='s',col=col1,lwd=plot.configs[["lwd"]], lty=plot.configs[["lty"]],
    main=stname, xlab=xlab,ylab=ylab,ylim=ylim,...)
  points(tt,cif0,type='s',col=col0,lwd=plot.configs[["lwd"]], lty=plot.configs[["lty"]])

  if (!is.null(conf.int)){
    z = -qnorm((1-conf.int)/2)
    se1 = fit$se1
    se0 = fit$se0
    points(tt,cif1+se1[ti]*z,type='s',lty=plot.configs[["ci.lty"]],lwd=plot.configs[["ci.lwd"]],col=col1)
    points(tt,cif1-se1[ti]*z,type='s',lty=plot.configs[["ci.lty"]],lwd=plot.configs[["ci.lwd"]],col=col1)
    points(tt,cif0+se0[ti]*z,type='s',lty=plot.configs[["ci.lty"]],lwd=plot.configs[["ci.lwd"]],col=col0)
    points(tt,cif0-se0[ti]*z,type='s',lty=plot.configs[["ci.lty"]],lwd=plot.configs[["ci.lwd"]],col=col0)
  }

  if(!plot.configs[["show.p.value"]]) legend(x,legend=plot.configs[["legend"]],cex=plot.configs[["legend.cex"]],col=c(col1,col0),lwd=rep(plot.configs[["lwd"]],2),lty=rep(plot.configs[["lty"]],2))
    else legend(x,legend=c(plot.configs[["legend"]], paste0('P = ', round(p,3))),cex=plot.configs[["legend.cex"]],
                col=c(col1,col0, NA),lwd=c(rep(plot.configs[["lwd"]],2,),NA),lty=c(rep(plot.configs[["lty"]],2),NA))
}

