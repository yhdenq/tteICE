#' @title
#' Plot estimated treatment effects
#'
#' @description This function plots the estimated treatment effect,
#' defined as the difference in potential cumulative incidences under treated and control groups,
#' along with pointwise confidence intervals.
#'
#' @param fit
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param decrease
#' A logical value indicating the type of curve difference to display.
#' If \code{decrease = FALSE} (default), the difference in cumulative
#' incidence functions (CIFs) is plotted. If \code{decrease = TRUE},
#' the difference in survival functions is plotted instead.
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
#' Defaults to \code{ylim = c(-1, 1)}.
#'
#' @param plot.configs
#' A named \code{list} of additional plot configurations. Common entries include:
#'
#' \itemize{
#'     \item \code{ylab}: character, label for the y-axis (default: \code{ylab=NULL}, use the default label).
#'     \item \code{main}: character, title for the plot (default: \code{main=NULL}, use the default label).
#'     \item \code{lty}: line type for effect curve (default: \code{lty=1}).
#'     \item \code{lwd}: line width for effect curve (default: \code{lwd=2}).
#'     \item \code{col}: line color for effect curve (default: \code{col="black"}).
#'     \item \code{add.null.line}: logical, whether to draw a horizontal line at 0 (default: \code{add.null.line=TRUE}, add the null line).
#'     \item \code{null.line.lty}: line type for horizontal line at 0 (default: \code{null.line.lty=2}.
#'     \item \code{ci.lty}: line type for confidence interval curves (default: \code{ci.lty=5}).
#'     \item \code{ci.lwd}: line width for confidence interval curves (default: \code{ci.lwd=1.5}).
#'     \item \code{ci.col}: line color for confidence interval curves (default: \code{ci.col="darkgrey"}).
#'   }
#'
#' @param ... Additional graphical arguments passed to function \code{\link{plot.default}} or function \code{\link{curve}}
#'
#' @importFrom graphics plot abline points legend
#' @importFrom stats qnorm
#' @importFrom stats binomial fitted glm.fit reformulate
#'
#' @seealso
#' \code{\link[graphics]{plot.default}},
#' \code{\link[graphics]{points}},
#' \code{\link[graphics]{curve}},
#' \code{\link[tteICE]{plot.tteICE}}
#'
#' @examples
#' ## Load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' bmt$A = A
#'
#' ## simple model fitting and plotting
#' library(survival)
#' fit = tteICE(Surv(t2,d4,type = "mstate")~A, data=bmt)
#' plot_ate(fit)
#'
#' ## model fitting using competing risk data
#' fit1 = surv.tteICE(A, bmt$t2, bmt$d4, 'composite')
#'
#' ## Plot asymptotic confidence intervals based on explicit formulas
#' plot_ate(fit1, ylim=c(-0.4,0.4))
#'
#' ## Plot bootstrap confidence intervals
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, 'natural', nboot=50) ## SE=0??
#' plot_ate(fit2, ylim=c(-0.4,0.4))
#'
#' ## Model with semicompeting risk data
#' fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#'
#' ## Plot asymptotic confidence intervals based on explicit formulas
#' plot_ate(fit3, ylim=c(-0.4,0.4),
#'          plot.configs=list(add.null.line=FALSE))
#'
#' ## Plot bootstrap confidence intervals
#' fit4 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2,
#'                   "composite", nboot=50)
#'
#' plot_ate(fit4, ylim=c(-0.4,0.4),
#'          plot.configs=list(add.null.line=FALSE, lty=2, main=""))
#'
#' @return Plot the average treatment effect (ATE) results from a tteICE object
#' @keywords internal

plot_ate <- function(fit,decrease=FALSE,conf.int=.95,
  xlab='Time',ylim=c(-1,1),xlim=NULL,
  plot.configs=list(ylab=NULL, main=NULL,
                    lty=1, lwd=2, col="black",
                    add.null.line=TRUE, null.line.lty=2,
                    ci.lty=5, ci.lwd=1.5, ci.col="darkgrey"),...){

  #---- validate input ----#
  adj <- .plot_ate_validate(fit, decrease, conf.int, xlab, xlim, ylim)
  conf.int <- adj$conf.int

  #---- set configurations ----#
  default.configs <- list(
    ylab=NULL,main=NULL,
    lty=1, lwd=2, col="black",
    add.null.line=TRUE, null.line.lty=2,
    ci.lty=5, ci.lwd=1.5, ci.col="darkgrey")
  plot.configs <- modifyList(default.configs, plot.configs)

  #---- plot ----#
  if(is.null(plot.configs[["main"]])){
  if (fit$strategy=='treatment') stname = 'Treatment policy strategy'
  if (fit$strategy=='composite') stname = 'Composite variable strategy'
  if (fit$strategy=='natural') stname = 'Hypothetical I (natural) strategy'
  if (fit$strategy=='removed') stname = 'Hypothetical II (removed) strategy'
  if (fit$strategy=='whileon') stname = 'While on treatment strategy'
  if (fit$strategy=='principal') stname = 'Principal stratum strategy'
} else {stname=plot.configs[["main"]]}

  tt = fit$time
  dcif = fit$ate
  se = fit$se
  ciu = dcif - qnorm((1-conf.int)/2)*se
  cil = dcif + qnorm((1-conf.int)/2)*se
  if (is.null(plot.configs[["ylab"]])) ylab = 'Difference in CIFs' else ylab=plot.configs[["ylab"]]

  if (decrease==TRUE){
    dcif = -dcif
    ciu = -ciu
    cil = -cil
    if(is.null(plot.configs[["ylab"]])) ylab = 'Difference in Survivals' else ylab=plot.configs[["ylab"]]
  }
  if (!is.null(xlim)) {
    ti = (tt>=xlim[1])&(tt<=xlim[2])
    tt = tt[ti]
    dcif = dcif[ti]
    ciu = ciu[ti]
    cil = cil[ti]
  }

  #---- main plot ----#
  plot(tt,dcif,type='s',main=stname,ylim=ylim,xlab=xlab,ylab=ylab,lty=plot.configs[["lty"]], lwd=plot.configs[["lwd"]],col=plot.configs[["col"]],...)
  points(tt,ciu,type='s',lty=plot.configs[["ci.lty"]],lwd=plot.configs[["ci.lwd"]],col=plot.configs[["ci.col"]])
  points(tt,cil,type='s',lty=plot.configs[["ci.lty"]],lwd=plot.configs[["ci.lwd"]],col=plot.configs[["ci.col"]])
  if(plot.configs[["add.null.line"]]) abline(h=0,lty=plot.configs[["null.line.lty"]])
}
