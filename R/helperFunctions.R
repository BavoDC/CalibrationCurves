#' Infix operator to run background jobs
#'
#' This infix operator can be used to create a background job in RStudio/Posit and, once completed, the value of rhs is
#' assigned to lhs.
#'
#' @param lhs the object that the rhs value is assigned to
#' @param rhs the value you want to assign to lhs
#'
#' @return prints the ID of the background job in the console and, once completed, the value of lhs is assigned to rhs
#'
#' @examples
#' # Can only be executed in Rstudio
#' \dontrun{x %<=% rnorm(1e7)}
`%<=%` <- function(lhs, rhs) {
  Pkgs = names(sessionInfo()$otherPkgs)
  lhs  = as.character(enquote(substitute(lhs))[2])
  rhs  = as.character(enquote(substitute(rhs))[2])
  Job  = c(if(is.null(Pkgs)) NULL else paste0("LibraryM(", paste0(Pkgs, collapse = ", "), ")"),
           paste0(lhs, " <- ", rhs))
  tmpR = tempfile(fileext = ".R")
  writeLines(Job, tmpR)
  if(!file.exists(tmpR))
    stop("Temporary R script not created")
  jobRunScript(tmpR, exportEnv = "R_GlobalEnv", importEnv = TRUE)
}

#' Infix operator to run background jobs
#'
#' This infix operator can be used to create a background job for a block of code in RStudio/Posit and, once completed,
#' all objects created in the block of code are imported into the global environment.
#'
#' @param lhs not used, see details and examples
#' @param rhs the block of code that you want to run
#'
#' @return prints the ID of the background job in the console and, once completed, the objects created in the block of code
#' are imported into the global environment
#'
#' @details
#' You can use this infix operator in two different ways. Either you set the left-hand side to \code{NULL} or you use the syntax
#' \code{`\%{}\%`  ({BlockOfCode})}
#'
#' @examples
#' # Can only be executed in Rstudio
#' \dontrun{
#' NULL %{}% {
#'  x = rnorm(1e7)
#'  y = rnorm(1e7)
#' }
#' `%{}%` ({
#'  x = rnorm(1e7)
#'  y = rnorm(1e7)
#' })
#' }
`%{}%` <- function(lhs, rhs) {
  Pkgs = names(sessionInfo()$otherPkgs)
  Job  = c(if(is.null(Pkgs)) NULL else paste0("LibraryM(", paste0(Pkgs, collapse = ", "), ")"),
           as.character(enquote(substitute(rhs))[2]))
  tmpR = tempfile(fileext = ".R")
  writeLines(Job, tmpR)
  if(!file.exists(tmpR))
    stop("Temporary R script not created")
  jobRunScript(tmpR, exportEnv = "R_GlobalEnv", importEnv = TRUE)
}


#' Function to load multiple packages at once
#'
#'
#' @param ... the packages that you want to load
#'
#' @return invisible \code{NULL}
#'
#' @examples
#' LibraryM(CalibrationCurves)
LibraryM <- function(...) {
  libs = as.list(substitute(list(...)))[-1L]
  if(!is.character(libs)) {
    if(length(libs[[1]]) > 1)
      libs = sapply(substitute(libs)[-1], as.character)
    else
      libs = as.character(substitute(libs))
  }
  lapply(libs, library, character.only = TRUE)
  invisible(NULL)
}

## loess.as function from the fANCOVA package
## loess with Automatic Smoothing Parameter Selection
loess.as <-
  function(x, y, degree=1, criterion=c("aicc", "gcv"),
           family = c("gaussian", "symmetric"), user.span=NULL, plot=FALSE, ...)
  {
    criterion <- match.arg(criterion)
    family <- match.arg(family)
    x <- as.matrix(x)

    if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
    if (!is.numeric(x)) stop("argument 'x' must be numeric!")
    if (!is.numeric(y)) stop("argument 'y' must be numeric!")
    if (any(is.na(x))) stop("'x' contains missing values!")
    if (any(is.na(y))) stop("'y' contains missing values!")
    if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span)))
      stop("argument 'user.span' must be a numerical number!")
    if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
    if(length(y) < 3) stop("not enough observations!")

    data.bind <- data.frame(x=x, y=y)
    if (ncol(x) == 1) {
      names(data.bind) <- c("x", "y")
    } else { names(data.bind) <- c("x1", "x2", "y") }

    opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){
      as.crit <- function (x) {
        span <- x$pars$span
        traceL <- x$trace.hat
        sigma2 <- sum(x$residuals^2 ) / (x$n-1)
        aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
        gcv <- x$n*sigma2 / (x$n-traceL)^2
        result <- list(span=span, aicc=aicc, gcv=gcv)
        return(result)
      }
      criterion <- match.arg(criterion)
      fn <- function(span) {
        mod <- update(model, span=span)
        as.crit(mod)[[criterion]]
      }
      result <- optimize(fn, span.range)
      return(list(span=result$minimum, criterion=result$objective))
    }

    if (ncol(x)==1) {
      if (is.null(user.span)) {
        fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
        span1 <- opt.span(fit0, criterion=criterion)$span
      } else {
        span1 <- user.span
      }
      fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
    } else {
      if (is.null(user.span)) {
        fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
        span1 <- opt.span(fit0, criterion=criterion)$span
      } else {
        span1 <- user.span
      }
      fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family = family, data=data.bind,...)
    }
    if (plot){
      if (ncol(x)==1) {
        m <- 100
        x.new <- seq(min(x), max(x), length.out=m)
        fit.new <- predict(fit, data.frame(x = x.new))
        plot(x, y, col="lightgrey", xlab="x", ylab="m(x)", ...)
        lines(x.new,fit.new, lwd=1.5, ...)
      } else {
        m <- 50
        x1 <- seq(min(data.bind$x1), max(data.bind$x1), len=m)
        x2 <- seq(min(data.bind$x2), max(data.bind$x2), len=m)
        x.new <- expand.grid(x1=x1, x2=x2)
        fit.new <- matrix(predict(fit, x.new), m, m)
        persp(x1, x2, fit.new, theta=40, phi=30, ticktype="detailed", xlab="x1", ylab="x2", zlab="y", col="lightblue", expand=0.6)
      }
    }
    return(fit)
  }


## Logit and inverse logit functions
Ilogit <- function(x) binomial()$linkinv(x)
Logit  <- function(x) binomial()$linkfun(x)
