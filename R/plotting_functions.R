#' Plot elicited data, fitted marginals or model output
#' 
#' @param Z list object that contains matrix \code{theta} of elicitations,
#'   character \code{link} and character \code{target}
#' @param X design matrix (can be \code{NULL}, unless modelled output is
#'   requested)
#' @param design.pt single integer that denotes design point of interest
#' @param elicited.fractiles logical, plot vertical lines for elicited 
#'   fractiles?
#' @param fitted.fractiles logical, plot vertical lines for fitted marginal 
#'   fractiles? Alternatively, a numeric vector of arbitrary fractiles to be
#'   plotted from the fitted elicitation distribution. If \code{TRUE} then the
#'   fractiles corresponding to the median, upper and lower level central CI
#'   are plotted
#' @param fitted.curve, logical plot fitted elicitation density?
#' @param CI.prob numeric scalar, locally specified probability assigned to the
#'   elicited central credible interval of the current design point. Defaults to
#'   \code{NULL} in which case the global value initially assigned by
#'   \code{designLink()} is used
#' @param estimated.probs numeric vector of values for which estimated
#'  probabilities are to be estimated from the fitted elicitation 
#'  distribution for the target theta. Default is \code{NULL}. 
#'  The result is output to the console.
#' @param modelled.fractiles logical, plot vertical lines for modelled 
#'   fractiles?
#' @param modelled.curve logical, plot modelled conditional mean prior density?
#' @param cumul.prob.bounds numeric vector of length two, giving plot bounds by 
#'   cumulative probability. This argument is ignored if there is not enough data
#'   to fit a parametric distribution or if \code{theta.bounds} is not \code{NULL}
#' @param theta.bounds numeric vector giving support of response for plotting
#'   purposes (can be \code{NULL}). This will overwrite \code{cumul.prob.bounds},
#'   if applicable
#' @param ylim.max numeric maximum value of y-axis (can be \code{NULL})
#' @param xlog logical log x-axis
#' @param design.table logical include design dataframe, elicited fractiles and 
#'   modelled or fitted fractiles
#' @param n.pts numeric giving number of point to evalate density curve (if 
#'   plotted)
#' @param fit.method character method used to fit conditional mean prior:
#'   \code{KL}, \code{moment}
#' @return a plot to the current device. See \code{dev.cur()} to check.
#' @examples 
#' # design matrix: two scenarios
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) 
#' rownames(X) <- c("scenario1", "scenario2")
#' colnames(X) <- c("covariate1", "covariate2")
#' 
#' # logit link
#' # central credible intervals with probability = 1/2
#' Z <- designLink(design = X, link = "logit", CI.prob = 0.5)
#' 
#' # 1st design point
#' # no elicited fractiles
#' indirect::plotDesignPoint(Z, design.pt = 1) 
#' # elicited median
#' Z <- indirect::elicitPt(Z, design.pt = 1, 
#'   lower.CI.bound = NA,
#'   median = 0.4,
#'   upper.CI.bound = NA,
#'   CI.prob = NULL)
#' indirect::plotDesignPoint(Z, design.pt = 1,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1))
#' # lower and upper quartiles and median
#' Z <- indirect::elicitPt(Z, design.pt = 1, 
#'   lower.CI.bound = 0.2,
#'   median = 0.4,
#'   upper.CI.bound = 0.6,
#'   comment = "Completed.")
#' indirect::plotDesignPoint(Z, design.pt = 1,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1),
#'   fitted.fractiles = TRUE, fitted.curve = TRUE)
#' indirect::plotDesignPoint(Z, design.pt = 1,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1),
#'   fitted.fractiles = c(1/10, 1/4, 1/2, 3/4, 9/10), 
#'   fitted.curve = TRUE)   
#'   
#' # second design point   
#' # central credible intervals with probability = 1/3 
#' # elicit upper and lower tertiles
#' Z <- elicitPt(Z, design.pt = 2,
#'   lower.CI.bound = 0.1,
#'   upper.CI.bound = 0.3,
#'   CI.prob = 1/3,
#'   comment = "Switched to tertiles.")
#' indirect::plotDesignPoint(Z, design.pt = 2,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1))   
#' indirect::plotDesignPoint(Z, design.pt = 2,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1),
#'   fitted.fractiles = TRUE, fitted.curve = TRUE)
#' indirect::plotDesignPoint(Z, design.pt = 1,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1),
#'   fitted.fractiles = c(1/10, 1/3, 1/2, 2/3, 9/10), 
#'   fitted.curve = TRUE) 
plotDesignPoint <- function(Z, X = NULL, design.pt = NULL,   
  elicited.fractiles = TRUE, fitted.fractiles = FALSE, fitted.curve = FALSE,
  CI.prob = NULL, estimated.probs = NULL,
  modelled.fractiles = FALSE, modelled.curve = FALSE, 
  cumul.prob.bounds = c(0.05, 0.95), theta.bounds = NULL, ylim.max = NULL, xlog = FALSE,
  design.table = TRUE, n.pts = 101, fit.method = "KL") {
  
  if (is.null(design.pt)) stop("specify design point")
  
  theta <- Z$theta
  link <- Z$link
  design <- Z$design
  
  # elicited data could be missing but design point must be specified
  dat <- theta[design.pt, 1:(ncol(theta) - 1)]
  
  # find fractiles
  ci.prob <- theta[design.pt, "CI_prob"]
  cumul.probs <- c((1 - ci.prob)/2, 0.5, 1 - (1 - ci.prob)/2)
  # if specific fitted fractiles aren't requested then use fractiles
  # for lower and upper central CI and median
  if (is.logical(fitted.fractiles)) {
    cumul.probs.fit <- cumul.probs
  } else if (is.numeric(fitted.fractiles)) {
    cumul.probs.fit <- fitted.fractiles
    fitted.fractiles <- TRUE
  } else {
    stop("fitted.fractiles should be either a logical value or numeric")
  }
  
  if (fitted.fractiles || fitted.curve || modelled.fractiles || modelled.curve) {
    
    x.axis.mf <- matrix(NA, nrow = 2, ncol = 2,
      dimnames = list(c("min", "max"), c("fitted", "modelled")))
    
    # transform elicitations to moments  
    mV.ls <- mV(Z, fit.method = fit.method)
    
    fitted.mean <- mV.ls$m[design.pt]
    fitted.sd <- sqrt(mV.ls$V[design.pt, design.pt])
    
    if (!is.null(cumul.prob.bounds)) {
      x.axis.mf[ , "fitted"] <- qnorm(cumul.prob.bounds, fitted.mean, fitted.sd)
    } else {
      x.axis.mf[ , "fitted"] <- qnorm(cumul.probs.fit, fitted.mean, fitted.sd)     
    }
    
    if (modelled.fractiles || modelled.curve) {
      
      if (is.null(X)) stop("design matrix X required")
      
      # estimate mu and Sigma for given elicitations
      muSigma.ls <- muSigma(Z, X, fit.method = fit.method, wls.method = "default")
      
      # predictions for g(theta(design.pt))
      modelled.mean <- X[design.pt, ]%*%muSigma.ls$mu
      modelled.sd <- sqrt(X[design.pt, , drop = FALSE]%*%muSigma.ls$Sigma%*%t(X[design.pt, , drop = FALSE]))
      
      if (!is.null(cumul.prob.bounds)) {
        x.axis.mf[ , "modelled"] <- qnorm(cumul.prob.bounds, modelled.mean, modelled.sd)
      } else {
        x.axis.mf[ , "modelled"] <- qnorm(cumul.probs.fit, modelled.mean, modelled.sd)    
      }
      
    }

    x.axis.mf <- c(min(x.axis.mf["min", ], na.rm = TRUE), max(x.axis.mf["max", ], na.rm = TRUE))

    # limits for plot's x axis
    # theta.bounds takes priority over cumul.prob.bounds if not NULL
    x.lim.cumul.prob <- FALSE
    if (!is.null(theta.bounds)) {
      x.axis.lim <- theta.bounds 
    } else if (!is.null(x.axis.mf)) {
      x.axis.lim <- x.axis.mf
      x.lim.cumul.prob <- TRUE
    } else {
      stop("not enough information provided for plotting fitted values")
    }
      
    # transform by link function
    switch(link,
      log = {
        if (x.lim.cumul.prob) x.axis.lim <- exp(x.axis.lim)
        if (fitted.fractiles) dat.fit <- exp(qnorm(cumul.probs.fit, fitted.mean, fitted.sd))
        if (modelled.fractiles) {
          dat.mod <- exp(qnorm(cumul.probs.fit, modelled.mean, modelled.sd))
        }
        if (fitted.curve || modelled.curve) {          
          if (xlog) {
            x <- seq(from = log(x.axis.lim[1]), to = log(x.axis.lim[2]), length.out = n.pts)
            x <- exp(x)
          } else {
            x <- seq(from = x.axis.lim[1], to = x.axis.lim[2], length.out = n.pts)
          }
          if (fitted.curve) dens.fit <- dlnorm(x, fitted.mean, fitted.sd)
          if (modelled.curve) {
            dens.mod <- dlnorm(x, modelled.mean, modelled.sd)
          }
        }
      },
      cloglog = {
        if (x.lim.cumul.prob) x.axis.lim <- 1 - exp(-exp(x.axis.lim))
        if (fitted.fractiles) dat.fit <- 1 - exp(-exp(qnorm(cumul.probs.fit, fitted.mean, fitted.sd)))
        if (modelled.fractiles) {
          dat.mod <- 1 - exp(-exp(qnorm(cumul.probs.fit, modelled.mean, modelled.sd)))
        }
        if (fitted.curve || modelled.curve) {          
          if (xlog) {
            x <- seq(from = log(x.axis.lim[1]), to = log(x.axis.lim[2]), length.out = n.pts)
            x <- exp(x)
          } else {
            x <- seq(from = x.axis.lim[1], to = x.axis.lim[2], length.out = n.pts)
          }
          if (fitted.curve) dens.fit <- dGompertzNorm(x, fitted.mean, fitted.sd)
          if (modelled.curve) {
            dens.mod <- dGompertzNorm(x, modelled.mean, modelled.sd)
          }
        }
      },
      logit = {
        if (x.lim.cumul.prob) x.axis.lim <- exp(x.axis.lim)/(1 + exp(x.axis.lim))
        if (fitted.fractiles) dat.fit <- exp(qnorm(cumul.probs.fit, fitted.mean, fitted.sd))/(1 + exp(qnorm(cumul.probs.fit, fitted.mean, fitted.sd)))
        if (modelled.fractiles) {
          dat.mod <- exp(qnorm(cumul.probs.fit, modelled.mean, modelled.sd))/(1 + exp(qnorm(cumul.probs.fit, modelled.mean, modelled.sd)))
        }    
        if (fitted.curve || modelled.curve) {          
          if (xlog) {
            x <- seq(from = log(x.axis.lim[1]), to = log(x.axis.lim[2]), length.out = n.pts)
            x <- exp(x)
          } else {
            x <- seq(from = x.axis.lim[1], to = x.axis.lim[2], length.out = n.pts)
          }
          if (fitted.curve) dens.fit <- dLogitNorm(x, fitted.mean, fitted.sd)
          if (modelled.curve) {
            dens.mod <- dLogitNorm(x, modelled.mean, modelled.sd)
          }
        }  
      },
      identity = {
        if (fitted.fractiles) dat.fit <- qnorm(cumul.probs.fit, fitted.mean, fitted.sd)
        if (modelled.fractiles) {
          dat.mod <- qnorm(cumul.probs.fit, modelled.mean, modelled.sd)
        }
        if (fitted.curve || modelled.curve) {          
          if (xlog) {
            x <- seq(from = log(x.axis.lim[1]), to = log(x.axis.lim[2]), length.out = n.pts)
            x <- exp(x)
          } else {
            x <- seq(from = x.axis.lim[1], to = x.axis.lim[2], length.out = n.pts)
          }
          if (fitted.curve) dens.fit <- dnorm(x, fitted.mean, fitted.sd)
          if (modelled.curve) {
            dens.mod <- dnorm(x, modelled.mean, modelled.sd)
          }    
        }
      },
      stop("link unspecified")
    )
    
    if (!fitted.fractiles) dat.fit <- rep(NA, 3)
    
    if (is.null(ylim.max)) {
      if (!modelled.curve && fitted.curve) max.y <- max(dens.fit, na.rm = TRUE)
      if (modelled.curve && !fitted.curve) max.y <- max(dens.mod, na.rm = TRUE)
      if (modelled.curve && fitted.curve) max.y <- max(dens.fit, dens.mod, na.rm = TRUE)
      if (!modelled.curve && !fitted.curve) max.y <- 1 # arbitrary value
      y.axis.lim <- c(0, max.y)
    } else {
      y.axis.lim <- c(0, ylim.max)
    }
    
  } else {
    
    # elicited data only - no fitting
    # info for design table, if requested: no additional fitted fractiles
    cumul.probs.fit <- cumul.probs
    dat.fit <- rep(NA, 3)
    
    # plotting limits
    if (!is.null(theta.bounds)) {
      x.axis.lim <- theta.bounds 
    } else if (all(!is.na(theta[design.pt, ]))) {
      x.axis.lim <- c(theta[design.pt, c("lower", "upper")])
      x.axis.lim[1] <- x.axis.lim[1] - x.axis.lim[1]/10
      x.axis.lim[2] <- x.axis.lim[2] + x.axis.lim[2]/10
    } else {
      x.axis.lim <- c(-10^6, 10^6)
    }
    if (is.null(ylim.max)) {
      y.axis.lim <- c(0, 1)
    } else {
      y.axis.lim <- c(0, ylim.max)
    }
    
  }
  
  # create plot
  if (design.table) layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2, byrow = TRUE))
  par(mar = c(4, 4, 2, 0) + 0.1)
  if (xlog) {
    xlog <- "x"
  } else {
    xlog <- ""
  }
  
  plot(x.axis.lim, y.axis.lim, type = "n", xlab = "", ylab = "Density", main = Z$target.caption, log = xlog, axes = FALSE)
  axis(1, cex.axis = 2)
  axis(2)
  box()
  if (elicited.fractiles) {
    abline(v = dat, lty = 2)
    for (i in 1:3) {
      f <- MASS::fractions(cumul.probs[i])
      f.num <- strsplit(as.character(f), split = "/")[[1]][1]
      f.den <- strsplit(as.character(f), split = "/")[[1]][2]
      text(dat[i], 0.98*max(y.axis.lim), substitute(f[z/zz], list(z = f.num, zz = f.den)), pos = 4)
    }
  }
  if (fitted.curve) lines(x, dens.fit, col = "blue")
  if (fitted.fractiles) {
    abline(v = dat.fit, lty = 2, col = "blue")
    for (i in 1:length(dat.fit)) {
      f <- MASS::fractions(cumul.probs.fit[i])
      f.num <- strsplit(as.character(f), split = "/")[[1]][1]
      f.den <- strsplit(as.character(f), split = "/")[[1]][2]
      text(dat.fit[i], 0, substitute(f[z/zz], list(z = f.num, zz = f.den)), pos = 4, col = "blue")
    }
  }
  if (modelled.curve) lines(x, dens.mod, col = "orange")
  if (modelled.fractiles) {
    abline(v = dat.mod, lty = 2, col = "orange")
    for (i in 1:length(dat.mod)) {
      f <- MASS::fractions(cumul.probs.fit[i])
      f.num <- strsplit(as.character(f), split = "/")[[1]][1]
      f.den <- strsplit(as.character(f), split = "/")[[1]][2]
      text(dat.mod[i], 0, substitute(f[z/zz], list(z = f.num, zz = f.den)), pos = 4, col = "orange")
    }
  } 
  
  # add design table if specified
  if (design.table) {
    designpoint <- matrix(signif(design[design.pt, ], 4), nrow = length(design[design.pt, ]), ncol = 1,
      dimnames = list(Covariate = colnames(design), "Value"))
    gplots::textplot(designpoint, valign = "top", mar = c(1, 1, 2, 1))
    title(paste("Scenario:", row.names(design)[design.pt]))
    box("figure")
    all.cumul.probs <- sort(unique(c(cumul.probs, cumul.probs.fit)))
    elicit.ind <- match(cumul.probs, all.cumul.probs)
    fitted.ind <- match(cumul.probs.fit, all.cumul.probs)
    fractiles <- matrix(NA, nrow = length(all.cumul.probs), ncol = 2)
    fractiles[elicit.ind, 1] <- signif(dat, 3)
    if (modelled.fractiles) {
      fractiles[fitted.ind, 2] <- signif(dat.mod, 3)
      dimnames(fractiles) <- list(Fractiles = as.character(MASS::fractions(all.cumul.probs)), c("Elicited", "Modelled"))
    } else {
      fractiles[fitted.ind, 2] <- signif(dat.fit, 3)
      dimnames(fractiles) <- list(Fractiles = as.character(MASS::fractions(all.cumul.probs)), c("Elicited", "Fitted"))
    }
    gplots::textplot(fractiles, valign = "top", mar = c(1, 1, 2, 1))
    title("Fractiles")
    box("figure")
  }
  
  if (!is.null(estimated.probs)) print(pdist(estimated.probs, Z, design.pt, fit.method = fit.method))
  
}

#' density for Gompertz transformed univariate Gaussian
#' 
#' @param x, numeric real
#' @param mu, numeric real
#' @param sigma, numeric real positive
#' @return tranformed density on support (0, 1)
#' @examples 
#' mu <- -1
#' sigma <- 1
#' z <- rnorm(10000, mu, sigma)
#' hist(1 - exp(-exp(z)), freq = FALSE)
#' curve(dGompertzNorm(x, mu = mu, sigma = sigma), col = 'red', add = TRUE)
#' integrate(function(x) dGompertzNorm(x, mu = mu, sigma = sigma), lower = 0, upper = 1) # equals 1
dGompertzNorm <- function(x, mu, sigma) {
  dnorm(log(-log(1 - x)), mu, sigma)/((x - 1)*log(1 - x))
}

#' density for logit transformed univariate Gaussian
#' 
#' @param x, numeric real
#' @param mu, numeric real
#' @param sigma, numeric real positive
#' @return tranformed density on support (0, 1)
#' @examples 
#' mu <- -1
#' sigma <- 1
#' z <- rnorm(10000, mu, sigma)
#' hist(exp(z)/(1 + exp(z)), freq = FALSE)
#' curve(dLogitNorm(x, mu = mu, sigma = sigma), col = 'red', add = TRUE)
#' integrate(function(x) dLogitNorm(x, mu = mu, sigma = sigma), lower = 0, upper = 1) # equals 1
dLogitNorm <- function(x, mu, sigma) {
  dnorm(log(x/(1 - x)), mu, sigma)/(x*(1 - x))
}

