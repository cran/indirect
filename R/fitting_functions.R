#' Helper function that translates elicited quantiles of target into independent
#' conditional means normal prior for a defined inverse link function.
#' 
#' The default for \code{fit.method} is option \code{KL}. This option uses an 
#' objective function that minimises a discretised directed divergence from a 
#' cumulative distribution implied by raw elicited fractiles to a normal 
#' conditional mean prior for the linear predictor. An alterative method 
#' \code{moment} assigns the location parameter of the normal conditional mean 
#' prior to the elicited median on the linear predictor scale. The variance 
#' parameter is estimated as \eqn{V = ((g(f_u) - g(f_l)/(qnorm(u) -
#' qnorm(l)))^2}, where \eqn{l} is the probability associated with the fractile
#' \eqn{f_l} that defines the lower bound for the central credible interval and 
#' \eqn{u} is the probability associated with the fractile \eqn{f_u} that
#' defines the upper bound for the central credible interval. This is also used
#' to initialise the optimisation for the \code{KL} method. Another optimsation 
#' method that minimises the sum of squares is also available as method 
#' \code{SS}. See the vignette for more details on the choice of objective 
#' function for \code{KL} and \code{SS}.
#' 
#' @param Z list object that contains matrix \code{theta} of elicitations and
#'   character \code{link}, see \code{\link{plotDesignPoint}}
#' @param fit.method character, \code{moment}, \code{KL}, \code{SS}. Default is
#'   \code{KL}.
#' @return A list with vector of means \code{m} and diagonal covariance matrix 
#'   \code{V}.
mV <- function(Z, fit.method = "KL") {
  theta <- Z$theta
  g.theta <- theta
  link <- Z$link
  lo.CI <- (1 - theta[ , "CI_prob"])/2
  hi.CI <- lo.CI + theta[ , "CI_prob"]
  
  # any misspecified link functions?
  if (!any(link%in%c("log", "cloglog", "logit", "identity"))) stop("unsupported link function")  
  
  switch(link, 
    log = g.theta[ , c("lower", "median", "upper")] <- log(theta[ , c("lower", "median", "upper")]),
    cloglog = g.theta[ , c("lower", "median", "upper")] <- log(-log(1 - theta[ , c("lower", "median", "upper")])),
    logit = g.theta[ , c("lower", "median", "upper")] <- log(theta[ , c("lower", "median", "upper")]) - log(1 - (theta[ , c("lower", "median", "upper")])),
    identity = g.theta[ , c("lower", "median", "upper")] <- theta[ , c("lower", "median", "upper")],
    stop("specify link function")
  )
  
  # initialise from "moments"
  m.init <- g.theta[ , 'median']
  V.init <- ((g.theta[ , "upper"] - g.theta[ , "lower"])/(qnorm(hi.CI) - qnorm(lo.CI)))^2
  m.init[is.na(m.init)] <- (g.theta[is.na(m.init), "upper"] + g.theta[is.na(m.init), "lower"])/2
  
  # allow for only partially specified bounds for central CI
  uppers.and.medians <- (!is.na(g.theta[ , "median"])) & (is.na(g.theta[ , "lower"]) & (!is.na(g.theta[ , "upper"])))
  lowers.and.medians <- (!is.na(g.theta[ , "median"])) & (is.na(g.theta[ , "upper"]) & (!is.na(g.theta[ , "lower"])))  
  if (any(uppers.and.medians)) {
    V.init[uppers.and.medians] <- ((g.theta[uppers.and.medians, "upper"] - g.theta[uppers.and.medians, "median"])/(qnorm(hi.CI)[uppers.and.medians] - qnorm(0.5)))^2
  }
  if (any(lowers.and.medians)) {
    V.init[lowers.and.medians] <- ((g.theta[lowers.and.medians, "median"] - g.theta[lowers.and.medians, "lower"])/(qnorm(0.5) - qnorm(lo.CI)[lowers.and.medians]))^2
  }
  
  # estimation
  if (fit.method == "moment") {
    m <- m.init
    V <- diag(V.init)
  } else if (fit.method == "KL" || fit.method == "SS") {
    # elicited probability intervals
    cumul.probs.elicit <- cbind(0, lo.CI, 0.5, hi.CI, 1)
    prob.intervals.elicit <- diff(t(cumul.probs.elicit))
    um.cumul.probs.elicit <- cbind(0, 0.5, hi.CI, 1)
    um.prob.intervals.elicit <- diff(t(um.cumul.probs.elicit))
    lm.cumul.probs.elicit <- cbind(0, lo.CI, 0.5, 1)
    lm.prob.intervals.elicit <- diff(t(lm.cumul.probs.elicit))
    lu.cumul.probs.elicit <- cbind(0, lo.CI, hi.CI, 1)
    lu.prob.intervals.elicit <- diff(t(lu.cumul.probs.elicit))
    
    lowers.and.uppers <- is.na(g.theta[ , "median"]) & (!is.na(g.theta[ , "lower"])) & (!is.na(g.theta[ , "upper"]))
    
    # fitted cumulative probabilities
    prob_fun <- function(x, elicited.fractiles) {
      xx <- c(x[1], sqrt(exp(x[2])))
      c(0, pnorm(elicited.fractiles, mean =  xx[1], sd = xx[2]), 1)
    }
    if (fit.method == "KL") {
      # objective function
      obj.fun <- function(x, elicited.fractiles, prob.intervals.elicit) {
        cumul.prob.fit <- prob_fun(x, elicited.fractiles)
        prob.intervals.fit <- diff(cumul.prob.fit)
        sum((log(prob.intervals.elicit) - log(prob.intervals.fit))*prob.intervals.elicit)
      }
    }
    if (fit.method == "SS") {
      obj.fun <- function(x, elicited.fractiles, prob.intervals.elicit) {
        cumul.prob.fit <- prob_fun(x, elicited.fractiles)
        prob.intervals.fit <- diff(cumul.prob.fit)
        sum(((prob.intervals.elicit) - (prob.intervals.fit))^2)
      }
    }
    m <- rep(NA, nrow(g.theta))
    V <- diag(m, nrow = nrow(g.theta))
    diag(V) <- NA
    # initialisation 
    pars <- cbind(m.init, log((V.init)))
    colnames(pars) <- c("m", "log.V")
    for (i in 1:nrow(g.theta)) {
      # optimisation initialised from moment method - no check for multimodality
      if (all(!is.na(g.theta[i, c("lower", "median", "upper")]))) {
        out <- optim(pars[i, ], function(x) obj.fun(x, g.theta[i, c("lower", "median", "upper")], 
                                                    prob.intervals.elicit[ , i]), method = "Nelder-Mead")
        m[i] <- out$par[1]
        V[i, i] <- exp(out$par[2])
      } else if (uppers.and.medians[i]) {
        out <- optim(pars[i, ], function(x) obj.fun(x, g.theta[i, c("median", "upper")], 
                                                    um.prob.intervals.elicit[ , i]), method = "Nelder-Mead")
        m[i] <- out$par[1]
        V[i, i] <- exp(out$par[2])
      } else if (lowers.and.medians[i]) {
        out <- optim(pars[i, ], function(x) obj.fun(x, g.theta[i, c("lower", "median")], 
                                                    lm.prob.intervals.elicit[ , i]), method = "Nelder-Mead")
        m[i] <- out$par[1]
        V[i, i] <- exp(out$par[2])
      } else if (lowers.and.uppers[i]) {
        out <- optim(pars[i, ], function(x) obj.fun(x, g.theta[i, c("lower", "upper")], 
                                                    lu.prob.intervals.elicit[ , i]), method = "Nelder-Mead")
        m[i] <- out$par[1]
        V[i, i] <- exp(out$par[2])
      }
    }
  } else {
    stop("specify fit.method argument")
  }
  
  list(m = m, V = V)

}

#' Function to estimate mean and covariance for unknown parameters 
#' \eqn{\beta}.
#' 
#' @param Z list of design points and link function that is an output of 
#'   function \code{designLink}
#' @param X model matrix for model formula and design points. The covariates 
#'   must correspond to the description of design points in \code{Z}, but can be
#'   transformed etc. If \code{NULL} then \code{X} will be coerced by applying 
#'   \code{as.matrix()} to \code{Z$design}. The matrix \code{X} should be full 
#'   rank when subsetted to the elicited design points. If a column of \code{X} 
#'   has the name \code{offset} then this column is treated as an offset during 
#'   estimation
#' @param fit.method character, \code{moment}, \code{KL}. See \code{\link{mV}}. Default
#'   is \code{KL}.
#' @param wls.method character giving the numerical solution method: \code{QR}, 
#'   using the QR decomposition, \code{SVD}, using the singular value 
#'   decomposition, or option \code{default} that uses \code{solve()}
#' @return list of \code{mu}, numeric vector of location parameters for the 
#'   normal prior; \code{Sigma}, the covariance matrix; and \code{log.like}, a 
#'   scalar
#' @examples 
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) # design
#' Z <- designLink(design = X)
#' Z <- elicitPt(Z, design.pt = 1,
#'   lower.CI.bound = -1,
#'   median = 0,
#'   upper.CI.bound = 1,
#'   comment = "The first completed elicitation scenario.")
#' Z <- elicitPt(Z, design.pt = 2,
#'   lower.CI.bound = -2,
#'   median = 1,
#'   upper.CI.bound = 2,
#'   comment = "The second completed elicitation scenario.")
#' prior <- muSigma(Z, X, fit.method = "KL")
#' prior$mu
#' prior$Sigma   
muSigma <- function(Z, X = NULL, fit.method = "KL", wls.method = "default") {
  
  theta <- Z$theta
  design <- Z$design
  if (is.null(X)) {
    warning("argument X unspecified: coercing design matrix from elicitation list object")
    X <- as.matrix(Z$design)
  }
  
  if (!all.equal(nrow(design), nrow(X), nrow(theta))) stop("at least one of design, X or theta has different dimensions")
  
  # subset to design points with at least two fractiles elicited
  elicited.design.points <- rowSums(!is.na(theta[ , c("lower", "median", "upper"), drop = FALSE])) >= 2

  if (sum(elicited.design.points) == 0) stop("no design points elicited")
  X.elicit <- X[elicited.design.points, , drop = FALSE]
  
  # conditional fits
  mV.ls <- mV(Z, fit.method = fit.method)
  
  # elicited moments
  m.elicit <- mV.ls$m[elicited.design.points]
  V.elicit <- mV.ls$V[elicited.design.points, elicited.design.points, drop = FALSE]
  
  # account for offset if present
  if (any(grepl("offset", colnames(X.elicit)))) {
    if (sum(grepl("offset", colnames(X.elicit))) > 1) stop("only one offset permitted")
    offset <- X.elicit[ , grep("offset", colnames(X.elicit))]
    m.elicit <- m.elicit - offset
    X.elicit <- X.elicit[ , -grep("offset", colnames(X.elicit))]
  }
  
  # check design matrix given elicited design points
  # assume full rank design matrix X
  checkX(X.elicit)

  # evaluate mu and Sigma
  if (wls.method == "QR") {
    # construct Xw and Yw (Seber and Lee 2003, Ch 11)
    DsqrtVe <- (1/sqrt(diag(V.elicit)))
    if (length(DsqrtVe) == 1) {
      DsqrtVe <- matrix(DsqrtVe)
    } else {
      DsqrtVe <- diag(DsqrtVe)
    }
    Xw <- DsqrtVe%*%X.elicit
    Yw <- DsqrtVe%*%m.elicit
    # QR decomposition
    decomp <- qr(Xw)
    mu <- backsolve(qr.R(decomp), qr.qty(decomp, Yw))
    # XX' = R'R Chambers (1971) Eq. 2.2
    # (XX')^-1 = Sigma
    Sigma <- chol2inv(qr.R(decomp))
  } else if (wls.method == "SVD") {
    # construct Xw and Yw (Seber and Lee 2003, Ch 11)
    DsqrtVe <- (1/sqrt(diag(V.elicit)))
    if (length(DsqrtVe) == 1) {
      DsqrtVe <- matrix(DsqrtVe)
    } else {
      DsqrtVe <- diag(DsqrtVe)
    }
    Xw <- DsqrtVe%*%X.elicit
    Yw <- DsqrtVe%*%m.elicit
    # now use SVD approach to OLS with full rank design matrix
    # Seber and Lee (2003) Ch 11
    # Xw = UDV' (V now referring to right singular vectors of Xw)
    # Note: U is n x p
    decomp <- svd(Xw, nu = ncol(Xw), nv = ncol(Xw))
    if (length(decomp$d) == 1) {
      Dm1 <- matrix((1/decomp$d))
      Dm2 <- matrix((1/decomp$d^2))
    } else {
      Dm1 <- diag(1/decomp$d)
      Dm2 <- diag(1/decomp$d^2)
    }
    mu <- decomp$v%*%Dm1%*%t(decomp$u)%*%Yw
    # X'X = VD^2V'
    # (X'X)^-1 = VD^-2V'
    Sigma <- decomp$v%*%Dm2%*%t(decomp$v)
  } else {
    Sigma <- solve(t(X.elicit)%*%solve(V.elicit)%*%X.elicit)
    mu <- Sigma%*%t(X.elicit)%*%solve(V.elicit)%*%m.elicit
  }

  # log like
  Xbeta <- X.elicit%*%mu
  ll <- -sum(elicited.design.points)/2*log(2*pi) - 0.5*(log(det(V.elicit))) - 0.5*t(m.elicit - Xbeta)%*%chol2inv(chol(V.elicit))%*%(m.elicit - Xbeta)  

  if (!is.null(colnames(X))) {
    offset <- FALSE
    if (any(grepl("offset", colnames(X)))) {
      beta.names <- colnames(X)[-grep("offset", colnames(X))]
      offset <- TRUE
    } else {
      beta.names <- colnames(X)
    }
    rownames(mu) <- beta.names
    rownames(Sigma) <- colnames(Sigma) <- beta.names
  } else {
    beta.names <- NULL
  }
    
  list(mu = mu, Sigma = Sigma, log.like = ll)
  
}


#' Helper function that checks for sensible covariate matrix.
#' 
#' @param X numeric matrix of covariates, \eqn{n} design points by \eqn{p}
#'   covariates, for a given model and design points.
#' @return throws an error if not full rank.
checkX <- function(X) {
  # check design matrix
  if (any(colSums(X != 0) == 0)) stop("each column must have at least one non-zero entry in X")
  qr.X <- qr(X)
  if (qr.X$rank < ncol(X)) stop("design matrix X columns are not linearly independent")
}

#' Helper function that gives the probability distribution function for design 
#' point.
#' 
#' @param x numeric: coordinate
#' @param Z list of design points and link function, see \code{\link{designLink}}
#' @param design.pt integer: design point
#' @param fit.method character: method for fit in \code{\link{mV}}, default is \code{KL}
#' @examples 
#' # design matrix: two scenarios
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) 
#' rownames(X) <- c("scenario1", "scenario2")
#' colnames(X) <- c("covariate1", "covariate2")
#' #' # logit link
#' # central credible intervals with probability = 1/2
#' Z <- designLink(design = X, link = "logit", CI.prob = 0.5)
#' #' # lower and upper quartiles and median
#' Z <- indirect::elicitPt(Z, design.pt = 1, 
#'   lower.CI.bound = 0.2,
#'   median = 0.4,
#'   upper.CI.bound = 0.6,
#'   comment = "Completed.")
#' indirect::plotDesignPoint(Z, design.pt = 1,   
#'   elicited.fractiles = TRUE, theta.bounds = c(0, 1),
#'   fitted.fractiles = TRUE, fitted.curve = TRUE)
#'   
#' # probability that target is below 0.1 and
#' # probability that target is below 0.9   
#' indirect::pdist(c(0.1, 0.9), Z, design.pt = 1)
pdist <- function(x, Z, design.pt = NULL, fit.method = "KL") {
  
      if (is.null(design.pt)) stop("specify design point")
    
      x <- sort(x)    
  
      # transform coordinate to linear predictor scale
      switch(Z$link, 
                       log = g.x <- log(x),
                       cloglog = g.x <- log(-log(1 - x)),
                       logit = g.x <- log(x) - log(1 - x),
                       identity = g.x <- x,
                       stop("specify link function")
                )
    
      f <- mV(Z, fit.method)
      vals <- pnorm(g.x, f$m[design.pt], sqrt(f$V[design.pt, design.pt]))
      names(vals) <- paste0("P(x <= ", x, ")")
      vals
      
}
