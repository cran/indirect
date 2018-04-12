#' Create list with information for the elicitation session
#' 
#' This builds the structure that will store elicited data. The linear predictor
#' has a normal prior \eqn{g(\theta) ~ N(m, V)}, \eqn{\theta} is the elicitation
#' target. Link functions \eqn{g(.)}: \code{logit}, \code{log}, \code{cloglog}, 
#' \code{identity}.
#' 
#' Assumption: at least two fractiles selected from the median, upper and lower 
#' bounds of hte central credible interval of probability \code{CI.prob} will be
#' elicited at each design point. The probabilities assigned to the central 
#' credible intervals  can vary across design points. The argument 
#' \code{CI.prob} can later be adjusted by design point during the elicitation 
#' exercise, see function \code{\link{elicitPt}}. In the first instance, it is
#' set to a global value specified by \code{CI.prob} in function
#' \code{\link{designLink}} with default value \eqn{0.5}.
#' @param design a dataframe with covariate values that will be displayed to the
#'   expert(s) during the elicitation session.
#' @param link character \code{logit}, \code{log}, \code{cloglog}, 
#'   \code{identity}
#' @param target character, name of target parameter of elicitation exercise
#' @param CI.prob numeric, a fraction between 0 and 1 that defines probability 
#'   attributed to central credible interval. For example, 1/2 for a central 
#'   credible interval of probability 0.5, or 1/3 for a central credible 
#'   interval of probablity 0.333... The default is probability 1/2.
#' @param expertID character, identifier for expert or group of experts
#' @param facilitator character, facilitator identifier
#' @param rapporteur character, rapporteur identifier. Default "none".
#' @param intro.comments character, text with any prefacing comments. This may 
#'   include, for example, the definition of the target parameter for the 
#'   elictation session. Beware of non-ASCII text and special characters, which 
#'   may affect the ability to save the elicitation record with function \code{\link{saveRecord}}
#'   or create a summary report with function \code{\link{makeSweave}}
#'   if called by the function \code{\link{makeSweave}} may affect ability to render  by
#'   means of \code{\link[utils]{Sweave}} or \code{knitr} etc.
#' @param fit.method character, method used to fit conditional means prior: 
#'   \code{KL} (default), \code{moment}, \code{SS} (see vignette and
#'   \code{\link{mV}} for more information on these options) 
#' @return list of \code{design} with entries: \code{theta}, a \eqn{n x 4} 
#'   matrix with columns that give lower, median and upper quantiles followed by
#'   \code{CI.prob} and \eqn{n} equal to the number of design points 
#'   (scenarios); \code{link}, the link function used; \code{target}; 
#'   \code{expert} \code{facilitator}; \code{rapporteur}; \code{date}; 
#'   \code{intro.comments}; \code{fit.method}.
#' @examples
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) # design
#' Z <- designLink(design = X, link = "logit", target = "target",
#'  CI.prob = 1/2, expertID = "Expert", facilitator = "facilitator")
designLink <- function(design, link = "identity", target = "Target", CI.prob = 1/2,
                       expertID = "Expert", facilitator = "Facilitator",
                       rapporteur = "none",  
                       intro.comments = "This is a record of the elicitation session.",
                       fit.method = "KL") {
  
  # append columns for elicited responses 
  theta <- matrix(NA, nrow = nrow(design), ncol = 4,
    dimnames = list(NULL, c("lower", "median", "upper", "CI_prob")))
  if (!is.null(CI.prob)) {
    if (length(CI.prob) != 1 || CI.prob < 0 || CI.prob > 1) stop("CI.prob is a scalar probability")
    theta[ , "CI_prob"] <- CI.prob
  }
  if (!is.null(link)) {
    # check that link function is supported
    if (!any(link%in%c("log", "cloglog", "logit", "identity"))) stop("unsupported link function")    
  }

  Z <- list(design = design, theta = theta, link = link, target = target,
            expertID = expertID, facilitator = facilitator, 
            rapporteur = rapporteur, intro.comments = intro.comments,
            comments = rep(" ", nrow(design)), fit.method = fit.method)
  
}

#' Function to create or update elicitation at a given design point.
#' 
#' @param Z list of \code{design} with entries: \code{theta}, a \eqn{n x 4} 
#'   matrix with columns that give lower, median and upper quantiles of the 
#'   central credible interval followed by the probability \code{CI.prob} 
#'   allocated to the interval; \code{link}, the link function used; and 
#'   \code{target}. This list object is created by \code{\link{designLink}}
#' @param design.pt single integer that denotes design point of interest
#' @param lower.CI.bound scalar that gives the lower bound of the central 
#'   credible interval, default \code{NA}.
#' @param median scalar value, default \code{NA}
#' @param upper.CI.bound scalar that gives the upper bound of the central 
#'   credible interval, default \code{NA}.
#' @param CI.prob numeric, a fraction between 0 and 1 that defines probability 
#'   attributed to central credible interval. For example, 1/2 for quartiles or
#'   1/3 for tertiles. Default \code{NULL} uses the initial \code{CI.prob} as
#'   defined by \code{\link{designLink}}.
#' @param comment character, ASCII text providing contributed commentary associated 
#' with elicitation design point. It is recommended to avoid special characters
#' such as quotation marks etc.
#' @return \code{Z}, a list of \code{design} with entries: \code{theta}, a 
#'   \eqn{n x 4} matrix with columns that give lower, median and upper quantiles
#'   followed by \code{CI.prob}  with updated entries for row specified by
#'   argument \code{design.pt}; \code{link}, the link function used; and
#'   \code{target}.
#' @examples 
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) # design
#' Z <- designLink(design = X)
#' Z <- elicitPt(Z, design.pt = 1,
#'   lower.CI.bound = -1,
#'   median = 0,
#'   upper.CI.bound = 1,
#'   comment = "A completed elicitation scenario.")
elicitPt <- function(Z, design.pt = NULL, 
                     lower.CI.bound = NA, median = NA, upper.CI.bound = NA, 
                     CI.prob = NULL, comment = " ") {
  
  if (is.null(design.pt)) stop("specify design point")
  
  CI.bounds <- c(lower.CI.bound, upper.CI.bound)
  
  if (all(!is.na(CI.bounds))) {
    if (CI.bounds[1] > CI.bounds[2]) {
      CI.bounds <- sort(CI.bounds)
      warning("resorted CI.bounds into ascending order")
    } else if (length(unique(CI.bounds)) == 1) {
      stop("lower and upper bounds of central credible interval are identical")
    }
  }
  
  if (!is.na(median) && !is.na(CI.bounds[1]) && median <= CI.bounds[1]) {
    stop("lower bound of central credible interval must be less than median")
  }
  if (!is.na(median) && !is.na(CI.bounds[2]) && median >= CI.bounds[2]) {
    stop("upper bound of central credible interval must be greater than median")
  }
  
  Z$theta[design.pt, c(1, 3)] <- CI.bounds
  Z$theta[design.pt, 2] <- median
  
  if (!is.null(CI.prob)) {
    Z$theta[design.pt, 4] <- CI.prob
  } 
  if (is.na(Z$theta[design.pt, 4]) || Z$theta[design.pt, 4] <= 0 || Z$theta[design.pt, 4] >= 1) {
    stop("CI.prob must be numeric scalar between (0, 1)")
  }
  
  if (comment != " ") {
    
    if (Z$comments[design.pt] != " ") warning("writing over previous comment")
    ascii.conv <- iconv(comment, "latin1", "ASCII")
    if (any(is.na(ascii.conv)) || any(ascii.conv != comment)) warning("non-ASCII text detected in commment")
    specials <- strsplit("#$%^&\"'", split = "")[[1]]
    lspec <- logical(length(specials))
    for (i in 1:length(lspec)) {
      lspec[i] <- grepl(specials[i], comment, fixed = TRUE)
    }
    if (any(lspec)) warning("at least one special character was detected in comment")
    
    Z$comments[design.pt] <- comment
    
  }

  Z
  
}

#' Function to check condition number diagnostic.
#' 
#' This function calculates the condition number of  the rescaled \eqn{n x
#' p} design matrix \eqn{X} such that each column has unit length.
#' 
#' @param X Design matrix
#' @return a scalar giving the condition number of the rescaled design matrix
#' @examples 
#' X <- matrix(rnorm(16), nrow = 4)
#' CNdiag(X)
CNdiag <- function(X) {
  indirect::checkX(X)
  Xs <- apply(X, 2, function(x) x/sqrt(sum(x^2)))
  kappa(Xs)
}

#' Function to save elicitation record.
#' 
#' @param designLink.obj list object initally created by function \code{\link{designLink}} 
#'   and subsequently updated by function \code{\link{elicitPt}}
#' @param conclusion.comments character, comments to conclude session. Beware of
#'   non-ASCII text and special characters, which may affect ability to save or
#'   generate a \code{Sweave} document by using \code{\link{makeSweave}} 
#' @param file character providing filename. 
#' @return an RDS file is created with filename \code{file}. A timestamp is
#'   added to \code{designLink.obj} using \code{Sys.time()}.
#' @examples
#' \dontrun{
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) # design
#' Z <- designLink(design = X)
#' tmp <- tempfile(pattern = "report", fileext =".rds")
#' saveRecord(Z, file = tmp)
#' }
saveRecord <- function(designLink.obj, conclusion.comments = "This concludes the elicitation record.",
                       file = "") {
  designLink.obj$conclusion.comments <- conclusion.comments
  designLink.obj$timestamp <- Sys.time()
  saveRDS(designLink.obj, file = file)
}

#' Function to create summary document from a saved elicitation record.
#' 
#' Creates a Sweave file that can be used to generate a pdf document of the 
#' summary report.
#' 
#' @param filename.rds character, filename of the record saved as an RDS object,
#'   see \code{?saveRDS}.
#' @param reportname character, filename without extension to be used for the 
#'   generated Sweave (\code{.Rnw}) file. The Sweave file supports the creation
#'   of report (\code{.pdf}) documentation and accompanying files such as the
#'   \code{.tex} file generated by using \code{\link[utils]{Sweave}} followed by
#'   \code{tools::texi2pdf()}.
#' @param title character, a title for the report
#' @param contact.details character, an email address or other mechanism by 
#'   which the expert may contact the facilitator or rapporteur
#' @param fitted.fractiles logical or numeric vector. A logical value of
#'   \code{FALSE} will not plot any fitted fractiles from the fitted subjective
#'   probability distribution. A logical value of \code{TRUE} will plot the
#'   fitted fractiles that correspond to the final iteration of the raw elicited
#'   fractiles. Alternatively, a numeric vector can specify arbitrary fractiles
#'   for plotting from the fitted distribution, e.g., \code{c(1/10, 1/4, 1/2,
#'   3/4, 9/10)}
#' @param cumul.prob.bounds numeric vector that specifies the upper and lower 
#'   plot bounds determined by this credible interval. The default is the 0.90 
#'   central credible interval, \code{c(0.05, 0.95)}
#' @examples 
#' \dontrun{ 
#' X <- matrix(c(1, 1, 0, 1), nrow = 2) # design
#' Z <- designLink(design = X)
#' Z <- elicitPt(Z, design.pt = 1,
#'   lower.CI.bound = -1,
#'   median = 0,
#'   upper.CI.bound = 1,
#'   comment = "A completed elicitation scenario.")
#' tmp.rds <- tempfile(pattern = "record", fileext =".rds")
#' saveRecord(Z, file = tmp.rds)
#' tmpReport <- tempfile(pattern = "report")
#' makeSweave(filename.rds = tmp.rds, reportname = tmpReport)
#' setwd(tempdir())
#' utils::Sweave(paste0(tmpReport, ".Rnw"))
#' tools::texi2pdf(paste0(tmpReport, ".tex")) 
#' }
makeSweave <- function(filename.rds = "", reportname = "",
                       title = "Elicitation record", contact.details = "none",
                       fitted.fractiles = TRUE, cumul.prob.bounds = c(0.05, 0.95)) {

  if (reportname == "") stop("'reportname' must be non-empty string")
  if (grepl("\\\\", filename.rds)) filename.rds <- gsub("\\\\", "/", filename.rds)
  if (grepl("\\\\", reportname)) gsub("\\\\", "/", reportname)  
  record <- readRDS(filename.rds)
  dpts <- 1:nrow(record$theta)
  printpdfs <- paste0(
    "\\includegraphics[width=\\textwidth]{designPt_", dpts, 
    ".pdf}\\\\ \\Sexpr{record$comments[", dpts, "]}\\\\ "
  )
cat("\\documentclass{article}
\\title{", title, "}
\\date{\\today}
\\begin{document}
\\maketitle
\\begin{description}
\\item[Session Completed:]",  as.character(strftime(record$timestamp, format = "%e %B %Y %H:%M"), usetz = TRUE), "
\\item[Target:]", record$target, "
\\item[Expert:]", record$expertID, "
\\item[Facilitator:]", record$facilitator, "
\\item[Rapporteur:]", record$rapporteur, "
\\item[Contact details for follow up:]", contact.details, "
\\end{description}
\\section{Introduction}
", record$intro.comments, "
<<load, echo = FALSE>>=
record <- readRDS(\"", filename.rds, "\")
@
\\section{Design}
All design points are printed below.
<<design, echo = FALSE>>=
record$design
@
\\section{Elicitations}
The fitted distribution represented by the density curve and fitted fractiles are the final elicited subjective 
probability distributions accepted into the elicitation record. The elicited subjective probability distribution 
were the best fits to the raw elicited fractiles. The raw fractiles are also reported here for each design point.
\\\\
<<plots,echo=FALSE>>=
for (i in 1:nrow(record$design)) {

  design.pt <- i

  if (sum(!is.na(record$theta[design.pt, 1:3])) < 2) {
    fpctls <- FALSE
    fcrv <- FALSE
    if (record$link == \"logit\" || record$link == \"cloglog\") tbnds <- c(0, 1)
    if (record$link == \"log\") tbnds <- c(0, 10^6)
    if (record$link == \"identity\") tnbds <- c(-10^6, 10^6)
  } else {
    if (is.numeric(c(", paste(fitted.fractiles, collapse = ","), "))) { 
      fpctls <- c(", paste(fitted.fractiles, collapse = ","), ")
    } else {
      fpctls <- TRUE
    }
    fcrv <- TRUE
  }

  fname <- paste0(\"designPt_\", i)
  pdf(file = paste0(fname, \".pdf\"), width = 7, height = 7)
  par(oma = rep(1, 4))
  indirect::plotDesignPoint(record, design.pt = design.pt, theta.bounds = NULL, xlog = FALSE,
                          elicited.fractiles = TRUE, fitted.fractiles = fpctls,
                          fitted.curve = fcrv, n.pts = 501,
                          cumul.prob.bounds = c(", paste(cumul.prob.bounds, collapse = ","), ")
  )
  dev.off()
}
@
", printpdfs, "
\\section{Conclusion}
", record$conclusion.comments, "
\\end{document}
", 
file = paste0(reportname, ".Rnw"), sep = ""
)
}

