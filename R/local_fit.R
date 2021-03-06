#' @title Local fit functions
#' @name local_fit
#' @aliases local_fit
#' @aliases local_fit_pls
#' @aliases local_fit_wapls
#' @aliases local_fit_gpr
#' @description
#' \loadmathjax
#' These functions define the way in which each local fit/prediction is done
#' within each iteration in the \code{\link{mbl}} function.
#' @usage
#' local_fit_pls(pls_c)
#'
#' local_fit_wapls(min_pls_c, max_pls_c)
#'
#' local_fit_gpr(noise_variance = 0.001)
#'
#' @param pls_c an integer indicating the number of pls components to be used in
#' the local regressions when the partial least squares (\code{local_fit_pls})
#' method is used.
#' @param min_pls_c an integer indicating the minimum number of pls components
#' to be used in the local regressions when the weighted average partial least
#' squares (\code{local_fit_wapls}) method is used. See details.
#' @param max_pls_c integer indicating the maximum number of pls components
#' to be used in the local regressions when the weighted average partial least
#' squares (\code{local_fit_wapls}) method is used. See details.
#' @param noise_variance a numeric value indicating the variance of the noise
#' for Gaussian process local regressions (\code{local_fit_gpr}). Default is
#' 0.001.
#'
#' @details
#' These functions are used to indicate how to fit
#' the regression models within the \code{\link{mbl}} function.
#'
#' There are three possible options for performing these regressions:
#' \itemize{
#'  \item{Partial least squares (pls, \code{local_fit_pls}):}{ It uses the
#'  orthogonal scores (non-linear iterative partial least squares, nipals)
#'  algorithm. The only parameter which needs to be optimized is the number of
#'  pls components.}
#'  \item{Weighted average pls (\code{local_fit_wapls}):}{ This method was
#'  developed by Shenk et al. (1997) and it used as the regression method in the
#'  widely known LOCAL algorithm. It uses multiple models generated by multiple
#'  pls components (i.e. between a minimum and a maximum number of pls
#'  components). At each local partition the final predicted value is a ensemble
#'  (weighted average) of all the predicted values generated by the multiple pls
#'  models. The weight for each component is calculated as follows:
#'
#'  \mjdeqn{w_{j}  =  \frac{1}{s_{1:j}\times g_{j}}}{w_j  =  1/(s_{1:j} xx g_{j})}
#'
#'  where \mjeqn{s_{1:j}}{s_{1:j}} is the root mean square of the spectral residuals of the
#'  unknown (or target) obasevation(s) when a total of \mjeqn{j}{j} pls components are
#'  used and \mjeqn{g_{j}}{g_{j}} is the root mean square of the regression coefficients
#'  corresponding to the \mjeqn{j}{j}th pls component (see Shenk et al., 1997 for
#'  more details).}
#'  \item{Gaussian process with dot product covariance (\code{local_fit_gpr):}{
#'  Gaussian process regression is a probabilistic and non-parametric Bayesian
#'  method. It is commonly described as a collection of random variables which
#'  have a joint Gaussian distribution and it is characterized by both a mean
#'  and a covariance function (Rasmussen and Williams, 2006). The covariance
#'  function used in the implemented method is the dot product. The only
#'  parameter to be taken into account in this method is the noise. In this
#'  method, the process for predicting the response variable of a new sample
#'  (\mjeqn{y_{u}}{y_{u}}) from its predictor variables
#'  (\mjeqn{x_{u}}{x_{u}}) is carried out first by computing a prediction
#'  vector (\mjeqn{A}{A}). It is derived from a reference/training observations
#'  congaing both a response vector (\mjeqn{Y}{Y}) and predictors (\mjeqn{X}{X}) as follows:
#'
#'  \mjdeqn{A = (X  X^{T} + \sigma^2 I)^{-1} Y}{A = (X  X^T + sigma^2 I)^{-1} Y}
#'
#'   where  \mjeqn{\sigma^2}{sigma^2} denotes the variance of the noise and \mjeqn{I}{I} the
#'   identity matrix (with dimensions equal to the number of observations in
#'   \mjeqn{X}{X}). The prediction of \mjeqn{y_{u}}{y_u} is then done as follows:
#'
#'   \mjdeqn{\hat{y}_{u} = (x_{u}x_{u}^{T}) A}{hat y_{u} = (x_{u} x_{u}^T) A}
#'
#'  }
#'  }
#'  }
#' @return An object of class \code{local_fit} mirroring the input arguments.
#' @author \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @references
#' Shenk, J., Westerhaus, M., and Berzaghi, P. 1997. Investigation of a LOCAL
#' calibration procedure for near infrared instruments. Journal of Near Infrared
#' Spectroscopy, 5, 223-232.
#'
#' Rasmussen, C.E., Williams, C.K. Gaussian Processes for Machine Learning.
#' Massachusetts Institute of Technology: MIT-Press, 2006.
#'
#' @seealso \code{\link{mbl}}
#' @examples
#' local_fit_wapls(min_pls_c = 3, max_pls_c = 12)
#' @export

## History:
## 28.05.2020 Leo     Hello world!
## 19.07.2020 Leo     arguments pls_max_iter and pls_tol were removed fas they
##                    are only required when modleing for more than one response
##                    variable, i.e. pls2 (which is not implemented for mbl)


local_fit_pls <- function(pls_c) {
  if (missing(pls_c)) {
    stop("'pls_c' must be specified")
  }

  if (length(pls_c) != 1 | !is.numeric(pls_c)) {
    stop(paste0(
      "'pls_c' must be a single numerical ",
      "value specifiying the maximum number of pls components to ",
      "be evaluated"
    ))
  }

  fit_type <- list(
    method = "pls",
    pls_c = pls_c
  )
  class(fit_type) <- c("local_fit", "list")
  fit_type
}

#' @aliases local_fit
#' @export local_fit_wapls
local_fit_wapls <- function(min_pls_c, max_pls_c) {
  if (missing(min_pls_c) | missing(max_pls_c)) {
    stop("Both 'min_pls_c' and 'max_pls_c' must be specified")
  }

  if (length(min_pls_c) != 1 | !is.numeric(min_pls_c)) {
    stop(paste0(
      "'min_pls_c' must be a single numerical ",
      "value specifiying the minimum number of pls components to ",
      "be evaluated"
    ))
  }

  if (length(max_pls_c) != 1 | !is.numeric(max_pls_c)) {
    stop(paste0(
      "'max_pls_c' must be a single numerical ",
      "value specifiying the maximum number of pls components to ",
      "be evaluated"
    ))
  }

  if (min_pls_c >= max_pls_c) {
    stop("min_pls_c must be smaller than max_pls_c")
  }

  fit_type <- list(
    method = "wapls",
    pls_c = c(min_pls_c = min_pls_c, max_pls_c = max_pls_c)
  )

  class(fit_type) <- c("local_fit", "list")
  fit_type
}

#' @aliases local_fit
#' @export local_fit_gpr
local_fit_gpr <- function(noise_variance = 0.001) {
  if (length(noise_variance) != 1 | !is.numeric(noise_variance)) {
    stop("'noise_variance' must be a single numerical value")
  }
  fit_type <- list(
    method = "gpr",
    noise_variance = noise_variance
  )
  class(fit_type) <- c("local_fit", "list")
  fit_type
}
