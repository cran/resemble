% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbl_control.R
\name{mbl_control}
\alias{mbl_control}
\title{A function that controls some few aspects of the memory-based learning
process in the \code{mbl} function}
\usage{
mbl_control(return_dissimilarity = FALSE,
            validation_type = c("NNv", "local_cv"),
            tune_locally = TRUE,
            number = 10,
            p = 0.75,
            range_prediction_limits = TRUE,
            allow_parallel = TRUE)
}
\arguments{
\item{return_dissimilarity}{a logical indicating if the dissimilarity matrix
between \code{Xr} and \code{Xu} must be returned.}

\item{validation_type}{a character vector which indicates the (internal) validation
method(s) to be used for assessing the global performance of the local models.
Possible options are: \code{"NNv"} and \code{"local_cv"}. Alternatively
\code{"none"} can be used when cross-validation is not required (see details
below).}

\item{tune_locally}{a logical. It only applies when
\code{validation_type = "local_cv"} and "pls" or "wapls" fitting algorithms are
used.  If \code{TRUE}, then the parameters of the local pls-based models
(i.e. pls factors for the "pls" method and minimum and maximum pls factors
for the "wapls" method) are optimized. Default is \code{TRUE}.}

\item{number}{an integer indicating the number of sampling iterations at
each local segment when \code{"local_cv"} is selected in the
\code{validation_type} argument. Default is 10.}

\item{p}{a numeric value indicating the percentage of calibration observations
to be retained at each sampling iteration at each local segment when \code{"local_cv"}
is selected in the \code{validation_type} argument. Default is 0.75 (i.e. 75 "\%").}

\item{range_prediction_limits}{a logical. It indicates whether the prediction
limits at each local regression are determined by the range of the response
variable within each neighborhood. When the predicted value is outside
this range, it will be automatically replaced with the value of the nearest
range value. If \code{FALSE}, no prediction limits are imposed.
Default is \code{TRUE}.}

\item{allow_parallel}{a logical indicating if parallel execution is allowed.
If \code{TRUE}, this parallelism is applied to the loop in \code{\link{mbl}}
in which each iteration takes care of a single observation in \code{Xu}. The
parallelization of this for loop is implemented using the
code{\link[foreach]{foreach}} function of the package \code{foreach}.
Default is \code{TRUE}.}
}
\value{
a \code{list} mirroring the specified parameters
}
\description{
\loadmathjax
This function is used to further control some aspects of the memory-based
learning process in the \code{mbl} function.
}
\details{
The validation methods available for assessing the predictive performance of
the memory-based learning method used are described as follows:
 \itemize{
 \item{Leave-nearest-neighbor-out cross-validation (\code{"NNv"}):}{ From
 the group of neighbors of each observation to be predicted, the nearest observation
 (i.e. the most similar observation) is excluded and then a local model is fitted
 using the remaining neighbors. This model is then used to predict the value
 of the response variable of the nearest observation. These predicted
 values are finally cross validated with the actual values (See Ramirez-Lopez
 et al. (2013a) for additional details). This method is faster than
 \code{"local_cv"}.}
 \item{Local leave-group-out cross-validation (\code{"local_cv"}):}{ The
 group of neighbors of each observation to be predicted is partitioned into
 different equal size subsets. Each partition is selected based on a
 stratified random sampling that uses the the distribution of
 the response variable in the corresponding set of neighbors. When
 \code{p} \mjeqn{>=}{\geqslant} 0.5 (i.e. the number of calibration
 observations to retain is larger than 50% of the total samples in the neighborhood),
 the sampling is conducted for selecting the validation samples, and when
 \code{p} < 0.5 the sampling is conducted for selecting the calibration
 samples (samples used for model fitting). The model fitted with the selected
 calibration samples is used to predict the response values of the local
 validation samples and the local root mean square error is computed.
 This process is repeated \mjeqn{m}{m} times and the final local
 error is computed as the average of the local root mean square errors
 obtained for all the \mjeqn{m}{m} iterations. In the \code{mbl_control} function
 \mjeqn{m}{m} is controlled by the \code{number} argument and the size of the
 subsets is controlled by the \code{p} argument which indicates the
 percentage of observations to be selected from the subset of nearest neighbors.
 The global error of the predictions is computed as the average of the local
 root mean square errors.}
 \item{No validation (\code{"none"}):}{ No validation is carried out.
 If \code{"none"} is seleceted along with \code{"NNv"} and/or
 \code{"local_cv"}, then it will be ignored and the respective
 validation(s) will be carried out.}
 }
}
\examples{
# A control list with the default parameters
mbl_control()
}
\references{
Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M.,
Scholten, T. 2013a. The spectrum-based learner: A new local approach for
modeling soil vis-NIR spectra of complex data sets. Geoderma 195-196, 268-279.

Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte,
J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for
use with soil vis-NIR spectra. Geoderma 199, 43-53.
}
\seealso{
\code{\link{f_diss}}, \code{\link{cor_diss}}, \code{\link{sid}},
\code{\link{ortho_diss}}, \code{\link{mbl}}
}
\author{
\href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez} and Antoine Stevens
}
