\name{optimiseBiomarker}
\alias{optimiseBiomarker}
\docType{package}
\title{
Estimates optimal number of biomarkers at a given error tolerance level
for various classification rules
}
\description{
 Using interactive control panel (see \code{\link[rpanel:rpanel.package]{rpanel}}) and 3D
 real-time rendering system (\code{\link[rgl:rgl-package]{rgl}}), this package provides a
 user friendly GUI for estimating the minimum  number of biomarkers
 (variables) needed to achieve  a given level of accuracy for two-group
 classification problems  based on microarray data. 
}

\usage{
optimiseBiomarker (error,
                   errorTol = 0.05,
                   method = "RF", nTrain = 100, 
                   sdB = 1.5,
                   sdW = 1,
                   foldAvg = 2.88,
                   nRep = 3) 

}

\arguments{
  \item{error}{The database of classification errors. See
    \code{errorDbase} for details.}
  \item{errorTol}{Error tolerance limit.}
  \item{method}{Classification method. Can be one of \code{"RF"},
  \code{"SVM"}, and \code{"KNN"} for Random Forest,
  Support Vector Machines, Linear Discriminant Analysis and k-Nearest
  Neighbour respectively.}
   \item{nTrain}{Training set size, i.e., the total number of biological
    samples in group 1 and group 2.}
  \item{sdB}{Biological variation (\eqn{\sigma_b}) of  data in
  log (base 2) scale.}
  \item{sdW}{Experimental (technical) variation (\eqn{\sigma_e}) of  data in
    log (base 2) scale.}
  \item{foldAvg}{Average fold change of the biomarkers.}
    \item{nRep}{Number of technical replications.}
 
%   \item{rho}{Common Pearson correlation between biomarkers. To ensure  positive definiteness of the correlation matrix, allowed values of \code{rho} are restricted between 0 and 0.95 inclusive.}
   }



\details{
The function \code{\link{optimiseBiomarker}} is a user friendly GUI for
interrogating the database of leave-one-out cross-validation errors,
\code{\link{errorDbase}}, to estimate optimal number of biomarkers for
microarray based classifications. The database is built on the basis of
simulated data using the \code{\link{classificationError}} function. The
function \code{\link{simData}} is used for simulating microarray data
for various combinations of factors such as the number of biomarkers,
training set size, biological variation, experimental variation, fold
change, replication, and correlation.
  
}
\author{
  Mizanur Khondoker, Till Bachmann, Peter Ghazal \cr

  Maintainer: Mizanur Khondoker \email{mizanur.khondoker@googlemail.com}.
}
\references{

Breiman, L. (2001). \emph{Random Forests}, Machine Learning
\bold{45}(1), 5--32.

Chang, Chih-Chung and Lin, Chih-Jen: \emph{LIBSVM: a library for Support
  Vector Machines}, \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}. 

Ripley, B. D. (1996). \emph{Pattern Recognition and Neural
  Networks}.Cambridge: Cambridge University Press.

Efron, B.  and Tibshirani, R. (1997). Improvements on Cross-Validation:
The .632+ Bootstrap Estimator. \emph{Journal of the American Statistical
  Association} \bold{92}(438), 548--560.

Bowman, A., Crawford, E., Alexander, G. and Bowman,
R. W. (2007). rpanel: Simple interactive controls for R functions using
the tcltk package. \emph{Journal of Statistical Software} \bold{17}(9). 
  
}


\seealso{
  \code{\link{simData}}
  \code{\link{classificationError}}

}
\examples{
if(interactive()){
data(errorDbase)
optimiseBiomarker(error=errorDbase)
}
}


\keyword{optimize}
