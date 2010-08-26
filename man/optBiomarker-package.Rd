\name{optBiomarker-package}
\Rdversion{1.1}
\alias{optBiomarker-package}
\alias{optBiomarker}
\docType{package}
\title{
R package for estimating optimal number of biomarkers at a given error tolerance level
for various classification rules
}
\description{
 Using interactive control panel  (\code{\link[rpanel:rpanel.package]{rpanel}})  and 3D
 real-time rendering system  (\code{\link[rgl:rgl-package]{rgl}}) , this package provides a
 user friendly GUI for estimating the minimum  number of biomarkers
 (variables) needed to achieve  a given level of accuracy for two-group
 classification problems based on microarray data. 
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

  Maintainer: Mizanur Khondoker \email{mizanur.khondoker@gmail.com}.
}
\references{
Khondoker, M. R., Till T. Bachmann, T. T., Mewissen, M., Dickinson, P. \emph{et al.}(2010).
Multi-factorial analysis of class prediction error: estimating optimal number of biomarkers for various classification rules.
\emph{Journal of Bioinformatics and Computational Biology}, \bold{8}, 945-965.

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
  \code{\link{optimiseBiomarker}}
}
\examples{
if(interactive()){
data(errorDbase)
optimiseBiomarker(error=errorDbase)
}
}

\keyword{package}
