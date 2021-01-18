\name{classificationError}
\alias{classificationError}
\alias{classificationError.data.frame}
\alias{classificationError.default}

\title{Estimation of misclassification errors (generalisation errors)
  based on statistical and various machine learning methods}

\description{
Estimates misclassification errors (generalisation errors), sensitivity and specificity using cross-validation, 
bootstrap and \code{632plus} bias corrected bootstrap methods based on Random Forest, 
Support Vector Machines, Linear Discriminant Analysis and k-Nearest Neighbour methods.
}
\usage{
\method{classificationError}{data.frame}(
          formula,
          data, 
          method=c("RF","SVM","LDA","KNN"), 
          errorType = c("cv", "boot", "six32plus"),
	  senSpec=TRUE,
          negLevLowest=TRUE,
	  na.action=na.omit, 
          control=control.errorest(k=NROW(na.action(data)),nboot=100),
          ...)
}

\arguments{
  \item{formula}{A formula of the form \code{lhs ~ rhs} relating response (class) 
                 variable and the explanatory variables. See \code{\link[stats]{lm}} for 
                 more detail.}

  \item{data}{A data frame containing the response (class membership) variable and 
              the explanatory variables in the formula.}

  \item{method}{A character vector of length \code{1} to \code{4} representing the classification 
                methods to be used. Can be one or more of \code{"RF"} (Random Forest), \code{"SVM"} 
                (Support Vector Machines), \code{"LDA"} (Linear Discriminant Analysis) and \code{"KNN"} 
                (k-Nearest Neighbour). Defaults to all four methods.}

  
  \item{errorType}{A character vector of length \code{1} to \code{3} representing the type of 
                  estimators  to be used for computing misclassification errors. 
                  Can be one or more of the \code{"cv"} (cross-validation), \code{"boot"} 
                (bootstrap) and \code{"632plus"} (632plus bias corrected bootstrap) estimators.
                 Defaults to all three estimators.}

  \item{senSpec}{Logical. Should sensitivity and specificity (for cross-validation estimator only)
                 be computed? Defaults to \code{TRUE}.}
  \item{negLevLowest}{Logical. Is the lowest of the ordered levels of the class variable represnts
                     the negative control? Defaults to \code{TRUE}.}

   \item{na.action}{Function which indicates what should happen when the data
                   contains \code{NA}'s, defaults to \code{\link{na.omit}}.} 
  \item{control}{Control parameters of the the function \code{\link[ipred]{errorest}}.}
  \item{\dots}{additional parameters to \code{method}.}
}
\details{
 In the current version of the package, estimation of sensitivity and
 specificity is limited to cross-validation estimator only. For LDA
 sample size must be greater than the number of explanatory variables to
 avoid singularity. The function \code{classificationError} does not
 check if this is satisfied, but the underlying function
 \code{\link[MASS]{lda}} produces warnings if this condition is violated.
}
\value{
Returns an object of class \code{classificationError} with components
  \item{call}{The call of the \code{classificationError} function.}

  \item{errorRate}{A \code{length(errorType)} by \code{length(method)}
  matrix  of classification errors.}
  
  \item{rocData}{A \code{2} by \code{length(method)}  matrix  of
  sensitivities (first row) and specificities (second row).}
 
}

\references{

Khondoker, M. R., Till T. Bachmann, T. T., Mewissen, M., Dickinson, P. \emph{et al.}(2010).
Multi-factorial analysis of class prediction error: estimating optimal number of biomarkers for various classification rules.
\emph{Journal of Bioinformatics and Computational Biology}, \bold{8}, 945-965.


Breiman, L. (2001). \emph{Random Forests}, Machine Learning
\bold{45}(1), 5--32.

Chang, Chih-Chung and Lin, Chih-Jen: \emph{LIBSVM: a library for Support
  Vector Machines}, \url{https://www.csie.ntu.edu.tw/~cjlin/libsvm/}. 

Ripley, B. D. (1996). \emph{Pattern Recognition and Neural
  Networks}.Cambridge: Cambridge University Press.

Efron, B.  and Tibshirani, R. (1997). Improvements on Cross-Validation:
The .632+ Bootstrap Estimator. \emph{Journal of the American Statistical
  Association} \bold{92}(438), 548--560.

}
\author{
Mizanur Khondoker, Till Bachmann, Peter Ghazal \cr
 Maintainer: Mizanur Khondoker \email{mizanur.khondoker@gmail.com}.}

\seealso{\code{\link{simData}}}

\examples{

\dontrun{
mydata<-simData(nTrain=30,nBiom=3)$data
classificationError(formula=class~., data=mydata)
}

}

\keyword{classif}


