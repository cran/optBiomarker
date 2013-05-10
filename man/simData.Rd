\name{simData}
\alias{simData}

\title{Simulation of microarray data}

\description{
The function simulates microarray data for two-group comparison with
user supplied parameters such as  number of biomarkers (genes or proteins),
sample size, biological and experimental (technical) variation,
replication, differential expression, and correlation between biomarkers. 
}
\usage{
simData(nTrain=100,
        nGr1=floor(nTrain/2),
        nBiom=50,nRep=3,
        sdW=1.0,
        sdB=1.0,
        rhoMax=NULL, rhoMin=NULL, nBlock=NULL,bsMin=3, bSizes=NULL, gamma=NULL,
        sigma=0.1,diffExpr=TRUE,
        foldMin=2,
        orderBiom=TRUE,
        baseExpr=NULL)
}

\arguments{

  \item{nTrain}{Training set size,.i.e., the total number of biological
    samples in group 1 (\code{nGr1}) and group 2.}
  
  \item{nGr1}{Size of group 1. Defaults to \code{floor(nTrain/2)}.}
  \item{nBiom}{Number of biomarkers (genes, probes or proteins).}
  \item{nRep}{Number of technical replications.}
  \item{sdW}{Experimental (technical) variation (\eqn{\sigma_e}) of  data in
  log (base 2) scale.}
  \item{sdB}{Biological variation (\eqn{\sigma_b}) of  data in
    log (base 2) scale.}
  \item{rhoMax}{Maximum Pearson's correlation coefficient between
  biomarkers. To ensure positive definiteness, allowed values are
  restricted between 0 and 0.95 inclusive.  If \code{NULL}, set to
  \code{runif(1,min=0.6,max=0.8)}.}

 \item{rhoMin}{Minimum Pearson's correlation coefficient between
  biomarkers. To ensure positive definiteness, allowed values are
  restricted between 0 and 0.95 inclusive. If \code{NULL}, set to
  \code{runif(1,min=0.2,max=0.4)}.}

\item{nBlock}{Number of blocks in the block diagonal (Hub-Toeplitz)
  correlation matrix.  If \code{NULL}, set to 1 for \code{nBiom<5} and
randomly selected from \code{c(1:floor(nBiom/bsMin))} for  \code{nBiom>=5}.}

\item{bsMin}{Minimum block size. \code{bsMin=3} by default.}
\item{bSizes}{A vector of length \code{nBlock} representing the block sizes
  (should sum to \code{nBlock}). If \code{NULL}, set to
  \code{c(bs+mod,rep(bs,nBlock-1)}, where \code{bs} is the integer
  part of \code{nBiom/nBlock} and \code{mod} is the remainder after
  integer division.}

\item{gamma}{Specifies a correlation structure. If \code{NULL}, assumes
  independence.\code{gamma=0} indicates a single block exchangeable
  correlation marix with constant correlation
  \code{rho=0.5*(rhoMin+rhoMax)}. A value greater than zero indicates
  block diagonal (Hub-Toeplitz) correlation matrix with decline rate
  determined by the value of \code{gamma}. Decline rate is linear for
  \code{gamma=1}.}
  
\item{sigma}{Standard deviation of the normal distribution (before
  truncation) where fold changes are generated from. See details.}
  \item{diffExpr}{Logical. Should systematic difference be introduced
  between the data of the two groups?}
 \item{foldMin}{Minimum value of fold changes. See details.}
 \item{orderBiom}{Logical. Should columns (biomarkers) be
  arranged in order of differential expression?}
  \item{baseExpr}{A vector of length \code{nBiom} to be used as base
  expressions \eqn{\mu}. See \code{realBiomarker} for details.}
}

\details{
Differential expressions are introduced by adding \eqn{z\delta} to the data
of group 2 where \eqn{\delta} values are generated from a truncated normal
distribution and \eqn{z} is randomly selected from \code{(-1,1)} to
characterise up- or down-regulation.

Assuming that \eqn{Y ~is~  N(\mu, \sigma^2)}, and \eqn{A=[a_1,a_2]}, a subset of
  \eqn{-Inf <y < Inf}, the conditional distribution of \eqn{Y} given \eqn{A}
  is called  truncated normal  distribution:
  

\deqn{f(y, \mu, \sigma)= (1/\sigma) \phi((y-\mu)/\sigma) / (\Phi((a2-\mu)/\sigma) -
\Phi((a_1-\mu)/\sigma))}

 for \eqn{a_1 <= y <= a_2}, and 0 otherwise,

where \eqn{\mu} is the mean of the original Normal distribution before truncation,
\eqn{\sigma} is the corresponding standard deviation,\eqn{a_2} is the upper truncation point,
\eqn{a_1} is the lower truncation point, \eqn{\phi(x)} is the density of the
standard normal distribution, and \eqn{\Phi(x)} is the distribution function
of the standard normal distribution. For \code{simData} function, we
consider \eqn{a_1=log_2(\code{foldMin})} and \eqn{a_2=Inf}. This ensures that the
biomarkers are differentially expressed by a fold change of
\code{foldMin} or more.
  
}

\value{
 A dataframe of dimension \code{nTrain} by \code{nBiom+1}. The first
 column is a factor (\code{class}) representing the group memberships of
 the samples.
}


\author{

 Mizanur Khondoker, Till Bachmann, Peter Ghazal \cr
  Maintainer: Mizanur Khondoker \email{mizanur.khondoker@gmail.com}.}

\references{
Khondoker, M. R., Till T. Bachmann, T. T., Mewissen, M., Dickinson, P. \emph{et al.}(2010).
Multi-factorial analysis of class prediction error: estimating optimal number of biomarkers for various classification rules.
\emph{Journal of Bioinformatics and Computational Biology}, \bold{8}, 945-965.
}

\seealso{\code{\link{classificationError}}}
\examples{

simData(nTrain=10,nBiom=3)
}

\keyword{datagen}
