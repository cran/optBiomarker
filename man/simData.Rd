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
        sdB=1.0,rho=0,
        sigma=1,diffExpr=TRUE,
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
   \item{rho}{Common Pearson correlation between biomarkers. To ensure
  positive definiteness, allowed values of \code{rho} are restricted between 0 and 0.95 inclusive.}
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
  Maintainer: Mizanur Khondoker \email{mizanur.khondoker@googlemail.com}.}


\seealso{\code{\link{classificationError}}}
\examples{

simData(nTrain=10,nBiom=3)
}

\keyword{datagen}
