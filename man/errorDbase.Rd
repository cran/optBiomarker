\name{errorDbase}
\docType{data}
\alias{errorDbase}
\title{
Database of leave-one-out cross validation errors for various
combinations of data characteristics
}
\description{
  This is a 7-dimensional array (database) of  leave-one-out cross
  validation errors for Random Forest, Support Vector Machines, Linear
  Discriminant Analysis and k-Nearest Neighbour classifiers. The
  database is the basis for estimating the optimal number of biomarkers
  at a given error tolerance level using \code{\link[optBiomarker]{optimiseBiomarker}}
  function. See \bold{Details} for more information.
  }

  \usage{data(errorDbase)}

  \details{
    The following table gives the dimension names, lengths and values/levels
    of the data object \code{errorDbase}.
    
    \tabular{lll}{
      \bold{Dimension name} \tab \bold{Length} \tab \bold{Values/Levels} \cr
      No. of biomarkers     \tab 14\tab (1-6, 7, 9, 11, 15, 20, 30, 40,
      50, 100) \cr
      Size of replication   \tab 5 \tab (1, 3, 5, 7, 10) \cr
      Biological variation (\eqn{\sigma_b})\tab 4\tab (0.5, 1.0, 1.5,
      2.5) \cr
      Experimental variation (\eqn{\sigma_e}) \tab 4 \tab (0.1, 0.5,
      1.0, 1.5) \cr
      Minimum (Average) fold change \tab 4 \tab (1 (1.73), 2(2.88), 3(4.03),
      5(6.33))\cr
      Training set size   \tab 5 \tab (10, 20, 50, 100, 250) \cr
      Classification method \tab 3 \tab (Random Forest, Support Vector
      Machine, k-Nearest Neighbour)\cr
      }
      
      We have a plan  to expand the database to a 8-dimensional one by
      adding  another dimension to store error rates at different level
      of correlation between  biomarkers. Length of each dimension will
      also be increased leading to a bigger database with a wider coverage
      of the parameter space. Current version of the database contain error rates
      for independent (correlation = 0) biomarkers only. Also, it does not
      contain error rates for Linear Discriminant Analysis, which we plan
      to implement in the next release of the package. With the current
      version of the database, optimal number of biomarkers can be
      estimated using the  \code{\link[optBiomarker]{optimiseBiomarker}}
      function for any intermediate values of the factors represented by
      the dimensions of the database. 
    }

    \format{7-dimensional numeric array.}

    \author{
Mizanur Khondoker, Till Bachmann, Peter Ghazal \cr
 Maintainer: Mizanur Khondoker \email{mizanur.khondoker@googlemail.com}.}

   
    \seealso{
      \code{\link[optBiomarker]{optimiseBiomarker}}
      }
   

    \keyword{datasets}
