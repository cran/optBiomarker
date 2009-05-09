\name{internal}
\alias{tnormAvg}
\alias{objFun}
\alias{avg2Lower}
\alias{plot3dFun}
\alias{dimSelect}
\alias{yapprox}
\alias{matapprox}
\alias{classpredict.lda}
\alias{classpredict.knn}


\title{Internal functions for optBiomarker package.}
\description{
 Internal functions for \code{optBiomarker} package.
}
\usage{
tnormAvg(a1=-Inf, a2=Inf,mu=0,sigma=1)
avg2Lower(truncAvg=log2(2), mu=0,sigma=1)
objFun(para, mu=0,sigma=1,truncAvg=log2(2))
plot3dFun(panel,...)
dimSelect(x, val)
yapprox(y, x, xout, method="linear", n=50,
          yleft, yright, rule = 1, f = 0, ties = mean)
matapprox<-function(mat, rx, cx, nr=100,nc=20,rout, cout)

classpredict.lda(object, newdata)
classpredict.knn(object, newdata) 

}

\details{
These functions are not to be called by the user.
}

\keyword{internal}


%###############################################################################
%## End of                           internal.Rd
%###############################################################################
