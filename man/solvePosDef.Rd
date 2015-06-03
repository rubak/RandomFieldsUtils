\name{solve}
\alias{solvePosDef}
\alias{solve}
\title{Solve a System of Equations for Positive Definite Matrices}
\description{
  This function solves the equality \eqn{a x = b} for \eqn{x}
  where \eqn{a} is a positive definite matrix and \eqn{b} is a vector or a
  matrix. It is slightly faster than the inversion by the
  \code{\link[base]{chol}}esky decomposition and clearly faster than
  \code{\link[base]{solve}}.
  It also returns the logarithm of the
  determinant at no additional computational costs.
}

\usage{
solvePosDef(a, b=NULL, logdeterminant=FALSE)
%, sparse=NA, method=-1)
}
\arguments{
  \item{a}{a square numeric or complex matrix containing the
          coefficients of the linear system.  Logical matrices are
          coerced to numeric.
  }
  \item{b}{
    a numeric or complex vector or matrix giving the right-hand
    side(s) of the linear system.  If missing, \code{b} is taken to be
    an identity matrix and \code{solvePosDef} will return the inverse of
    \code{a}.
  }
  \item{logdeterminant}{logical.
    whether the logarithm of the determinant should also be returned
  }
 % \item{sparse}{logical or \code{NA}.
%    If \code{NA} the function determines whether a sparse
%    matrix algorithm of the package \pkg{spam} should be used.
%  }
%  \item{method}{integer vector.
%    If the sparse matrix algorithm is not used, \code{method}
%    determines the alternative algorithm. See Details.
%  }  
}
\value{
  If \code{logdeterminant=FALSE} the function returns a vector or a matrix,
  depending on \code{b} which is the solution to the linear equation.
  Else the function returns a list containing both
  the solution to the linear equation and
  the logarithm of the
  determinant of \code{a}.
}
\details{
%  The values of \code{method} could be:
%  \itemize{
%    \item \code{<0} :
If the matrix is diagonal direct calculations are performed.

Else if the matrix is sparse the package \pkg{spam} is used.

Else the Cholesky decomposition is tried.
If it fails, the SVD is tried.
%    \item \code{0} : Cholesky
%    \item \code{1} : QR (not programmed yet)
%    \item \code{2} : SVD
%    \item \code{3} : LU
%    \item \code{4} : sparse matrix algorithm by \pkg{spam}
% }
}

\references{
  See \link[spam]{chol.spam} of the package \pkg{spam}
}
\seealso{
  \link[spam]{chol.spam} in the package \pkg{spam}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}, and
  Reinhard Furrer
}
\keyword{math}
\examples{

## This examples shows that 'solvePosDef' can be much faster
## than 'solve'

## creating a covariance matrix for a temporal process
covmatrix <- function(model, x) {
  x <- as.matrix(dist(x))
  return(eval(substitute(model)))
}

size <- 600
x <- runif(size, 0, size / 10)
M <- covmatrix((1 - x) * (x < 1) , x) ## Askey's model of covariance
b <- seq(0, 1, len=size)
unix.time(inv2 <- solve(M, b))
unix.time(inv1 <- solvePosDef(M, b)) ## much faster in this case
range(inv2 - inv1)

}

