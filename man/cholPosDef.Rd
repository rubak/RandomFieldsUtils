\name{chol}
\alias{cholPosDef}
\alias{chol}
\title{Cholesky Decomposition of Positive Definite Matrices}
\description{
  This function calculates the Choleskey decomposition of a matrix.
}

\usage{
cholPosDef(a)
%, sparse=NA, method=-1)
}
\arguments{
  \item{a}{a square real-valued positive definite matrix
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
  a matrix containing the Choleskey decomposition (in its upper part)
}
\details{
If the matrix is diagonal direct calculations are performed.

%Else if the matrix is sparse the package \pkg{spam} is used.

Else the Cholesky decomposition is tried.

}

% \references{ See \link[spam]{chol.spam} of the package \pkg{spam} }

\seealso{
  \link[spam]{chol.spam} in the package \pkg{spam}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{math}
\examples{

## This examples shows that 'cholPosDef' can be much faster
## than 'chol'

## creating a covariance matrix for a temporal process
covmatrix <- function(model, x) {
  x <- as.matrix(dist(x))
  return(eval(substitute(model)))
}

size <- 600
x <- runif(size, 0, size / 10)
M <- covmatrix((1 - x) * (x < 1) , x) ## Askey's model of covariance
b <- seq(0, 1, len=size)
unix.time(C2 <- chol(M))
unix.time(C1 <- cholPosDef(M))
range(C2 - C1)
stopifnot(all(abs(C2 - C1) < 10^{-9}))

} 

