\name{Instruction Set}
\alias{uses.simd.instruction}
\alias{misses.simd.instruction}
\title{
  CPU instruction set
}
\description{
   The function checks whether a certain instruction is used (missed) under
   the current compilation of a package.
}
\usage{
uses.simd.instruction(which=NULL, pkgs=NULL)
misses.simd.instruction(which=NULL, pkgs=NULL)
}
\arguments{ 
  \item{which}{character vector with values in \code{"SSE2", "SSSE3", "AVX", "AVX2", "CUDA"}}
  \item{pkgs}{character vector or missing.}
}

%\details{}
 
\value{
  logical vector  of length \code{which}
  or matrix with number of rows equal to
  the length of \code{which}. An element is \code{TRUE} if
  the instruction set is used (missed) by the package.

  If an arguments is \code{NULL} all available information is given.
}


\me

\examples{ %    library(RandomFieldsUtils)
  uses.simd.instruction()
  misses.simd.instruction()
}

%library(gtools); keywords()
\keyword{sysdata}
