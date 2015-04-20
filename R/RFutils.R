
.struve <- function(x, nu, sign, expon.scaled) {
  storage.mode(x) <- "double"
  storage.mode(nu) <- "double"
  storage.mode(expon.scaled) <- "logical"
  storage.mode(sign) <- "double"
#  res <- double(max(length(x), length(nu)))
 .Call("struve", x, nu, sign, expon.scaled)
}


struveH <- function(x, nu)  .struve(x, nu, -1, FALSE)
struveL <- function(x, nu, expon.scaled=FALSE)  .struve(x, nu, 1, expon.scaled)
I0L0 <- function(x) {
  storage.mode(x) <- "double"
#  res <- double(length(x))
  .Call("I0ML0", x)
}


solvePosDef <- function(a, b, logdeterminant=FALSE) {  
  if (!is.double(a)) storage.mode(a) <- "double"
  if (ncol(a) != nrow(a)) stop("'a' is not a square matrix")

  if (missing(b)) {
     b <- NULL
  } else {
    if (!is.double(b)) storage.mode(b) <- "double"
    if (ncol(a) != if (is.vector(b)) length(b) else nrow(b))
      stop("'a' and 'b' do not match in size")
  }
  
  if (logdeterminant) {
    logdet <- double(1)
    res <- .Call("solvePosDef", a, b, logdet)
    return(list(inv=res, logdet=logdet))
  } else {
    .Call("solvePosDef", a, b, double(0))
  }
}


Print <- function(..., digits=6, empty.lines=2) { #
  max.elements <- 99
  l <- list(...)
  n <- as.character(match.call())[-1]
  cat(paste(rep("\n", empty.lines), collapse="")) #
  for (i in 1:length(l)) {
    cat(n[i]) #
    if (!is.list(l[[i]]) && is.vector(l[[i]])) {
      L <- length(l[[i]])
      if (L==0) cat(" = <zero>")#
      else {
        cat(" [", L, "] = ", sep="")
        cat(if (is.numeric(l[[i]]))
            round(l[[i]][1:min(L , max.elements)], digits=digits)#
            else l[[i]][1:min(L , max.elements)]) #
        if (max.elements < L) cat(" ...")
      }
    } else {
       if (is.list(l[[i]])) {
        cat(" =  ") #
        str(l[[i]], digits.d=digits) #
      } else {
        cat(" =")
        if (length(l[[i]]) <= 100 && FALSE) {
          print(if (is.numeric(l[[i]])) round(l[[i]], digits=digits)#
                else l[[i]])
        } else {
          if (length(l[[i]]) > 1 && !is.vector(l[[i]]) && !is.matrix(l[[i]])
              && !is.array(l[[i]])) cat("\n")
          str(l[[i]]) #
        }
      }
    }
    cat("\n")
  }
}

