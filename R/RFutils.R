## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  


sleep.milli <- function(n) {
  .C(C_sleepMilli, as.integer(n))
  invisible(NULL)
}

sleep.micro <- function(n) {
  .C(C_sleepMicro, as.integer(n))
  invisible(NULL)
}
	
hostname<-function(){.C(C_hostname, h=paste(seq(0,0,l=100), collapse=""),
                        as.integer(100))$h}

pid <- function() {.C(C_pid, i=integer(1))$i}


FileExists <- function(file, printlevel=RFoptions()$basic$printlevel) {
    ## for parallel simulation studies: the same data output file should not
  ## be created twice. So:
  ## 1. if file exists then assume another process has done the work already
  ## 2. if file.lock existss then assume another process is doing the work
  ## 3.a. otherwise create file.lock to show other processes that the process
  ##      will do the work
  ## 3.b. check if another process has started with the same work at the same
  ##      time it may happen that in case of simulatenous creation of file.lock
  ##      no process will do the work...(then the lock file will rest.)
  PL_ERRORS <- 6
  lock.ext <- ".lock";
  if (file.exists(file)) { #1.
    if (printlevel>=PL_ERRORS ) cat("'", file, "' already exists.\n");
    return(1)
  } else { 
    LockFile <- paste(file, lock.ext, sep="")
    if (file.exists(LockFile)) { #2.
      if (printlevel>=PL_ERRORS ) cat("'",file,"' is locked.\n");
      return(2);
    }
    PID <- pid();
    write(file=LockFile,c(PID,hostname()),ncolumns=2,append=TRUE); #3.a.
    Pid <- matrix(scan(LockFile,what=character(0), quiet=TRUE),nrow=2)
    if ((sum(Pid[1,]==PID)!=1) || (sum(Pid[1,]>PID)>0)){ #3.b.
      if (printlevel>PL_ERRORS )
        cat("Lock file of '", file, "' is knocked out.\n");
      return(3);
    }
  }
  return(0);
}

LockRemove <- function(file) {
  ## removes auxiliary files created by FileExists
  lock.ext <- ".lock";
  file.remove(paste(file, lock.ext, sep=""))
}



Print <- function(..., digits=6, empty.lines=2) { #
  ## ?"..1"
#  print(..1)
#  print(substitute(..1))
#   print(missing(..100))
   
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



cholPosDef <- function(a) {
 # return(.Call("CholPosDef", a, PACKAGE="RandomFieldsUtils"))
  .Call(C_CholPosDef, a)
 }


solvePosDef <- function(a, b=NULL, logdeterminant=FALSE) {
  if (logdeterminant) {
    logdet <- double(1)
    res <- .Call(C_SolvePosDef, a, b, logdet)
    return(list(inv=res, logdet=logdet))
  } else {
    .Call(C_SolvePosDef, a, b, double(0))
  }
}

sortx <- function(x, from=1, to=length(x),
                        decreasing=FALSE, na.last = NA) {
  n <- length(x)
 if (n <= 4000 || (to - from) < (0.35 + is.double(x) * 0.15) * n) {   
    if (decreasing) {
      x <- -x
      if (!is.na(na.last)) na.last <- !na.last
    }
    ans <- .Call(C_sortX, x, as.integer(from), as.integer(to),
                 as.logical(na.last))
    return(if (decreasing) -ans else ans)
 } else {
    return(if (from==1 && to==n)
           sort(x, decreasing=decreasing, na.last=na.last) else
           sort(x, decreasing=decreasing, na.last=na.last)[from:to])
  }
}


orderx <- function(x, from=1, to=length(x),
                         decreasing=FALSE, na.last = NA) {
#  cat(length(x) / (0.35 + 0.14 * log(length(x)))^2, "")
  if ((to - from) * (0.35 + 0.14 * log(length(x)))^2 > length(x)) { #10^2:1, 10^3:1.5, 10^4:3 10^5:5 10^6:5, 10^7: 8, 10^8:10,
    #cat("old\n");
    return(if (from==1 && to==length(x))
           order(x, decreasing=decreasing, na.last=na.last) else
           order(x, decreasing=decreasing, na.last=na.last)[from:to])
  }
  if (decreasing) {
    x <- -x
    if (!is.na(na.last)) na.last <- !na.last
  }

  .Call(C_orderX, x, as.integer(from), as.integer(to),
        as.logical(na.last))
}
