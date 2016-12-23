

summary.RFopt <- function(object, ...) {
  object <- lapply(object, function(z) z[order(names(z))])
  object <- object[c(1, 1 + order(names(object[-1])))]
  class(object) <- "summary.RFopt"
  object
}


print.summary.RFopt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFopt <- function(x, ...) {
  print.summary.RFopt(summary.RFopt(x, ...)) #
  invisible(x)
}

summary.RFoptElmnt <- function(object, ...) {
  object <- object[order(names(object))]
  class(object) <- "summary.RFoptElmt"
  object
}

print.summary.RFoptElmnt <- function(x, ...) {
  str(x, give.attr=FALSE, ...) #
  invisible(x)
}

print.RFoptElmnt <- function(x, ...) {
  print.summary.RFoptElmnt(summary.RFoptElmnt(x, ...)) #
  invisible(x)
}

RFoptions <- function(..., no.readonly=TRUE) {
##  on.exit(.C("RelaxUnknownRFoption", FALSE))
  ##  .C("RelaxUnknownRFoption", TRUE)
  opt <- lapply(.External(C_RFoptions, ...),
                function(x) {
                  class(x) <- "RFoptElmnt"
                  x
                })
  if (length(opt)!=0) {      
    class(opt) <-  "RFopt"
    if (!no.readonly) {
      opt$readonly <- list()
    }
  }
  if (length(opt)==0) {
 #   O <- opt
#    names(O) <- NULL
#    opt <- c(opt, unlist(O))
    invisible(opt)
  } else opt
}
