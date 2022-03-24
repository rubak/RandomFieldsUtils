
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
## Copyright (C) 2018 -- 2021 Martin Schlather
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


##.C <- function(...) { cat("\ndebugging version: "); x <- ...elt(1); print(if (is.list(x)) x[1] else x); .Call(C_DebugCall); z <- .C(...);  cat("\nend call: ");  .Call(C_DebugCall); z}

##.Call <- function(...) { cat("\ndebugging version: "); x <- ...elt(1); print(if (is.list(x)) x[1] else x); .Call(C_DebugCall); z <- .Call(...); cat("\nend call: "); .Call(C_DebugCall); z}



.onLoad <- function(lib, pkg) {
  .C("loadoptions")
}


.onAttach <- function (lib, pkg) {
  n <- as.integer(0.75 * parallel::detectCores(logical=FALSE) +
                  0.25 * parallel::detectCores(logical=TRUE))
  if (!is.na(n)) .C("setCPUs", n) # or 1 if not fully installed
  n <- .C("recompilationNeeded", n)
  if (n[[1]] > 0) packageStartupMessage(paste("Consider the use of one of \n\tRFoptions(install.control=list(force=FALSE))\n\tRFoptions(install.control=list(repos=NULL))\n\tRFoptions(install.control=NULL)\nIf you use it, the package might run faster the next time it is loaded.\nRecompilation is at your own risk. There is currently no guarantee that it will work.\nIn case of problems:", CONTACT))
  else if (n[[1]] < 0)
    packageStartupMessage(paste(if (n[[1]] == -1) "R" else "Auto-r",
                                "ecompilation didn't (fully) succeed, see 'RFoptions(install.control=list(pkgs=NULL))'.\n"))
 }

.onDetach <- function(lib) {
## do not use the following commmands in .onDetach!
}

.onUnload <- function(lib, pkg){
#   RFoptions(storing=FALSE) ## delete everything
  .C("detachoptions")
}
