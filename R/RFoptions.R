
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2021 Martin Schlather
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

detach_packages <- function(pkgs) {
  for (pkg in pkgs) {
    pkg <- paste0("package:", pkg)
   while(pkg %in% search()) detach(pkg, unload = TRUE, character.only=TRUE)
  }
}
libraries <- function(pkgs, control, verbose=FALSE) {
  if (length(control) > 0) {
    idx <- pmatch(names(control), names(as.list(args(library))))
    control <- control[idx[!is.na(idx)]]
  }
  for (pkg in pkgs) do.call("library", c(list(pkg), control))
  if (verbose) message("libraries attached.")
}

OneTo <- function(n)
  return(if (length(n) > 1) stop("invalid end of loop") else if (n < 1)
                                                              NULL else 1:n)
S <- function(x) if (length(x) > 1) "s" else ""
ARE <- function(x) if (length(x) > 1) "are" else "is"
HAVE <- function(x) if (length(x) > 1) "have" else "has"


sources <- function(pkgs, raw=FALSE, repos=NULL) {
  gitrepos <- "schlather/PACKAGES"
  gitinfo <- "https://github.com/"
  gitdownload <- "https://raw.githubusercontent.com/"


  debug <- FALSE
  ip <- installed.packages()[pkgs, "Version"] # OK
  names(ip) <- pkgs
  s <- c("local", "cran", "github")
  found <- matrix(FALSE, nrow=length(pkgs), ncol=length(s))
  V <- where <- matrix("", nrow=length(pkgs), ncol=length(s))
  dimnames(V) <- dimnames(where) <- dimnames(found) <-list(pkgs, s)

  for (frm in c("local0", s)) {
    from <- frm
    if (from == "local0") {
      from <- "local"
      url <- ""
    } else if (from == "local") url <- getwd()
    else {
       if (from == "cran") {     
        type <- "source"
        if (length(repos) == 0) repos <- getOption("repos")
        if (debug) print(repos) ## OK
        cran <- NULL
        url <- try(contrib.url(repos=repos, type="source"))
        if (!is(url, "try-error")) {
          cran <- try(available.packages(contriburl = url)[pkgs, "Version"])
        if (is(cran, "try-error") || length(cran) == 0) next
        }
        if (length(cran) == 0) next
      } else if (from == "github") {    
        url <- paste0(gitinfo, gitrepos)
        github <- try(grep("tar.gz", fixed=TRUE, readLines(url), value = TRUE))
        if (is(github, "try-error") || length(github) == 0) next
      } else stop("BUG")
    }

    for (i in 1:length(pkgs)) {
      if (from == "cran") {
        versions <- cran[i] ## length 1
      } else {
        if (from == "local") {
          if (url == "") f <- dir(pattern=paste0(pkgs[i], "_.*\\.tar\\.gz"))
          else f <- dir(pattern=paste0(pkgs[i], "_.*\\.tar\\.gz"), path=url)
        } else {
           f <- grep(paste0(pkgs[i],"_"), github, value = TRUE)
        }
        if (length(f) > 0) {
          pkg <- paste0(pkgs[i],"_")
          versions <- sapply(strsplit(f, "\\.tar\\.gz"), function(x) {
             s <- strsplit(x[1], pkg)[[1]]
            s[length(s)]
            })
         } else versions <- NULL
      }
      old.version <- ip[i]
      where[i, from] <- url
      for (j in OneTo(length(versions))) {        
        cmp <- compareVersion(versions[j], ip[i])
        if (cmp >= 0) {
          found[i, from] <- TRUE
          if (compareVersion(versions[j], old.version)) {
            old.version <- versions[j]
            V[i, from] <- versions[j]
          }
        }
      }
    }

    if (frm == "local") { ## NOT 'from '
      if (all(anyfound <- apply(found, 1, any))) break; ## all found locally
    }
  }

  
  if (debug) Print(list(where=where, found=found, newer.version=V, ip=ip)) ## OK
  if (raw) return(list(where=where, found=found, newer.version=V, ip=ip))

  failed <- !apply(found, 1, any)
  if (any(failed)) {
    if (all(failed)) return(list(what=NULL, failed=failed))
    where <- where[!failed, , drop=FALSE]
    found <- found[!failed, , drop=FALSE]
    V <- V[!failed, , drop=FALSE]
    ip <- ip[!failed]
    pkgs <- pkgs[!failed]
  }

  what <- matrix("", nrow=length(ip), ncol=4)
  dimnames(what) <- list(names(ip), c("how", "where", "version", "call"))
  method <- colnames(V)

  if (all(apply(V == "", 1, any, na.rm=TRUE))) {## take
    ##  all current iff all current are available. This is the safest.
    found[V != ""] <- FALSE 
    dim(found) <- dim(where)
   } else if (all(what[, "cran"] != "")) found[,  method != "cran"] <- FALSE ## take
  ## cran versions if all cran vesions are available; second safest since this necessitates
  ## that R version is recent enough
  ## Otherwise try the best, i.e. take always the most recent ones -- this reduced
  ## probability of version incompatibilities
  for (i in 1:length(ip)) {
    if (length(f <- which(found[i,])) == 0) next    
    newest <- f[1]
    for (j in f[-1]) if (compareVersion(V[j], V[newest]) > 0) newest <- j
    what[i, 1:3] <- c(method[newest], where[i, newest],
                      if (V[i, newest] == "") ip[i] else V[i, newest])
  }
  idx <- what[, "how"] == "local"
  path <- what[idx, "where"]
  add <- path != "" & substring(path, nchar(path)) != .Platform$file.sep
  path[add] <- paste0(path, .Platform$file.sep)
  what[idx, "call"] <- paste0(path, pkgs[idx], "_", what[idx, "version"], ".tar.gz")
  idx <- what[, "how"] == "github"
  what[idx, "call"] <- paste0(gitdownload, gitrepos, "/main/", pkgs[idx], "_",
                              what[idx, "version"], ".tar.gz")
  idx <- what[, "how"] == "cran"
  what[idx, "call"] <- pkgs[idx]

  if (debug) Print(t(what), failed) ## OK
  return(list(what=what, failed=failed))
}

#    pkgs <- c("RandomFieldsUtils", "miraculix", "RandomFields");print("XX");  print(s <- sources(pkgs));  tmp <- apply(found, 1, any)
# https://raw.githubusercontent.com/schlather/PACKAGES/main/miraculix_1.0.2.tar.gz

reinstallPackages <- function(ic, basic, install.control) {
  install <- basic$install
##  Print(basic)
  verbose <- FALSE
  force <- quiet <- SERVER <- pkgs.given <- path.given <- FALSE
  repos <- path <- pkgs <- NULL
  if (ic) {
    N <- names(install.control)
    if ("pkg" %in% N)
      stop("'pkg' is an invalid option for 'install.control'. Did you mean 'pkgs'?")
    pkgs.given <- "pkgs" %in% N
    path.given <- "path" %in% N
    path <- install.control$path
    delete <- c("repos", "path", "force", "pkgs", "SERVER")
    for (arg in c(delete, "verbose", "quiet"))
      if (length(install.control[[arg]]) > 0) {
        assign(arg, install.control[[arg]])
        if (arg %in% delete) install.control[[arg]] <- NULL
      }

    if (length(install.control$force) > 0 && !force) install <- "ask"
    else if (length(install) > 0 && install %in% c("ask", "none"))
      install <- "install"
  }

  if (!pkgs.given) pkgs <- .Call(C_getPackagesToBeInstalled, force) 

  verbose <- verbose && !quiet
  if (length(pkgs) == 0) {
    .Call(C_AVXmessages, "all")
    if (!quiet)
      message(if (!pkgs.given) "No packages found to be installed.",
              if (!path.given && !pkgs.given)
                " Consider setting, in 'install.control', a path to a local directory.",
              if (verbose) " This happens particularly if the the installation process was interrupted. Try it again in the next session or use 'RFoptions(install.control=list(force=TRUE))' for instance.")
    return()
  }
  
  if (install == "ask") {
    if (!quiet)
      cat("The package", S(pkgs), " ", paste0("'", pkgs, "'", collapse=", "),
          " ", HAVE(pkgs), " been compiled without appropriate SIMD/AVX2 flags. So, calculations can be slow. If the package",
          S(pkgs), " ", ARE(pkgs),
          " recompiled with the necessary flags, the calculations might be faster.\nR should be restarted after re-compiling. The argument 'install.control' might be used to run the re-compilation without asking and to pass further arguments to 'install.packages', e.g., 'RFoptions(install.control=list(verbose=TRUE))'\nTo avoid this feedback, set 'RFoptions(install=\"none\")' or 'RFoptions(install=\"install\")' before calling any other function of '",
          pkgs[length(pkgs)],"'.\n\n", sep="")

    omp <- .Call(C_AVXmessages, pkgs)
  }

  ## pkgs <- c("RandomFieldsUtils", "miraculix", "RandomFields");print("XX")
  if (!quiet) cat("Searching for tar balls... ")
  s <- sources(pkgs,repos=repos)
  cat("\n")
  if (all(s$failed)) {
    if (!quiet) cat("Not a single source found for re-installation.\n")
    return()
  }

  tell.which <- function(s, verbose) {
    cat("The following package", S(!s$failed), " will be re-installed:\n", sep="",
        paste0(if (!verbose) "\t",
               rownames(s$what), "_", s$what[, "version"],
               " from ", s$what[, "how"],
               if (verbose) ", ",  if (verbose) s$what[, "where"], "\n")
        )
    if (any(s$failed)) {
      cat("No recent tar ball found for ",
          paste0("'", names(s$failed)[s$failed], "'", collapse=", ", sep=""),
          ". ", sep="")
      if (verbose) 
        cat("Consider calling\n\t'RFoptions(install.control=list(path=\"<local directory>\",\n\t\t\tverbose=TRUE))'")
      cat("\n")
    }
  }
##  tell.which(s, verbose)
  
  if (install == "ask") {
    if (!quiet) tell.which(s, verbose)
    repeat {
      txt <- paste0("Shall '", rownames(s$what)[1],
                    "' and all further packages based on 'RandomFieldsUtils' be recompiled (Y/n/h/s)erver/<args>) ? ")
     install.control <- readline(txt)
      if (install.control %in% c("h", "H")) {
        cat("\nHelp info\n=========\n")
        cat("Y : installation from \n")
        cat("n : interruption.\n     No further re-installation in this session possible\n")
        cat("s : SERVER=\"avx\"\n     This option guarantees downwards compatibility to avx. See ?RFoptions for details.\n")
        cat("<args>: any arguments for 'install.packages',\n    e.g. 'lib = \"~\", quite=TRUE'\n")
        cat("\n")
      } else break
    }
     
    install <- if (install.control %in% c("n", "N")) "none" else "install"
    path <- NULL
    if (install.control %in% c("s", "S")) SERVER <- "avx"
    if (nchar(install.control) <= 3)  install.control <-""
    if (verbose) {
      if (install == "none") {
        cat("If you have stopped the re-compilation since does not work, consider one of the following possiblities:")
        .Call(C_AVXmessages, NULL)
        cat("\nIf all fails, call 'RFoptions(install=\"none\")' after any loading of the package.\nOtherwise you will be bothered with the above question again and again.\n")
      } else {
        S <- "\t*************************************************\n"
        cat("\n", S, "\t***         Do not forget to restart R.       ***\n",S)
        sleep.milli(1500)
      }
    }
  } else {
    omp <- .Call(C_AVXmessages, "OMP")
    if (!quiet) tell.which(s, verbose)
  }


  if (install != "none") {
    if (is.character(install.control)) 
      install.control <- eval(parse(text=paste("list(", install.control, ")")))
    SIMD_FLAGS <- CXX_FLAGS <- args <- ""
    if (length(install.control$configue.args) > 0) {
      args <- install.control$configue.args
      install.control$configue.args <- NULL
    }
    if (length(install.control$CXX_FLAGS) > 0) {
      CXX_FLAGS <- install.control$CXX_FLAGS
      install.control$CXX_FLAGS <- omp <- NULL
    }
     if (length(install.control$CXX_FLAGS) > 0) {
      SIMD_FLAGS <- install.control$SIMD_FLAGS
      install.control$SIMD_FLAGS <- NULL
    }
    if (length(install.control$USE_GPU) > 0) {
      usegpu <- install.control$USE_GPU
      install.control$USE_GPU <- NULL
    } else usegpu <- "\"try\""
    
    idx <- pmatch(names(install.control),names(as.list(args(install.packages))))
    install.control <- install.control[which(!is.na(idx))]
        
    args <- paste0(args,
                   " SERVER=", SERVER, 
                   if (length(SIMD_FLAGS) > 0)
                     paste0(" SIMD_FLAGS='", SIMD_FLAGS, "'"),
                   if (length(CXX_FLAGS) + length(omp) > 0)
                     paste0(" CXX_FLAGS='", CXX_FLAGS, " ", omp, "'")
                   )
    if (verbose) Print(install.control, args) ## OK

    how <- s$what[, "how"]
    pkgs <- s$what[, "call"]
    for (p in 1:nrow(s$what)) {
       z <- Try(do.call("install.packages",
                        c(list(pkgs=pkgs[p], type="source",
                               repos =  if (how[p] == "cran")
                                          s$what[p, "where"] else NULL),
                          install.control,
                          configure.args=args)))
      if (is(z, "try-error")) print(z) ## OK
    }
    ## on.exit({detach_packages(rev(pkgs)); libraries(pkgs)}, add=TRUE)
  }
  cat("\n\n")
}


RFoptions <- function(..., no.readonly=TRUE, no.class=FALSE, install.control=NULL) {
  opt <- .External(C_RFoptions, ...)  
  ##  if (is.list(opt)) Print(basic) else Print(opt)
  ic <- hasArg("install.control")
##   print(opt)
  ## print(ic)
   if (ic || (length(opt) > 0 && is.list(opt) && is.list(opt$basic) &&
             opt$basic$installPackages && interactive())) {
     reinstallPackages(ic=ic, basic=opt$basic, install.control=install.control)
    if (ic) return(invisible(NULL))
  }
  if (length(opt) == 0 || no.class) return(invisible(opt))
  if (is.list(opt[[1]])) {
    opt <- lapply(opt,
		  function(x) {
		    class(x) <- "RFoptElmnt"
		    x
		})
    class(opt) <-  "RFopt"
  } else class(opt) <-  "RFoptElmnt"
  if (!no.readonly) {
    opt$readonly <- list()
  }
  opt
}
