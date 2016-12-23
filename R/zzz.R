
.onLoad <- function(lib, pkg) {
  .Call("attachRFoptionsUtils")
}

.onAttach <- function (lib, pkg) {
#   packageStartupMessage("This is RandomFieldsUtils Version: 0.3.8");
}

.onDetach <- function(lib) {
# .Call("detachRFoptionsUtils")
}

.onUnload <- function(lib, pkg){
  .Call("detachRFoptionsUtils")
}
