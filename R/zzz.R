
.onLoad <- function(lib, pkg) {
  .Call("attachRFoptionsUtils")
}

.onAttach <- function (lib, pkg) {
   packageStartupMessage("This is RandomFieldsUtils Version: 0.3.3");
}

.onDetach <- function(lib) {
}

.onUnload <- function(lib, pkg){
    .Call("detachRFoptionsUtils")
}
