
.onLoad <- function(lib, pkg) {
  .Call("attachRandomFieldsUtils", interactive())
}

.onAttach <- function (lib, pkg) {
#   packageStartupMessage("This is RandomFieldsUtils Version: 0.5");
}

.onDetach <- function(lib) {
# .Call("detachRanodmFieldsUtils")
}

.onUnload <- function(lib, pkg){
  .Call("detachRandomFieldsUtils")
}
