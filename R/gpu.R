
gpu_info <- function(devices = NULL){
    # Query function for device info
    # devices should be a vector of integers with the correct device number
    # if no device argument is given, the standard value is taken from RFoptions
  
  if (length(devices) == 0) dev <- RFoptions()$installNrun$gpuDevices
  else if (!is.vector(devices) || any(devices != (dev <- as.integer(devices))))
    stop("Devices have to be a vector of integers")

  gpu_list <- .Call("gpu_info", as.integer(dev))
  class(gpu_list) <- "gpu_list"

  return(gpu_list)
}

print.gpu_list <- function(x, ...){
  ## pretty print function for a gpu_list instance returned by gpu_info
  if(!is(x, "gpu_list")) stop("Wrong argument type")
  df <- do.call(rbind.data.frame, x)
  rownames(df) <- NULL
  print(df)  # OK
}
