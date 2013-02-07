setClass("ALKr",
  representation(
    alk = "matrix",
    N = "matrix",
    method = "character",
    parameters = "list",
    name = "character",
    description = "character"),
  validity = function(object) {
    if (nrow(object@alk) != length(object@length_classes))
      return("Incorrect number of length classes")
    if (ncol(object@alk) != length(object@age_classes))
      return("Incorrect number of age classes")
    if (!identical(dim(alk), dim(N)))
      return("alk and N matrices have different dimensions")
    return(TRUE)
  }
)

setMethod("initialize", "ALKr",
  function(.Object, alk, N, method, parameters, name = "", description = "") {
    .Object@alk <- alk
    .Object@N <- N
    .Object@method <- method
    .Object@parameters <- parameters
    .Object@name <- name
    .Object@description <- description
    .Object
  }
)

setMethod("show",
  signature(object = "ALKr"),
  function(object) {
    if (object@name != "") cat(paste(object@name,"\n"))
    if (object@description != "") cat(paste(object@description,"\n\n"))
    print(object@alk)
    cat(paste("\nMethod:", object@method,"\n"))
    if (length(object@parameters) > 0)
      print(matrix(object@parameters, dimnames = list(names(object@parameters), "Value")))
  }
)