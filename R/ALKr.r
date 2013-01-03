setClass("ALKr",
  representation(
    alk = "matrix",
    N = "matrix",
    method = "character",
    parameters = "list"),
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
  function(.Object, alk, N, method, parameters) {
    .Object@alk <- alk
    .Object@N <- N
    .Object@method <- method
    .Object@parameters <- parameters
    .Object
  }
)

setMethod("show",
  signature(object = "ALKr"),
  function(object) {
    print(object@alk)
    cat(paste("\nMethod:", object@method,"\n"))
    print(matrix(object@parameters, dimnames = list(names(object@parameters), "Value")))
  }
)

summary.ALKr <-
  function(object, length_classes = as.numeric(rownames(object@alk))) {
    
    nj <- colSums(object@N)
    i <- colSums(length_classes * object@N)
    lj <- i / nj
      
    vlj <- colSums(object@N * i^2) / (nj - 1) - 2 * lj * i / (nj - 1) + lj^2 * nj / (nj - 1)
    
    result <- list(pj = nj / sum(nj), mean_lj = lj, var_lj = vlj)
    
    cat("Proportion of age-class:")
    print(result$pj)
    cat("\nMean length at age:")
    print(result$mean_lj)
    cat("\nVariance of length at age:")
    print(result$var_lj)
    cat("\n")
    cat(paste("\nMethod:", object@method,"\n"))
    print(matrix(object@parameters, dimnames = list(names(object@parameters), "Value")))
    
    invisible(result)
  }