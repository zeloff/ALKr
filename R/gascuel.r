#' @rdname inverse_ALK
#' @export
gascuel <- function(x, fi1, fi2, initial_values,
                    threshold = 0.0001, maxiter = 2000,
                    age_classes = colnames(x), length_classes = rownames(x)) {
  
  converged <- FALSE
  
  innerGascuel <- function(params, optimizing = TRUE, threshold) {
    sigma <- params[1] + params[2] * lj + params[3] * c(lj[1], diff(lj))
    
    # Penalize impossible sigma values
    if(any(sigma <= 0, na.rm = TRUE)) return(Inf)
    
    invALK <- outer(li, lj, dnorm,
                    sd = matrix(rep(sigma, length(li)),
                                nrow = length(li),
                                ncol = length(lj),
                                byrow = TRUE))
    
    pj <- rep(1 / length(lj), length(lj))
    pj.prev <- pj + threshold * 2
    
    iter <- 0
    
    while(sum(abs(pj.prev - pj)) > threshold & iter < maxiter) {
      pj.prev <- pj
      
      pij <- t(t(invALK) * pj)
      pij[is.na(pij)] <- 0
      
      pj <- colSums(pij * (pi_ / rowSums(pij)))
      iter <- iter + 1
    }
    
    if (iter < maxiter) converged <<- TRUE
    
    if (optimizing) return(sum((pj.prev - pj)^2))
    else {
      result <- t(t(invALK) * pj) * sum(fi2)
      result[is.na(result)] <- 0
      return(result)
    }
  }
  
  nij <- fi1 * calc_ALK(x)
  li <- as.numeric(rownames(nij))
  lj <- colSums(nij * li) / colSums(nij)
  pi_ <- fi1 / sum(fi1)
  
  optimal <- optim(initial_values, fn = innerGascuel, threshold = threshold)
  
  if (optimal$convergence == 1) warning("Parameter optimization exceeded maxiter")
  if (optimal$convergence == 10) warning("Degeneracy of the Nelder-Mead simplex optimization")
  
  result <- calc_ALK(innerGascuel(optimal$par, optimizing = FALSE, threshold = threshold))
  rownames(result) <- rownames(nij)
  
  new("ALKr",
      alk = result,
      N = nij,
      age_classes = age_classes,
      length_classes = length_classes,
      method = "Gascuel",
      parameters = list(
        ConvergenceThreshold = threshold,
        alpha = optimal$par[1],
        beta = optimal$par[2],
        gamma = optimal$par[3],
        Converged = optimal$convergence == 0)
    )
}
