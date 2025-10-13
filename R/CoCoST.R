# Save this code in a file named "CoCoST_function.R"

library(Matrix)
library(RSpectra)

CoCoST <- function(data1, W1, data2, W2, para = 0.2, Dim = 2) {
  # Input validation
  if (missing(data1)) stop("Foreground data not provided, please provide data")
  if (missing(data2)) stop("Needs to provide the background data")
  if (missing(W1)) stop("Affinity matrix for foreground data not provided!")
  if (missing(W2)) stop("Affinity matrix for background data not provided!")
  
  message(sprintf("Extracting %d contrastive components with parameter: %f", Dim, para))
  
  # Use sparse matrices for efficiency
  W1 <- as(W1, "sparseMatrix")
  W2 <- as(W2, "sparseMatrix")
  
  # Construct normalized graph Laplacians
  D1 <- Diagonal(x = rowSums(W1))
  D1inv_sqrt <- Diagonal(x = 1 / sqrt(rowSums(W1)))
  L1 <- D1 - W1
  normL1 <- D1inv_sqrt %*% L1 %*% D1inv_sqrt
  
  D2 <- Diagonal(x = rowSums(W2))
  D2inv_sqrt <- Diagonal(x = 1 / sqrt(rowSums(W2)))
  L2 <- D2 - W2
  normL2 <- D2inv_sqrt %*% L2 %*% D2inv_sqrt
  
  # Adjust Laplacians
  L11 <- Diagonal(nrow(D1)) - (0.2) * normL1
  L22 <- Diagonal(nrow(D2)) - (0.2) * normL2
  
  # Compute covariance matrices
  cov1 <- crossprod(data1, L11 %*% data1)
  cov2 <- crossprod(data2, L22 %*% data2)
  
  # Compute the projection matrix W
  RR <- cov1 - para * cov2
  eig <- eigs_sym(RR, k = Dim)  # RSpectra for faster eigen decomposition
  
  projMatrix <- eig$vectors
  compNames <- paste("CoCoST", 1:Dim, sep = " ")
  colnames(projMatrix) <- compNames
  
  # Extract contrastive principal components
  fgComponents <- data1 %*% projMatrix
  colnames(fgComponents) <- compNames
  
  bgComponents <- data2 %*% projMatrix
  colnames(bgComponents) <- compNames
  
  # Return the model
  CoCoModel <- list(
    projMatrix = projMatrix,
    fgComponents = fgComponents,
    bgComponents = bgComponents
  )
  
  return(CoCoModel)
}
