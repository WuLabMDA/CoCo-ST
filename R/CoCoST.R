#' CoCoST: Compare-and-Contrast Spatial Transcriptomics
#'
#' CoCoST identifies spatial structures by comparing two datasets
#' (foreground and background) using normalized graph Laplacians and
#' contrastive component analysis.
#'
#' @param data1 Foreground data matrix (e.g., tumor sample).
#' @param W1 Affinity matrix for the foreground data.
#' @param data2 Background data matrix (e.g., normal sample).
#' @param W2 Affinity matrix for the background data.
#' @param para Contrastive parameter controlling the subtraction weight (default = 0.2).
#' @param Dim Number of components to extract (default = 2).
#'
#' @return A list containing:
#' \item{projMatrix}{Projection matrix.}
#' \item{fgComponents}{Foreground contrastive components.}
#' \item{bgComponents}{Background contrastive components.}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' data1 <- matrix(rnorm(1000), 100, 10)
#' data2 <- matrix(rnorm(1000), 100, 10)
#' library(kernlab)
#' rbf <- laplacedot(sigma = 0.5)
#' W1 <- kernelMatrix(rbf, t(data1))@.Data
#' W2 <- kernelMatrix(rbf, t(data2))@.Data
#' result <- CoCoST(t(data1), W1, t(data2), W2)
#' }
#'
#' @export
CoCoST <- function(data1, W1, data2, W2, para = 0.2, Dim = 2) {
  if (missing(data1)) stop("Foreground data not provided, please provide data")
  if (missing(data2)) stop("Needs to provide the background data")
  if (missing(W1)) stop("Affinity matrix for foreground data not provided!")
  if (missing(W2)) stop("Affinity matrix for background data not provided!")
  
  message(sprintf("Extracting %d contrastive components with parameter: %f", Dim, para))
  
  # Dependencies
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Please install 'Matrix'")
  if (!requireNamespace("RSpectra", quietly = TRUE)) stop("Please install 'RSpectra'")
  library(Matrix)
  library(RSpectra)
  
  # Sparse affinity matrices
  W1 <- as(W1, "sparseMatrix")
  W2 <- as(W2, "sparseMatrix")
  
  # Normalized Laplacians
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
  
  # Covariances and projection
  cov1 <- crossprod(data1, L11 %*% data1)
  cov2 <- crossprod(data2, L22 %*% data2)
  
  RR <- cov1 - para * cov2
  eig <- eigs_sym(RR, k = Dim)
  
  projMatrix <- eig$vectors
  compNames <- paste("CoCoST", 1:Dim, sep = " ")
  colnames(projMatrix) <- compNames
  
  fgComponents <- data1 %*% projMatrix
  bgComponents <- data2 %*% projMatrix
  colnames(fgComponents) <- colnames(bgComponents) <- compNames
  
  CoCoModel <- list(
    projMatrix = projMatrix,
    fgComponents = fgComponents,
    bgComponents = bgComponents
  )
  
  return(CoCoModel)
}
