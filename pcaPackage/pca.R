#' @title Principal Component Analysis
#' @description This allows for the exploration of the Principle component analysis of mass spectometry proteomics data
#' @param Data Load in data set with data frame 
#' @param numerics Filter our columns with no numerical values 
#' @param NA's  replace Inf/-Inf with NA
#' @param filtering remove columns that are all NA and omitting columns with 0
#' @param checks Check to make sure remaining data is usable 
#' @param Transformation transforms the data set 
#' @param covariance generates matrixs 
#' @param Eigen eigen vectors and eigen values
#' @keywords PCA
#' @export
#' @examples 
#' Principal component analysis (PCA) is used to reduce the dimensionality of high-throughput proteomics data by summarizing variation in protein abundance across conditions. This allows for easier identification of

pca <- function(data, eigen_opt = NA, threshold = 85) { # conducts principal component analysis
  data <- as.data.frame(data)

  # keep only numeric columns
  data <- data[, sapply(data, is.numeric), drop = FALSE]

  # replace Inf/-Inf with NA
  data[] <- lapply(data, function(x) {
    x[is.infinite(x)] <- NA
    x
  })

  # remove columns that are all NA
  data <- data[, colSums(!is.na(data)) > 0, drop = FALSE]

  # remove columns with sd = 0 or non-finite sd
  data <- data[, sapply(data, function(x) {
    s <- sd(x, na.rm = TRUE)
    is.finite(s) && s > 0
  }), drop = FALSE]

  # remove rows with missing values
  data <- na.omit(data)

  # safety checks
  if (ncol(data) < 2) {
    stop("Not enough usable numeric columns remain after cleaning.")
  }

  if (nrow(data) < 2) {
    stop("Not enough usable rows remain after cleaning.")
  }

  std <- data[0, , drop = FALSE] # exclude non-numeric properties from analysis
  for(i in colnames(std)) { # standardized values for each factor
    std[1, i] <- mean(as.numeric(data[, i]), na.rm = TRUE)
    std[2, i] <- sd(as.numeric(data[, i]), na.rm = TRUE)
  }
  rownames(std) <- c("mean", "sd")

  trns <- data # transform (AKA standardize) data so that each factor has a mean of 0 and a variance of 1
  for(i in colnames(std)){
    trns[, i] <- (data[, i] - std[1, i]) / std[2, i]
  }

  cov_mat <- data[0, , drop = FALSE] # generate a covariance matrix for each factor
  for(i in colnames(cov_mat)) { # i is column and j is row, so should go [j, i]
    for(j in colnames(cov_mat)) {
      cov_mat[j, i] <- stats::cov(trns[, j], trns[, i])
    }
  }

  eigen <- eigen(data.matrix(cov_mat)) # generate eigenvectors and eigenvalues from the covariance matrix

  if(is.na(eigen_opt)) { # if the desired number of principal components is not specified, determine from threshold
    eigen_opt <- scree_plot(eigen$values, threshold = threshold)$optimal_number_of_components
  }

  pca <- as.matrix(trns) %*% eigen$vectors[, 1:eigen_opt, drop = FALSE] # conduct matrix multiplication

  columns <- c()
  for(i in 1:ncol(pca)) { # name columns to their appropriate principal components
    columns <- append(columns, paste0("PC", i))
  }
  colnames(pca) <- columns

  output <- list(pca, eigen_opt, eigen$values, eigen$vectors, as.matrix(cov_mat), as.matrix(trns), as.matrix(std))
  names(output) <- c("pca", "optimal_number_of_principal_components", "eigen_values", "eigen_vectors", "covariance_matrix", "transformed_data", "standardized_values")
  return(output)
}
