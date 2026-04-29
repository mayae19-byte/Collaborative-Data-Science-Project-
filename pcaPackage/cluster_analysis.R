#' @title k-means Clustering Analysis of Principal Components
#' @description This function performs k-means clustering on PCA results to group similar observations. It also determines the optimal number of clusters using the elbow method and generates plots.
#' @param pca matrix of principal component values for each observation
#' @param num_clusters number of clusters to use for analysis; if left blank, the function will use the optimal number of clusters determined by the elbow method
#' @param max_clusters maximum number of clusters to analyze for determining the optimal number of clusters
#' @param x the principal component to be plotted on the x axis of the scatterplot
#' @param y the principal component to be plotted on the y axis of the scatterplot
#' @param nstart number of times to run k-means clustering
#' @export
#' @examples
#' cluster_analysis(pca(iris)$pca, num_clusers = NA, max_clusters = 10, x = "PC1", y = "PC2", nstart = 10)


cluster_analysis <- function(pca, num_clusters = NA, max_clusters = 10, x = "PC1", y = "PC2", nstart = 10) {
  
  # Compute TWSS
  twss <- c()
  for(i in 1:max_clusters) {
    twss <- append(twss, stats::kmeans(pca$pca, centers = i, nstart = nstart)$tot.withinss)
  }
  
  # Data for elbow plot
  clust_coords <- data.frame(k = 1:max_clusters, twss = twss)
  
  # Find optimal cluster number (elbow)
  clust_opt <- which.min(diff(twss)) + 1
  
  if(is.na(num_clusters)) {
    num_clusters <- clust_opt
  }
  
  # Final clustering
  kmeans_result <- stats::kmeans(pca$pca, centers = num_clusters, nstart = nstart)
  
  cluster_values <- kmeans_result$cluster
  qual <- 100 * kmeans_result$betweenss / kmeans_result$totss
  
  # Elbow plot
  cluster_determine <- ggplot2::ggplot(clust_coords, ggplot2::aes(x = k, y = twss)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = clust_opt, col = "red") +
    ggplot2::xlab("Number of Clusters") +
    ggplot2::ylab("Total Within Sum of Squares") +
    ggplot2::theme_minimal()
  
  # PCA scatter with clusters
  cluster_scatter <- ggplot2::ggplot(
    data = as.data.frame(pca$pca),
    ggplot2::aes(x = .data[[x]], y = .data[[y]], color = as.factor(cluster_values))
  ) +
    ggplot2::geom_point() +
    ggplot2::xlab(x) +
    ggplot2::ylab(y) +
    ggplot2::theme_minimal()
  
  # Output
  output <- list(
    cluster_values = cluster_values,
    optimal_number_of_clusters = clust_opt,
    quality = qual,
    optimal_cluster_plot = cluster_determine,
    pca_scatterplot = cluster_scatter
  )
  
  return(output)
}