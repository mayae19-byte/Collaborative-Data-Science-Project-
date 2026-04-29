#A scree plot displays the contribution of each principal component to the total variance in the dataset. By examining the point where the curve begins to level off, users can determine the number of components that capture most of the meaningful variation.
#' @title Scree Plot
#' @description Generates a scree plot, determines the optimal number of principal components, and calculates the total proportion of variance explained
#' @param eigen_values Matrix of eigenvalues
#' @param threshold percent of variance the minimum number of principal components must explain
#' @param y_axis y axis for the scree plot; use "Eigenvalue" for eigenvalues to be plotted and "Percentage" for percentage of explained variance
#' @keywords pca principal component analysis eigen eigenvalue eigenvalues
#' @export
#' @examples
#' scree_plot(pca(iris)$eigen_values, threshold = 85, y_axis = "Eigenvalue")

#A scree plot displays the contribution of each principal component to the total variance in the dataset. By examining the point where the curve begins to level off, users can determine the number of components that capture most of the meaningful variation.

scree_plot <- function(eigen_values, threshold = 85, y_axis = "Eigenvalue") { # generates a scree plot and determines how many principal components are optimal for analysis
  
  # Initialize number of components and explained variance 
  num_comp <- 0
  exp_var <- 0
  while(exp_var < threshold) { # determines the minimum number of principal components that explains the minimal threshold percent of variance (default is 85 percent)
    num_comp <- num_comp + 1
    exp_var <- exp_var + 100 * eigen_values[num_comp]/sum(eigen_values)
  }
  
  # Calculate percentage of variance explained by each component
  percentage <- t(as.matrix(100 * eigen_values/sum(eigen_values))) # determines how much of the dataset's variance can be explained by the given eigenvalues
  columns <- c()
  for(i in 1:ncol(percentage)) { # name columns to their appropriate principal components
    columns <- append(columns, paste0("PC", i))
  }
  
  # Create column names (PC1, PC2, etc.)
  colnames(percentage) <- columns
  if(y_axis == "Percentage") { # y axis is the actual eigenvalue of each number of components, but can be modified to show percentage of explained variance instead
    y_axis <- "Percentage of Variance Explained"
    eigen_values <- percentage[1,]
  }
  
  # Create data frame for plotting
  value_coords <- data.frame(1:length(eigen_values), eigen_values)
  
  # Generate scree plot
  plot <- ggplot2::ggplot(data = value_coords, aes(x = value_coords[, 1], y = value_coords[, 2])) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::xlab("Principal Component") +
    ggplot2::ylab(y_axis) +
    ggplot2::geom_vline(xintercept = num_comp, col = "red") +
    ggplot2::theme(panel.border = element_rect(colour = "black"))
  
  # Store outputs in a list
  output <- list(plot, num_comp, percentage)
  
  # Name output elements
  names(output) <- c("plot", "optimal_number_of_components", "percentage_of_variance_explained")
  return(output)
}
