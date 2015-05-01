# Author: Adam Gower

fastLM <- function (dependent, mod) {
  # Function to quickly compute a set of linear models that all have the same right-hand side (rhs)
  # (essentially a wrapper for the lmFit() function in the limma package)
  # Adam Gower, 2009
  #
  # INPUT
  # dependent         a matrix that contains instances of the dependent variable (lhs of the linear model),
  #                   (e.g., genes), one in each row, with observations (e.g., samples) by columns
  # independents      a matrix or data frame that contains the independent variables (rhs of the linear model),
  #                   (e.g., demographic covariates), one in each column, with observations (e.g., samples) by rows
  #
  # OUTPUT
  # A list with five elements:
  # coefficients      a matrix of the linear model coefficients, one dependent variable instance per row
  # stdev.unscaled    a matrix of the unscaled standard deviation of the linear model coefficients,
  #                   one dependent variable instance per row
  # sigma             a vector of the residual standard deviation for each dependent variable instance
  #                   NOTE: standard error = stdev.unscaled * sigma
  # t                 a matrix of t statistics of the linear model coefficients, one dependent variable instance per row
  # p                 a matrix of p values of the linear model coefficients, one dependent variable instance per row
  
  require(limma);
  
  model <- lmFit(dependent,mod);
  
  result <- list();
  result$coefficients <- model$coefficients;
  result$stdev.unscaled <- model$stdev.unscaled;
  result$sigma <- model$sigma;
  # t statistic = coefficient / standard error = coefficient / (unscaled standard deviation * sigma)
  result$t <- model$coefficients / (model$stdev.unscaled * model$sigma);
  # p value is computed using degrees of freedom based on number of observations and number of independent variables
  result$p <- 2 * pt(-abs(result$t), df=ncol(dependent)-ncol(mod)-1);
  
  return(result);
}