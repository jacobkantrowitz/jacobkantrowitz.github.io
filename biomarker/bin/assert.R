# Author: Joseph Perez-Rogers
# Date: 2014-06-12
# Purpose: Function to enforce the validity of expressions and return error messages
# Obtained from: http://stackoverflow.com/questions/8343509/better-error-message-for-stopifnot
# Input: An evaluatable expression and an error message
# Output: Error if expr is FALSE
# Usage: assert(expression,"error message")

assert <- function (expr, error) {
  if (! expr) stop(error, call. = FALSE)
}