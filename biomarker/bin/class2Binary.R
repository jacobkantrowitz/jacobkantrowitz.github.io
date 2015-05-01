# Author: Joe Perez-Rogers
# Date: 2014-11-28
# Purpose: Convert a two level factor or string into a binary (0/1) factor after defining the case/control levels

class2Binary <- function(x,case,control){
  to.case <- x==case
  to.control <- x==control
  x[to.case] <- 1
  x[to.control] <- 0
  return(as.factor(x))
}