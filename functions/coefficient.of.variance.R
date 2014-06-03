# calculates the coefficient of variance
coefficient.of.variance <- function(x)
{
	# coefficient.of.variance is a simple function that returns the coefficient of variance
	# for a variable X (e.g. a data frame)

	sd(x) / abs(mean(x))
}
