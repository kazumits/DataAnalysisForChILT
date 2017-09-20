# geometric normalization between samples (columns)
gscale <- function(X)
{
  if(any(X<0)) stop("All values should be non-negative.")
  scale(X,
    scale = apply(X/apply(X, 1, function(x)
      exp(mean(log(x)))), 2, function(x)
        median(x,na.rm=TRUE)
      ),
    center=FALSE
  )
}