# Truncated Normal Distribution in C++
dtnorm and rtnorm functions in C++ codes.

## Getting Started

```
require(tnorm)

## Plot regular normal and trunicated normal distributions
x <- seq(-2.5, 2.5, length.out=1e2)
y <- dnorm(x)
plot(x, y, type="l", lty="dashed", main = "Normal Distribution", ylim=c(0,0.45), ylab="Density")
lines(x, dtn(x, 0, 1, -2, 2), col="tomato", lwd=2)
mtext("dnorm(x)", adj = 0)
mtext("dtnorm(x)", col="tomato", adj = 1)

```

## Installation 

```
## From github
devtools::install_github("TasCL/tnorm")
## From source: 
install.packages("tnormc_0.2.0.0.tar.gz", repos = NULL, type="source")

```

## Prerequisities

 - R (>= 3.0.2)
 - Rtools
 - Rcpp package

## References

 - Robert, C. P. (1995). Simulation of truncated normal variables, Satatistics and 
Computing, 5, 121--125. http://dx.doi.org/10.1007/BF00143942
 - RcppTN 0.1-8 (https://github.com/olmjo/RcppTN) 
 - msm package (https://cran.r-project.org/web/packages/msm/index.html). 
