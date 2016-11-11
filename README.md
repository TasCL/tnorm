# Truncated Normal Distribution in C++
Implement dtnorm and rtnorm functions use C++ codes.


## Getting Started

```
require(tnorm)

## Plot regular normal and trunicated normal distributions
plot(function(x) dnorm(x, log = FALSE), -2.5, 2.5,
     main = "Normal Distribution", ylim=c(0,0.45), ylab="Density")
curve(dtnorm(x, lower=-2, upper=2), add=TRUE, col="tomato", lwd=2)
mtext("dnorm(x)", adj = 0)
mtext("dtnorm(x)", col = "tomato", adj = 1)

```

## Installation 

```

## From github
devtools::install_github("TasCL/tnorm")
## From source: 
install.packages("tnormc_0.1.0.1.tar.gz", repos = NULL, type="source")

```

## Prerequisities
 - R (>= 3.0.2)
 - Rtools
