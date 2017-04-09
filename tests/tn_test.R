rm(list=ls())
require(tnorm)
require(microbenchmark)
x <- seq(-50, 10, length.out=1e2)
mean <- 0
sd <- 1
lower <- -Inf
upper <- Inf

dat1 <- ptn(x, mean, sd, lower, upper)
dat2 <- msm::ptnorm(x, mean, sd, lower, upper)

plot(x, log(dat1[,1]))
lines(x, log(dat2), col="red", lwd=2)
mtext("pnorm(x, log=TRUE)", adj = 0)
mtext("log(pnorm(x))", col = "red", adj = 1)

all.equal(dat1[,1], dat2)
all(dat1[,1] == dat2)

res <- microbenchmark(
    ptn(x, mean, sd, lower, upper),
    msm::ptnorm(x, mean, sd, lower, upper))
print(res)


x <- seq(-50, 10, length.out=1e3)
mean <- 0
sd <- 1
lower <- 0
upper <- 5

dat1 <- ptn(x, mean, sd, lower, upper, lp=TRUE)
dat2 <- msm::ptnorm(x, mean, sd, lower, upper, log.p=TRUE)
all.equal(dat1[,1], dat2)
all(dat1[,1] == dat2)

dat1
dat2
plot(x, log(dat1[,1]))
lines(x, log(dat2), col="red", lwd=2)
mtext("pnorm(x, log=TRUE)", adj = 0)
mtext("log(pnorm(x))", col = "red", adj = 1)

res <- microbenchmark(
    ptn(x, mean, sd, lower, upper, lp=T),
    msm::ptnorm(x, mean, sd, lower, upper, log.p=T))
print(res)


x <- seq(-5, 5, length.out=1e3)
mean <- 0
sd <- 1
lower <- -Inf
upper <- Inf

dat1 <- dtn(x, mean, sd, lower, upper, 0)
dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
plot(x, dat1, type="l", lwd=2)
lines(x, dat2, lwd=2, lty="dashed", col="red")

res <- microbenchmark(
    dtn(x, mean, sd, lower, upper, 0),
    msm::dtnorm(x, mean, sd, lower, upper, 0))
print(res)
boxplot(res)
ggplot2::autoplot(res)
## Unit: microseconds
##                                       expr    min      lq     mean median
##          dtn(x, mean, sd, lower, upper, 0) 43.163 43.9665 46.33943 46.061
##  msm::dtnorm(x, mean, sd, lower, upper, 0) 76.268 78.6080 96.26420 81.296

x <- seq(-5, 5, length.out=1e3)
mean <- 0
sd <- 1
lower <- 0
upper <- Inf

dat1 <- dtn(x, mean, sd, lower, upper, 0)
dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
plot(x, dat1, type="l", lwd=2)
lines(x, dat2, lwd=2, lty="dashed", col="red")

res <- microbenchmark(
    dtn(x, mean, sd, lower, upper, 0),
    msm::dtnorm(x, mean, sd, lower, upper, 0))
print(res)

x <- seq(-5, 5, length.out=1e3)
mean <- 0
sd <- 1
lower <- -Inf
upper <- 2

dat1 <- dtn(x, mean, sd, lower, upper, 0)
dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
plot(x, dat1, type="l", lwd=2)
lines(x, dat2, lwd=2, lty="dashed", col="red")

res <- microbenchmark(
    dtn(x, mean, sd, lower, upper, 0),
    msm::dtnorm(x, mean, sd, lower, upper, 0))
print(res)

x <- seq(-5, 5, length.out=1e3)
mean <- 0
sd <- 1
lower <- -2
upper <- 2

dat1 <- dtn(x, mean, sd, lower, upper, 0)
dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
plot(x, dat1, type="l", lwd=2)
lines(x, dat2, lwd=2, lty="dashed", col="red")

res <- microbenchmark(
    dtn(x, mean, sd, lower, upper, 0),
    msm::dtnorm(x, mean, sd, lower, upper, 0))
print(res)


n <- 1e5
mean <- 0
sd <- 1
lower <- -Inf
upper <- Inf

dat1 <- rtn(n, mean, sd, lower, upper)
dat2 <- msm::rtnorm(n, mean, sd, lower, upper)
den2 <- density(dat2)
hist(dat1, breaks="fd", freq=F)
lines(den2$x, den2$y, lwd=2.5)

res <- microbenchmark(
    rtn(n, mean, sd, lower, upper),
    msm::rtnorm(n, mean, sd, lower, upper))

print(res)
boxplot(res)
ggplot2::autoplot(res)


n <- 1e5
mean <- 0
sd <- 1
lower <- 0
upper <- Inf

dat1 <- rtn(n, mean, sd, lower, upper)
dat2 <- msm::rtnorm(n, mean, sd, lower, upper)
den2 <- density(dat2)
hist(dat1, breaks="fd", freq=F)
lines(den2$x, den2$y, lwd=2.5)

res <- microbenchmark(
    rtn(n, mean, sd, lower, upper),
    msm::rtnorm(n, mean, sd, lower, upper))
print(res)


n <- 1e5
mean <- 0
sd <- 1
lower <- 0
upper <- 3

dat1 <- rtn(n, mean, sd, lower, upper)
dat2 <- msm::rtnorm(n, mean, sd, lower, upper)
den2 <- density(dat2)
hist(dat1, breaks="fd", freq=F)
lines(den2$x, den2$y, lwd=2.5)

res <- microbenchmark(
    rtn(n, mean, sd, lower, upper),
    msm::rtnorm(n, mean, sd, lower, upper))
print(res)

