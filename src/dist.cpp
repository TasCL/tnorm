#include <tnorm.hpp>

#define SQRT_2PI   2.5066282746310007e+0 /* sqrt(2 x pi) */
#define M_E	   2.7182818284590452354	/* e */

// Accept-Reject Algorithm 0; Naive method A-R method
inline double rtnorm0(const double lower, const double upper)
{
  int valid = 0;
  double z  = 0.0;

  while (valid == 0) {
    z = R::rnorm(0.0, 1.0);
    if (z <= upper && z >= lower) { valid = 1; }
  }
  return z ;
}

// Algorithm 1; 'expl'; use when lower > mean; upper = INFINITY; p 122, right
inline double rtnorm1(const double lower, const double upper)
{
  int valid = 0;
  double z = 0, u = 0, r = 0; // a stands for alphaStar in Robert (1995)
  double a = 0.5 * (sqrt(lower * lower + 4.0) + lower);

  while (valid == 0) {
    z = R::rexp(a) + lower; // control lower boundary
    u = R::runif(0, 1);
    r = exp(-0.5 * (z - a) * (z - a));
    if (u <= r && z <= upper) { valid = 1; }
  }
  return z ;
}

// Algorithm 2; 'expu'; use when upper < mean; lower = -INFINITY.
inline double rtnorm2(const double lower, const double upper)
{
  int valid  = 0;
  double z = 0, u = 0, r = 0; // a stands for alphaStar in Robert (1995)
  double a = 0.5 * (sqrt(upper * upper + 4.0) - upper);

  while (valid == 0) {
    z = R::rexp(a) - upper;
    u = R::runif(0, 1);
    r = exp(-0.5 * (z - a) * (z - a));
    if (u <= r && z <= -lower) { valid = 1 ; }
  }
  return -z;  // note the negative
}

// Algorithm 3; u;  page 123. 2.2. Two-sided truncated normal dist.
inline double rtnorm3(const double lower, const double upper)
{
  int valid  = 0;
  double z = 0, u = 0, r = 0; // a stands for alphaStar in Robert (1995)

  while (valid == 0) {
    z = R::runif(lower, upper) ;

    if (lower > 0) {
      r = exp( 0.5 * (lower*lower - z*z) );
    } else if (upper < 0) {
      r = exp( 0.5 * (upper*upper - z*z) );
    } else  {
      r = exp( -0.5 * z * z ) ;
    }

    u = R::runif(0, 1) ;
    if (u <= r) { valid = 1 ; }
  }
  return z ;
}

inline double rtn_scalar(const double mean,  const double sd, const double lower,
  const double upper)
{
  double z, l, u, eq_a1, eq_a2; // Standardised lower and upper
  bool a0, a1, a2, a3;
  l = (lower - mean) / sd; // l == stdlower, u == stdupper
  u = (upper - mean) / sd;

  // Accept-Reject Algorithm 0;
  a0 = (l < 0 && u == INFINITY) || (l == -INFINITY && u > 0) ||
    (std::isfinite(l) && std::isfinite(u) && l < 0 && u > 0 && (u - l) > SQRT_2PI);

  // Algorithm (1): Use Proposition 2.3 with only lower truncation. upper==INFINITY
  // rejection sampling with exponential proposal. Use if lower > mean
  eq_a1 = l + (2.0 * std::sqrt(M_E) / (l + std::sqrt(l * l + 4.0))) *
    (std::exp( 0.25 * (2.0 * l - l * std::sqrt(l * l + 4.0))));
  a1 = (l >= 0) && (u > eq_a1);

  // Algorithm (2): Use -x ~ N_+ (-mu, -mu^+, sigma^2) on page 123. lower==-INFINITY
  // rejection sampling with exponential proposal. Use if upper < mean.
  eq_a2 = -u + (2.0 * std::sqrt(M_E) / (-u + std::sqrt(u*u + 4.0))) *
    (std::exp( 0.25 * (2.0 * u + u * std::sqrt(u*u + 4.0))));
  a2 = (u <= 0) && (-l > eq_a2);

  if (a0) {
    z = rtnorm0(l, u);
  } else if (a1) {
    z = rtnorm1(l, u);
  } else if (a2) {
    z = rtnorm2(l, u);
  } else {              // rejection sampling with uniform proposal.
    z = rtnorm3(l, u);  // Use if bounds are narrow and central.
  }
  return z*sd + mean;
}

inline double dtn_scalar(const double x, const double mean, const double sd,
  const double lower, const double upper, int lp)
{
  double out, numer, denom;
  if ((x >= lower) && (x <= upper)) {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0);
    numer = R::dnorm(x, mean, sd, lp);
    out = lp ? (numer - std::log(denom)) : (numer/denom);
  } else {
    out = lp ? -INFINITY : 0;
  }
  return(out);
}

double ptn_scalar(const double q, const double mean, const double sd,
                  const double lower, const double upper, int lt, int lp) {
    double out, numer, denom, qtmp;
    if (lt) {
        out = (q < lower) ? 0 : 1;
    } else {
        out = (q < lower) ? 1 : 0;
    }
    if ((q >= lower) && (q <= upper)) {
        denom = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0);
        qtmp  = lt ? (R::pnorm(q, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0)) :
            (R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(q, mean, sd, 1, 0));
        out  = lp ? (std::log(qtmp)-std::log(denom)) : (qtmp/denom);
    }
    return(out);
}


//' Truncated Normal Distribution
//'
//' Three functions implementing random number generator, density function and
//' cumulative distribution function (\code{rtn}, \code{dtn}, and \code{ptn})
//' for truncation normal distribution.
//'
//' @param x,n,q x in \code{dtn} is a vector of quantiles; n in \code{rtn}
//' is number of observations. n must be a scalar. q is a vector of quantiles.
//' @param mean mean (must be scalar).
//' @param sd standard deviation (must be scalar).
//' @param lower lower truncation value (must be scalar).
//' @param upper upper truncation value (must be scalar).
//' @param lt lower.tail. a boolean switch; if TRUE (default) probabilities are
//' \code{P[X <= x]}, otherwise, \code{P[X > x]}.
//' @param lp log.p. a boolean switch; if TRUE (default is FALSE) probabilities p
//' are given as \code{log(p)}.
//' @return a column vector.
//' @examples
//' ## rtn example
//' dat1 <- rtn(1e5, 0, 1, 0, Inf)
//' ## dat2 <- msm::rtnorm(n, mean, sd, lower, upper)
//' ## den2 <- density(dat2)
//' hist(dat1, breaks="fd", freq=F)
//' ## lines(den2$x, den2$y, lwd=2.5)
//' ## res <- microbenchmark(
//' ##     rtn(n, mean, sd, lower, upper),
//' ##     msm::rtnorm(n, mean, sd, lower, upper))
//' ## print(res)
//'
//' ## dtn example
//' x <- seq(-5, 5, length.out=1e3)
//' dat1 <- dtn(x, 0, 1, -2, 2, 0)
//' ## dat2 <- msm::dtnorm(x, mean, sd, lower, upper, 0)
//' plot(x, dat1, type="l", lwd=2)
//' ## lines(x, dat2, lwd=2, lty="dashed", col="red")
//'
//' ## res <- microbenchmark(
//' ##     dtn(x, mean, sd, lower, upper, 0),
//' ##     msm::dtnorm(x, mean, sd, lower, upper, 0))
//' ## print(res)
//'
//' ## ptn example
//' x <- seq(-50, 10, length.out=1e3)
//' mean <- 0
//' sd <- 1
//' lower <- 0
//' upper <- 5
//' dat1 <- ptn(x, 0, 1, 0, 5, lp=TRUE)
//' ## dat2 <- msm::ptnorm(x, mean, sd, lower, upper, log.p=TRUE)
//' ## all(dat1[,1] == dat2)
//'
//' plot(x, log(dat1[,1]))
//' ## lines(x, log(dat2), col="red", lwd=2)
//' ## mtext("pnorm(x, log=TRUE)", adj = 0)
//' ## mtext("log(pnorm(x))", col = "red", adj = 1)
//' @export
// [[Rcpp::export]]
arma::vec dtn(arma::vec x, double mean, double sd, double lower, double upper,
              int lp=false) {
    if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
    if (sd < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
    if (sd==R_NegInf   || sd==R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
    if (mean==R_NegInf || mean==R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}

    int idx;
    arma::vec out(x.n_elem);

    for(arma::vec::iterator i = x.begin(); i < x.end(); i++) {
        idx = std::distance(x.begin(), i);
        out[idx] = dtn_scalar(*i, mean, sd, lower, upper, lp);
    }
    return out;
}

//' @rdname dtn
//' @export
// [[Rcpp::export]]
arma::vec rtn(int n, double mean, double sd, double lower, double upper) {
    if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
    if (sd < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
    if (sd==R_NegInf   || sd==R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
    if (mean==R_NegInf || mean==R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}
    arma::vec out(n);
    for(arma::vec::iterator i = out.begin(); i < out.end(); i++) {
        *i = rtn_scalar(mean, sd, lower, upper);
    }
    return out;
}

//' @rdname dtn
//' @export
// [[Rcpp::export]]
arma::vec ptn(arma::vec q, double mean, double sd, double lower, double upper,
              int lt=true, int lp=false) {
    if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
    if (sd < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
    if (sd==R_NegInf   || sd==R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
    if (mean==R_NegInf || mean==R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}

    int idx;
    arma::vec out(q.n_elem);

    for(arma::vec::iterator i = q.begin(); i < q.end(); i++) {
        idx = std::distance(q.begin(), i);
        out[idx] = ptn_scalar(*i, mean, sd, lower, upper, lt, lp);
    }
    return out;

}
