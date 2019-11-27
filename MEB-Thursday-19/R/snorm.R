# Computes the density of a skew-normal distribution. 
#
# Args:
#   x: vector of quantiles.
#   mu: location parameter(s). 
#   sd: scale parameter(s). 
#   delta: skew parameter(s).
#   log: logical for whether the log density should be returned.
#
# Returns:
#   log density or density. 
dsnorm <- function(x, mu = 0, sd = 1, delta = 0, log = FALSE){
  x <- x - mu
  is_log <- log
  log <- base::log
  
  if(length(is_log) == 1 && length(x) > 1)
    is_log <- rep(is_log, length(x))
  
  log_d <- log(2) + dnorm(x, sd = sd, log = TRUE) + 
    pnorm(delta * x, log.p = TRUE)
  
  ifelse(is_log, log_d, exp(log_d))
}

# Computes the density of a multivariate skew-normal distribution. 
#
# Args:
#   x: matrix of quantiles. Each rows is a quantile.
#   mu: vector with location parameters.
#   Sig: matrix with scale parameters.
#   delta: vector with shape/skewness parameters.
dmsnorm <- function(x, mu, Sig, delta, log = FALSE){
  if(!require(mvtnorm))
    stop(sprintf("Requires %s", sQuote("mvtnorm")))
  x <- t(t(x) - mu)
  is_log <- log
  log <- base::log
  
  log_d <- log(2) + dmvnorm(x, sigma = Sig, log = TRUE) + 
    pnorm(drop(x %*% delta), log.p = TRUE)
  
  if(is_log) log_d else exp(log_d)
}

# tests
local({
  # log and non-log gives consistent results
  test_vals <- expand.grid(
    x = seq(-4, 4), mu = seq(-4, 4), sd = 2^(0:4), delta = seq(-2, 2))
  
  with(test_vals, stopifnot(isTRUE(all.equal(
        dsnorm(x = x, mu = mu, sd = sd, delta = delta), 
    exp(dsnorm(x = x, mu = mu, sd = sd, delta = delta, log = TRUE))
  ))))
  
  # integrate to one
  test_vals <- expand.grid(
    mu = seq(-4, 4), sd = 2^(0:4), delta = seq(-2, 2))
  consts <- with(test_vals, mapply(function(mu, sd, delta){
    integrate(dsnorm, -Inf, Inf, mu = mu, sd = sd, delta = delta, 
              rel.tol = 1e-10)$value
  }, mu = mu, sd = sd, delta = delta))
  
  stopifnot(isTRUE(all.equal(rep(1, length(consts)), consts)))
  
  xy <- as.matrix(expand.grid(-2:2, -2:2))
  Sig <- matrix(c(1, 1, 1, 3), 2)
  rho <- c(-1, 1)
  mu  <- c(1, 0) 
  d1 <- dmsnorm(xy, mu, Sig, rho, log = FALSE)
  d2 <- dmsnorm(xy, mu, Sig, rho, log = TRUE)
  stopifnot(isTRUE(all.equal(d1, exp(d2))))
  
  const <- integrate(function(x)
    sapply(x, function(xx)
      integrate(function(y) dmsnorm(cbind(xx, y), mu, Sig, rho, FALSE), 
                    -100, 100)$value), -100, 100) 
  
  stopifnot(isTRUE(all.equal(const$value, 1, 
                             tolerance = .Machine$double.eps^(1/5))))

})
