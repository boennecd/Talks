# Simulate potentially right censored outcomes from a generalized survival
# model with a given link function and a mixed weibull distribution as the 
# baseline. 
#
# Args:
#   b: coefficients for the fixed effects. 
#   X: design matrix for the fixed effects.
#   Z: design matrix for the random effects. 
#   cvmat: covariance matrix for the random effects. 
#   grp: integer or factor vector with group identifiers. 
#   ps: mixture probabilities. 
#   gs: mixture shape parameters. 
#   ls: mixture scale parameters.
#   cens_func: function to simulate censoring time given integer with number
#              of observations.
#   link: character specifying the link function.
# 
# Returns: 
#   The observed observation times, the event indicators, and the drawn  
#   random effects.
sim_wei_mix_base <- function(
  b, X, Z, cvmat, grp, ps = c(.5, .5), gs = c(1.5, .5), ls = c(1, 1), 
  cens_func = function(n_obs) runif(n_obs, 1, 5), 
  link = c("PH", "PO", "probit")){
  n_obs <- NROW(X) 
  n_fix <- NCOL(X)
  n_rng <- NCOL(Z)
  n_mix <- length(gs)
  link <- link[1]
  
  # check args
  chol_cv <- chol(cvmat)
  stopifnot(
    length(ps) == n_mix, all(ps > 0),
                         all(gs > 0), 
    length(ls) == n_mix, all(ls > 0), 
    length(b) == n_fix, is.numeric(b),
    all(dim(as.matrix(cvmat)) == c(n_rng, n_rng)), 
    length(grp) == n_obs, is.integer(grp) || is.factor(grp),
    NROW(Z) == n_obs    , is.numeric(Z), 
                          is.numeric(X), 
    is.character(link), link %in% c("PH", "PO", "probit"))
  
  # make transformations and assign variables used later
  if(is.factor(grp))
    grp <- as.integer(grp)
  n_grp <- length(unique(grp))
  
  # assign survial function 
  ps <- ps / sum(ps)
  S0 <- function(ti)
    colSums(ps * exp(-ls * t(outer(ti, gs, `^`))))
  
  # simulate linear predictor
  rngs <- matrix(rnorm(n_grp * n_rng), nc = n_rng) %*% chol_cv
  colnames(rngs) <- colnames(Z)
  lp <- drop(X %*% b + rowSums(Z * rngs[grp, ]))
  
  # simulate survival time and censoring time
  # TODO: can we just simulate from a mixture?
  unifs <- runif(n_obs, 0, 1)
  root_func <- switch(
    link,
    PH = function(ti, u, lp) S0(ti) - u^exp(-lp), 
    PO  = function(ti, u, lp) {
      sb <- S0(ti)
      f  <- exp(lp / 2)
      (1 - sb) / (sb + .Machine$double.eps) * f - (1 - u) / u / f
    }, 
    probit = function(ti, u, lp){
      sb <- S0(ti)
      eta_min <- qnorm(.Machine$double.eps)
      trunc <- function(eta)
        pmax(eta_min, pmin(-eta_min, eta))
      lp_2 <- lp / 2
      trunc(-qnorm(sb) + lp_2) + trunc(lp_2 + qnorm(u))
    })
  y <- mapply(
    uniroot, MoreArgs = list(
      f = root_func, interval = c(0, .Machine$double.xmax^(1/16)), 
      tol = .Machine$double.eps^(3/4)), u = unifs, lp = lp)
  y <- unlist(y["root", ])
  y <- pmax(.Machine$double.eps^(1/2), y)
  
  cens <- cens_func(n_obs)
  stopifnot(length(cens) == n_obs, all(cens > 0))
  event <- y < cens
  
  y <- pmin(y, cens)
  
  list(y = y, event = event, rngs = rngs)  
}

# test the function
if(TRUE)
  (function(){
    # reset seed after test
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    old_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- old_seed)
    
    #####
    # simulate data
    get_dat <- function(link){
      dat <- within(list(), {
        n_p_grp <- 2L
        n_grp <- 5L
        n_obs <- n_p_grp * n_grp
        grp <- as.integer(gl(n_grp, n_p_grp))
        
        X <- cbind(`(Intercept)` = 1, treatment = runif(n_obs) > .5, 
                   x = rnorm(n_grp)[grp])
        Z <- cbind(`(Intercept)` = 1, x = X[, "x"])
        beta <- c(.5, 1/4, 1/4)
        
        stds <- c(.5, 1/4)
        rho <- .5 * prod(stds)
        cvmat <- matrix(c(stds[1]^2, rho, rho, stds[2]^2), nc = 2)
      })
      
      with(dat, sim_wei_mix_base(
        b = beta, X = X, Z = Z, cvmat = cvmat, grp = grp))
    }
    
    set.seed(1)
    sim_dat <- get_dat("PH")
    
    stopifnot(isTRUE(all.equal(
      sim_dat, 
      list(y = c(
        0.0261795159458619, 0.113360951088831, 0.0277792241363675, 
        0.138872464541846, 0.369658829436196, 0.0531238909295669, 2.265086828731, 
        0.969466868448416, 0.0391259321937663, 0.0843747507945069), 
        event = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), 
        rngs = structure(c(
          0.755890584225424, 0.194921618205716, -0.310620290270902, -1.10734994358875, 0.562465459071554, 
          0.179244234333629, 0.0452251097670186, 0.126691461298552, -0.0990378816268971, 
          0.269199772646765),  .Dim = c(5L, 2L), 
          .Dimnames = list(NULL, c("(Intercept)", "x")))))))
    
    set.seed(1)
    sim_dat <- get_dat("PO")
    
    stopifnot(isTRUE(all.equal(
      sim_dat, 
      list(
        y = c(0.0261795159458619, 0.113360951088831, 0.0277792241363675, 
              0.138872464541846, 0.369658829436196, 0.0531238909295669, 2.265086828731, 
              0.969466868448416, 0.0391259321937663, 0.0843747507945069), 
        event = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), 
        rngs = structure(
          c(0.755890584225424, 0.194921618205716, -0.310620290270902, -1.10734994358875, 0.562465459071554, 
            0.179244234333629, 0.0452251097670186, 0.126691461298552, -0.0990378816268971, 0.269199772646765), 
          .Dim = c(5L, 2L), .Dimnames = list(NULL, c("(Intercept)", "x")))))))
    
    set.seed(1)
    sim_dat <- get_dat("probit")
    
    stopifnot(isTRUE(all.equal(
      sim_dat, 
      list(
        y = c(0.0261795159458619, 0.113360951088831, 0.0277792241363675, 
              0.138872464541846, 0.369658829436196, 0.0531238909295669, 2.265086828731, 
              0.969466868448416, 0.0391259321937663, 0.0843747507945069), 
        event = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), 
        rngs = structure(
          c(0.755890584225424, 0.194921618205716, -0.310620290270902, -1.10734994358875, 0.562465459071554, 
            0.179244234333629, 0.0452251097670186, 0.126691461298552, -0.0990378816268971, 
            0.269199772646765), 
          .Dim = c(5L, 2L), .Dimnames = list(NULL, c("(Intercept)", "x")))))))
  })()

invisible()