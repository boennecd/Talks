#####
# build and load the dll/dso
.tmb_lib_f <- (function(){
  if(!require(TMB))
    stop(sprintf("Requires %s", sQuote("TMB")))
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  
  # set working dir to the source dir
  f_sep <- .Platform$file.sep
  root_dir_regx <- paste0("^.+", f_sep, "MEB-Thursday-19")
  stopifnot(grepl(root_dir_regx, getwd()))
  root_dir <- gsub(sprintf("(%s)(.+$)", root_dir_regx), "\\1", getwd())
  setwd(paste0(root_dir, f_sep, "src"))
  
  # build cpp file
  f <- "implementation.cpp"
  f_strip <- gsub("(.+)(\\.cpp$)", "\\1", f)
  stopifnot(file.exists(f))
  lib_file_path <- file.path(getwd(), f_strip)
  lib_file <- dynlib(lib_file_path)
  if(lib_file %in% unlist(sapply(getLoadedDLLs(), `[`, "path")))
    # notice the comments in `help('dyn.unload')` though...
    # In particular, see 
    #   https://github.com/kaskr/adcomp/issues/27#issuecomment-161990918
    dyn.unload(lib_file)
  
  stopifnot(
    compile(f, PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)") == 0)
  dyn.load(lib_file)

  f_strip
})()

.laplace_char <- "Laplace"
.gva_char     <- "GVA"
.snva_char    <- "SNVA"

.alpha_to_gamma <- function(alpha){
  mu <- sqrt(2 / pi) * alpha / sqrt(1 + alpha^2)
  (4 - pi) / 2 * mu^3 / (1 - mu^2)^(3/2)
}
.gamma_to_alpha <- function(gamma){
  cf <- 2 * gamma / (4 - pi)
  v1 <- local({
    cf <- ( cf)^(1/3)
    cf / sqrt(1 + cf^(2))
  })
  v2 <- local({
    cf <- (-cf)^(1/3)
    -cf / sqrt(1 + cf^(2))
  })
  mu <- ifelse(gamma > 0, v1, v2)
  o <- sqrt(pi) * mu 
  o / sqrt(2 - pi * mu^2)
}

#####
# Define function to get objective function.
# 
# Args:
#   formula: two-sided formula where the left-handside is a `Surv` object  
#            and the right-handside is the fixed effects. 
#   data: `data.frame` with variables used in the model. 
#   df: integer scalar with the degrees of freedom used for the baseline 
#       spline. 
#   Z: one-sided formula where the right-handside are the random effects.
#   cluster: vector with integers or factors for group identifiers.
#   do_setup: character vector indicating which approximation to setup. 
#             It is included to test computation time.
#   n_nodes: number of nodes to use in (adaptive) Gauss–-Hermite quadrature.
#   param_type: characters for the parameterization used with the SNVA.
#   skew_start: starting value for the Pearson's moment coefficient of 
#               skewness parameter when a SNVA is used. 
#   link: character specifying the link function.
#   theta: starting values for covariance matrix. 
#   beta: starting values for fixed effect coefficients.
# 
# Returns: 
#   TODO: what?
get_gsm_obj_func <- function(
  formula, data, df, Z, cluster, do_setup = c("Laplace", "GVA", "SNVA"), 
  n_nodes = 20L, param_type = c("DP", "CP_trans", "CP"), 
  skew_start = .alpha_to_gamma(-1), link = c("PH", "PO", "probit"), 
  theta = NULL, beta = NULL){
  if(!require(rstpm2))
    stop(sprintf("Requires %s", sQuote("rstpm2")))
  if(!require(TMB))
    stop(sprintf("Requires %s", sQuote("TMB")))
  link <- link[1]
  param_type <- param_type[1]
  skew_boundary <- 0.99527
  stopifnot(
    is.integer(df), df > 0L, inherits(formula, "formula"), 
    inherits(Z, "formula"), is.data.frame(data), 
    !missing(cluster), 
    all(do_setup %in% c(.laplace_char, .gva_char, .snva_char)), 
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L, 
    is.numeric(skew_start), length(skew_start) == 1L, 
    link %in% c("PH", "PO", "probit"), 
    param_type %in% c("DP", "CP_trans", "CP"))
  eval(bquote(stopifnot(
    .(-skew_boundary) < skew_start && skew_start < .(skew_boundary))))
  
  #####
  # get design matrices and outcome
  mf_X <- model.frame(formula, data = data)
  mt_X <- terms(mf_X)
  X <- model.matrix(mt_X, mf_X)
  n_x_fix <- NCOL(X)
  
  y <- model.response(mf_X)
  stopifnot(inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"))
  event <- y[, 2]
  tobs  <- y[, 1]
    
  formula_Z <- Z
  mf_Z <- model.frame(formula_Z, data = data)
  mt_Z <- terms(mf_Z)
  Z <- model.matrix(mt_Z, mf_Z)
  n_rng <- NCOL(Z)
  stopifnot(NCOL(X) > 0, NCOL(Z) > 0, NROW(X) == NROW(Z), 
            length(y) == NROW(X))
  
  # get the cluster variable
  cluster <- substitute(cluster)
  grp <- eval(cluster, data)
  stopifnot(is.factor(grp) || is.integer(grp) || is.character(grp), 
            isTRUE(length(grp) == NROW(X)))
  if(!is.factor(grp))
    grp <- as.factor(grp)
  grp <- as.integer(grp)
  n_grp <- length(unique(grp))
  
  #####
  # add time-varying baseline
  time_var <- formula[[2L]][[2L]]
  formula_b <- eval(
    bquote(~ nsx(log(.(time_var)), df = .(df), intercept = FALSE) - 1))
  mt_b <- terms(model.frame(formula_b, data = data[event > 0, ]))
  X <- cbind(X, model.matrix(mt_b, data))
  
  #####
  # approximate derivatives w/ finite difference
  is_fix <- 1:n_x_fix 
  XD <- X
  XD[,  is_fix] <- 0.
  XD[, -is_fix] <- local({
    dt <- .Machine$double.eps^(1/3)
    dat_m1 <- dat_p1 <- eval(bquote(data.frame(.(time_var))), data)
    x <- dat_m1[[1L]]
    dat_p1[[1L]] <- x * exp( dt)
    dat_m1[[1L]] <- x * exp(-dt)
    (model.matrix(mt_b, dat_p1) - model.matrix(mt_b, dat_m1)) / 2 / dt / x
  })
  
  #####
  # get starting values w/ coxph and then a lm fit
  stopifnot(isTRUE(attr(mt_Z, "intercept") == 1L),  
            colnames(Z)[1] == "(Intercept)")
  inits <- local({
    link_hat <- local({
      cox_fit <- coxph(formula, data, model = TRUE)
      S_hat <- rstpm2:::Shat(cox_fit) # was lazy...
      if(link == "PH")
        pmax(log(.Machine$double.eps) / 4, log(-log(S_hat)))
      else if(link == "PO")
        log((1 - S_hat) / S_hat)
      else if(link == "probit")
        -qnorm(S_hat)
      else
        stop(sprintf("%s not implemented", sQuote(link)))
    })
    
    cox_frm <- eval(bquote(
      update(formula, . ~ . + frailty(.(cluster), 
                                      distribution = "gaussian"))))
    cox_fit <- coxph(cox_frm, data)
    
    list(theta = cox_fit$history[[1L]]$theta, link_Hat = link_hat)
  })
  
  inits$coef <- local({
    keep <- event > 0
    lm.fit(x = X[keep, , drop = FALSE], 
           y = inits$link_Hat[keep])$coefficients
  })
  
  #####
  # setup ADFun object for the Laplace approximation
  data_ad_func <- list(
    tobs = tobs, event = event, X = X, XD = XD, Z = Z, grp = grp - 1L,
    link = link)
  theta_start <- if(n_rng == 1L)
    log(inits$theta) / 2 else local({
      out <- matrix(0., n_rng, n_rng)
      out[1L, 1L] <- sqrt(inits$theta)
      diag(out)[-1] <- 1
      keep <- upper.tri(out, TRUE)
      nam <- paste0("R", outer(1:n_rng, 1:n_rng, paste, sep = "."))
      structure(out[keep], names = nam[keep])
    })
  
  # the user may have provided values
  theta <- if(!is.null(theta)){
    stopifnot(all(dim(theta) == c(n_rng, n_rng)))
    structure(
      if(n_rng > 1L) chol(theta)[upper.tri(theta, TRUE)] else 
        log(theta) / 2, 
      names = names(theta_start))
  } else theta_start
  beta <- if(!is.null(beta)){
    stopifnot(length(beta) == length(inits$coef))
    structure(beta, names = names(inits$coef))
  } else inits$coef
  
  # assign parameter list
  params = list(
    eps = .Machine$double.eps^(1/2), kappa = 1e8, b = beta, 
    theta = theta)
  
  if(.laplace_char %in% do_setup){ 
  # get Laplace AD function
  adfunc_laplace <- local({
    data_ad_func <- c(
      list(app_type = .laplace_char), data_ad_func)
    # TODO: initialize random effects in a smarter way
    params$u <- matrix(0., n_rng, n_grp)
    MakeADFun(
      data = data_ad_func, parameters = params, DLL = .tmb_lib_f,
      silent = TRUE, random = "u")
  })
  
  # we make a wrapper object to account for the eps and kappa and allow the 
  # user to change these
  laplace_out <- with(new.env(), {
    eps <- adfunc_laplace$par["eps"]
    kappa <- adfunc_laplace$par["kappa"]
    fn <- adfunc_laplace$fn
    gr <- adfunc_laplace$gr
    he <- adfunc_laplace$he
    get_x <- function(x)
      c(eps = eps, kappa = kappa, x)
    
    out <- adfunc_laplace[
      !names(adfunc_laplace) %in% c("par", "fn", "gr", "he")]
    
    par <- adfunc_laplace$par[-(1:2)]
    names(par)[seq_along(inits$coef)] <- names(inits$coef)
    
    c(list(
      par = par,  
      fn = function(x, ...){ fn(get_x(x))                               }, 
      gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
      he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] }, 
      # function to set penalty parameters
      update_pen = function(eps, kappa){
        p_env <- parent.env(environment())
        if(!missing(eps))
          assign("eps"  , eps  , p_env)
        if(!missing(kappa))
          assign("kappa", kappa, p_env)
        invisible(with(p_env, c(eps = eps, kappa = kappa)))
      }
    ), out)
  })
  } else 
    laplace_out <- NULL
  
  #####
  # setup ADFun object for the GVA
  if(.gva_char %in% do_setup){
  # set the initial values
  n_mu     <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)
  
  # set means to zero
  idx_mean <- sapply(1:n_grp - 1L, function(i) i * n_p_grp + 1:n_mu)
  theta_VA[idx_mean] <- 0.
  # set parameters for Lambda = L^\top L
  idx_Lambda <- sapply(1:n_grp - 1L, function(i) 
    i * n_p_grp + n_mu + 1:n_Lambda)
  theta_VA[idx_Lambda] <- params$theta
  # set names
  theta_VA_names <- if(n_rng == 1L)
    c("mu1", "log_sd_va1") else
      c(paste0("mu", 1:n_rng), gsub("^R", "L", names(params$theta)))
  theta_VA_names <- c(outer(theta_VA_names, 1:n_grp, paste, sep = ":"))
  names(theta_VA) <- theta_VA_names
  
  adfunc_VA <- local({
    data_ad_func <- c(
      list(app_type = .gva_char), data_ad_func, 
      list(n_nodes = n_nodes))
    params$theta_VA <- theta_VA
      
    MakeADFun(
      data = data_ad_func, parameters = params, DLL = .tmb_lib_f,
      silent = TRUE)
  })
  
  # we make a wrapper object to account for the eps and kappa and allow the 
  # user to change these
  gva_out <- with(new.env(), {
    eps <- adfunc_VA$par["eps"]
    kappa <- adfunc_VA$par["kappa"]
    fn <- adfunc_VA$fn
    gr <- adfunc_VA$gr
    he <- adfunc_VA$he
    get_x <- function(x)
      c(eps = eps, kappa = kappa, x)
    
    out <- adfunc_VA[
      !names(adfunc_VA) %in% c("par", "fn", "gr", "he")]
    
    par <- adfunc_VA$par[-(1:2)]
    names(par)[seq_along(inits$coef)] <- names(inits$coef)
    idx_va <- (length(par) - length(theta_VA_names) + 1):length(par)
    names(par)[idx_va] <-
      theta_VA_names
    
    c(list(
      par = par,  
      fn = function(x, ...){ fn(get_x(x))                               }, 
      gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
      he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] }, 
      # function to set penalty parameters
      update_pen = function(eps, kappa){
        p_env <- parent.env(environment())
        if(!missing(eps))
          assign("eps"  , eps  , p_env)
        if(!missing(kappa))
          assign("kappa", kappa, p_env)
        invisible(with(p_env, c(eps = eps, kappa = kappa)))
      },
      # function to get parameters
      get_params = function(x)
        x[-idx_va], 
      control = list(maxit = 1000L)
    ), out)
  })
  } else 
    gva_out <- NULL
  
  #####
  # setup ADFun object for the SNVA
  if(.snva_char %in% do_setup){
  # set the initial values
  n_mu     <- n_rng
  n_rho    <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda + n_rho
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)
  
  # set the means
  theta_VA <- 
    rep(NA_real_, (n_rng * 2 + (n_rng * (n_rng + 1L)) / 2L) * n_grp)
  
  idx_mu     <- 
    sapply(1:n_grp - 1L, function(i) i * n_p_grp + 1:n_mu)
  idx_lambda <- 
    sapply(1:n_grp - 1L, function(i) i * n_p_grp + 1:n_Lambda + n_mu)
  idx_rho    <- 
    sapply(1:n_grp - 1L, function(i) i * n_p_grp + 1:n_rho + n_mu + n_Lambda)
  
  # TODO: do something smarter...
  if(param_type != "DP"){
    if(n_rng == 1L){
      stopifnot(n_rng == 1L)
      
      # set the mean to zero
      theta_VA[idx_mu] <- 0
      # set the log variance to match that of the random effect distribution
      theta_VA[idx_lambda] <- params$theta
      # set the Pearson's moment coefficient of skewness
      do_trans <- param_type == "CP_trans"
      theta_VA[idx_rho] <- 
        if(do_trans)
          log((skew_boundary + skew_start) / (skew_boundary - skew_start))
        else 
          skew_start
      # set the names
      theta_VA_names <- c("mu1", "log_sd_va1", 
                          if(do_trans) "trans_skew1" else "skew1")
      theta_VA_names <- c(outer(theta_VA_names, 1:n_grp, paste, sep = ":"))
      names(theta_VA) <- theta_VA_names
      
    } else {
      # set the means to zero
      theta_VA[idx_mu] <- 0
      # add upper triangle of matrix C such that Cov = C^\top C
      theta_VA[idx_lambda] <- params$theta
      # set the Pearson's moment coefficient of skewness
      do_trans <- param_type == "CP_trans"
      theta_VA[idx_rho] <- 
        if(do_trans)
          log((skew_boundary + skew_start) / (skew_boundary - skew_start))
        else 
          skew_start
      # set the names
      L <- matrix(paste0("L", outer(1:n_rng, 1:n_rng, paste, sep = ".")), 
                  n_rng)
      theta_VA_names <- c(
        paste0("mu", 1:n_rng), L[upper.tri(L, TRUE)], 
        paste0(if(do_trans) "trans_skew" else "skew", 1:n_rng))
      theta_VA_names <- c(outer(theta_VA_names, 1:n_grp, paste, sep = ":"))
      names(theta_VA) <- theta_VA_names
      
    }
    
  } else {
    if(n_rng == 1L){
      alpha <- .gamma_to_alpha(skew_start)
      nu    <- sqrt(2 / pi) * alpha / sqrt(1 + alpha^2) 
      omega <- exp(params$theta) / sqrt(1 - nu^2)
      rho   <- alpha / omega
      xi    <- -omega * nu 
      # set the log scale parameters
      theta_VA[idx_mu] <- log(omega)
      # set the skew/shape parameters
      theta_VA[idx_rho] <- rho
      # set the location parameters such that the mean is zero
      theta_VA[idx_lambda] <- xi
      # set names
      theta_VA_names <- c("mu1", "log_sd_va1", "rho1")
      theta_VA_names <- c(outer(theta_VA_names, 1:n_grp, paste, sep = ":"))
      names(theta_VA) <- theta_VA_names

    } else {
      alpha <- .gamma_to_alpha(skew_start)
      nu    <- sqrt(2 / pi) * alpha / sqrt(1 + alpha^2)
      
      # get the starting value for the covariance matrices
      R_start <- matrix(0., n_rng, n_rng)
      R_start[upper.tri(R_start, TRUE)] <- params$theta
      R_start <- crossprod(R_start)
      omega <- sqrt(diag(R_start) / (1 - nu^2))
      
      # map to Lambda
      dnu <- omega * nu
      R_start <- R_start + outer(dnu, dnu)
      keep <- upper.tri(R_start, TRUE)
      R_start <- chol(R_start)[keep]
      
      # insert values 
      theta_VA[idx_mu]     <- -dnu
      theta_VA[idx_lambda] <- R_start
      theta_VA[idx_rho]    <- alpha / omega
      
      # set names
      theta_VA_names <- c(
        paste0("mu" , 1:n_mu), 
        paste0("L", outer(1:n_rng, 1:n_rng, paste, sep = "."))[keep], 
        paste0("rho", 1:n_rho))
      theta_VA_names <- c(outer(theta_VA_names, 1:n_grp, paste, sep = ":"))
      names(theta_VA) <- theta_VA_names
    }
  }
  
  adfunc_VA <- local({
    data_ad_func <- c(
      list(app_type = .snva_char), data_ad_func, 
      list(n_nodes = n_nodes, param_type = param_type))
    params$theta_VA <- theta_VA
    
    MakeADFun(
      data = data_ad_func, parameters = params, DLL = .tmb_lib_f,
      silent = TRUE)
  })
  
  # we make a wrapper object to account for the eps and kappa and allow the 
  # user to change these
  snva_out <- with(new.env(), {
    eps <- adfunc_VA$par["eps"]
    kappa <- adfunc_VA$par["kappa"]
    fn <- adfunc_VA$fn
    gr <- adfunc_VA$gr
    he <- adfunc_VA$he
    get_x <- function(x)
      c(eps = eps, kappa = kappa, x)
    
    out <- adfunc_VA[
      !names(adfunc_VA) %in% c("par", "fn", "gr", "he")]
    
    par <- adfunc_VA$par[-(1:2)]
    names(par)[seq_along(inits$coef)] <- names(inits$coef)
    idx_va <- (length(par) - length(theta_VA_names) + 1):length(par)
    names(par)[idx_va] <-
      theta_VA_names
    
    c(list(
      par = par,  
      fn = function(x, ...){ fn(get_x(x))                               }, 
      gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
      he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] }, 
      # function to set penalty parameters
      update_pen = function(eps, kappa){
        p_env <- parent.env(environment())
        if(!missing(eps))
          assign("eps"  , eps  , p_env)
        if(!missing(kappa))
          assign("kappa", kappa, p_env)
        invisible(with(p_env, c(eps = eps, kappa = kappa)))
      },
      # function to get parameters
      get_params = function(x)
        x[-idx_va], 
      control = list(maxit = 1000L)
    ), out)
  })
  } else 
    snva_out <- NULL
  
  list(laplace = laplace_out, gva = gva_out, snva = snva_out, y = y, 
       event = event, X = X, XD = XD, Z = Z, grp = grp, terms = list(
         X = mt_X, Z = mt_Z, baseline = mt_b))
}

# run tests 
local({
  func <- MakeADFun(
    data = list(app_type = "run tests"), parameters = list(tobs = 0), 
    DLL = .tmb_lib_f, silent = TRUE)
})

invisible()