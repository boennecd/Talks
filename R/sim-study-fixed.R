library(pedmod)

# the seeds we will use
seeds <- c(53236043L, 47565336L, 54716919L, 52688711L, 48994978L, 6229287L, 11669871L, 5399628L, 78071011L, 87439610L, 12972561L, 53161981L, 33971704L, 3969741L, 22596857L, 33375033L, 9874108L, 82974789L, 74279235L, 65729841L, 49504044L, 66004537L, 41917002L, 47290080L, 7606093L, 50492586L, 1632442L, 9652979L, 37762574L, 57034155L, 44965518L, 39814183L, 28347229L, 87472255L, 82416354L, 72685955L, 25187024L, 21176973L, 95538326L, 53922491L, 78608450L, 65020476L, 82682654L, 34389520L, 59008433L, 56776517L, 20767476L, 33745720L, 39107668L, 48055072L, 86806447L, 9133178L, 85628899L, 33707378L, 76349038L, 46513772L, 70394680L, 54804785L, 41060870L, 51220936L, 47148984L, 87842321L, 13879797L, 41503242L, 41088927L, 63924334L, 82560061L, 33620057L, 16944874L, 18707922L, 45964288L, 64592760L, 47959719L, 58529746L, 90825649L, 89114177L, 99916412L, 66589690L, 30833593L, 38718675L, 46217135L, 95830312L, 7031732L, 24918922L, 93634584L, 42506661L, 67635440L, 81179108L, 24920900L, 92604592L, 95409602L, 41703707L, 43978274L, 83385074L, 57230964L, 93605861L, 15407548L, 66466998L, 69007420L, 24482521L)

# get the data for the family structure we will use
base_fam <- local({
  # assign the family data
  fam_dat <- data.frame(
    id = 1:16, 
    sex = c(rep(1:2, 2), 1L, 2L, 2L, 1:2, 1, rep(1:2, 3)), 
    dadid = c(rep(NA, 4), NA, 1L, 1L, 3L, 3L, NA, 5L, 5L, 8L, 8L, 10L, 10L), 
    momid = c(rep(NA, 4), NA, 2L, 2L, 4L, 4L, NA, 6L, 6L, 7L, 7L,  9L,  9L))

  # get the scale matrix for the genetic effect
  library(kinship2)
  ped_dat <- with(fam_dat, pedigree(
    id = id, dadid = dadid, momid = momid, sex = sex, 
    famid = rep(1, NROW(fam_dat))))["1"]
  rel_mat <- 2 * kinship(ped_dat)

  # get the scale matrix for the maternal effect
  met_mat <- diag(NROW(fam_dat))
  for(i in seq_len(NROW(fam_dat)))
    for(j in seq_len(NROW(fam_dat))){
      mi <- fam_dat$momid[i]
      mj <- fam_dat$momid[j]
      if(is.na(mi) || is.na(mj)) 
        next
      
      same_grand <- fam_dat$momid[mi] == fam_dat$momid[mj]
      if(mi == mj)
        met_mat[i, j] <- 1
      else if(!is.na(same_grand) && same_grand)
        met_mat[i, j] <- .5
    }
  
  # fixed effect design matrix
  X <- cbind(
    `(Intercept)` = 1, X = seq(-2, 2, length.out = NROW(fam_dat)), 
    Dummy = seq_len(NROW(fam_dat)) %% 2 == 0)
  
  list(rel_mat = rel_mat, met_mat = met_mat, X = X, ped_dat = ped_dat)
})

# parameters
b <- c(-2, .5, 1)
sig <- c("Genetic" = .5, "Maternal" = .33)
n_families <- 10000L

# run the simulation study
library(microbenchmark)
fix_sim_res <- lapply(seeds, function(s){
  f <- file.path(
    "cache", "fixed", sprintf("%d-%d.RDS", n_families, s))
  if(!file.exists(f)){
    # simulate the data
    set.seed(s)
    Sig <- 
      diag(NROW(base_fam$met_mat)) + base_fam$rel_mat * sig["Genetic"] +
      base_fam$met_mat * sig["Maternal"]
    Sig_chol <- chol(Sig)
    eta <- drop(base_fam$X %*% b)
    
    dat <- replicate(n_families, {
      y <- drop(eta + rnorm(NROW(Sig)) %*% Sig_chol > 0)
      
      with(base_fam, list(
        y = y, X = X, 
        scale_mats = list(rel_mat = rel_mat, met_mat = met_mat)))
    }, simplify = FALSE)
    
    # get starting values
    ti <- get_nanotime() # to track the computation time
    y <- unlist(lapply(dat, `[[`, "y"))
    X <- do.call(rbind, lapply(dat, `[[`, "X"))
    start_fit <-  glm.fit(X, y, family = binomial("probit"))
    beta <- start_fit$coefficients
    rm(y, X, start_fit)
    
    # starting values for the scale parameters
    sc <- setNames(rep(log(.2), 2), names(sig))
    
    # create a C++ object
    ll_terms <- get_pedigree_ll_terms(dat, max_threads = 4L)
    
    # assign a function to approximate the log likelihood and the gradient
    fn <- function(par, seed = 1L, rel_eps = 1e-2, use_aprx = TRUE, 
                   n_threads = 4L, indices = NULL, maxvls = 10000L, 
                   do_print = FALSE){
      set.seed(seed)
      out <- -eval_pedigree_ll(
        ll_terms, par = par, maxvls = maxvls, abs_eps = 0, rel_eps = rel_eps, 
        minvls = 1000L, use_aprx = use_aprx, n_threads = n_threads, 
        indices = indices)
      if(do_print){
        cat(sprintf("Negative log likelihood and parameters: %15.2f\n", 
                    out))
        print(par)
      }
      out
    }
    gr <- function(par, seed = 1L, rel_eps = 1e-2, use_aprx = TRUE, 
                   n_threads = 4L, indices = NULL, maxvls = 10000L, 
                   do_print = FALSE){
      set.seed(seed)
      out <- -eval_pedigree_grad(
        ll_terms, par = par, maxvls = maxvls, abs_eps = 0, rel_eps = rel_eps, 
        minvls = 1000L, use_aprx = use_aprx, n_threads = n_threads, 
        indices = indices)
      structure(c(out), value = -attr(out, "logLik"), 
                n_fails = attr(out, "n_fails"))
    }
    ti_setup <- (get_nanotime() - ti) / 1e9
    
    # optimize
    ti <- get_nanotime()
    opt <- optim(c(beta, sc), fn, gr, method = "BFGS", maxvls = 5000L)
    ti_opt <- (get_nanotime() - ti) / 1e9 + ti_setup
    
    # start w/ higher precision
    ti <- get_nanotime()
    opt_xtr <- optim(opt$par, fn, gr, method = "BFGS", maxvls = 25000L)
    ti_opt_xtr <- (get_nanotime() - ti) / 1e9 + ti_setup + ti_opt
    
    # save the output 
    opt$time <- ti_opt
    opt_xtr$time <- ti_opt_xtr
    out <- list(optim = opt, optim_xtr = opt_xtr, seed = s)
    
    message(sprintf("Times (w/o w/ extra samples):      %10.2f %10.2f", 
                    out$optim$time, out$optim_xtr$time))
    message(sprintf("Max log ML (w/o w/ extra samples): %10.2f %10.2f", 
                    -out$optim$value, -out$optim_xtr$value))
    message(paste0(capture.output(rbind(
      truth                    = c(b, log(sig)),
      optim                    = out$optim$par, 
      `optim w/ extra samples` = out$optim_xtr$par)), 
      collapse = "\n"))
    message("")
    
    # save
    saveRDS(out, f)
  } else
    message(sprintf("Loading results with seed %d", s))
    
  readRDS(f)
})


## ----fixed_bias_ests-------------------------------------------------------------------
# get the estimates
estimates <- sapply(fix_sim_res, function(x) 
  cbind(`w/o extra` = x$optim$par, `w/ extra` = x$optim_xtr$par), 
  simplify = "array")

# transform back the scale parameters
estimates[-seq_along(b), , ] <- 
  exp(estimates[-seq_along(b), , ])
dimnames(estimates)[[1]] <- c(colnames(base_fam$X), names(sig))

# compute the errors
errors <- estimates - c(b, sig)

# compute bias estimates and standard errors
comp_se <- function(x)
  apply(x, 1, sd) / sqrt(NCOL(x))
bias <- rbind(`Bias w/o extra` = rowMeans(errors[, "w/o extra", ]),
              `SE w/o extra`   = comp_se (errors[, "w/o extra", ]), 
              `Bias w/ extra`  = rowMeans(errors[, "w/ extra", ]),
              `SE w/ extra`    = comp_se (errors[, "w/ extra", ]))

# gather the computation times
times <- sapply(fix_sim_res, function(x) 
  c(`W/o extra` = x$optim$time, `W/ extra` = x$optim_xtr$time))

# compute summary statistics
time_stats <- t(apply(times, 1, function(x)
  c(Mean = mean(x), Median = median(x), Sd = sd(x))))
