#####
# define function to winsorize
wz <- function(x, probs = c(.01, .99), lb = NULL, ub = NULL, do_scale = FALSE, 
               scale = NULL, do_center = FALSE, mu = NULL){
  if(is.null(lb))
    lb <- unname(quantile(x, probs[1]))
  if(is.null(ub))
    ub <- unname(quantile(x, probs[2]))
  
  x <- pmin(pmax(x, lb), ub)
  if(is.null(scale))
    scale <- sd(x)
  if(is.null(mu))
    mu <- mean(x)
  if(do_center)
    x <- x - mu
  if(do_scale)
    x <- x / scale
  
  structure(x, lb = lb, ub = ub, class = "wz", scale = scale, 
            do_scale = do_scale, do_center = do_center, mu = mu)
}

makepredictcall.wz <- function(var, call){
  if(as.character(call)[1L] != "wz") 
    return(call)
  at <- attributes(var)[c("lb", "ub", "scale", "do_scale", "do_center", 
                          "mu")]
  xxx <- call[1L:2L]
  xxx[names(at)] <- at
  xxx
}

local({
  set.seed(48350025)
  df <- data.frame(y = rnorm(100), x = rnorm(100, 1, sd = 2))
  f <- lm(y ~ wz(x), df)
  mm <- model.matrix(f)
  stopifnot(min(mm[, 2]) == quantile(df$x, .01), 
            max(mm[, 2]) == quantile(df$x, .99))
  stopifnot(
    predict(f, newdata = df)[1:10] == predict(f, newdata = df[1:10, ]))
  
  f <- lm(y ~ wz(x, do_scale = TRUE), df)
  mm <- model.matrix(f)
  stopifnot(
    isTRUE(all.equal(sd(mm[, 2]), 1)),
    predict(f, newdata = df)[1:10] == predict(f, newdata = df[1:10, ]))
  
  f <- lm(y ~ wz(x, do_scale = TRUE, do_center = TRUE), df)
  mm <- model.matrix(f)
  stopifnot(
    isTRUE(all.equal(sd(mm[, 2]), 1)),
    isTRUE(all.equal(mean(mm[, 2]), 0)),
    predict(f, newdata = df)[1:10] == predict(f, newdata = df[1:10, ]))
})

#####
# function to winsorize and then make B-spline basis matrix for a 3. order
# polynomial spline with a sum-to-zero constraint and which is orthogonal to 
# 1. order term. The latter is to make it easy to test the significance w/
# `drop1`
sp_w_c <- function(
  x, df = NULL, lb = NULL, ub = NULL, probs = c(.01, .99), 
  Boundary.knots = NULL, knots = NULL, Z = NULL, do_excl_slope = TRUE, 
  do_center = FALSE, mu = NULL){
  # first winsorize
  x <- wz(x, lb = lb, ub = ub, do_center =  do_center, mu = mu)
  lb <- attr(x, "lb")
  ub <- attr(x, "ub")
  do_center <- attr(x, "do_center")
  mu <- attr(x, "mu")
  
  if(is.null(Boundary.knots))
    Boundary.knots <- c(lb, ub) - do_center * mu
  
  # apply a sum-to-zero constraint
  require(splines)
  X <- bs(x, df = df + 1L, intercept = TRUE, knots = knots, 
          Boundary.knots = Boundary.knots)
  knots <- attr(X, "knots")
  
  if(is.null(Z)){
    if(do_excl_slope){
      C <- crossprod(cbind(rep(1, nrow(X)), x), X)
      qrc <- qr(t(C))
      Z <- qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C)]
      
    } else {
      C <- rep(1, nrow(X)) %*% X
      qrc <- qr(t(C))
      Z <- qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C)]
      
    }
  }
  
  structure(
    X %*% Z, lb = lb, ub = ub, knots = knots, Boundary.knots = Boundary.knots,
    do_excl_slope = do_excl_slope, do_center = do_center, mu = mu, 
    Z = Z, class = "sp_w_c")
}

makepredictcall.sp_w_c <- function(var, call){
  if(as.character(call)[1L] != "sp_w_c") 
    return(call)
  at <- attributes(var)[c("lb", "ub", "knots", "Boundary.knots", "Z", 
                          "do_excl_slope", "do_center", "mu")]
  xxx <- call[1L:2L]
  xxx[names(at)] <- at
  xxx
}

local({
  set.seed(48350025)
  df <- data.frame(y = rnorm(100), x = rnorm(100, 1, sd = 2))
  f <- lm(y ~ wz(x) + sp_w_c(x, 4L), df)
  mm <- model.matrix(f)
  stopifnot(
    ncol(mm) == 5L,
    sum(attr(mm, "assign") == 2) == 3,
    predict(f, newdata = df)[1:10] == predict(f, newdata = df[1:10, ]), 
    # are orthogonal
    isTRUE(all.equal(c(mm[, 1] %*% mm[, -(1:2)]), rep(0, 3))),
    isTRUE(all.equal(c(mm[, 2] %*% mm[, -(1:2)]), rep(0, 3))))
  
  f <- lm(y ~ wz(x, do_center = TRUE) + sp_w_c(x, 4L, do_center = TRUE), df)
  mm <- model.matrix(f)
  stopifnot(
    ncol(mm) == 5L,
    sum(attr(mm, "assign") == 2) == 3,
    predict(f, newdata = df)[1:10] == predict(f, newdata = df[1:10, ]), 
    # are orthogonal
    isTRUE(all.equal(c(mm[, 1] %*% mm[, -(1:2)]), rep(0, 3))),
    isTRUE(all.equal(c(mm[, 2] %*% mm[, -(1:2)]), rep(0, 3))))
  
  
  # we get the same if we call `bs` and omit the first order term
  f2 <- lm(y ~ bs(wz(x), 4L), df) 
  stopifnot(isTRUE(all.equal(predict(f), predict(f2))))
  
  f3 <- lm(y ~ sp_w_c(x, 4L, do_excl_slope = FALSE), df)
  stopifnot(
    isTRUE(all.equal(predict(f), predict(f3))), 
    sum(attr(model.matrix(f3), "assign") == 1) == 4)
  
  f4 <- lm(y ~ wz(x) + sp_w_c(x, 4L, do_center = TRUE), df)
  stopifnot(isTRUE(all.equal(predict(f), predict(f4))))
})
