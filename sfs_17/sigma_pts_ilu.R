library(mvtnorm)
x.points <- seq(-3,3,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mu <- c(0,0)
sigma <- matrix(c(2,1,1,1),nrow=2)
for (i in 1:100) {for (j in 1:100) {
  z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                    mean=mu,sigma=sigma)
}}

plot(range(x.points), range(y.points), type = "n", xlab = "", ylab = "", 
     frame = FALSE, axes = FALSE)

lvls <- seq(0, .17, by = .01)
.filled.contour(x.points,y.points, z, levels = lvls, col = rev(heat.colors(length(lvls), alpha = .9)))

decomp <- chol(sigma)

q <- 2
lambda <- .1

pts <- rbind(
  c(0, 0), 
  sqrt(q + lambda) * decomp,
  - sqrt(q + lambda) * decomp)

points(pts[, 1], pts[, 2], pch = 16, cex = par()$cex * .75)