### Final SE as of Dec 20, 2018

score.i <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  expt <- dBvZINB5.Expt(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, debug = FALSE)
  score <- 
    c(expt["log.R0.E"] - log(b1) - digamma(a0),
      expt["log.R1.E"] - log(b1) - digamma(a1),
      expt["log.R2.E"] - log(b1) - digamma(a2),
      - (expt["v.E"] + a0 + a1 + a2)/b1 + (expt["R0.E"] + expt["R2.E"]) * (b2 + 1)/b1^2 + 
        expt["R1.E"]/b1^2,
      expt["v.E"]/b2 - (expt["R0.E"] + expt["R2.E"]) /b1,
      expt["E1.E"]/p1 - expt["E4.E"]/p4,
      expt["E2.E"]/p2 - expt["E4.E"]/p4,
      expt["E3.E"]/p3 - expt["E4.E"]/p4)
  names(score) <- c(paste0("a", 0:2), paste0("b", 1:2), paste0("p", 1:3))
  score
}
score <- Vectorize(score.i, vectorize.args = c("x", "y"))

info <- function(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, exact = FALSE, ...) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  
  s <- score(xy.reduced$x, xy.reduced$y, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...)
  result <- s %*% (t(s) * xy.reduced$freq)
  
  if (exact) {
    s.sum <- s %*% xy.reduced$freq
    result <- result - s.sum %*% t(s.sum)/n
  }
  result
}
if (FALSE) {
  do.call(info, list(3,4,1,2,3,1,2,.3,.4,.1,.2))
  do.call(info, list(c(3,4), c(4,5),1,2,3,1,2,.3,.4,.1,.2))
  }

BZINB5.se <- function(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...) {
  if (any(!is.finite(c(a0,a1,a2,b1,b2,p1,p2,p3,p4)))) {return(rep(NA, 10))}
  info.mat <- info(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...)
  cov.mat <- try(solve(info.mat))
  if (class(cov.mat) == "try-error") {
    se <- rep(NA, 10)
    names(se) <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4", "rho")
    return(se)
  }
  
  # variance of rho hat
  rho <- a0/sqrt((a0 + a1) * (a0 + a2)) *sqrt(b1 *b2 /(b1 + 1) /(b2 + 1))
  d.g <- rho * c(1/a0 - 1/{2*(a0 + a1)} - 1/{2*(a0 + a2)}, - 1/{2*(a0 + a1)}, - 1/{2*(a0 + a2)},
           1/{2 *b1 *(b1 + 1)}, 1/{2 *b2 *(b2 + 1)})
  var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g
  
  # variance of p4 hat
  d.p4 <- -c(1, 1, 1)
  var.p4 <- t(d.p4) %*% cov.mat[6:8, 6:8] %*% d.p4
  
  se <- sqrt(c(diag(cov.mat), p4 = var.p4, rho=var.rho))
  se
}

if (FALSE) {
  ### test!!
  library(dplyr)
  set.seed(2)
  xy <- rBvZINB5(500, 1,10,11,.5,.5,.4,.2,.2,.2)
  x = as.numeric(xy[,1])
  y = as.numeric(xy[,2])
  info(x, y, 1,2,3,1,1,.5,.15,.25,.1) %>% round(0)
  info(x, y, 1,2,3,1,1,.5,.15,.25,.1) %>% solve
  tmp.parm <- ML.BvZINB5(x, y, maxiter = 50)
  do.call(BZINB5.se, c(list(x,y), as.list(tmp.parm[1:9])))
  BZINB5.se(x, y, 1,1,1,1,1,.5,.2,.2,.1) %>% round(1)
  
  tt(1)
  do.call(BZINB5.se, c(list(x,y), as.list(tmp.parm[1:9])))
  tt(2)
  
}
