### Final SE as of Dec 20, 2018

score.i <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  expt <- dBvZINB5.Expt.se(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4)
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


info <- function(xvec, yvec, freq = rep(1, length(xvec)), a0, a1, a2, b1, b2, p1, p2, p3, p4, exact = FALSE, ...) {
  
  s <- score(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...)
print(class(s))
head(s)
  result <- s %*% (t(s) * freq)
  
  if (exact) {
    n <- sum(xy.reduced$freq)
    s.sum <- s %*% freq
    result <- result - s.sum %*% t(s.sum)/n
  }

  result
}

if (FALSE) {
  do.call(info, list(3,4,1,2,3,1,2,.3,.4,.1,.2))
  do.call(info, list(c(3,4), c(4,5),1,2,3,1,2,.3,.4,.1,.2))
}

#' @useDynLib bzinb
#' @export
BZINB5.se <- function(xvec, yvec, freq = rep(1, length(xvec)), a0, a1, a2, b1, b2, p1, p2, p3, p4, ...) {
  if (any(!is.finite(c(a0,a1,a2,b1,b2,p1,p2,p3,p4)))) {return(rep(NA, 10))}
  info.mat <- info(xvec, yvec, freq, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...)
  
  # singular matrix
  if (qr(info.mat)$rank < 8) {
    warnings("The information matrix is not full rank, and thus is not invertible.")
    se <- rep(NA, 11)
    names(se) <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4", "rho", "logit.rho")
    return(se)
  } 
  
  cov.mat <- try(solve(info.mat))
  if (class(cov.mat) == "try-error") {
    se <- rep(NA, 10)
    names(se) <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4", "rho", "logit.rho")
    return(se)
  }

  # variance of p4 hat
  d.p4 <- -c(1, 1, 1)
  var.p4 <- t(d.p4) %*% cov.mat[6:8, 6:8] %*% d.p4

  # variance of rho hat
  rho <- a0/sqrt((a0 + a1) * (a0 + a2)) *sqrt(b1 *b2 /(b1 + 1) /(b2 + 1))
  d.g <- rho * c(1/a0 - 1/{2*(a0 + a1)} - 1/{2*(a0 + a2)}, - 1/{2*(a0 + a1)}, - 1/{2*(a0 + a2)},
                 1/{2 *b1 *(b1 + 1)}, 1/{2 *b2 *(b2 + 1)})
  var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g

  # variance of logit(rho hat)
  var.logit.rho <- var.rho / rho / (1-rho)

  se <- sqrt(c(diag(cov.mat), p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
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

#'
#'
#' @examples
#' set.seed(2)
#' xy <- rBvZINB5(500, 1,10,11,.5,.5,.4,.2,.2,.2)
#' x = as.numeric(xy[,1])
#' y = as.numeric(xy[,2])
#' ci.rho(x, y, 1.3, 7.4, 6.4, 0.6, 0.8, 0.4, 0.2, 0.2, 0.2)
#'
#' @useDynLib bzinb
#' @export
ci.rho <- function(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, logit = TRUE, alpha = 0.05, ...) {
  if (alpha > 0.5) {stop("alpha is too large (> 0.5).")}
  rho <- a0/sqrt((a0 + a1) * (a0 + a2)) *sqrt(b1 *b2 /(b1 + 1) /(b2 + 1))
  se <- BZINB5.se(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, ...)

  # logit or not
  if (logit) {
    m <- qlogis(rho)
    s <- se["logit.rho"]
  } else {
    m <- rho
    s <- se["rho"]
  }

  # CI
  ci <- m + c(lb = -1, ub = 1) * qnorm(1 - alpha/2) * s
  if (logit) ci <- plogis(ci)

  return(c(rho = rho, ci))
}
