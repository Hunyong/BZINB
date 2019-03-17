##########################################################################################
## 1. Bivariate Poisson
##########################################################################################

#' Density of bivariate Poisson distribution (bp)
#'
#' This function calculates the density of bp model
#' at a data point (x, y) under a set of parameters
#' (m0, m1, m2).
#'
#' @param x a number. The first element of the data point.
#' @param y a number. The second element of the data point.
#'
#' @param m0 a positive number. The shared mean parameter.
#' @param m1 a positive number. The second mean parameter.
#' @param m2 a positive number. The third mean parameter.
#'
#' @param log logical. If TRUE, the log-density is returned.
#' @param max a positive integer. If at least one of \code{x}
#' and \code{y} is greater than this value, the density of
#' zero is returned.
#'
#' @return the bp density at the data point
#'
#' @examples
#'
#' dbp(2, 2, 0, 1, 1)
#' dbp(100, 100, 0, 1, 1)
#'
#' # Sum of joint density should equal the marginal density
#' sum(sapply(0:100, function(s) dbp(s,2,1,2,2)))
#' dpois(2, 3)
#'
#' # joint(bp) > pois^2 (around diagonal)
#' dbp(2,2,1,2,2)
#' dpois(2,3)^2
#'
#' # joint(bp) < pois*pois (off-diagonal)
#' dbp(1,3,1,2,2)
#' dpois(1,3)*dpois(3,3)
#'
dbp <- function(x, y, m0, m1, m2, log = FALSE, max = 500) {
  # max = 500: when counts exceed the max, p = 0
  if (m1 * m2 == 0) {  # previous code doesn't consider when mu1, mu2, or both is 0
    if ((x - y) %*% (m1 - m2) >= 0) {
      m <- min(x,y); s <- x + y - 2*m; mm <- max(m1, m2)
      result <- dpois (m, m0) * dpois (s, mm)
      # print(c("m", m, "s", s, "mm", mm, "result", result))
    } else { result <- 0 }
    if (log) {result <- log(result)}
    return(result)
  }
  f1 <- dpois(x, m1, log = TRUE)
  f2 <- dpois(y, m2, log = TRUE)
  m <- min(x,y); p <- m0/m1/m2; if (m == 0) {p <- 0} # in order not to make p NaN
  fun.a <- function(x, y, s, p, adj) {
    ifelse (p == 0, ifelse(s==0, 1, 0), exp(lchoose(x,s)+lchoose(y,s)+lfactorial(s)+ s*log(p) - adj))
  }
  if (max(x,y) > 100) {adj = 300} else {adj = 0}  # Handle numerical error for large numbers
  if (max(x,y) > max) {result <- ifelse(log, -Inf, 0)} else {
    f3 <- log(sum(sapply(0:m, function(s) fun.a(x=x, y=y, s=s, p=p, adj=adj)))) + adj
    result <- f1 + f2 - m0 + f3
    if (!log) {result <- exp(result)}
  }
  if (!is.finite(ifelse(log,exp(result),result))) {result <- 0}
  return(result)
}
dbp.vec <- Vectorize(dbp)


bp <- function(xvec, yvec, tol = 1e-6) {
  .check.input()
  
  xvec <- as.numeric(xvec); yvec <- as.numeric(yvec)
  len <- length(xvec)
  # xbar <- mean(xvec)
  # ybar <- mean(yvec)
  vec <- cbind(xvec,yvec)
  bar <- apply(vec,2,mean)
  m <- min(bar)
  lik.bp <- function(a) (-sum(apply(vec, 1, dbp, m0 = a, m1 = bar[1] - a, m2 = bar[2] - a, log = TRUE)))
  # result <- nlminb(m, lik, lower = 0)
  if (m == 0) {  # when either is all zero, then mu0 is automatically 0.
    result = data.frame(par = 0, value = NA)  # lik not actually NA but can be obtained if necessary
  } else {
    result <- optim(m, lik.bp, lower = 0, upper = m, method="Brent")
  }
  result <- data.frame(mu0 = result$par, mu1 = bar[1] - result$par, mu2 = bar[2] - result$par, likelihood = -result$value)
  rownames(result) <- NULL
  return(result)
}
if (FALSE) { # example
  bp(x,y)
  bp(x,z)
  set.seed(1000); a1 <- rpois(20,1); a2 <- rpois(20,2); a3 <- rpois(20,3)
  bp(a1+a2, a1+a3)
}


##########################################################################################
## 2. Bivariate ZIP
##########################################################################################

## 2.1 basic functions for BvZIP
dbzip <- function(x, y = NULL, pp, m0, m1, m2, log = FALSE) {
  fxy <- (1-pp) * dbp (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) + ifelse((x == 0 & y == 0), pp,0)
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dbzip.vec <- Vectorize(dbzip)
lik.bzip <- function(x, y, param, pp = NULL, m0  = NULL, m1  = NULL, m2  = NULL){
  if (is.null(pp)|is.null(m0)|is.null(m1)|is.null(m2)) {
    pp = param[1]; m0 = param[2]; m1 = param[3]; m2 = param[4]
  }
  sum(dbzip.vec(x, y, pp = pp, m0 = m0, m1 = m1, m2 = m2, log=TRUE))
}
rbzip <- function(n, pp, m0, m1, m2) {
  z <- rbinom(n,1,pp)
  common.tmp <- rpois(n,m0)
  x.tmp <- rpois(n,m1)
  y.tmp <- rpois(n,m2)
  x <- (common.tmp + x.tmp)*(1-z)
  y <- (common.tmp + y.tmp)*(1-z)
  return(data.frame(x=x, y=y))
}
rbzip.vec <- Vectorize(rbzip)

## 2.2 param estimation
# formal EM algorithm
bzip <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  .check.input()
  
  # counter, likelihood, param.prev for recursive function
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(pp = 0, mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0

  # E-step
  fun.cond.exp <- function(x, y, pp, m0, m1, m2) {
    if (x+y == 0) {
      cond.prob <- pp/(pp+(1-pp)*exp(-m0-m1-m2))
      cond.expt <- c(1, m0, m1, m2) * cond.prob
    } else {
      m = min(x,y)
      prob <- sapply(0:m, function(u) {dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)})
      # print(prob); print(prob==0) #debug
      if (sum(prob) == 0) {
        # using the ratio of two probs
        if (m==0) {prob <- 1} else {
          prob <- sapply(0:(m-1), function(u) (m0/m1/m2*(x-u)*(y-u)/(u+1)))
          prob <- cumprod(c(1,prob))
          #prob <- sapply(0:m, function(u) {dbp(x = u, y=0, m0 = 0, m1 = m0, m2 = 0) * dbp(x = x-u, y=0, m0 = 0, m1 = m1, m2 = 0) *
          #    dbp(x = 0, y= y-u, m0 = 0, m1 = 0, m2 = m2)*10^100})
          # used dbp instead of simply dpois to handle large numbers (when mu=1, x = 172, dpois=0 not small number)
          # and adjust by 10^100 (num/den cancel out adj.)
        }
        #print(c(1,prob)) # debug
        #print(prob)
      }
      eu <- sum((0:m)*prob)/sum(prob)
      cond.expt <- c(0, 0, x, y) + eu *c(0, 1, -1, -1)
    }
    return(cond.expt)
  }
  fun.cond.exp <- Vectorize(fun.cond.exp)

  # M-step
  param.update <- function(x, y, pp, m0, m1, m2) {
    result <- fun.cond.exp(x = x, y = y, pp = pp, m0 = m0, m1 = m1, m2 = m2)
    return(apply(result,1,mean))
  }

  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- c(0.5,1,1,1)
  }

  # Repeat
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", "pi", paste0("mu",0:2)))}
  repeat {
    iter = iter + 1
    #print(c(param))
    #print(lik(vec, pp=param[1], m0=param[2], m1=param[3], m2=param[4])) # debug
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, pp = param[1], m0 = param[2], m1 = param[3], m2 = param[4])
    if (showFlag) {print(c(iter, round(param,5), lik.bzip(xvec, yvec, param)))}
    if (max(abs(param - param.old)) <= tol) {
      param <- data.frame(matrix(param,1,4))
      names(param) <- c("pp", "mu0", "mu1", "mu2")
      return(param)
      break
    }
  }
  return(param)
}

# bzip examples
if (FALSE) {
  bzip(vec[,1], vec[,2], showFlag=TRUE, initial = rep(0,4), tol = 1e-10)
  bzip(vec[,1], vec[,2], showFlag=TRUE, tol = 1e-8)
  bzip(extractor(1), extractor(2), initial = rep(0,4), showFlag=TRUE)
  bzip(rep(0,100), rep(0,100), initial = rep(0,4), showFlag=TRUE)
  # For all 0 pairs, identifiability issue: same likelihood for (pi=1) and (pi=0)
  bzip(c(1,rep(0,99)), rep(0,100), initial = rep(0,4), showFlag=TRUE)
  bzip(c(1,rep(0,799)), rep(0,800), initial = rep(0,4), showFlag=TRUE)
  bzip(c(1,rep(0,799)), rep(0,800), initial = c(799/800,0,0,0), showFlag=TRUE) # not converging. initial should be all zero
}


moment.bzip <- function(pp, m0, m1, m2) {
  mean.bp <- (m0 + c(m1, m2))
  var.bp <- diag(mean.bp); var.bp[c(2,3)] <- m0
  mean <- (1-pp)*mean.bp
  var <- (1-pp)*var.bp + pp*(1-pp) * mean.bp %o% mean.bp
  cor <- var[2] / sqrt(prod(diag(var)))
  return(list(mean = mean, var = var, cor = cor))
}


##########################################################################################
## 3. Bivariate zip.b: General bzip with marginal ZIP condition (6params)
##########################################################################################

## 3.1 basic functions for bzip.b
dbzip.b <- function(x, y = NULL, p1, p2, p3, p4 = 1 - p1 - p2 - p3, m0, m1, m2, log = FALSE) {
  fxy <- p1 * dbp (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) +
    p2 * {if (y == 0) dbp (x=x, y=y, m0 = 0, m1 = m0 + m1, m2 = 0) else 0} +
    p3 * {if (x == 0) dbp (x=x, y=y, m0 = 0, m1 = 0, m2 = m0 + m2) else 0} +
    p4 * {if (x + y == 0) 1 else 0}
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dbzip.b.vec <- Vectorize(dbzip.b)
lik.bzip.b <- function(x, y, param, p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, m0  = NULL, m1  = NULL, m2  = NULL){
  if (is.null(p1)|is.null(p2)|is.null(p3)|is.null(p4)|is.null(m0)|is.null(m1)|is.null(m2)) {
    p1 = param[1]; p2 = param[2]; p3 = param[3]; p4 = param[4];
    m0 = param[5]; m1 = param[6]; m2 = param[7]
  }
  sum(dbzip.b.vec(x, y, p1 = p1, p2 = p2, p3 = p3, p4 = p4,
                  m0 = m0, m1 = m1, m2 = m2, log=TRUE))
}
rbzip.b <- function(n, p1, p2, p3, p4, m0, m1, m2) {
  z <- rmultinom(n,1,c(p1,p2,p3,p4))
  u.tmp <- rpois(n,m0)
  v.tmp <- rpois(n,m1)
  w.tmp <- rpois(n,m2)
  #(Xi, Yi) = ((Ei1 + Ei2)(Ui + Vi), (Ei1 + Ei3)(Ui + Wi)).
  x <- (z[1,] + z[2,]) * (u.tmp + v.tmp)
  y <- (z[1,] + z[3,]) * (u.tmp + w.tmp)
  return(data.frame(x=x, y=y))
}

## 3.2 parameter estimation
# Method of moment, for starting values of MLE
mme.bzip.b <- function(xvec, yvec) {
  .check.input()
  
  m.x <- mean(xvec); m.x2 <- mean(xvec^2)
  m.y <- mean(yvec); m.y2 <- mean(yvec^2)
  m.x2.y <- mean(xvec^2*yvec); m.y2.x <- mean(yvec^2*xvec)
  m.x.y <- mean(xvec*yvec)
  w <- (m.x2.y + m.y2.x) / (2 * m.x.y)
  lambda1 <- m.x2 / m.x - 1
  lambda2 <- m.y2 / m.y - 1
  mu0 <- (lambda1*lambda2) * (w - (2 + lambda1 + lambda2)/2) / (lambda1 + lambda2 + 1 - w)
  #print(c(w,lambda1, lambda2, mu0))
  mu0 <- max(1e-10, min(mu0, m.x, m.y))
  mu1 <- m.x - mu0
  mu2 <- m.y - mu0

  pi1 <- m.x.y / (lambda1*lambda2 + mu0)
  pi2 <- m.x / lambda1 - pi1
  pi3 <- m.y / lambda2 - pi1
  pi4 <- 1 - pi1 - pi2 - pi3
  pi <- c(pi1, pi2, pi3, pi4)
  pi <- pmax(0, pmin(pi,1))
  pi <- pi/sum(pi)
  return(data.frame(pi1 = pi[1], pi2 = pi[2], pi3 = pi[3], pi4 = pi[4], mu0 = mu0, mu1 = mu1, mu2 = mu2))
}

bzip.b <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE, maxiter = 200) {
  .check.input()
  
  # counter, likelihood, param.prev for recursive function
  # initial guess manipulation (for small counts pi4=0)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 1,  mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0

  # E-step
  fun.cond.exp.a <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    pra <- p1 * dbp(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    prb <- p2 * (y==0) * dbp(x=x, y=0, m0 =0, m1 = m0 + m1, m2 = 0)
    prc <- p3 * (x==0) * dbp(x=0, y=y, m0 =0, m1 = 0, m2 = m0 + m2)
    prd <- p4 * (x + y ==0)
    pragr <- c(pra, prb, prc, prd)
    prsum <- sum(pragr)

    ### Conditional expectations (EE) for all profile cases of (x, y)
    EE <- (pragr/prsum)

    ### Conditional expectations (EkU) for all profile cases of (x, y)
    m = min(x,y)
    num <- sum(sapply(0:m, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)}))
    den <- dbp(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    EE1U <- EE[1] * num/den

    num <- sum(sapply(0:x, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1)}))
    den <- dbp(x=x, y=0, m0 = 0, m1 = m0 + m1, m2 = 0)
    EE2U <- EE[2] * num/den

    num <- sum(sapply(0:y, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = y - u, lambda = m2)}))
    den <- dbp(x=0, y=y, m0 = 0, m1 = 0, m2 = m0 + m2)
    EE3U <- EE[3] * num/den

    EE4U <- EE[4] * m0

    m0.den <- 1
    m0.num <- (EE1U + EE2U + EE3U + EE4U)
    m1.den <- sum(EE[1:2])
    m1.num <- m1.den * x - EE1U - EE2U
    m2.den <- sum(EE[c(1,3)])
    m2.num <- m2.den * y - EE1U - EE3U

    #print(c(ee,eu,ev,ew))   ####debug
    return(c(ee, m0.den, m0.num, m1.den, m1.num, m2.den, m2.num))
  }
  fun.cond.exp <- Vectorize(fun.cond.exp.a)

  # M-step
  param.update <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    result <- fun.cond.exp(x, y, p1, p2, p3, p4, m0, m1, m2)
    #print(result)  ####debug
    result[result < 0] <- 0
    result <- apply(result,1,mean)
    result2 <- result[1:4]
    result2[5] <- result[6]/result[5]
    result2[6] <- result[8]/result[7]
    result2[7] <- result[10]/result[9]
    result2[is.na(result2)] <- as.numeric(c(p1, p2, p3, p4, m0, m1, m2))[is.na(result2)]
    return(result2)
  }

  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- mme.bzip.b(xvec=xvec, yvec=yvec)
    if (sum(is.na(initial))>0) {
      initial <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
      initial <- initial/sum(initial)      # relative freq
      initial[5:7] <- c(0.001, sum.x.y/len/(1-initial[4]))
    }
    if (min(sum.x.y) < 5 & min(sum.x.y) > 0) {
      initial[2:3] <- (initial[4] - 1e-10) * sum.x.y/sum(sum.x.y)
      initial[4] <- 1e-10}
    #print(initial)
  }
  initial <- pmax(initial, rep(1e-10, 7))
  # cat("initial = ", paste0(initial, collapse = ", "), "\n")

  # Repeat
  iter = 0
  param = initial
  if (showFlag) {print(c("iter", paste0("p",1:4), paste0("mu",0:2)))}
  repeat {
    if (iter >= maxiter) { warning("EM exceeded maximum number of iterations")
      param <- data.frame(matrix(NA,1,7))
      names(param) <- c("p1", "p2", "p3", "p4", "mu0", "mu1", "mu2")
      return(param)}

    iter = iter + 1
    # print(lik(vec, pp=param[1:4], m0=param[5], m1=param[6], m2=param[7])) # debug
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, p1 = param[1], p2 = param[2], p3 = param[3], p4 = param[4], m0=param[5], m1=param[6], m2=param[7])
    # cat("param = ", paste0(param, collapse = ", "), "\n param.old = ", paste0(param.old, collapse = ", "), "\n")
    if (showFlag) {print(c(iter, round(param,5), lik.bzip.b(xvec,yvec,param=param) ))} #lik.bzip(xvec, yvec, param)
    if (max(abs(param - param.old)) <= tol) {
      param <- data.frame(matrix(param,1,7))
      names(param) <- c("p1", "p2", "p3", "p4", "mu0", "mu1", "mu2")
      return(param)
      break
    }
  }
  return(param)
}

# Comparing algorithms
if (FALSE) {
  tt(1)
  bzip.b(extractor(1),extractor(2), showFlag=TRUE)   #EM (bzip.b)
  # p1         p2          p3        p4       mu0     mu1     mu2
  # 0.002569217 0.03868079 0.006423043 0.952327 3.272441e-09 16.72727 3.61411
  # 1285 iterations, 2.67 mins (lik=-487.8919)
  # lik.bzip.b(extractor(1),extractor(2), c(0.002569217, 0.03868079, 0.006423043, 0.952327, 3.272441e-09, 16.72727, 3.61411)) #-487.8919
  tt(2)

  tt(1)
  bzip.b(extractor(1),extractor(2), showFlag=TRUE)   #EM2 (bzip.b)
  # p1         p2          p3        p4       mu0     mu1     mu2
  # 0.002569209 0.03868079 0.006423021 0.952327 7.816525e-09 16.72727 3.614232
  # 10 iterations, 2.45 mins  (lik=-487.8919)
  #lik.bzip.b(extractor(1),extractor(2), c(0.002569209, 0.03868079, 0.006423021, 0.952327, 7.816525e-09, 16.72727, 3.614232)) #-487.8919
  tt(2)

  tt(1)
  bzip(extractor(1),extractor(2), showFlag=TRUE)      #EM (bzip)
  #   pp          mu0      mu1       mu2
  # 0.9525 0.0001888937 14.52612 0.6840215
  # 1122 iterations, 58 sec
  tt(2)

  tt(1)
  bzip.old(extractor(1),extractor(2), showFlag=TRUE)   #Direct maximization
  # pp          mu0      mu1       mu2
  # 0.9525 7.699164e-09 14.52631 0.6842103
  # 3 iterations, 0.53 sec
  tt(2)
}


### xy.reduced function
#
# reduce <- function(xvec, yvec) {
#   xy.reduced <- as.data.frame(table(xvec,yvec))
#   names(xy.reduced) <- c("x", "y","freq")
#   xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
#   xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
#   xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
#   xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
#   return(xy.reduced)
# }
