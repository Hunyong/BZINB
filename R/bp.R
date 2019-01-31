##########################################################################################
## 1. Bivariate Poisson
##########################################################################################

#' Density of bivariate Poisson distribution (BP)
#'
#' This function calculates the density of BP model
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
#' @return the BP density at the data point
#'
#' @examples
#'
#' dBP(2, 2, 0, 1, 1)
#' dBP(100, 100, 0, 1, 1)
#'
#' # Sum of joint density should equal the marginal density
#' sum(sapply(0:100, function(s) dBP(s,2,1,2,2)))
#' dpois(2, 3)
#'
#' # joint(BP) > pois^2 (around diagonal)
#' dBP(2,2,1,2,2)
#' dpois(2,3)^2
#'
#' # joint(BP) < pois*pois (off-diagonal)
#' dBP(1,3,1,2,2)
#' dpois(1,3)*dpois(3,3)
#'
dBP <- function(x, y, m0, m1, m2, log = FALSE, max = 500) {
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
dBP.vec <- Vectorize(dBP)


ML.BP <- function(xvec, yvec, tol = 1e-6) {
  xvec <- as.numeric(xvec); yvec <- as.numeric(yvec)
  len <- length(xvec)
  # xbar <- mean(xvec)
  # ybar <- mean(yvec)
  vec <- cbind(xvec,yvec)
  bar <- apply(vec,2,mean)
  m <- min(bar)
  lik.BP <- function(a) (-sum(apply(vec, 1, dBP, m0 = a, m1 = bar[1] - a, m2 = bar[2] - a, log = TRUE)))
  # result <- nlminb(m, lik, lower = 0)
  if (m == 0) {  # when either is all zero, then mu0 is automatically 0.
    result = data.frame(par = 0, value = NA)  # lik not actually NA but can be obtained if necessary
  } else {
    result <- optim(m, lik.BP, lower = 0, upper = m, method="Brent")
  }
  result <- data.frame(mu0 = result$par, mu1 = bar[1] - result$par, mu2 = bar[2] - result$par, likelihood = -result$value)
  rownames(result) <- NULL
  return(result)
}
if (FALSE) { # example
  ML.BP(x,y)
  ML.BP(x,z)
  set.seed(1000); a1 <- rpois(20,1); a2 <- rpois(20,2); a3 <- rpois(20,3)
  ML.BP(a1+a2, a1+a3)
}


##########################################################################################
## 2. Bivariate ZIP
##########################################################################################

## 2.1 basic functions for BvZIP
dBZIP <- function(x, y = NULL, pp, m0, m1, m2, log = FALSE) {
  fxy <- (1-pp) * dBP (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) + ifelse((x == 0 & y == 0), pp,0)
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dBZIP.vec <- Vectorize(dBZIP)
lik.BZIP <- function(x, y, param, pp = NULL, m0  = NULL, m1  = NULL, m2  = NULL){
  if (is.null(pp)|is.null(m0)|is.null(m1)|is.null(m2)) {
    pp = param[1]; m0 = param[2]; m1 = param[3]; m2 = param[4]
  }
  sum(dBZIP.vec(x, y, pp = pp, m0 = m0, m1 = m1, m2 = m2, log=TRUE))
}
rBZIP <- function(n, pp, m0, m1, m2) {
  z <- rbinom(n,1,pp)
  common.tmp <- rpois(n,m0)
  x.tmp <- rpois(n,m1)
  y.tmp <- rpois(n,m2)
  x <- (common.tmp + x.tmp)*(1-z)
  y <- (common.tmp + y.tmp)*(1-z)
  return(data.frame(x=x, y=y))
}
rBZIP.vec <- Vectorize(rBZIP)

## 2.2 param estimation
# formal EM algorithm
ML.BZIP <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
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
          #prob <- sapply(0:m, function(u) {dBP(x = u, y=0, m0 = 0, m1 = m0, m2 = 0) * dBP(x = x-u, y=0, m0 = 0, m1 = m1, m2 = 0) *
          #    dBP(x = 0, y= y-u, m0 = 0, m1 = 0, m2 = m2)*10^100})
          # used dBP instead of simply dpois to handle large numbers (when mu=1, x = 172, dpois=0 not small number)
          # and adjust by 10^100 (num/den cancel out adj.)
        }
        #print(c(1,prob)) # debug
        #print(prob)
      }
      EU <- sum((0:m)*prob)/sum(prob)
      cond.expt <- c(0, 0, x, y) + EU *c(0, 1, -1, -1)
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
    if (showFlag) {print(c(iter, round(param,5), lik.BZIP(xvec, yvec, param)))}
    if (max(abs(param - param.old)) <= tol) {
      param <- data.frame(matrix(param,1,4))
      names(param) <- c("pp", "mu0", "mu1", "mu2")
      return(param)
      break
    }
  }
  return(param)
}

# ML.BZIP examples
if (FALSE) {
  ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, initial = rep(0,4), tol = 1e-10)
  ML.BZIP(vec[,1], vec[,2], showFlag=TRUE, tol = 1e-8)
  ML.BZIP(extractor(1), extractor(2), initial = rep(0,4), showFlag=TRUE)
  ML.BZIP(rep(0,100), rep(0,100), initial = rep(0,4), showFlag=TRUE)
  # For all 0 pairs, identifiability issue: same likelihood for (pi=1) and (pi=0)
  ML.BZIP(c(1,rep(0,99)), rep(0,100), initial = rep(0,4), showFlag=TRUE)
  ML.BZIP(c(1,rep(0,799)), rep(0,800), initial = rep(0,4), showFlag=TRUE)
  ML.BZIP(c(1,rep(0,799)), rep(0,800), initial = c(799/800,0,0,0), showFlag=TRUE) # not converging. initial should be all zero
}


moment.BZIP <- function(pp, m0, m1, m2) {
  MEAN.BP <- (m0 + c(m1, m2))
  VAR.BP <- diag(MEAN.BP); VAR.BP[c(2,3)] <- m0
  MEAN <- (1-pp)*MEAN.BP
  VAR <- (1-pp)*VAR.BP + pp*(1-pp) * MEAN.BP %o% MEAN.BP
  COR <- VAR[2] / sqrt(prod(diag(VAR)))
  return(list(mean = MEAN, var = VAR, cor = COR))
}


##########################################################################################
## 3. Bivariate ZIP.B: General BZIP1 with marginal ZIP condition (6params)
##########################################################################################

## 3.1 basic functions for BvZIP.B
dBZIP.B <- function(x, y = NULL, p1, p2, p3, p4 = 1 - p1 - p2 - p3, m0, m1, m2, log = FALSE) {
  fxy <- p1 * dBP (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) +
    p2 * {if (y == 0) dBP (x=x, y=y, m0 = 0, m1 = m0 + m1, m2 = 0) else 0} +
    p3 * {if (x == 0) dBP (x=x, y=y, m0 = 0, m1 = 0, m2 = m0 + m2) else 0} +
    p4 * {if (x + y == 0) 1 else 0}
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dBZIP.B.vec <- Vectorize(dBZIP.B)
lik.BZIP.B <- function(x, y, param, p1 = NULL, p2 = NULL, p3 = NULL, p4 = NULL, m0  = NULL, m1  = NULL, m2  = NULL){
  if (is.null(p1)|is.null(p2)|is.null(p3)|is.null(p4)|is.null(m0)|is.null(m1)|is.null(m2)) {
    p1 = param[1]; p2 = param[2]; p3 = param[3]; p4 = param[4];
    m0 = param[5]; m1 = param[6]; m2 = param[7]
  }
  sum(dBZIP.B.vec(x, y, p1 = p1, p2 = p2, p3 = p3, p4 = p4,
                  m0 = m0, m1 = m1, m2 = m2, log=TRUE))
}
rBZIP.B <- function(n, p1, p2, p3, p4, m0, m1, m2) {
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
MME.BZIP.B <- function(xvec, yvec) {
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

ML.BZIP.B <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE, maxiter = 200) {
  # counter, likelihood, param.prev for recursive function
  # initial guess manipulation (for small counts pi4=0)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(data.frame(p1 = 0, p2 = 0, p3 = 0, p4 = 1,  mu0 = 0, mu1 = 0, mu2 = 0))} # everything is zero ==> nonestimable, set pi = 0

  # E-step
  fun.cond.exp.a <- function(x, y, p1, p2, p3, p4, m0, m1, m2) {
    PrA <- p1 * dBP(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    PrB <- p2 * (y==0) * dBP(x=x, y=0, m0 =0, m1 = m0 + m1, m2 = 0)
    PrC <- p3 * (x==0) * dBP(x=0, y=y, m0 =0, m1 = 0, m2 = m0 + m2)
    PrD <- p4 * (x + y ==0)
    PrAgr <- c(PrA, PrB, PrC, PrD)
    PrSum <- sum(PrAgr)

    ### Conditional expectations (EE) for all profile cases of (x, y)
    EE <- (PrAgr/PrSum)

    ### Conditional expectations (EkU) for all profile cases of (x, y)
    m = min(x,y)
    num <- sum(sapply(0:m, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)}))
    den <- dBP(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    EE1U <- EE[1] * num/den

    num <- sum(sapply(0:x, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1)}))
    den <- dBP(x=x, y=0, m0 = 0, m1 = m0 + m1, m2 = 0)
    EE2U <- EE[2] * num/den

    num <- sum(sapply(0:y, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = y - u, lambda = m2)}))
    den <- dBP(x=0, y=y, m0 = 0, m1 = 0, m2 = m0 + m2)
    EE3U <- EE[3] * num/den

    EE4U <- EE[4] * m0

    m0.den <- 1
    m0.num <- (EE1U + EE2U + EE3U + EE4U)
    m1.den <- sum(EE[1:2])
    m1.num <- m1.den * x - EE1U - EE2U
    m2.den <- sum(EE[c(1,3)])
    m2.num <- m2.den * y - EE1U - EE3U

    #print(c(EE,EU,EV,EW))   ####debug
    return(c(EE, m0.den, m0.num, m1.den, m1.num, m2.den, m2.num))
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
    initial <- MME.BZIP.B(xvec=xvec, yvec=yvec)
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
    if (showFlag) {print(c(iter, round(param,5), lik.BZIP.B(xvec,yvec,param=param) ))} #lik.BZIP(xvec, yvec, param)
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
  ML.BZIP.B(extractor(1),extractor(2), showFlag=TRUE)   #EM (BZIP.B)
  # p1         p2          p3        p4       mu0     mu1     mu2
  # 0.002569217 0.03868079 0.006423043 0.952327 3.272441e-09 16.72727 3.61411
  # 1285 iterations, 2.67 mins (lik=-487.8919)
  # lik.BZIP.B(extractor(1),extractor(2), c(0.002569217, 0.03868079, 0.006423043, 0.952327, 3.272441e-09, 16.72727, 3.61411)) #-487.8919
  tt(2)

  tt(1)
  ML.BZIP.B2(extractor(1),extractor(2), showFlag=TRUE)   #EM2 (BZIP.B)
  # p1         p2          p3        p4       mu0     mu1     mu2
  # 0.002569209 0.03868079 0.006423021 0.952327 7.816525e-09 16.72727 3.614232
  # 10 iterations, 2.45 mins  (lik=-487.8919)
  #lik.BZIP.B(extractor(1),extractor(2), c(0.002569209, 0.03868079, 0.006423021, 0.952327, 7.816525e-09, 16.72727, 3.614232)) #-487.8919
  tt(2)

  tt(1)
  ML.BZIP(extractor(1),extractor(2), showFlag=TRUE)      #EM (BZIP)
  #   pp          mu0      mu1       mu2
  # 0.9525 0.0001888937 14.52612 0.6840215
  # 1122 iterations, 58 sec
  tt(2)

  tt(1)
  ML.BZIP.old5(extractor(1),extractor(2), showFlag=TRUE)   #Direct maximization
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
