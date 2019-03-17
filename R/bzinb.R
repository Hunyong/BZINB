#' @import Rcpp BH


dbzinb <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE) {
  dxy <- dbnb(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a2, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(ifelse(log, log(result), result))
}
dbzinb.vec <- Vectorize(dbzinb)

lik.bzinb <- function(x, y, param) sum(log(do.call(dbzinb.vec, c(list(x, y), as.list(param)))))

#' @useDynLib bzinb
#' @export
rbzinb <- function(n, a0, a1, a2, b1, b2, p1, p2, p3, p4, param=NULL) {
  if (!is.null(param)) {a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
                        p1 = param[6]; p2 = param[7]; p3 = param[8]; p4 = param[9]
  }

  rmat <- matrix(rgamma(n*3, shape = c(a0, a1, a2), rate = 1/b1), n, 3, byrow=TRUE)
  rmat2 <- rmat
  rmat2[,3] <- rmat2[,1] + rmat2[,3]
  rmat2[,2] <- rmat2[,1] + rmat2[,2]
  rmat2 <- rmat2[,2:3]
  rmat2[,2] <- rmat2[,2]*b2/b1
  uv <- matrix(rpois(n*2, rmat2), n, 2)

  E <- t(rmultinom(n, 1, c(p1, p2, p3, p4)))
  z <- cbind(E[,1]+E[,2], E[,1]+E[,3])

  xy <- uv * z
  colnames(xy) <- c("x", "y")

  return(xy)
}

trueCor <- function(a0, a1, a2, b1, b2) a0 * sqrt(b1 * b2 /(b1+1) /(b2+1) /(a0 + a1) /(a0 + a2))

if (FALSE) {
  dbzinb(1,1,1,1,1,1,.5,.25,.25,.25,.25)
  tmp <- sapply(0:50, function(r) sapply (0:50, function(s) dbzinb(s,r,1,1,1,1,.5,.25,.25,.25,.25)))
  sum(tmp)
}

### 2.em

#' @useDynLib bzinb
#' @export
dbzinb.expt.vec <- function(xvec, yvec, freq, n = length(freq), a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  .check.input()
  
  result <- rep(0, 12)
  dbzinb_expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  result
}

# step-by-step expt calculation for information matrix
dbzinb.expt.se <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  result <- rep(0, 12)
  dbzinb_expt(x, y, freq = 1, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  result
}
dbzinb.expt.se.mat <- Vectorize(dbzinb.expt.se)

dbzinb.expt.mat <- function(xvec, yvec, freq, n = sum(freq), a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  .check.input()
  
  result <- rep(0, 12)
  dbzinb_expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  result
}

if (FALSE) {
  tmp <- dbzinb.expt.vec(c(1,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
  tmp <- dbzinb.expt.vec(c(0,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
  tmp <- dbzinb.expt.vec(extractor(1),extractor(2),1,1,1,1,2,.25,.25,.25,.25)
  t(tmp)[21:40,]
  dbzinb.expt.vec(c(10,1,2),c(10,1,1), 1.193013282, 0.003336139, 0.002745513, 3.618842924, 3.341625901, .25,.25,.25,.25)
}

abp.names <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4") # global variable
expt.names <- c("lik", "ER0", "ER1", "ER2", "ElogR0", "ElogR1", "ElogR2", "EE1", "EE2", "EE3", "EE4", "EV")
# EM with booster
# maxiter control added, output =param + lik + #iter
# Mar 15, 2018: Print pureCor instead of cor
bzinb.base <- function (xvec, yvec, initial = NULL, tol = 1e-8, maxiter=50000, showFlag=FALSE,
                             debug = FALSE, SE = TRUE, vcov = FALSE) {
  .check.input()
  if (!SE & vcov) {warning("To get covariance matrix (vcov), SE should be TRUE. The covariance matrix will not be obtained.")}
  require(rootSolve)
  if (debug) {showFlag=TRUE}
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n.reduced <- as.integer(length(xy.reduced$freq))
  # XMAX <- max(xy.reduced$x)
  # YMAX <- max(xy.reduced$y)
  # XYMAX(XMAX, YMAX)
  
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0, if (SE) {rep(NA, 11)}))} # 9 params, lik, iter, pureCor, and 11 SE's
  #print(xy.reduced)
  info <- if (SE) {matrix(0, ncol = 8, nrow = 8, dimnames = list(abp.names[-9], abp.names[-9]))} else {0}
#   if (SE) { #internal SE function (only use after param is obtained.)
#   
#     .se <- function(param) {
# # print(names(c(list(xy.reduced$x, xy.reduced$y, xy.reduced$freq), as.list(param))))
# # print(c("freq", xy.reduced$freq))
# # print(c("param", param))
# # print(as.list(param))
#       # se <- bzinb.se(xvec = xy.reduced$x, yvec = xy.reduced$y, freq = xy.reduced$freq, 
#       #                 a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4],  b2 = param[5], 
#       #                 p1 = param[6], p2 = param[7], p3 = param[8], p4 = param[9])
#       se <- do.call(bzinb.se, c(list(xy.reduced$x, xy.reduced$y, freq = xy.reduced$freq), as.list(param)))
#       names(se) <- paste0("se.", c(abp.names, "rho", "logit.rho"))
#       se
#   }
# }

  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)|is.na(s2.y)) {s2.x <- s2.y <- 1}
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n

    initial <- rep(NA,9)
    names(initial) <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4")
    initial[4] <- s2.x /ifelse(xbar==0,1e-4, xbar) #%>% print
    initial[5] <- s2.y /ifelse(ybar==0,1e-4, ybar) #%>% print
    initial[2:3] <- c(xbar,ybar)/pmax(initial[4:5], c(0.1,0.1)) #%>% print
    initial[1] <- min(initial[2:3]) * abs(cor.xy) #%>% print
    initial[2:3] <-  initial[2:3] - initial[1] #%>% print

    initial[6:9] <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
    initial[6:9] <- initial[6:9]/sum(initial[6:9])      # relative freq
    initial <- pmax(initial, 1e-5)
    if(is.na(sum(initial))) { initial[is.na(initial)] <- 1}
  # print(initial) ###
  } else {
    names(initial) <- abp.names
  }
  
  iter = 0L
  param = initial
  lik = -Inf
  expt = setNames(as.double(rep(0, 12)), expt.names)
  # print(c(param, showFlag, iter))
  em(param = param, xvec = xy.reduced$x, yvec = xy.reduced$y, 
     freq = xy.reduced$freq, n = n.reduced, expt = expt, info = info,
     se = as.integer(SE), iter = as.integer(iter), 
     maxiter = as.integer(maxiter), tol = as.double(tol), 
     showFlag = as.integer(showFlag))

  # underlying correlation (rho)
  rho <- param[1]/sqrt((param[1] + param[2]) * (param[1] + param[3])) *
    sqrt(param[4] * param[5] /(param[4] + 1) /(param[5] + 1))
  logit.rho <- qlogis(rho)
  
  if (SE) {
    qr.info <- try(qr(info))
    if (class(qr.info) == "try-error") {
      warning ("The information matrix has NA/NaN/Inf and thus the standard error is not properly estimatd.")
      std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
      cov.mat <- NA
    } else if (qr(info)$rank < 8) {
      warning ("The information matrix is (essentially) not full rank, and thus the standard error is not reliable.")
      std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
      cov.mat <- NA
    } else {
      cov.mat <- try(solve(info))
      if (class(cov.mat) == "try-error") {
        std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
        cov.mat <- NA
      } else {
        # variance of p4 hat
        var.p4 <- sum (cov.mat[6:8, 6:8]) # = sum_i,j cov(pi, pj)
        
        # variance of rho hat
        # rho <- a0/sqrt((a0 + a1) * (a0 + a2)) *sqrt(b1 *b2 /(b1 + 1) /(b2 + 1))
        
        # d.g <- rho * c(1/a0 - 1/{2*(a0 + a1)} - 1/{2*(a0 + a2)}, - 1/{2*(a0 + a1)}, - 1/{2*(a0 + a2)},
        #                1/{2 *b1 *(b1 + 1)}, 1/{2 *b2 *(b2 + 1)})
        d.g <- rho * c(1/param[1] - 1/{2*(param[1] + param[2])} - 1/{2*(param[1] + param[3])}, 
                       - 1/{2*(param[1] + param[2])}, 
                       - 1/{2*(param[1] + param[3])},
                       1/{2 *param[4] *(param[4] + 1)}, 
                       1/{2 *param[5] *(param[5] + 1)})
        
        var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g
        
        # variance of logit(rho hat)
        var.logit.rho <- var.rho / rho^2 / (1-rho)^2
        # std.param = sqrt(c(setNames(diag(cov.mat), abp.names[1:8]), 
        #                    p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
        std.param = sqrt(c(diag(cov.mat), 
                           p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
      }
    } 
    
  }

  # print(c(param, showFlag, iter))  
  result <- list(call = call,
                 rho = matrix(c(rho, logit.rho, if(SE) std.param[c("rho", "logit.rho")] else rep(NA, 2)),
                              ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                 coefficients = matrix(c(param, if(SE) std.param[1:9] else rep(NA, 9)),
                                       ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                 lik = expt[1],
                 iter = iter)
  if (SE & vcov) {
    result$info = info
    result$vcov = cov.mat
  }
  return(result)
}
#' @useDynLib bzinb
#' @export
bzinb <- function(xvec, yvec, ...) {
  .check.input()
  
  # if (!is.integer(xvec) | !is.integer(yvec)) stop("xvec and yvec should be integers.")
  # nonnegative
  # len(xvec) == len(yvec)
  # any(is.na(xvec))
  xvec = as.integer(round(xvec, digits = 0))
  yvec = as.integer(round(yvec, digits = 0))
  call <- match.call()
  result <- try(bzinb.base(xvec,yvec, ...))
  if (class(result)=="try-error") {
    result <- list(rho = matrix(rep(NA, 4),
                                ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                   coefficients = matrix(rep(NA, 9),
                                         ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                   lik = NA,
                   iter = NA)
    if (SE & vcov) {
      result$info = NA
      result$vcov = NA
    }
  }
  return(result)
}

if (FALSE) {
  # param for pair 1 and 2
  param <- c(1,1,1,1,1,  9.486775e-01,  1.893068e-02,  1.847954e-02,  1.391224e-02)
  set.seed(1)
  tmp <- rbzinb (800, param=param)
  table(tmp[,1], tmp[,2])

  microbenchmark::microbenchmark(
    bzinb(tmp[,1], tmp[,2], maxiter=1000, SE = FALSE), 
    bzinb.old(tmp[,1], tmp[,2], maxiter=1000, SE = FALSE),
    time = 5)
  bzinb(tmp[,1], tmp[,2], maxiter=100)
  
  bzinb(tmp[,1], tmp[,2], maxiter=20000, showFlag=F, SE = TRUE, initial = initialize(tmp[,1], tmp[,2]))
  
  library(devtools); load_all()
  set.seed(11)
  tmp <- rbzinb (800, param=param)
  # bzinb.old(tmp[,1], tmp[,2], maxiter=10, showFlag=T)
  # bzinb(tmp[,1], tmp[,2], maxiter=10, showFlag=T)
  # iter = 2059, likelihood = -2987.08
  bzinb(tmp[,1], tmp[,2], maxiter=2, showFlag=F, initial = c(a0 = 0.865977, a1 = 1.00028, a2 = 1.11896, b1 = 1.11955, b2 = 1.01823, p1 = 0.962907, p2 = 0.00529951, p3 = 1.62319e-27, p4 = 0.0317939))
  
  maxiter = 60000; a1 <- Sys.time();( b1 <- bzinb.old(tmp[,1], tmp[,2], maxiter=maxiter, showFlag=F)); a2 <- Sys.time(); cat("New method\n"); (b2 <- bzinb(tmp[,1], tmp[,2], maxiter=maxiter, showFlag=F)); a3 <- Sys.time();  a2-a1; a3-a2; cat("diff(lik) = ", b2[10] - b1[10], ", iters (old, new) ", b2[11], b1[11])
}

#' @useDynLib bzinb
#' @export
idgamma <- function(y) {
  x = 1
  inv_digamma(x = x, y)
  c(x, y)
}

#' @useDynLib bzinb
#' @export
opt <- function(b1, expt, a) {
  opt_b1(b1, expt, a)
  return(c(a, b1))
}


#' @useDynLib bzinb
#' @export
bzinb.se <- function(xvec, yvec, param = NULL, ...) {
  .check.input()
  if (any(!is.finite(param))) {return(rep(NA, 10))}
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n.reduced <- as.integer(length(xy.reduced$freq))
  
  if (is.null(param)) {
    warning("param was not provided. It will be estimated.")
    est <- bzinb(xvec, yvec, ...)
    param <- unlist(est$coefficients[,1])
    iter <- est$iter
  } else {
    iter <- NA
  }
  
  if (any(is.na(param))) {
    warning("params include NA.")
    result <- list(rho = matrix(rep(NA, 4),
                                ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                   coefficients = matrix(rep(NA, 18),
                                         ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                   lik = NA, iter = NA, info = NA, vcov = NA)
    return(result)
  }
  
  rho <- param[1]/sqrt((param[1] + param[2]) * (param[1] + param[3])) * 
    sqrt(param[4] *param[5] /(param[4] + 1) /(param[5] + 1))
  logit.rho <- qlogis(rho)
  
  expt = setNames(as.double(rep(0, 12)), expt.names)
  s_i = setNames(as.double(rep(0, 8)), abp.names[-9])
  info <- matrix(0, ncol = 8, nrow = 8, dimnames = list(abp.names[-9], abp.names[-9]))
  dbzinb_expt_vec (xvec = xy.reduced$x, yvec = xy.reduced$y, freq = xy.reduced$freq, n = n.reduced, 
                    a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], 
                    p1 = param[6], p2 = param[7], p3 = param[8], p4 = param[9], 
                    expt = expt, s_i = s_i, info = info, se = 1)
  
  
  # inverse of info
  qr.info <- try(qr(info))
  if (class(qr.info) == "try-error") {
    warning ("The information matrix has NA/NaN/Inf and thus the standard error is not properly estimatd.")
    std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
    cov.mat <- NA
  } else if (qr(info)$rank < 8) {
    warning ("The information matrix is (essentially) not full rank, and thus the standard error is not reliable.")
    std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
    cov.mat <- NA
  } else {
    cov.mat <- try(solve(info))
    if (class(cov.mat) == "try-error") {
      std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
      cov.mat <- NA
    } else {
      # variance of p4 hat
      var.p4 <- sum (cov.mat[6:8, 6:8]) # = sum_i,j cov(pi, pj)
      
      # variance of rho hat
      d.g <- rho * c(1/param[1] - 1/{2*(param[1] + param[2])} - 1/{2*(param[1] + param[3])}, 
                     - 1/{2*(param[1] + param[2])}, 
                     - 1/{2*(param[1] + param[3])},
                     1/{2 *param[4] *(param[4] + 1)}, 
                     1/{2 *param[5] *(param[5] + 1)})
      
      var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g
      
      # variance of logit(rho hat)
      var.logit.rho <- var.rho / rho^2 / (1-rho)^2
      # std.param = sqrt(c(setNames(diag(cov.mat), abp.names[1:8]), 
      #                    p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
      std.param = sqrt(c(diag(cov.mat), 
                         p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
    }
  } 
  
  result <- list(rho = matrix(c(rho, logit.rho, std.param[c("rho", "logit.rho")]),
                              ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                 coefficients = matrix(c(param, std.param[1:9]),
                                       ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                 lik = expt[1],
                 iter = iter,
                 info = info,
                 vcov = cov.mat)
  return(result)
}
