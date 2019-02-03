# BvZINB4: BvZINB3 + varying zero inflation parameters
#' library(rootSolve)

dBvZINB5 <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE) {
  dxy <- dBvNB3(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a2, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(ifelse(log, log(result), result))
}
dBvZINB5.vec <- Vectorize(dBvZINB5)

lik.BvZINB5 <- function(x, y, param) sum(log(do.call(dBvZINB5.vec, c(list(x, y), as.list(param)))))

#' @useDynLib bzinb
#' @export
rBvZINB5 <- function(n, a0, a1, a2, b1, b2, p1, p2, p3, p4, param=NULL) {
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
  dBvZINB5(1,1,1,1,1,1,.5,.25,.25,.25,.25)
  tmp <- sapply(0:50, function(r) sapply (0:50, function(s) dBvZINB5(s,r,1,1,1,1,.5,.25,.25,.25,.25)))
  sum(tmp)
}

### 2.EM
# `dBvZINB5.Expt()` written in R was replaced by `dBvZINB5.Expt.cpp()`

#' @useDynLib bzinb
#' @export
dBvZINB5.Expt.vec <- function(xvec, yvec, freq, n = length(freq), a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  result <- rep(0, 12)
  dBvZINB_Expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  result
}

# step-by-step Expt calculation for information matrix
dBvZINB5.Expt.se <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  result <- rep(0, 12)
  dBvZINB_Expt(x, y, freq = 1, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  result
}
dBvZINB5.Expt.se.mat <- Vectorize(dBvZINB5.Expt.se)

dBvZINB5.Expt.mat <- function(xvec, yvec, freq, n = sum(freq), a0, a1, a2, b1, b2, p1, p2, p3, p4) {
  result <- rep(0, 12)
  dBvZINB_Expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  result
}

if (FALSE) {
  tmp <- dBvZINB5.Expt.vec(c(1,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
  tmp <- dBvZINB5.Expt.vec(c(0,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
  tmp <- dBvZINB5.Expt.vec(extractor(1),extractor(2),1,1,1,1,2,.25,.25,.25,.25)
  t(tmp)[21:40,]
  dBvZINB5.Expt.vec(c(10,1,2),c(10,1,1), 1.193013282, 0.003336139, 0.002745513, 3.618842924, 3.341625901, .25,.25,.25,.25)
}

# EM with booster
# maxiter control added, output =param + lik + #iter
# Mar 15, 2018: Print pureCor instead of cor
ML.BvZINB5.base <- function (xvec, yvec, initial = NULL, tol=1e-8, maxiter=200, showFlag=FALSE,
                             showPlot=FALSE, cor.conv = FALSE, boosting=FALSE, debug = FALSE, SE = TRUE) {
  require(rootSolve)
  if (debug) {showFlag=TRUE}
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n.reduced <- length(xy.reduced$freq)
  XMAX <- max(xy.reduced$x)
  YMAX <- max(xy.reduced$y)
  
  
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0, if (SE) {rep(NA, 10)}))} # 9 params, lik, iter, pureCor, and 10 SE's
  #print(xy.reduced)
  if (SE) { #internal SE function (only use after param is obtained.)
    .se <- function(param) {
      se <- do.call(BZINB5.se, c(list(xy.reduced$x, xy.reduced$y, xy.reduced$freq), as.list(param)))
      names(se) <- paste0("se.", c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4", "rho", "logit.rho"))
      se
    }
  }

  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)|is.na(s2.y)) {s2.x <- s2.y <- 1}
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n

    initial <- rep(NA,9)
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
  }
# print(initial)
  booster <- function (param.matrix, xvec, yvec, n.cand = 10) {
    param.matrix[,6:9] <- qlogis(param.matrix[,6:9])  # logit transformation for probs
    param.matrix[,1:5] <- log(param.matrix[,1:5])  # log transformation for positives
    a <- param.matrix[1,]
    b <- param.matrix[5,]
    candidate <- matrix(b, byrow=TRUE, ncol=9, nrow = n.cand)
    index <- which((abs(b-a) > 1e-5) & is.finite(b) & is.finite(a))  # target param for grid search

    for (s in 1:n.cand) {
      candidate[s,index] <- b[index] + (b[index] - a[index]) * 3^(s-1)
    }
    candidate[,6:9] <- plogis(candidate[,6:9])  # back-transformation
    candidate[,6:9] <- candidate[,6:9]/ apply(candidate[,6:9],1,sum) # normalize
    candidate[,1:5] <- exp(candidate[,1:5])  # back-transformation for probs
    #print(candidate[,1:4])  #debug

    lik <- sapply(1:n.cand, function(s) {lik.BvZINB4(xvec, yvec, candidate[s,])})
    lik <- ifelse(is.infinite(lik), -Inf, lik)  # sometimes likelihood is inf which is nonsense. force it to -Inf
    if (sum(!is.finite(lik)) > 0) {
      return(cbind(candidate,lik)[1:max(min(which(!is.finite(lik)))-1,1),])
    } else {return(cbind(candidate,lik))}
  }

  # cor.trace <<- data.frame(iter=1, pureCor=1)
  iter = 0
  param = initial
  lik = Inf
  pureCor = 0
  boost = 0
  index = 1 # previous boosting index
  if (showPlot) {
    par(mfrow=c(2,1))
    par(mar=c(2,4,1,4))
  }
# cat(442)
  repeat {
    iter = iter + 1
    param.old <- param # saving old parameters
    if (debug) {lik.old <- lik} #debugging
    pureCor.old <- pureCor
    # updating
# cat(449)
    
    expt <- do.call(dBvZINB5.Expt.vec, c(list(xy.reduced$x, xy.reduced$y, xy.reduced$freq, n.reduced), as.list(param)))
    
tmp.expt <<- expt
# cat(453)
# cat("sum(expt): ", expt, "\n")
    # loglik = expt[1] * n
    delta <- expt[12] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    param[6:9] = expt[8:11]
# pi = E(Z)
# cat(457)
    opt.vec <- function(par.ab) {
      par.ab <- exp(par.ab)
      r1 <- sum(expt[2:4]) - sum(par.ab[1:3]) * par.ab[4]
      r2 <- expt[5:7] - digamma(par.ab[1:3]) - log(par.ab[4])
      # print(c(r1,r2)) ###
      return(c(r1,r2))
    }

    param.l <- log(param)
#expt %>% print
#param.l %>% print

    result <- try(multiroot(opt.vec, start=param.l[1:4])$root, silent=TRUE)
    if (class(result)=="try-error") {
      initial = rep(1,4)
      result <- multiroot(opt.vec, start = initial[1:4], rtol=1e-20)$root
    }
    param[1:4] <- exp(result)
    param[5]   <- param[4] * delta                                       # b2
    pureCor <- do.call(trueCor, as.list(param[1:5]))
# cat("param: ", param, "\n")
    if (debug) {
      lik <- lik.BvZINB5(xvec, yvec, param = param)  #debugging
      if (lik < lik.old) warnings("likelihood decreased!")
    }
    # cor.trace[iter,] <<- c(iter,pureCor)
    if (showPlot & (iter %% 20 == 0)) {
      span <- min(max(iter-200+1,1),101):iter
      span2 <- max(iter-100+1,1):iter
      yspan <- c(min(0.2, min(cor.trace[span,2]-0.05)),max (max(cor.trace[span,2])+0.05,0.4))
      yspan2 <- c(min(max(cor.trace[span2,2]) - 0.001, min(cor.trace[span2,2]-0.001)),max (max(cor.trace[span2,2])+0.001,0.4))
      plot(cor.trace[span,"iter"], cor.trace[span,"pureCor"], xlab="iteration", ylab="pureCorrelation", pch=".", col="blue", ylim = yspan)
      plot(cor.trace[span2,"iter"], cor.trace[span2,"pureCor"], xlab="iteration", ylab="pureCorrelation", pch=20, col="red")
    }

    # boosting
    if (boosting) {
      if (iter == 6 + boost*5) {  # Creating an empty matrix
        param.boost <- matrix(NA, nrow = 5, ncol = 9)
      }
      if (iter >= 6 + boost*5 & iter <= 10 + boost*5 ) {  # Storing last ten params
        param.boost[iter - (5 + boost*5),] <- param
      }
      if (iter == 10 + boost*5) {
        param.boost <- booster(param.boost, xvec, yvec, n.cand = min(max(5, index * 2),20))
        # print(dim(param.boost)); print(length(param.boost))
        if (showFlag) {print(param.boost)}

        if (is.null (dim(param.boost))) {
          param <- param.boost[1:9]
        } else {
          index <- which.max(param.boost[,10])
          param <- param.boost[index,1:9]
          if (showFlag) {print(paste0("Jump to the ",index, "th parameter"))}
        }
        boost <- boost + 1
      }
    }

    #print (expt) #####
    if (showFlag) {cat("iter ", iter, "parm:", round(param,4), if (debug) {c("D.lik=", round(lik - lik.old, 2))},
                       "lik=", expt[1] * n, "p.Cor=", pureCor, "\n")} #lik: lik of previous iteration
    if (maxiter <= iter) {
      lik <- lik.BvZINB5(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "rho")
      if (SE) result <- c(result, .se(param))
      return(result)
      }
    if (max(abs(param - param.old)) <= tol) {
      lik <- lik.BvZINB5(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "rho")
      if (SE) result <- c(result, .se(param))
      return(result)
    }
    if (cor.conv & abs(pureCor - pureCor.old) <= tol) {  # if pureCor is converged, then done!
      lik <- lik.BvZINB5(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "rho")
      if (SE) result <- c(result, .se(param))
      return(result)
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}
#' @useDynLib bzinb
#' @export
ML.BvZINB5 <- function(xvec, yvec, ...) {
  if (!is.integer(xvec) | !is.integer(yvec)) stop("xvec and yvec should be integers.")
  # nonnegative
  # len(xvec) == len(yvec)
  # any(is.na(xvec))
  
  result <- try(ML.BvZINB5.base(xvec,yvec, ...))
  if (class(result)=="try-error") {
    result <- c(rep(NA,5+4), NA, 0)
  }
  return(result)
}

if (FALSE) {
  # param for pair 1 and 2
  param <- c(1,1,1,1,1,  9.486775e-01,  1.893068e-02,  1.847954e-02,  1.391224e-02)
  set.seed(1)
  tmp <- rBvZINB5 (800, param=param)
  table(tmp[,1], tmp[,2])

  ML.BvZINB5(tmp[,1], tmp[,2], maxiter=20, showFlag=TRUE, SE = TRUE)
}
