# BvZINB4: BvZINB3 + varying zero inflation parameters
library(rootSolve)
source("R/F0001-base.R")
source("R/F0101-BvNB3.R")
source("R/F0102-BvZINB5-SE.R")

dBvZINB5 <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE) {
  dxy <- dBvNB3(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a2, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(ifelse(log, log(result), result))
}
dBvZINB5.vec <- Vectorize(dBvZINB5)

lik.BvZINB5 <- function(x, y, param) sum(log(do.call(dBvZINB5.vec, c(list(x, y), as.list(param)))))

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

### updated in 2018.12. Decreasing Likelihood problem solved.
dBvZINB5.Expt <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, debug = FALSE) {
  # Base density
  t1 = (b1 + b2 + 1) /(b1 + 1); t2 = (b1 + b2 + 1) /(b2 + 1)
  adj.A <- adj.B1 <- adj.C <- adj.sum <- 0
  l1 <- function(k, m, adjj=0) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) 
                                   + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) - adjj)
  l1.C <- function(k, m, adjj=0) exp(k *log(t1) + m *log(t2) - adjj)
  l1.B <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
  # l1.B to be updated several lines later depending on l2.B ~ l4.B
  l2.B <- exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj.B1) * p2 * ifelse(y==0, 1, 0)
  l3.B <- exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj.B1) * p3 * ifelse(x==0, 1, 0)
  l4.B <- p4 * ifelse(x + y == 0, 1, 0) * exp(adj.B1)
  
  #l2.A, l3.A added.
  l2.A <- function(k, adjj=0) exp( lgamma(x +a0 -k) + lgamma(k + a1) - lgamma(a0) - lgamma(x-k+1) - lgamma(a1) - lgamma(k+1) - adjj)
  l3.A <- function(m, adjj=0) exp( lgamma(y +a0 -m) + lgamma(m + a2) - lgamma(a0) - lgamma(y-m+1) - lgamma(a2) - lgamma(m+1) - adjj)
  
  
  # l1.AC For numerical stability use only.
  l1.AC <- function(k, m, adjj=0) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) 
                                      + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + k *log(t1) + m *log(t2) - adjj)
  
  
  # cat("l1.B ", l1.B,"\n")  
  if (l1.B < - 200 & log(l2.B + l3.B + l4.B) < 0) {
    if (debug) cat("adjustment activated for l1.B\n")     
    adj.B1 = ((-l1.B - 200) %/% 100) * 100 # prevent exp(l1.B) from being 0
    l1.B = l1.B + adj.B1
  }
  l1.B <- exp(l1.B) * p1
  if (debug)  cat("l1.B ", l1.B,"\n")  
  
  
  l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
  l2.A.mat <- sapply(0:x, l2.A, adjj = adj.A)  # %>% print
  l3.A.mat <- sapply(0:y, l3.A, adjj = adj.A)  # %>% print
  
  l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
  while (log(sum( l.A.mat)) > 250) {
    ### may have to be updated for l2.A.mat and l3.A.mat ###
    if (debug)  cat("adjustment activated for A.mat\n")     
    adj.A = adj.A + 200
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
  }
  while (log(sum( l.C.mat)) > 250) {
    if (debug)  cat("adjustment activated for C.mat\n")     
    adj.C = adj.C + 200
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
  }
  
  # print(l.C.mat)
  
  #   if (is.infinite(sum( l.A.mat))) {
  # cat("activated once")    
  #     adj.A = 200
  #     l.A.mat <- sapply(0:x, function(k) sapply(0:y, function(m) {l1(k =k, m = m) *exp(-adj.A)}))
  #     if (is.infinite(sum( l.A.mat))) { ## added for further adjustment
  # cat("activated twice")      
  #       adj.A = 500
  #       l.A.mat <- sapply(0:x, function(k) sapply(0:y, function(m) {l1(k =k, m = m) *exp(-adj.A)}))
  #     }  
  #   } #%>%print
  #adjustment is cancelled out for each Expectation, so can be ignored. But for the final likelihood it should be adjusted at the end.
  sum.AC <- sum(l.A.mat * l.C.mat)
  
  if (is.infinite(sum.AC)| log(sum.AC) > 200) {
    if (debug) cat("adjustment activated for AC.mat (too large)\n")
    adj.A = adj.A + 100
    adj.C = adj.C + 100
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
    sum.AC <- sum(l.A.mat * l.C.mat)
  } else if (log(sum.AC) < - 100) {
    if (debug) cat("adjustment activated for AC.mat (too small)\n")
    adj.A = adj.A - 200  # floor(log(sum(l.A.mat)/x/y)*2/3)
    adj.C = adj.C - 200
    l.A.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
    l.C.mat <- sapply(0:x, function(k) sapply(0:y, l1.C, k = k, adjj = adj.C))  # %>% print
    l.AC.mat <- sapply(0:x, function(k) sapply(0:y, l1.AC, k = k, adjj = adj.C + adj.A))
    sum.AC <- sum(l.AC.mat)
  }
  sum.A <- sum(l.A.mat)
  l.sum <- sum.AC * l1.B + sum.A * sum (l2.B +  l3.B +  l4.B) * exp(-adj.C)
  if (l.sum == 0) {
    adj.sum = -floor(log(sum.AC)*2/3 + log(l1.B)*2/3)
    if (debug) cat("adjustment activated for l.sum (adj = ", adj.sum, ")\n")
    l.sum <- sum.AC * exp(adj.sum) * l1.B + sum.A * (exp(adj.sum) * sum (l2.B +  l3.B +  l4.B)) * exp(-adj.C)
    ## paranthesis matters. sum.A = some number, exp(adj.sum) = almost inf, sum(l2.B +  l3.B +  l4.B) = 0, ...
    # Then without paranthesis, Inf * 0 = NaN,
    # But with paranthesis, c * (large number * 0) = c * 0 = 0
  }
  
  if (debug) {
    cat("sum.AC", sum.AC,"\n\n")
    cat("sum.A", sum.A,"\n\n")
    cat("sum(l.C.mat)", sum(l.C.mat),"\n\n")
    cat("l1.B", l1.B,"\n\n")
    cat("l2.B", l2.B,"\n\n")
    cat("l3.B", l3.B,"\n\n")
    cat("l4.B", l4.B,"\n\n")
    cat("l.sum ", l.sum, "\n")
  }
  # print(c(l.sum, log(l.sum))); print(l.A.mat); print(l.C.mat); print(c(l1.B, l2.B, l3.B, l4.B, adj.A)) #####
  # expectation components
  R0.E1 <- function(k, m) {x - k + y - m + a0}
  log.R0.E1 <- function(k, m) {digamma(x - k + y - m + a0)}
  log.R0.E2 <- function(k) {digamma(x - k + a0)}
  log.R0.E3 <- function(m) {digamma(y - m + a0)}
  R0.E1.B <- b1/(1 + b1 + b2)
  R0.E2.B <- b1/(1 + b1)
  R0.E3.B <- b1/(1 + b2)
  R0.E4.B <- b1
  
  R1.E1 <- function(k) {k + a1}
  log.R1.E1 <- function(k) {digamma(k + a1)}
  log.R1.E2 <- function(k) {digamma(k + a1)}
  R1.E1.B <- b1/(1 + b1)
  R1.E2.B <- b1/(1 + b1)
  R1.E3.B <- b1
  R1.E4.B <- b1
  
  R2.E1 <- function(m) {m + a2}
  log.R2.E1 <- function(m) {digamma(m + a2)}
  log.R2.E3 <- function(m) {digamma(m + a2)}
  R2.E1.B <- b1/(1 + b2)
  R2.E2.B <- b1
  R2.E3.B <- b1/(1 + b2)
  R2.E4.B <- b1
  
  R0.mat <- sapply (0:x, function(k) sapply(0:y, R0.E1, k=k))
  R0.mat <- R0.mat * l.A.mat
  R0.E <- sum(R0.mat * l.C.mat * exp(adj.sum) * l1.B * R0.E1.B)  / l.sum + 
    sum(R0.mat*(  l2.B * R0.E2.B
                  + l3.B * R0.E3.B
                  + l4.B * R0.E4.B)*exp(-adj.C + adj.sum)) / l.sum
  
  # cat("R0.E ", R0.E, "\n")  
  R1.mat <- t(matrix(sapply(0:x, R1.E1), x+1, y+1))
  R1.mat <- R1.mat * l.A.mat
  R1.E <- sum(R1.mat * l.C.mat * exp(adj.sum) * l1.B * R1.E1.B) / l.sum +  
    sum(R1.mat*(  l2.B * R1.E2.B 
                  + l3.B * R1.E3.B 
                  + l4.B * R1.E4.B)*exp(-adj.C + adj.sum)) / l.sum
  # cat("R1.E ", R1.E, "\n")  
  
  R2.mat <- matrix(sapply(0:y, R2.E1), y+1, x+1) #%>% print
  R2.mat <- R2.mat * l.A.mat
  R2.E <- sum(R2.mat * l.C.mat * exp(adj.sum) * l1.B * R2.E1.B) / l.sum +  
    sum(R2.mat*(  l2.B * R2.E2.B 
                  + l3.B * R2.E3.B 
                  + l4.B * R2.E4.B)*exp(-adj.C + adj.sum)) / l.sum
  # cat("R2.E ", R2.E, "\n")  
  
  log.R0.mat <- sapply(0:x, function(k) sapply(0:y, log.R0.E1, k=k))
  log.R0.mat <- l.A.mat * (log.R0.mat + log (R0.E1.B))
  log.R0.mat2 <- sapply(0:x, log.R0.E2)
  log.R0.mat2 <- l2.A.mat * (log.R0.mat2 + log (R0.E2.B))
  log.R0.mat3 <- sapply(0:y, log.R0.E3)
  log.R0.mat3 <- l3.A.mat * (log.R0.mat3 + log (R0.E3.B))
  
  log.R0.E <- 
    sum(log.R0.mat * l.C.mat) * exp(adj.sum - adj.C) * l1.B + 
    sum(log.R0.mat2 * l2.B)  * exp(adj.sum) + 
    sum(log.R0.mat3 * l3.B)  * exp(adj.sum) + 
    (digamma(a0) + log(b1)) * exp(adj.sum) * l4.B
  log.R0.E <- log.R0.E / l.sum
  
  log.R1.mat <- sapply(0:x, log.R1.E1)
  log.R1.mat2 <- log.R1.mat  # saving a vector form
  log.R1.mat2 <- l2.A.mat * (log.R1.mat2 + log (R1.E2.B))
  log.R1.mat <- t(matrix(log.R1.mat, x+1, y+1))
  log.R1.mat <- l.A.mat * (log.R1.mat + log (R1.E1.B))
  log.R1.mat3 <- l3.A.mat * (digamma(a1) + log(R1.E3.B))
  
  log.R1.E <- 
    sum(log.R1.mat * l.C.mat) * exp(adj.sum - adj.C) * l1.B + 
    sum(log.R1.mat2 * l2.B)  * exp(adj.sum) + 
    sum(log.R1.mat3 * l3.B)  * exp(adj.sum) + 
    (digamma(a1) + log(b1)) * exp(adj.sum) * l4.B
  log.R1.E <- log.R1.E / l.sum
  
  
  log.R2.mat <- sapply(0:y, log.R2.E1)
  log.R2.mat3 <- log.R2.mat  # saving a vector form
  log.R2.mat3 <- l3.A.mat * (log.R2.mat3 + log (R2.E3.B))
  log.R2.mat <- matrix(log.R2.mat, y+1, x+1)
  log.R2.mat <- l.A.mat * (log.R2.mat + log (R2.E1.B))
  log.R2.mat2 <- l2.A.mat * (digamma(a2) + log(R2.E2.B))
  
  log.R2.E <- 
    sum(log.R2.mat * l.C.mat) * exp(adj.sum - adj.C) * l1.B + 
    sum(log.R2.mat2 * l2.B)  * exp(adj.sum) + 
    sum(log.R2.mat3 * l3.B)  * exp(adj.sum) + 
    (digamma(a2) + log(b1)) * exp(adj.sum) * l4.B
  log.R2.E <- log.R2.E / l.sum
  
  
  # cat("log.R2.E ", log.R2.E, "\n")  
  E.E <- c(sum.AC * exp(adj.sum) * l1.B, sum.A * c(l2.B, l3.B, l4.B)*exp(-adj.C + adj.sum))
  E.E <- E.E/sum(E.E)
  
  # cat("E.E ", E.E, "\n")  
  v.E <- ifelse(y == 0, 0, y) + (a0 + a2) * b2 * sum(E.E[c(2,4)])
  # v.E <- (sum.AC * exp(adj.sum) * l1.B * y + 
  #           sum.A * l2.B * a2 * b2*exp(-adj.C + adj.sum) +
  #           dnbinom(x, a0 + a1 + 1, b1/(1+b1)) * exp(-adj.A - adj.C + adj.sum) * a0 * b2 * p2 * ifelse(y==0, 1, 0) +
  #           sum.A * l3.B * y *exp(-adj.C + adj.sum) +
  #           sum.A * l4.B * (a0 + a2) * b2 *exp(-adj.C + adj.sum)) / l.sum
  
  result <- c(log(l.sum) + adj.A -adj.B1 + adj.C - adj.sum, R0.E, R1.E, R2.E, log.R0.E, log.R1.E, log.R2.E, E.E, v.E) #%>%print
  names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
  return(result)
}

dBvZINB5.Expt.vec <- Vectorize(dBvZINB5.Expt)

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
  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0, if (SE) {rep(NA, 10)}))} # 9 params, lik, iter, pureCor, and 10 SE's
  #print(xy.reduced)
  if (SE) { #internal SE function (only use after param is obtained.)
    .se <- function() {
      se <- do.call(BZINB5.se, c(list(xvec, yvec), as.list(param[1:9])))
      names(se) <- paste0("se.", c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4", "rho"))
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
    expt <- do.call(dBvZINB5.Expt.vec, c(list(xy.reduced$x, xy.reduced$y), as.list(param)))
    expt <- as.vector(expt %*% xy.reduced$freq / n)
# cat(453)
    # loglik = expt[1] * n
    delta <- expt[12] / (expt[2] + expt[4])                   # delta = E(V) / (E(xi0 + xi2))
    param[6:9] = expt[8:11]                                                    # pi = E(Z)
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
      if (SE) result <- c(result, .se())
      return(result)
      }
    if (max(abs(param - param.old)) <= tol) {
      lik <- lik.BvZINB5(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "rho")
      if (SE) result <- c(result, .se())
      return(result)
    }
    if (cor.conv & abs(pureCor - pureCor.old) <= tol) {  # if pureCor is converged, then done!
      lik <- lik.BvZINB5(xvec, yvec, param = param)
      result <- c(param, lik, iter, pureCor)
      names(result) <- c("a0", "a1", "a2", "b1", "b2", paste0("p",1:4), "lik","iter", "rho")
      if (SE) result <- c(result, .se())
      return(result)
    }
  }
  #result <- data.frame(a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], pi = param[6])
  #return(result)
}
ML.BvZINB5 <- function(xvec, yvec, ...) {
  result <- try(ML.BvZINB5.base(xvec,yvec,initial = c(as.numeric(ML.BvNB3(xvec,yvec)), .94,.02,.02,.02), ...))
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
