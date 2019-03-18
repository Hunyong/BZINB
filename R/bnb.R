### 1. Density, likelihood, gradient :  This function affects bzinb3, bzinb4
dbnb <- function(x, y, a0, a1, a2, b1, b2, log=FALSE) {
  p1 = (b1 + b2 + 1) /(b1 + 1); p2 = (b1 + b2 + 1) /(b2 + 1)
  adj = 0
  l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p1))
  l2 <- function(m) sum(sapply(0:x, l1, m = m))
  #print(sapply(0:x, l1, m = 0)); print(sum(sapply(0:x, l1, m = 0))) ##
  l3 <- function(m) log(l2(m)) + lgamma(m +a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p2)
  l4 <- sum(sapply(0:y, function(m) exp(l3(m)))) #%>%print
  if (is.infinite(l4)) {adj = 200; l4 <- sum(sapply(0:y, function(m) exp(l3(m) - adj)))} #%>%print
  l4 <- log(l4) - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1*log(1+b1) - a2*log(1+b2)  + adj
  return(ifelse(log, l4, exp(l4)))
}
if (FALSE) {
  dbnb(1,1,1,1,1,1,2) ; dbnb(1,1,1,1,1,1,1); dbnb(1,1,1,1,1,1)
  tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dbnb(s,r,1,1,1,1,.5)))
  sum(tmp) #1
}
dbnb.vec <- Vectorize(dbnb)
lik.bnb <- function(xvec, yvec, a0, a1, a2, b1, b2, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 5) stop("length(param) must be 5.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
  }
  sum(log(dbnb.vec(xvec, yvec, a0, a1, a2, b1, b2)))
}

#' @export
#' @rdname bzinb
rbzinb <- function(n, a0, a1, a2, b1, b2, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 5) stop("length(param) must be 5.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
  }
  if (length(n) != 1) {stop("length(n) must be 1.")}
  .check.ab(c(a0, a1, a2, b1, b2))
  
  rmat <- matrix(rgamma(n*3, shape = c(a0, a1, a2), rate = 1/b1), n, 3, byrow=TRUE)
  rmat2 <- rmat
  rmat2[,3] <- rmat[,1] + rmat[,3]
  rmat2[,2] <- rmat[,1] + rmat[,2]
  rmat2 <- rmat2[,2:3]
  rmat2[,2] <- rmat2[,2]*b2/b1
  xy <- matrix(rpois(n*2, rmat2), n, 2)
  colnames(xy) <- c("x", "y")
  
  return(xy)
}

dbnb.gr <- function(x, y, a0, a1, a2, b1, b2) {
  p1 = (b1 + b2 + 1) /(b1 + 1); p2 = (b1 + b2 + 1) /(b2 + 1)
  gr.b1.1 <- x/b1 - (x + y + a0) /(b1 + b2 + 1) - a1 / (b1 + 1)
  gr.b1 <- function(k, m) {(k + m) / (b1 + b2 + 1) - k / (b1 + 1) + gr.b1.1}
  gr.b2.1 <- y/b2 - (x + y + a0) /(b1 + b2 + 1) - a2 / (b2 + 1)
  gr.b2 <- function(k, m) {(k + m) / (b1 + b2 + 1) - m / (b2 + 1) + gr.b2.1}
  gr.a0.1 <- - digamma(a0) - log(1 + b1 + b2)
  gr.a0 <- function(k, m) {if (x==0 & y==0) {- log(1 + b1 + b2)} else  {digamma(x +y +a0 -k -m) + gr.a0.1}}

  gr.a1.1 <- - digamma(a1) - log(1 + b1)
  gr.a1 <- function(k) {if (k==0) {- log(1 + b1)} else  {digamma(a1 + k) + gr.a1.1}}
  gr.a2.1 <- - digamma(a2) - log(1 + b2)
  gr.a2 <- function(m) {if (m==0) {- log(1 + b2)} else  {digamma(a2 + m) + gr.a2.1}}
  
  l1 <- function(k, m) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p1)
                           + lgamma(m +a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p2))
  l2 <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
  l2 <- exp(l2)
  l1.mat <- sapply(0:x, function(k) sapply(0:y, l1, k=k))
  l1.mat <- (l1.mat * l2) #%>%print
  l1.sum <- sum(l1.mat) 
  
  gr.b1.mat <- sapply(0:x, function(k) sapply(0:y, gr.b1, k=k))
  gr.b1.mat <- l1.mat * gr.b1.mat
  gr.b1.sum <- sum(gr.b1.mat)/sum(l1.mat)
  
  gr.b2.mat <- sapply(0:x, function(k) sapply(0:y, gr.b2, k=k))
  gr.b2.mat <- l1.mat * gr.b2.mat
  gr.b2.sum <- sum(gr.b2.mat)/sum(l1.mat)
    
  gr.a0.mat <- sapply(0:x, function(k) sapply(0:y, gr.a0, k=k))
  gr.a0.mat <- l1.mat * gr.a0.mat
  gr.a0.sum <- sum(gr.a0.mat)/sum(l1.mat)
  
  gr.a1.mat <- matrix(sapply(0:x, gr.a1),x+1, y+1) #%>% print
  gr.a1.mat <- l1.mat * t(gr.a1.mat)
  gr.a1.sum <- sum(gr.a1.mat)/sum(l1.mat)
  
  gr.a2.mat <- matrix(sapply(0:y, gr.a2), y+1, x+1) #%>% print
  gr.a2.mat <- l1.mat * gr.a2.mat
  gr.a2.sum <- sum(gr.a2.mat)/sum(l1.mat)
  result <- list(logdensity=log(l1.sum), gradient = c(gr.a0.sum, gr.a1.sum, gr.a2.sum, gr.b1.sum, gr.b2.sum))
  return(result)
}
dbnb.gr.vec <- Vectorize(dbnb.gr)

if (FALSE) {
  tmp <- dbnb.gr.vec(c(1,1,1),c(0,1,2),1,1,1,1,2)[2,]
  sapply(tmp, cbind)
}

### 2. MLE
bnb <- function (xvec, yvec, tol = 1e-8, method="BFGS", showFlag=FALSE) {
  .check.input(xvec, yvec)
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  if (max(xvec)==0 & max(yvec)==0) {return(rep(1e-10,5))}
  #print(xy.reduced)

  fn.1 = function (param) {
    val <- dbnb.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5])
    lik <- sapply(val[1,],cbind) %*% xy.reduced$freq
    return(lik)
  }
  gr.1 = function (param) {
    val <- dbnb.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5])
    lik <- sapply(val[1,],cbind) %*% xy.reduced$freq
    gr <- sapply(val[2,],cbind) %*% xy.reduced$freq
    # print(c(param,gr)) ###
    return(gr)
  }

  
  #log-scaled params: param.l
  fn.log = function (param.l) { fn.1 (exp(param.l))}
  gr.log = function (param.l) { 
    if (showFlag) {print(exp(param.l))}
    (as.vector(gr.1 (exp(param.l))) * exp(param.l)) 
  }
  
  # initial guess
  xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
  s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
  cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
  initial <- rep(NA,5)
  initial[4] <- s2.x /xbar
  initial[5] <- s2.y /ybar
  initial[2:3] <- c(xbar,ybar)/initial[4:5]
  initial[1] <- min(initial[2:3]) * abs(cor.xy)
  initial[2:3] <-  initial[2:3] - initial[1]
  initial <- pmax(initial, 1e-5)
  #print(initial) ###
  
  # result <- multiroot(f = gr.log, start = c(1,1,1,1), positive=TRUE)$root
  # result <- exp(optim(par = log(initial), fn = fn.log, control=list(fnscale=-1), method = "Nelder-Mead")$par)
  result <- try(exp(optim(par = log(initial), fn = fn.log, gr = gr.log, control=list(fnscale=-1, abstol = tol), method = method)$par), silent=TRUE)
  if (class(result)=="try-error") {
    initial = rep(1,5)
    result <- exp(optim(par = log(initial), fn = fn.log, gr = gr.log, control=list(fnscale=-1, abstol = tol), method = method)$par)
  }
  result <- data.frame(a0 = result[1], a1 = result[2], a2 = result[3], b1 = result[4], b2 = result[5])
  return(result)
}
# simple tests
if (FALSE) {
  bnb(c(10,1,1),c(10,1,2)); bnb(c(10,1,1),c(10,1,2))
  tt(1)
  a <- bnb(extractor(1), extractor(38), tol=1e-15)
  tt(2) # 31secs
  
  bnb(extractor(2), extractor(5), method="BFGS", tol=1e-30)
  bnb(extractor(2), extractor(5), method="Nelder-Mead", tol=1e-30)
  sum(sapply(dbnb.gr.vec(extractor(2), extractor(5), 0.0004670043, 0.002597285, 0.01597736, 10.48014, 35.76876)[1,],cbind))
  sum(sapply(dbnb.gr.vec(extractor(2), extractor(5), 0.0004686926, 0.002604057, 0.01598007, 10.45238, 35.75869)[1,],cbind))
  
  tt(1)
  bnb(extractor(1), extractor(3), method="BFGS")
  tt(2) #8sec
  tt(1)
  bnb(extractor(1), extractor(38), method="BFGS", showFlag=TRUE)
  bnb(extractor(1), extractor(38), method="Nelder-Mead", showFlag=TRUE)
  tt(2) #31sec
  #lik.bnb(extractor(1), extractor(38), param = c(5.790158e-03, 4.300688e-03, 7.836757e-02, 7.586956e+01, 1.015767e+02))
  
}
# some deviance tests
if (FALSE) {
  tt(1)
  a <- bnb(extractor(11), extractor(16))
  tt(2) #1.8 sec!
  lik.bnb(extractor(11), extractor(16), param = a) #lik = -805
  
  tt(1)
  a <- bnb(extractor(11), extractor(16))
  tt(2) #1.9 sec!
  lik.bnb(extractor(11), extractor(16), param = a) #lik = -805
  
  tt(1)
  a <- bnb(extractor(30), extractor(50))
  tt(2) #0.2 sec!
  lik.bnb(extractor(30), extractor(50), param = a) #lik = -113
  
  tt(1)
  a <- bnb(extractor(8), extractor(11))
  tt(2) #17 sec!
  lik.bnb(extractor(8), extractor(11), param = a) #lik = -1938
  
  tt(1)
  a <- bnb(extractor(1), extractor(3))
  tt(2) #17 sec!
  lik.bnb(extractor(1), extractor(3), param = a) #lik = -1458
  
}

if (FALSE) {
  tt(1)
  a <- pairwise.MLE(data=data[iset.Mm.c2[[1]],], ML.fun = bnb, showFlag=TRUE)
  tt(2)
  bnb(extractor(1), extractor(12), showFlag=TRUE)
  bnb(extractor(51), extractor(12), showFlag=TRUE)
  
  tt(1)
  MLE.Geneset1$bnb <- pairwise.MLE(data[iset.Mm.c2[[1]],], ML.fun = bnb, showFlag=TRUE)  ## 1.8hrs
  tt(2)
  
  bnb(0,155)
  
}



### 3. Deviance
dev.bnb <- function(xvec, yvec, param = NULL, a0 = NULL, a1 = NULL, a2= NULL, b1 = NULL, b2 = NULL) {
  .check.input(xvec, yvec)
  
  # If params = NULL, apply bnb. if else, apply those params
  if (is.null (param)) { 
    if (is.null (a0) | is.null (a1) | is.null (a2) | is.null (b1) | is.null (b2)) {
      param = bnb (xvec = xvec, yvec = yvec)
    }
    else { param = c(a0, a1, a2, b1, b2)}
  }
  # Log-likelihood of the bnb model
  lik.model <- lik.bnb (xvec = xvec, yvec = yvec, param = param)
  
  # Reduced calculation
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  
  # Saturated model bzip params
  bnb.vec <- Vectorize(bnb)
  param.sat <- t(bnb.vec(xy.reduced$x, xy.reduced$y)) # %>%print
  param.sat <- do.call(rbind, param.sat)  #"matrix of lists" into "a matrix"
  lik.sat   <- sum(dbnb.vec(x= xy.reduced$x, y = xy.reduced$y, a0 = param.sat[,1], a1 = param.sat[,2], a2 = param.sat[,3], b1 = param.sat[,4], b2 = param.sat[,5],  log = TRUE) * xy.reduced$freq)
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}
if (FALSE) {
  dev.bnb(extractor(1), extractor(2), param=as.vector(a[,2]))  #455  vs bzip.B:3150
  dev.bnb(extractor(1), extractor(5), param=as.vector(a[,5]))  #849  vs bzip.B:2877
  dev.bnb(extractor(1), extractor(3))  #2028 vs bnb: 2030
  dev.bnb(extractor(1), extractor(8))  #2066 vs bnb: 2077
  
  ## heatmap
  a <- bnb(extractor(1), extractor(3), method="BFGS")
  tmp <- sapply(0:50, function(r) sapply(0:50, function(s) dbnb(s,r, a[1,1],a[1,2],a[1,3],a[1,4],a[1,5])))
  plot_ly(z = tmp[1:2,2:10], type = "heatmap")
  table(extractor(1), extractor(3))
  
  MLE.Geneset1$bnb ########
  a <- cbind(1,2,dev.bnb(extractor(1), extractor(2), param = c(1,1,1,1,1,1,1)))
  tt(1)
  for (i in 1:dim(MLE.Geneset1$bnb)[1]) {
    #if (i ==300 | i == 349 ) {a[i,] <- c(a1, a2, NA,NA,NA)}
    #if (i >= 350) {
    tmp <- MLE.Geneset1$bnb[i,]
    a1 <- tmp[1]; a2 <- tmp[2]
    a3 <- tmp[4:7] #params
    a[i,] <- c(a1,a2,dev.bnb(extractor(as.numeric(a1),1), extractor(as.numeric(a2),1), param = a3))
    print(a[i,])  
  }
  tt(2) # 5.8 hours
  
  MLE.Geneset1$bnb$dev <- a$deviance
  
  sum(MLE.Geneset1$bnb$dev>863); sum(MLE.Geneset1$bnb$dev<863)
  plot3.2.1 <-  ggplot(MLE.Geneset1$bnb, aes(non0.min, dev)) +
    geom_point(aes(color=1-non0.min)) + ggtitle(paste0("deviance of bnb model - Geneset 1")) +
    ylim(c(0,25000)) +
    xlab("nonzero-count proprotion(min)") + ylab("Deviance of bnb") + geom_abline(intercept=863, slope=0, linetype, color="red", size=0.1)
  png2(plot3.2.1, width = 720, height = 360)
}

