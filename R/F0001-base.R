## 1.1 pairwise association macro
pairwise.dep <- function(data, measure, name, nondata.col = 1, p.nonzero = TRUE, rm.id.vec = TRUE, log.data = FALSE, add.log = 1, p.value = FALSE, ...) {
  # data: a dataframe of which pairs are to be analyzed
  # measure can be either a function or a list of functions
  # nondata.col: the column number to be ignored for measurement eg. first col (gene names)
  # name: list of names of the association measures (need to find some code to extract from the list of functions)
  # log, add.log: log tansformation, addition before log transformation
  # p.nonzero: show proportion of non-zero values
  # p.value: Calculate the p-values by permutation test
  data <- as.matrix(data[,-nondata.col])
  if (log.data) {data <- log(data + add.log)}
  dim.p <- dim(data)[1]
  comb <- expand.grid(1:dim.p,1:dim.p)
  if (rm.id.vec)  {rw <- which(comb[,1] > comb[,2])} else {rw <- which(comb[,1] >= comb[,2])}
  comb <- data.frame(comb[rw,c(2,1)])
  # [lower.tri(matrix(1:dim.p^2,dim.p),TRUE),][,c(2,1)]
  # comb <- data.frame(t(combn(1:dim.p,2)))
  comb$pair <- apply(comb,1,function(x) paste(x,collapse="-"))
  
  # function of a measurement
  fun.measure <- function(measure0, comb0 = comb, data0 = data, p.value0 = p.value, ...) {
    a <- sapply(1:dim(comb0)[1], function(i) {
      x1 <- data0[comb0[i,1],]
      x2 <- data0[comb0[i,2],]
      point.est <- try(measure0(x1, x2, ...),silent=TRUE)
      point.est <- ifelse(class(point.est)=="try-error",NA, point.est)
      if (p.value0 == TRUE) {
        if (is.na(point.est) == FALSE) {
          p.val <- p.permutation(x1, x2, measure = measure0, point.est = point.est, ...)$p.val
        } else {p.val <- NA}
        return(data.frame(point.est=point.est, p.val=p.val))
      } else {return(data.frame(point.est=point.est, p.val=NA))}
    }
    )
    return(as.data.frame(t(a)))
  }
  if (class(measure) != "list") {measure <- list(measure)}
  ##sapply(measure, fun.measure)
  
  storage <- list(); j = 0
  for (x in measure) {
    j = j + 1
    storage[[j]] <- fun.measure(x)
    names(storage)[j] <- name[j]
  }
  for (k in 1:length(measure)) {comb <- cbind(comb, unlist(storage[[k]]$point.est))}
  for (k in 1:length(measure)) {comb <- cbind(comb, unlist(storage[[k]]$p.val))}
  names(comb) <- c(1,2,"pair", name, paste0(name,".p.val"))
  
  if (p.nonzero == TRUE) {
    n <- dim(data)[2]
    comb$non0.1 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,1],] != 0) / n})
    comb$non0.2 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,2],] != 0) / n})
    comb$non0.min <- pmin(comb$non0.1,comb$non0.2)
  }
  return(comb)
}

if (FALSE) {
  # example
  pairwise.dep(matrix(1:100,5), mi.empirical.2, name = "mi", p.value=TRUE,n.permutation=20)
  pairwise.dep(matrix(1:100,5), cor, name = "cor")
  pairwise.dep(matrix(0,5,20), mi.empirical.2, name = "mi")  #all zero data -> all NA
  pairwise.dep(matrix(1:100,5), list(cor, mi.empirical.2), name = c("cor","MI"), log = TRUE)
  pairwise.dep(matrix(1:100,5), list(cor, mi.empirical.2), name = c("cor","MI"), p.value=TRUE, n.permutation=20)
}


# remove global assign later !!!! <<- 
pairwise.MLE <- function(data, ML.fun, nondata.col = 1, p.nonzero = TRUE, rm.id.vec = TRUE, print.each.pair=TRUE, ...) {
  # data: a dataframe of which pairs are to be analyzed
  # nondata.col: the column number to be ignored for measurement eg. first col (gene names)
  # name: short name for the distribution
  # p.nonzero: include nonzero proportion in the result
  data <- as.matrix(data[,-nondata.col])
  dim.p <- dim(data)[1]
  comb <- expand.grid(1:dim.p,1:dim.p)
  if (rm.id.vec)  {rw <- which(comb[,1] > comb[,2])} else {rw <- which(comb[,1] >= comb[,2])}
  comb <- data.frame(comb[rw,c(2,1)])
  # [lower.tri(matrix(1:dim.p^2,dim.p),TRUE),][,c(2,1)]
  # comb <- data.frame(t(combn(1:dim.p,2)))
  comb$pair <- apply(comb,1,function(x) paste(x,collapse="-"))
  
  # empty.result <- ML.fun(xvec = rep(0,3), yvec = rep(0,3))
  # comb <- cbind(comb, matrix(ncol = dim(empty.result)[2]))
  # names(comb) <- c(1,2,"pair", names(empty.result))
  names(comb) <- c(1,2,"pair")
  MLE <- apply(comb[,1:2], 1, function(s) {
    # if (s[1] <= 6 | s[2] <=23) {return(data.frame(matrix(NA,1,4)))} # debug #7,24 has problem
    x <- data[s[1],]
    y <- data[s[2],]
    tt(1)
    print(c(Estimated.pair = s))
    tmp.sss <<- s   # debug
    tmp.aaa <<- cbind(x,y)   # debug
    # result <- ML.fun(xvec = x, yvec = y, ...)
    #result <- ifelse(class(result)=="try-error",rep(NA,length(ML.fun(xvec = 0, yvec = 0, ...))), result)   %>% print
    result <- try(ML.fun(xvec = x, yvec = y, ...),silent=TRUE)
    if (class(result)=="try-error") {
      result <- rep(NA,length(ML.fun(xvec = 0, yvec = 0, ...)))
    }
    if (print.each.pair) {print(result); print(tt(2))}
    return(result)
  })
  tmp.bbbb <<- MLE
  if (class(MLE) == "list") {MLE <- do.call(rbind, MLE)}
  else if (class(MLE) == "matrix") { MLE <- t(MLE)}
  comb <- cbind(comb, MLE)
  
  if (p.nonzero == TRUE) {
    n <- dim(data)[2]
    comb$non0.1 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,1],] != 0) / n})
    comb$non0.2 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,2],] != 0) / n})
    comb$non0.min <- pmin(comb$non0.1,comb$non0.2)
  }
  
  return(comb)
}
if (FALSE) { # example
  head(pairwise.MLE(data=data[iset.Mm.c2[[1]],], ML.fun = ML.BvNB1))
}

extractor <- function(gene.no , geneset.no = 1, dat=data, nondata.col = 1, genesetlist = iset.Mm.c2) {
  extracted <- dat[genesetlist[[geneset.no]][gene.no],-nondata.col]
  extracted <- as.numeric(extracted)
  return(extracted)
}
if (FALSE) {extractor(1)}

# time calculator
tt <- function(s){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    return(data.frame(begin = time.tmp, end = Sys.time(), elapsed = Sys.time() - time.tmp))
  }
}

tt.em <- function(s, step){
  if (s==0) { # initialize
    time.E.step <<- 0 # record time. initial
    time.M.step <<- 0 # record time. initial
    print("Time for E- and M-steps initialized.")
  } else if (s==1) { # begin time
    time.begin.tmp <<- Sys.time()
  } else if (s==2) { # measure elapsed time
    time.elapsed = Sys.time() - time.begin.tmp
    if (step == "E") {
      time.E.step <<- time.E.step + time.elapsed
    } else if (step == "M") {
      time.M.step <<- time.M.step + time.elapsed
    } else warning("Wrong step was given. E or M should be given.")
  } else if (s==3) {
    print(c("E-step" = time.E.step, "M-step" = time.M.step))
  }
}
if (FALSE) {
  tt.em(0) # initialize
  tt.em(1) # start measuring
  tt.em(2, "E") # stack up times for E-step
  tt.em(1)
  tt.em(2, "M")
  tt.em(3) # summary
  
}


nonzero <- function(x, cut = .1) {
  len = length(x); above = names(which(x > cut)); len.a = length(above)
  return(list(stat = c(n.total = len, n.above = len.a, p.above.cut = round(len.a/len,2)),
              which = above))
}
screen.zero <- function(data, geneset, cut = .1, output="which", exclude.col = 1) {
  # output = either "which" or "stat"
  data[geneset,-1] %>% apply(1,function(x) mean(x!=0)) %>% nonzero(cut=cut) %>%"[["(output)
}
screened.length <- function(data, geneset, cut = .1, exclude.col = 1) {
  screened.set = screen.zero(data = data, geneset = geneset, cut = cut, output="which", exclude.col = exclude.col)
  return(length(screened.set))
}
