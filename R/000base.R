#' @useDynLib bzinb
#' @export
#' @examples 
#' # A sample data with 4 genes and sample size of 30
#' set.seed(1)
#' tmp = t(sapply(1:4, rBvZINB5, n = 30, param = c(1, 1, 1, 1, 1, .5, .2, .2, .1)))
#' # Calculating all six pairwise gene-gene correlations
#' pairwise.bzinb(tmp, showFlag = TRUE)
pairwise.bzinb <- function(data, nonzero.prop = TRUE, fullParam = FALSE, 
                         showFlag = FALSE, ...) {
  # data: a dataframe of which pairs are to be analyzed. rows = genes, cols = sample
  # nonzero.prop: include nonzero proportion in the result
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  } else if (!is.matrix(data)) {
    stop ("data is neither matrix nor data.frame.")
  }
  
  dim.p <- dim(data)[1]
  comb <- expand.grid(1:dim.p, 1:dim.p)
  rw <- which(comb[, 1] > comb[, 2])
  comb <- data.frame(comb[rw, c(2,1)])
  comb$pair <- apply(comb, 1, function(x) paste(x, collapse = "-"))
  names(comb) <- c(1, 2, "pair")
  
  MLE <- apply(comb[,1:2], 1, function(s) {
    # if (s[1] <= 6 | s[2] <=23) {return(data.frame(matrix(NA,1,4)))} # debug #7,24 has problem
    x <- data[s[1], ]
    y <- data[s[2], ]
    result <- ML.BvZINB5(xvec = x, yvec = y, ...)
    
    if (showFlag) cat("pair ", s[1], "-", s[2], ": ", "(rho, se.rho) = (", 
                      formatC(result$rho[1,1], digits = 5, format = "f", flag = "0"), ", ", 
                      formatC(result$rho[1,2], digits = 5, format = "f", flag = "0"), ")\n")
    result2 <- c(rho = result$rho["rho", "Estimate"], 
                 se.rho = result$rho["rho", "Std.err"],
                 if (fullParam) result$coefficients[, "Estimate"],  
                 if (fullParam) result$coefficients[, "Std.err"], 
                 if (fullParam) c(result$lik, result$iter))
    if (fullParam) {
      names(result2) <- c("rho", "se.rho", abp.names, 
                          paste0("se.", abp.names), "logLik", "num.iter")
    }
    return(result2)
  })
  comb <- cbind(comb, t(MLE))
  
  if (nonzero.prop == TRUE) {
    n <- dim(data)[2]
    comb$nonzero.1 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,1],] != 0) / n})
    comb$nonzero.2 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,2],] != 0) / n})
    comb$nonzero.min <- pmin(comb$nonzero.1, comb$nonzero.2)
  }
  
  return(comb)
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


# binary profile function: for ML.BZIP.B
bin.profile <- function(xvec, yvec) {
  xvec[xvec != 0] = 1
  yvec[yvec != 0] = 1
  vec <- cbind(xvec,yvec)
  a <- rep(0,4)
  a[1] <- sum(apply(vec,1,prod)) # 1,1
  a[2] <- sum(xvec) - a[1]       # 1,0
  a[3] <- sum(yvec) - a[1]       # 0,1
  a[4] <- length(xvec) - sum(a)  # 0,0
  return(a)
}

.check.input <- function(xvec, yvec) {
  
}
