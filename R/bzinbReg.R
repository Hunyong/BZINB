##########################################################################################
## 4. Bivariate Zero-Inflated Negative Binomial Regression (BZINB-Reg)
##########################################################################################

#' @rdname bzinbReg
#' @name bzinbReg
#' @aliases bzinbReg
#' 
#' @title The bivariate zero-inflated negative binomial regression.
#' 
#' @description the bivariate zero-inflated negative regression.
#' 
#' @param y1,y2 a pair of bzinb random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param n number of observations.
#' @param initial starting value of param for EM algorithm, a vector of nine values.
#' @param tol tolerance for judging convergence. \code{tol = 1e-8} by default.
#' @param maxiter maximum number of iterations allowed. \code{tol = 50000} by default.
#' @param showFlag if \code{TRUE}, the updated parameter estimates for each iteration 
#'   are printed out. If a positive integer, the updated parameter estimates for 
#'   iterations greater than \code{showFlag} are printed out.
#' @param vcov if \code{TRUE}, the variance-covariance matrix and information matrix 
#'   are returned.
#' 
#' @return 
#'  \itemize{
#'    \item the regression coefficients by the MLE of the BZINB regression model.
#'      \itemize{
#'        \item \code{coefficients} estimate and standard error of the BZINB parameters
#'        \item \code{lik} log-likelihood of the maximum likelihood estimate
#'        \item \code{iter} total number of iterations
#'        \item \code{info} information matrix. Provided when \code{vcov} is \code{TRUE}.
#'        \item \code{vcov} variance-covariance matrix. Provided when \code{vcov} is \code{TRUE}.
#'      }
#'    \item \code{lik.bzinb} gives the log-likelihood of a set of parameters for a BZINB pair.
#'
#'  }
#'  
#' @details
#'  
#' 
#' @examples
#' bzinb(y1 = data1[,1], y2 = data1[,2], showFlag = FALSE)
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Preisser, J., Liu, C., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#'  Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from 
#'  incomplete data via the EM algorithm. Journal of the Royal Statistical Society: 
#'  Series B (Methodological), 39(1), 1-22.
#'  
#' @import Rcpp
#' @export
#' @useDynLib bzinb
#' 


abp.names <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4") # global variable
expt.names <- c("lik", "ER0", "ER1", "ER2", "ElogR0", "ElogR1", "ElogR2", "EE1", "EE2", "EE3", "EE4", "EV")

bzinb.base.reg <- 
  function (y1, y2, Z, W, 
            tol = 1e-8, maxiter = 50000, showFlag = FALSE, 
            vcov = FALSE, bnb = 0, initial = NULL) {
  se = TRUE  # se is estimated by default
  n <- length(y1)
  pZ <- dim(Z)[2]
  pW <- dim(W)[2]
  dim.param <- 3 + 2 * pZ + 3 * pW
  
  if (max(y2)==0) {
    stop("TBD!!!!!!!!!!")
    # return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0, if (se) {rep(NA, 11)}))
  } # 9 params, lik, iter, pureCor, and 11 se's
  
  # info <- if (se) {matrix(0, ncol = 8, nrow = 8, dimnames = list(abp.names[-9], abp.names[-9]))} else {0}
  ae.names  <- c("a0", "a1", "a2", 
                 paste0("eta1_", colnames(Z)),
                 paste0("eta2_", colnames(Z)))
  aeg.names <- c("a0", "a1", "a2", 
                 paste0("eta1_", colnames(Z)),
                 paste0("eta2_", colnames(Z)),
                 paste0("gamma1_", colnames(W)),
                 paste0("gamma2_", colnames(W)),
                 paste0("gamma3_", colnames(W)))
  # initial guess
  if (is.null(initial)) {
    y1bar <- mean(y1); y2bar <- mean(y2); xybar <- mean(c(y1bar, y2bar))
    s2.y1 <- var(y1); s2.y2 <- var(y2); if(is.na(s2.y1)|is.na(s2.y2)) {s2.y1 <- s2.y2 <- 1}
    cor.y <- cor(y1,y2); if (is.na(cor.y)) {cor.y <- 0}
    zero <- sum(y1 == 0 & y2 == 0) / n

    initial <- rep(0.1, dim.param)
    names(initial) <- aeg.names
    # initial[4] <- s2.x /ifelse(y1bar==0,1e-4, y1bar)
    # initial[5] <- s2.y /ifelse(y2bar==0,1e-4, y2bar)
    # initial[2:3] <- c(y1bar,y2bar)/pmax(initial[4:5], c(0.1,0.1))
    # initial[1] <- min(initial[2:3]) * abs(cor.y)
    # initial[2:3] <-  initial[2:3] - initial[1]
    # 
    # initial[6:9] <- bin.profile(y1, y2)   # freq of each zero-nonzero profile
    # initial[6:9] <- initial[6:9]/sum(initial[6:9])      # relative freq
    # initial <- pmax(initial, 1e-5)
    # if(is.na(sum(initial))) { initial[is.na(initial)] <- 1}
    
  } else {
    names(initial) <- aeg.names
  }
  ## tmp.initial <<- initial
  param = initial
  lik = -Inf
  lik.vec = rep(0, maxiter + 1)
  nonconv = 0L
  Z[] = as.double(Z) 
  W[] = as.double(W)
  
# print(param)
# print(y1)
# print(y2)
# print(Z)
# print(W)
# print(pZ)
# print(pW)
# print(n)

  em.out <- emReg(param2 = param, xvec = y1, yvec = y2, 
                  ZZ = Z, WW = W,
                  pZ = pZ, pW = pW,
                  n = n, se = as.integer(se), 
                  maxiter = as.integer(maxiter), tol = as.double(tol), 
                  showFlag = as.integer(showFlag), bnb = bnb)
  
  em.out[[3]] <- matrix(em.out[[3]], dim.param, dim.param)
  names(em.out) <- c("param2",              # "y1", "y2", "freq", "n",
                     "expt", "info",        # "se",
                     "iter",                # "maxiter", "tol", "showFlag",
                     "nonconv", "trajectory") # , "bnb")

  # overwriting the original object
  param = setNames(em.out$param2, aeg.names)
  expt  = setNames(em.out$expt, expt.names)
  info = if (se) {matrix(em.out$info, dim.param, dim.param,
                         dimnames = list(aeg.names, aeg.names))} else 0
  iter  = em.out$iter
  nonconv = em.out$nonconv
  trajectory = em.out$trajectory

  if (nonconv == 1) warning("The iteration exited before reaching convergence.")

  # tmp.expt <<- expt
  # tmp.traj <<- lik.vec
  # underlying correlation (rho)
  # rho <- param[1]/sqrt((param[1] + param[2]) * (param[1] + param[3])) *
  #   sqrt(param[4] * param[5] /(param[4] + 1) /(param[5] + 1))
  # logit.rho <- qlogis(rho)

  # if (bnb) {info = info[1:5, 1:5]}
  # 
  if (se) {
    par.names <- if (bnb) ae.names else aeg.names
    dim.par <- length(par.names)

    qr.info <- try(qr(info))
    if (class(qr.info)[1] == "try-error") {
      warning ("The information matrix has NA/NaN/Inf and thus the standard error is not properly estimatd.")
      std.param = setNames(rep(NA, dim.par), par.names)
      cov.mat <- NA
    } else if (qr(info)$rank < dim.par) {
      warning ("The information matrix is (essentially) not full rank, and thus the standard error is not reliable.")
      std.param = setNames(rep(NA, dim.par), par.names)
      cov.mat <- NA
    } else {
      cov.mat <- try(solve(info))
      if (class(cov.mat)[1] == "try-error") {
        std.param = setNames(rep(NA, dim.par), par.names)
        cov.mat <- NA
      } else {
        std.param = sqrt(diag(cov.mat))
      }
    }

  }

  result <- list(coefficients = 
                   matrix(c(param, if (se) std.param else rep(NA, dim.param)),
                          ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))),
                 lik = expt[1],
                 iter = iter)
  if (se & vcov) {
    result$info = info
    result$vcov = cov.mat
  }
  
  return(result)
}

#' @export
#' @rdname bzinbReg
#' @example 
#' set.seed(1)
#' dat <- rBzinbData()
#' bzinb.reg(cbind(y1, y2) ~ ., ~ X1, data = dat)
bzinb.reg <- function(mu.formula,       # mu.formula = cbind(Sepal.Length, Sepal.Width) ~ Species + Petal.Length
                      nu.formula = ~ 1, # nu.formula = ~ Species + Petal.Width
                      data,
                      tol = 1e-8, maxiter = 50000, showFlag = FALSE,
                      vcov = FALSE, initial = NULL) {
  
  formula = mu.formula
  formula[[3]] = rlang::expr(!!mu.formula[[3]] + !!nu.formula[[2]])
  if (missing(data)) data <- environment(formula)
  mf = model.frame(formula, data = data)
  mt <- attr(mf, "terms")
  mtNB <- terms(mu.formula, data = data)
  Z <- model.matrix(mtNB, mf)   # covariates for the NB model
  mtZI <- terms(nu.formula, data = data)
  W <- model.matrix(mtZI, mf)   # covariates for the zero model
  y <- model.response(mf, "numeric")
  if (dim(y)[2] != 2) stop("The dimension of the outcome variable (in mu.formula) must be two.")
  if (any(y < 0)) stop("Negative outcome values are not allowed.")
  y[] <- as.integer(round(y, 0))  # y values in integers. 

  ### from here!!!
  result <- try(bzinb.base.reg(y1 = y[, 1], y2 = y[, 2], 
                               Z = Z, W = W, tol = tol, maxiter = maxiter, 
                               showFlag = showFlag, vcov = vcov, initial = initial))
  if (class(result)[1] == "try-error") {
    result <- list(coefficients = matrix(rep(NA, 18),
                                         ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))),
                   lik = NA,
                   iter = NA)
    if (vcov) {
      result$info = NA
      result$vcov = NA
    }
  }
  return(result)
}

