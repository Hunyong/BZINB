#' rBzinbData
#' BZINB data generation with random normal covariates where the covariates corresponding to eta 
#' and gamma coincide up to the smaller dimension.
#' 
#' install.packages("MASS")
#' @examples
#' rBzinbData()
rBzinbData <- function(n = 200, alpha = c(1, 1, 1), 
                       eta1 = c(-3, 0.5, 0.5, 0.5, 0.5), eta2 = c(-2, 0, 1, 0, 1),
                       gamma1 = c(1,1), gamma2 = c(-1,1), gamma3 = c(-1,1), rho = 0.3) {
  if (length(alpha) != 3) stop("The length of alpha should be 3.")
  if (any(alpha<= 0 )) stop("alpha should be all positive.")
  dEta = length(eta1)
  if (dEta != length(eta2)) stop("The length of eta1 and eta2 should be equal.")
  if (dEta < 1) stop("The length of eta1 should be at least one.")
  dGamma = length(gamma1)
  if (dGamma != length(gamma2) | dGamma != length(gamma3)) stop("The length of gamma1, gamma2, and gamma3 should be equal.")
  if (dGamma < 1) stop("The length of gamma1 should be at least one.")
  d = max(dEta, dGamma)
  if (d > 1) {
    Sigma = diag(rep(1-rho, d - 1)) + rho
    covariate <- cbind(1, MASS::mvrnorm(n = n, mu = rep(0, d-1), Sigma = Sigma)) # adding the intercept term.
  } else {
    covariate <- matrix(1, n, 1)
  }
  Z <- covariate[, 1:dEta]
  W <- covariate[, 1:dGamma]
  b12 = exp(Z %*% cbind(eta1, eta2 - eta1))  # b1 and b2/b1
  p123 = exp(W %*% cbind(gamma1, gamma2, gamma3, 0))
  pSum = apply(p123, 1, sum)
  p123 = p123/pSum
  
  rmat <- matrix(rgamma(n*3, shape = rep(alpha, each = n), rate = 1/b12[, 1]), n, 3)
  rmat2 <- rmat
  rmat2[,3] <- rmat[,1] + rmat[,3]
  rmat2[,2] <- rmat[,1] + rmat[,2]
  rmat2 <- rmat2[,2:3]
  rmat2[,2] <- rmat2[,2]*b12[,2]
  uv <- matrix(rpois(n*2, rmat2), n, 2)
  
  E <- t(sapply(1:n, function(i) rmultinom(1, 1, p123[i, ])))
  z <- cbind(E[,1]+E[,2], E[,1]+E[,3])
  
  xy <- uv * z
  colnames(xy) <- c("y1", "y2")
  colnames(covariate) <- paste0("X", 1:d - 1)
  
  return(as.data.frame(cbind(xy, covariate[, -1, drop = FALSE])))
}

