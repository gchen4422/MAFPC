# -*- coding: utf-8 -*-
estimate_s_rss_joint <- function(z1, z2,
                                 R1, R2,
                                 n1 = NULL, n2 = NULL,
                                 r_tol = 1e-8) {
  #-- 1) finite-sample shrinkage on each z
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be >1")
      sigma2 <- (n-1)/(z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)
  
  #-- 2) stack z’s and build block-diag R
  z_joint <- c(z1, z2)
  # use Matrix::bdiag for efficiency if you like
  R_joint <- as.matrix(Matrix::bdiag(R1, R2))
  
  #-- 3) eigen-decomp (cached if already present)
  if (is.null(attr(R_joint, "eigen")))
    attr(R_joint, "eigen") <- eigen(R_joint, symmetric = TRUE)
  eigJ <- attr(R_joint, "eigen")
  # floor tiny negatives
  if (any(eigJ$values < -r_tol))
    warning("R_joint not PSD: flooring small negative eigenvalues")
  d <- eigJ$values
  d[d < r_tol] <- 0
  
  #-- 4) define joint neg-log-likelihood in s
  ztv <- crossprod(z_joint, eigJ$vectors)  # V^T z
  neglog <- function(s) {
    tmp <- (1 - s)*d + s
    0.5 * sum(log(tmp)) +
      0.5 * sum((ztv^2)/tmp)
  }
  
  #-- 5) optimize on [0,1]
  out <- optim(0.5, fn = neglog, method = "Brent", lower = 0, upper = 1)
  s_hat <- out$par
  s_hat
}


# ’ @title Joint kriging_rss for two ancestries
# ’
# ’ @description
# ’ Computes the conditional distribution of each variant’s z-score
# ’ given all others, and flags potential allele‐flip SNPs, for
# ’ two ancestries simultaneously.  We stack (z1,z2) and LD = blockdiag(R1,R2),
# ’ estimate a single s via “null-mle”, then compute Ω = ((1−s)LD + sI)^{-1},
# ’ conditional means/variances, fit a mixture‐of‐normals to residuals,
# ’ and compute likelihood ratios.
# ’
# ’ @param z1 Numeric vector of z‐scores from ancestry 1.
# ’ @param z2 Numeric vector of z‐scores from ancestry 2.
# ’ @param R1 Correlation matrix for ancestry 1.
# ’ @param R2 Correlation matrix for ancestry 2.
# ’ @param n1 Sample size for ancestry 1 (optional).
# ’ @param n2 Sample size for ancestry 2 (optional).
# ’ @param r_tol Eigenvalue tolerance for PSD enforcement (default 1e-8).
# ’ @param s Regularization parameter; by default estimated jointly.
# ’
# ’ @return A list with components:
# ’   * `plot`: a ggplot scatter of observed vs expected z’s (colored by ancestry),
# ’     with allele‐flip candidates in red.
# ’   * `conditional_dist`: a data.frame with columns
# ’     `ancestry, z, condmean, condvar, z_std_diff, logLR`.
# ’
# ’ @importFrom Matrix bdiag
# ’ @importFrom stats optim dnorm
# ’ @importFrom ggplot2 ggplot aes geom_point geom_abline theme_bw labs
# ’ @importFrom mixsqp mixsqp
# ’ @export


kriging_rss_joint<- function(z1, z2, R1, R2,
                                  n1 = NULL, n2 = NULL,
                                  r_tol = 1e-8,
                                  s = NULL) {
  #-- 1) finite‐sample shrinkage & NA‐impute
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be > 1")
      sigma2 <- (n - 1)/(z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)
  
  #-- 2) stack z and block‐diag LD
  z  <- c(z1, z2)
  R  <- as.matrix(Matrix::bdiag(R1, R2))
  
  #-- 3) estimate s if missing
  if (is.null(s)) {
    s <- estimate_s_rss_joint(z1, z2, R1, R2, n1, n2, r_tol)
  }
  force(s)
  if (s < 0) stop("s must be non‐negative")
  if (s > 1) {
    warning("s > 1; clamping to 0.8")
    s <- 0.8
  }
  
  #-- 4) eigen‐decompose & enforce PSD
  if (is.null(attr(R, "eigen")))
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  eig <- attr(R, "eigen")
  d   <- eig$values
  if (any(d < -r_tol))
    warning("R not PSD: flooring small negative eigenvalues")
  d[d < r_tol] <- 0
  
  #-- 5) build precision matrix Ω = [(1-s)R + sI]^{-1}
  dinv      <- 1/((1 - s)*d + s)
  dinv[is.infinite(dinv)] <- 0 ## Here we are using the psedo-inverse (if d_i = 0, then 1/d_i = 0, if d_i > 0, then 1/d_i)
  precision <- eig$vectors %*% (t(eig$vectors) * dinv)
  
  #-- 6) compute conditional means & variances
  p_total <- length(z)
  condmean <- numeric(p_total)
  condvar  <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- -(1/precision[i,i]) * precision[i,-i] %*% z[-i]
    condvar[i]  <- 1/precision[i,i]
  }
  z_std_diff <- (z - condmean)/sqrt(condvar)
  
  #-- 7) fit mixture‐of‐normals on residuals
  a_min <- 0.8
  a_max <- if (max(z_std_diff^2) < 1) 2 else 2*sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max/a_min)/log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max
  
  sd_mtx      <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  w           <- mixsqp::mixsqp(matrix_llik, log = TRUE,
                                control = list(verbose = FALSE))$x
  
  #-- 8) compute log‐likelihood ratios for allele‐flip
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  # rebuild for flipped sign
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix    <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix    <- logl1mix - logl0mix
  
  
  # 7) split into LR1, LR2, form joint
  p1 <- length(z1)
  logLR1       <- logLRmix[1:p1]
  logLR2       <- logLRmix[(p1+1):length(z)]
  logLR_joint  <- logLR1 + logLR2
  
  
  #-- 9) assemble results & plot
  ancestry <- rep(c("anc1","anc2"), times = c(length(z1), length(z2)))
  res      <- data.frame(
    ancestry   = ancestry,
    z          = z,
    condmean   = condmean,
    condvar    = condvar,
    z_std_diff = z_std_diff,
    logLR      = logLRmix,
    logLR_joint= logLR_joint
  )
  
  p <- ggplot2::ggplot(res, ggplot2::aes(x = condmean, y = z, color = ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected z", y = "Observed z")
  
  # — now flag per‐ancestry
  p1 <- length(z1)
  # ancestry1 positions 1…p1
  idx_flag1 <- which(logLRmix[1:p1] > 2 & abs(z_std_diff[1:(p1)])>4)
  
  # ancestry2 positions (p1+1)…end
  idx_flag2  <- which(logLRmix[(p1+1):length(z)] > 2 &
                        abs(z_std_diff[(p1+1):length(z)]) > 4)
  idx_flag_joint  <- which(logLR_joint > 2 & abs(z_std_diff) > 4)
  
  
  # highlight in red
  if (length(idx_flag1) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[idx_flag1, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  if (length(idx_flag2) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[p1+idx_flag2, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  if (length(idx_flag_joint) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[idx_flag_joint, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  list(plot = p, conditional_dist = res)
}









kriging_rss_joint_ori <- function(z1, z2, R1, R2,
                              n1 = NULL, n2 = NULL,
                              r_tol = 1e-8,
                              s = NULL) {
  #-- 1) finite‐sample shrinkage & NA‐impute
  shrink <- function(z, n) {
    if (!is.null(n)) {
      if (n <= 1) stop("n must be > 1")
      sigma2 <- (n - 1)/(z^2 + n - 2)
      z <- sqrt(sigma2) * z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)
  
  #-- 2) stack z and block‐diag LD
  z  <- c(z1, z2)
  R  <- as.matrix(Matrix::bdiag(R1, R2))
  
  #-- 3) estimate s if missing
  if (is.null(s)) {
    s <- estimate_s_rss_joint(z1, z2, R1, R2, n1, n2, r_tol)
  }
  force(s)
  if (s < 0) stop("s must be non‐negative")
  if (s > 1) {
    warning("s > 1; clamping to 0.8")
    s <- 0.8
  }
  
  #-- 4) eigen‐decompose & enforce PSD
  if (is.null(attr(R, "eigen")))
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
  eig <- attr(R, "eigen")
  d   <- eig$values
  if (any(d < -r_tol))
    warning("R not PSD: flooring small negative eigenvalues")
  d[d < r_tol] <- 0
  
  #-- 5) build precision matrix Ω = [(1-s)R + sI]^{-1}
  dinv      <- 1/((1 - s)*d + s)
  dinv[is.infinite(dinv)] <- 0 ## Here we are using the psedo-inverse (if d_i = 0, then 1/d_i = 0, if d_i > 0, then 1/d_i)
  precision <- eig$vectors %*% (t(eig$vectors) * dinv)
  
  #-- 6) compute conditional means & variances
  p_total <- length(z)
  condmean <- numeric(p_total)
  condvar  <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- -(1/precision[i,i]) * precision[i,-i] %*% z[-i]
    condvar[i]  <- 1/precision[i,i]
  }
  z_std_diff <- (z - condmean)/sqrt(condvar)
  
  #-- 7) fit mixture‐of‐normals on residuals
  a_min <- 0.8
  a_max <- if (max(z_std_diff^2) < 1) 2 else 2*sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max/a_min)/log2(1.05))
  a_grid <- 1.05^(seq(-npoint, 0)) * a_max
  
  sd_mtx      <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  w           <- mixsqp::mixsqp(matrix_llik, log = TRUE,
                                control = list(verbose = FALSE))$x
  
  #-- 8) compute log‐likelihood ratios for allele‐flip
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  # rebuild for flipped sign
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik, 1, max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix    <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix    <- logl1mix - logl0mix
  
  #-- 9) assemble results & plot
  ancestry <- rep(c("anc1","anc2"), times = c(length(z1), length(z2)))
  res      <- data.frame(
    ancestry   = ancestry,
    z          = z,
    condmean   = condmean,
    condvar    = condvar,
    z_std_diff = z_std_diff,
    logLR      = logLRmix
  )
  
  p <- ggplot2::ggplot(res, ggplot2::aes(x = condmean, y = z, color = ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected z", y = "Observed z")
  
  # — now flag per‐ancestry
  p1 <- length(z1)
  # ancestry1 positions 1…p1
  idx_flag1 <- which(logLRmix[1:p1] > 2 & abs(z[1:p1]) > 2)
  # ancestry2 positions (p1+1)…end
  idx_flag2  <- which(logLRmix[(p1+1):length(z)] > 2 &
                       abs(z[(p1+1):length(z)]) > 2)
  
  # highlight in red
  if (length(idx_flag1) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[idx_flag1, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  if (length(idx_flag2) > 0) {
    p <- p + ggplot2::geom_point(
      data = res[p1+idx_flag2, ],
      ggplot2::aes(x = condmean, y = z),
      color = "red"
    )
  }
  list(plot = p, conditional_dist = res)
}




estimate_s_rss_two <- function(z1, z2,
                               R1, R2,
                               n1 = NULL, n2 = NULL,
                               r_tol = 1e-8) {
  # 1) finite-sample shrinkage helper
  shrink <- function(z,n) {
    if (!is.null(n)) {
      if(n<=1) stop("n must be >1")
      sigma2 <- (n-1)/(z^2 + n - 2)
      z <- sqrt(sigma2)*z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1,n1)
  z2 <- shrink(z2,n2)
  
  # 2) eigen-decomps (cached)
  prep <- function(R) {
    if(is.null(attr(R,"eigen")))
      attr(R,"eigen") <- eigen(R, symmetric=TRUE)
    eig <- attr(R,"eigen")
    d <- eig$values
    if(any(d < -r_tol)) warning("flooring small negative eigs")
    d[d < r_tol] <- 0
    list(V = eig$vectors, d = d)
  }
  e1 <- prep(R1)
  e2 <- prep(R2)
  
  # 3) precompute V^T z
  y1 <- crossprod(z1, e1$V)
  y2 <- crossprod(z2, e2$V)
  
  # 4) joint neg-log-lik
  joint_neglog <- function(par) {
    λ1 <- par[1]; λ2 <- par[2]
    tmp1 <- pmax((1-λ1)*e1$d + λ1, 1e-6)
    tmp2 <- pmax((1-λ2)*e2$d + λ2, 1e-6)
    0.5*(sum(log(tmp1)) + sum((y1^2)/tmp1)
         + sum(log(tmp2)) + sum((y2^2)/tmp2))
  }
  
  # 5) optimize over [0,1]^2
  out <- optim(c(0.5,0.5), joint_neglog,
               method="L-BFGS-B",
               lower = c(1e-6, 1e-6),
               upper = c(1, 1))
  
  
  c(lambda1 = out$par[1], lambda2 = out$par[2])
}



estimate_s_rss_two_penalty <- function(z1, z2,
                               R1, R2,
                               n1 = NULL, n2 = NULL,
                               rho = 0,
                               alpha = 1,
                               r_tol = 1e-8) {
  # finite-sample shrinkage & NA-impute
  shrink <- function(z,n) {
    if(!is.null(n)){
      if(n<=1) stop("n must be >1")
      sigma2 <- (n-1)/(z^2 + n - 2)
      z <- sqrt(sigma2)*z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1,n1); z2 <- shrink(z2,n2)
  
  # prepare eigen-decomps
  prep <- function(R) {
    if(is.null(attr(R,"eigen"))) attr(R,"eigen") <- eigen(R, symmetric=TRUE)
    eig <- attr(R,"eigen")
    d <- eig$values
    if(any(d< -r_tol)) warning("flooring small negative eigenvalues")
    d[d<r_tol] <- 0
    list(V=eig$vectors, d=d)
  }
  e1 <- prep(R1); e2 <- prep(R2)
  
  # precompute projections
  y1 <- crossprod(z1, e1$V)
  y2 <- crossprod(z2, e2$V)
  
  # negative log-likelihood sum
  neglik <- function(lam1, lam2) {
    tmp1 <- (1-lam1)*e1$d + lam1
    tmp2 <- (1-lam2)*e2$d + lam2
    0.5*( sum(log(tmp1)) + sum((y1^2)/tmp1)
          + sum(log(tmp2)) + sum((y2^2)/tmp2) )
  }
  
  # negative log-prior for Beta(alpha,alpha)
  neglog_prior <- function(lam) {
    if(alpha==1) return(0)
    -(alpha-1)*(log(lam)+log(1-lam))
  }
  
  # joint objective: neg-loglik + penalty + neg-log-prior
  obj <- function(par) {
    lam1 <- par[1]; lam2 <- par[2]
    val <- neglik(lam1, lam2)
    # penalty
    if(rho>0) val <- val + rho*(lam1-lam2)^2
    # prior penalties
    val + neglog_prior(lam1) + neglog_prior(lam2)
  }
  
  # optimize on (0,1) with bounds away from edges
  lower <- rep(1e-6,2)
  upper <- rep(1-1e-6,2)
  out <- optim(c(0.5,0.5), fn=obj,
               method="L-BFGS-B",
               lower=lower, upper=upper)
  c(lambda1=out$par[1], lambda2=out$par[2])
}



#’ @title Joint kriging_rss with separate shrinkage λ₁,λ₂
#’
#’ @description
#’ Runs conditional‐mean/variance diagnostics for two ancestries,
#’ but first estimates (λ₁,λ₂) with a penalty and Beta(α,α) prior.
#’
#’ @param z1,z2      Numeric vectors of z‐scores for ancestry 1 and 2.
#’ @param R1,R2      Correlation matrices for the two ancestries.
#’ @param n1,n2      Sample sizes (optional) for finite‐sample z‐shrinkage.
#’ @param rho        Penalty strength on (λ₁−λ₂)² (default 0 = independent).
#’ @param alpha      Beta(α,α) prior shape (default α=1 = flat).
#’ @param r_tol      Eigenvalue tolerance (default 1e-8).
#’
#’ @return A list with  
#’   - `plot`: ggplot of observed vs expected z’s colored by ancestry,  
#’   - `conditional_dist`: data.frame of ancestry, z, condmean, condvar, z_std_diff, logLR.  
#’
#’ @importFrom Matrix bdiag
#’ @importFrom stats optim dnorm
#’ @importFrom ggplot2 ggplot aes geom_point geom_abline theme_bw labs
#’ @importFrom mixsqp mixsqp
#’ @export
kriging_rss_two <- function(z1, z2, R1, R2,
                              n1 = NULL, n2 = NULL,
                              rho   = 0,
                              alpha = 1,
                              r_tol = 1e-8) {
  # 1) finite-sample shrinkage
  shrink <- function(z,n) {
    if(!is.null(n)) {
      if(n<=1) stop("n must be > 1")
      sigma2 <- (n-1)/(z^2 + n - 2)
      z <- sqrt(sigma2)*z
    }
    z[is.na(z)] <- 0
    z
  }
  z1 <- shrink(z1, n1)
  z2 <- shrink(z2, n2)
  
  # 2) estimate (λ1,λ2) with penalty+prior
  lambdas <- estimate_s_rss_two_penalty(
    z1, z2, R1, R2, n1, n2,
    rho   = rho,
    alpha = alpha,
    r_tol = r_tol
  )
  λ1 <- lambdas["lambda1"]
  λ2 <- lambdas["lambda2"]
  
  # 3) build per‐ancestry Σ and Ω
  p1 <- length(z1); p2 <- length(z2)
  I1 <- diag(p1); I2 <- diag(p2)
  
  Σ1 <- (1-λ1)*R1 + λ1*I1
  Σ2 <- (1-λ2)*R2 + λ2*I2
  
  # enforce PSD if needed (small eigenvalues)
  fixPSD <- function(M) {
    eig <- eigen(M, symmetric=TRUE)
    d   <- eig$values
    d[d < r_tol] <- 0
    eig$vectors %*% diag(d) %*% t(eig$vectors)
  }
  Σ1 <- fixPSD(Σ1)
  Σ2 <- fixPSD(Σ2)
  
  Ω1 <- solve(Σ1)
  Ω2 <- solve(Σ2)
  
  # 4) stack Ω and z
  Ω  <- as.matrix(Matrix::bdiag(Ω1, Ω2))
  z  <- c(z1, z2)
  
  # 5) compute conditional means & variances
  p_total  <- p1 + p2
  condmean <- numeric(p_total)
  condvar  <- numeric(p_total)
  for (i in seq_len(p_total)) {
    condmean[i] <- - (1/Ω[i,i]) * (Ω[i,-i] %*% z[-i])
    condvar[i]  <-  1/Ω[i,i]
  }
  z_std_diff <- (z - condmean)/sqrt(condvar)
  
  # 6) fit mixture‐of‐normals on residuals (no change)
  a_min <- 0.8
  a_max <- if(max(z_std_diff^2)<1) 2 else 2*sqrt(max(z_std_diff^2))
  npoint <- ceiling(log2(a_max/a_min)/log2(1.05))
  a_grid <- 1.05^(seq(-npoint,0))*a_max
  
  sd_mtx      <- outer(sqrt(condvar), a_grid)
  matrix_llik <- stats::dnorm(z - condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik,1,max)
  matrix_llik <- matrix_llik - lfactors
  w           <- mixsqp::mixsqp(matrix_llik, log=TRUE,
                                control=list(verbose=FALSE))$x
  
  # 7) log‐likelihood‐ratio for allele‐flip
  logl0mix <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  matrix_llik <- stats::dnorm(z + condmean, sd = sd_mtx, log = TRUE)
  lfactors    <- apply(matrix_llik,1,max)
  matrix_llik <- matrix_llik - lfactors
  logl1mix    <- drop(log(exp(matrix_llik) %*% (w + 1e-15))) + lfactors
  logLRmix    <- logl1mix - logl0mix
  
  # 8) assemble results & plot
  ancestry <- rep(c("anc1","anc2"), times = c(p1,p2))
  res <- data.frame(
    ancestry   = ancestry,
    z          = z,
    condmean   = condmean,
    condvar    = condvar,
    z_std_diff = z_std_diff,
    logLR      = logLRmix
  )
  p <- ggplot2::ggplot(res, ggplot2::aes(x=condmean,y=z,color=ancestry)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept=0,slope=1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Expected z",y="Observed z")
  
  # 9) highlight per‐ancestry flips
  idx1 <- which(logLRmix[1:p1] > 2 & abs(z[1:p1]) > 2)
  idx2 <- which(logLRmix[(p1+1):p_total] > 2 &
                  abs(z[(p1+1):p_total]) > 2) #+ p1
  
  if(length(idx1)) p <- p + ggplot2::geom_point(
    data = res[idx1, ], ggplot2::aes(x=condmean,y=z), color="red")
  if(length(idx2)) p <- p + ggplot2::geom_point(
    data = res[idx2, ], ggplot2::aes(x=condmean,y=z), color="red")
  
  list(plot = p, conditional_dist = res)
}


