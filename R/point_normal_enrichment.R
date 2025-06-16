# The point-normal family uses the ashr class normalmix.
#
pne_checkg <- function(g_init, fix_g, mode, scale, pointmass, call) {
  check_g_init(g_init = g_init,
               fix_g = fix_g,
               mode = mode,
               scale = scale,
               pointmass = pointmass,
               call = call,
               class_name = "normalmix",
               scale_name = "sd")
}


# Point-normal parameters are alpha = logit(pi0), beta = log(s2), and mu.
#
pne_initpar <- function(g_init, mode, scale, pointmass, x, s) {
  if (!is.null(g_init) && length(g_init$pi) == 1) {
    par <- list(alpha = -Inf,
                beta = 2 * log(g_init$sd),
                mu = g_init$mean)
  } else if (!is.null(g_init) && length(g_init$pi) == 2) {
    par <- list(alpha = -log(1 / g_init$pi[1] - 1),
                beta = 2 * log(g_init$sd[2]),
                mu = g_init$mean[1])
  } else {
    par <- list()
    if (!pointmass) {
      par$alpha <- -Inf
    } else {
      par$alpha <- 0 # default
    }
    if (!identical(scale, "estimate")) {
      if (length(scale) != 1) {
        stop("Argument 'scale' must be either 'estimate' or a scalar.")
      }
      par$beta <- 2 * log(scale)
    } else {
      par$beta <- log(mean(x^2)) # default
    }
    if (!identical(mode, "estimate")) {
      par$mu <- mode
    } else {
      par$mu <- mean(x) # default
    }
  }

  return(par)
}


pne_scalepar <- function(par, scale_factor) {
  if (!is.null(par$beta)) {
    par$beta <- par$beta + 2 * log(scale_factor)
  }
  if (!is.null(par$mu)) {
    par$mu <- scale_factor * par$mu
  }

  return(par)
}


# Precomputations.
#
pne_precomp <- function(x, s, annot, par_init, fix_par) {
  fix_mu  <- fix_par[3]

  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  if (any(s == 0)) {
    which_s0 <- which(s == 0)
    which_x_nz <- which(x[which_s0] != par_init$mu)
    n0 <- length(which_s0) - length(which_x_nz)
    n1 <- length(which_x_nz)
    sum1 <- sum((x[which_s0[which_x_nz]] - par_init$mu)^2)
    x <- x[-which_s0]
    s <- s[-which_s0]
  } else {
    n0 <- 0
    n1 <- 0
    sum1 <- 0
  }
  n2 <- length(x)

  s2 <- s^2

  if (fix_mu) {
    z <- (x - par_init$mu)^2 / s2
    sum_z <- sum(z)
  } else {
    z <- NULL
    sum_z <- NULL
  }

  # add intercept for annotation matrix
  annot <- cbind(rep(1, length = n2),
                 annot)


  return(list(n0 = n0, n1 = n1, sum1 = sum1, n2 = n2, s2 = s2, z = z, sum_z = sum_z, annot = annot))
}


# The objective function (the negative log likelihood, minus a constant).
# k annotations including intercept
# (b_0, b_1, etc) first k  parameters
# log(sigma^2)
# mu
# par_init <- list(bs, beta, mu)
# bs has length k
pne_nllik <- function(par, x, s, par_init, fix_par,
                     n0, n1, sum1, n2, s2, z, sum_z, annot,
                     calc_grad, calc_hess) {

  k <- ncol(annot)

  fix_bs <- fix_par[1:k]
  fix_pi0 <- all(fix_bs)

  fix_s2  <- fix_par[k + 1]
  fix_mu  <- fix_par[k + 2]

  i <- 1
  while(i <= k){
    if(fix_bs[i]){
      bs[i] <- par_init$bs[i]
    } else {
      bs[i] <- par[i]
    }
    i <- i + 1
  }
  alpha <- rowSums(t(t(annot)*bs)) # logit(pi0)

  if (fix_s2) {
    beta <- par_init$beta
  } else {
    beta <- par[i]
    i <- i + 1
  }
  if (fix_mu) {
    mu <- par_init$mu
  } else {
    mu <- par[i]
    z <- (x - mu)^2 / s2
    sum_z <- sum(z)
  }

  logist.alpha  <- 1 / (1 + exp(-alpha)) # vector length n2
  logist.nalpha <- 1 / (1 + exp(alpha))

  logist.beta   <- 1 / (1 + s2 * exp(-beta)) # scalar or vector
  logist.nbeta  <- 1 / (1 + exp(beta) / s2)

  y <- 0.5 * (z * logist.beta + log(logist.nbeta)) # vector

  # Negative log likelihood.
  C <- pmax(y, alpha)
  if (n0 == 0 || logist.alpha == 0) {
    nllik <- 0
  } else {
    nllik <- -n0 * log(logist.alpha)
  }
  nllik <- nllik - sum(log(logist.nalpha))
  if (n1 > 0) {
    nllik <- nllik + 0.5 * n1 * beta
  }
  if (sum1 > 0) {
    nllik <- nllik + 0.5 * sum1 * exp(-beta)
  }
  nllik <- nllik + 0.5 * sum_z - sum(log(exp(y - C) + exp(alpha - C)) + C)

  if (calc_grad || calc_hess) {
    dlogist.beta  <- logist.beta * logist.nbeta

    logist.y  <- 1 / (1 + exp(alpha - y)) # vector
    logist.ny <- 1 / (1 + exp(y - alpha)) # also a vector

    # Gradient.
    grad <- numeric(length(par))
    i <- 1
    if (!fix_pi0) {
      grad_pi0 <- (sum(logist.alpha) - sum(logist.ny))*colSums(annot) # vector length k
      if(any(fix_bs)){
        grad_pi0[fix_bs == TRUE] <- 0
      }
      grad <- grad_pi0
      i <- i + k
    }
    if (!fix_s2) {
      dy.dbeta <- 0.5 * (z * dlogist.beta - logist.beta)
      grad[i] <- 0.5 * (n1 - sum1 * exp(-beta)) - sum(logist.y * dy.dbeta)
      i <- i + 1
    }
    if (!fix_mu) {
      dy.dmu <- (mu - x) * logist.beta / s2
      grad[i] <- sum((mu - x) / s2) - sum(logist.y * dy.dmu)
    }
    attr(nllik, "gradient") <- grad
  }

  if (calc_hess) {
    dlogist.alpha <- logist.alpha * logist.nalpha # vector
    dlogist.y <- logist.y * logist.ny

    # Hessian.
    hess <- matrix(nrow = length(par), ncol = length(par))
    i <- 1
    ### need to update. not ready.
    if (!fix_pi0) {
      tmp <- sum(dlogist.alpha)
      hess[i, i] <- tmp - sum(dlogist.y)
      j <- i + 1
      if (!fix_s2) {
        hess[i, j] <- hess[j, i] <- sum(dlogist.y * dy.dbeta)
        j <- j + 1
      }
      if (!fix_mu) {
        hess[i, j] <- hess[j, i] <- sum(dlogist.y * dy.dmu)
      }
      i <- i + 1
    }
    if (!fix_s2) {
      d2y.dbeta2 <- 0.5 * ((z * (logist.nbeta - logist.beta) - 1) * dlogist.beta)
      tmp <- 0.5 * sum1 * exp(-beta) - sum(dlogist.y * dy.dbeta^2)
      hess[i, i] <- tmp - sum(logist.y * d2y.dbeta2)
      j <- i + 1
      if (!fix_mu) {
        d2y.dbetadmu <- (mu - x) * dlogist.beta / s2
        tmp <- -sum(dlogist.y * dy.dbeta * dy.dmu)
        hess[i, j] <- hess[j, i] <- tmp - sum(logist.y * d2y.dbetadmu)
      }
      i <- i + 1
    }
    if (!fix_mu) {
      tmp <- sum(1 / s2 - dlogist.y * dy.dmu^2)
      hess[i, i] <- tmp - sum(logist.y * logist.beta / s2)
    }
    attr(nllik, "hessian") <- hess
  }

  return(nllik)
}


# Postcomputations. A constant was subtracted from the log likelihood and needs
#   to be added back in. We also check boundary solutions here.
#
pne_postcomp <- function(optpar, optval, x, s, par_init, fix_par, scale_factor,
                        n0, n1, sum1, n2, s2, z, sum_z) {
  llik <- pn_llik_from_optval(optval, n1, n2, s2)
  retlist <- list(par = optpar, val = llik)

  # Check the solution pi0 = 1.
  fix_pi0 <- fix_par[1]
  fix_mu  <- fix_par[3]
  if (!fix_pi0 && fix_mu) {
    pi0_llik <- sum(-0.5 * log(2 * pi * s^2) - 0.5 * (x - par_init$mu)^2 / s^2)
    pi0_llik <- pi0_llik + sum(is.finite(x)) * log(scale_factor)
    if (pi0_llik > llik) {
      retlist$par$alpha <- Inf
      retlist$par$beta <- 0
      retlist$val <- pi0_llik
    }
  }

  return(retlist)
}

pn_llik_from_optval <- function(optval, n1, n2, s2) {
  if (length(s2) == 1) {
    sum.log.s2 <- n2 * log(s2)
  } else {
    sum.log.s2 <- sum(log(s2))
  }

  return(-optval - 0.5 * ((n1 + n2) * log(2 * pi) + sum.log.s2))
}


# Summary results.
#
pne_summres <- function(x, s, optpar, output) {
  w  <- 1 - 1 / (exp(-optpar$alpha) + 1)
  a  <- exp(-optpar$beta)
  mu <- optpar$mu

  return(pn_summres_untransformed(x, s, w, a, mu, output))
}

pn_summres_untransformed <- function(x, s, w, a, mu, output) {
  wpost <- wpost_normal(x, s, w, a, mu)
  pmean_cond <- pmean_cond_normal(x, s, a, mu)
  pvar_cond <- pvar_cond_normal(s, a)

  posterior <- list()

  if (result_in_output(output)) {
    posterior$mean  <- wpost * pmean_cond + (1 - wpost) * mu
    posterior$mean2 <- wpost * (pmean_cond^2 + pvar_cond) + (1 - wpost) * (mu^2)
    posterior$mean2 <- pmax(posterior$mean2, posterior$mean^2)
    posterior$sd    <- sqrt(pmax(0, posterior$mean2 - posterior$mean^2))
  }

  if (lfsr_in_output(output)) {
    posterior$lfsr  <- (1 - wpost) + wpost * pnorm(0, abs(pmean_cond), sqrt(pvar_cond))
  }

  return(posterior)
}


