### test likelihood calculations
set.seed(2025)
theta <- rnorm(n = 10)
x <- rnorm(n = 10, mean = theta)
s <- rep(1, 10)

## params:
# alpha = logit(pi0)
# beta = log(sigma^2)
# mu
pn_likelihood_by_hand <- function(params, x, s){
  pi0 <- plogis(params[1])
  sigma_2 <- exp(params[2])
  mu <- params[3]
  n <- length(x)

  lik_part1 <- pi0*dnorm(x = x, mean = mu, sd = s)
  lik_part2 <- (1-pi0)*dnorm(x = x, mean = mu, sd = sqrt(s^2 + sigma_2))
  ll <- sum(log(lik_part1 + lik_part2))
  return(ll)
}


logit <- function(x){log(x)- log(1-x)}
params <- c(logit(0.2), log(0.5), 0)
pn_likelihood_by_hand(params = params,
                      x = x, s = s)

par_init <- list(alpha = 0, beta = 0, mu = 0)
precomp <- pn_precomp(x, s, par_init = par_init, fix_par = rep(FALSE, 3))


ll_first_part <- pn_nllik(par = params, x, s, par_init = par_init,
                          fix_par = rep(FALSE, 3),
                          n0 = precomp$n0,
                          n1 = precomp$n1,
                          sum1 = precomp$sum1,
                          n2 = precomp$n2,
                          s2= precomp$s2,
                          z = precomp$z,
                          sum_z = precomp$sum_z,
                          calc_grad = FALSE, calc_hess = FALSE)



pn_llik_from_optval(ll_first_part, precomp$n1, precomp$n2, precomp$s2)


myfunc <- function(params){
  -1*pn_likelihood_by_hand(params, x, s)
}
numDeriv::grad(myfunc, x = params)

ll_first_part <- pn_nllik(par = params, x, s, par_init = par_init,
                          fix_par = rep(FALSE, 3),
                          n0 = precomp$n0,
                          n1 = precomp$n1,
                          sum1 = precomp$sum1,
                          n2 = precomp$n2,
                          s2= precomp$s2,
                          z = precomp$z,
                          sum_z = precomp$sum_z,
                          calc_grad = TRUE, calc_hess = FALSE)



###### Test PNE


set.seed(2025)
theta <- rnorm(n = 10)
x <- rnorm(n = 10, mean = theta)
s <- rep(1, 10)
A <- rbinom(n= 10, size = 1, prob = 0.5)


## params:
# b0
# b1
# beta = log(sigma^2)
# mu
pne_likelihood_by_hand <- function(params, x, s, annot){
  b0 <- params[1]
  b1 <- params[2]
  pi0 <- plogis(b0 + b1*annot)
  sigma_2 <- exp(params[3])
  mu <- params[4]
  n <- length(x)

  lik_part1 <- pi0*dnorm(x = x, mean = 0, sd = s)
  lik_part2 <- (1-pi0)*dnorm(x = x, mean = mu, sd = sqrt(s^2 + sigma_2))
  ll <- sum(log(lik_part1 + lik_part2))
  return(ll)
}


logit <- function(x){log(x)- log(1-x)}
params <- c(0.11, 0.4, log(0.5), 0)
pne_likelihood_by_hand(params = params,
                      x = x, s = s, annot = A)

par_init <- list(bs = c(0, 0), beta = 0, mu = 0)
precomp <- pne_precomp(x, s, par_init = par_init, fix_par = rep(FALSE, 4),
                       annot = matrix(A, ncol = 1))


ll_first_part <- pne_nllik(par = params, x, s, par_init = par_init,
                          fix_par = rep(FALSE, 4),
                          n0 = precomp$n0,
                          n1 = precomp$n1,
                          sum1 = precomp$sum1,
                          n2 = precomp$n2,
                          s2= precomp$s2,
                          z = precomp$z,
                          sum_z = precomp$sum_z,
                          annot = precomp$annot,
                          calc_grad = FALSE, calc_hess = FALSE)

pn_llik_from_optval(ll_first_part, precomp$n1, precomp$n2, precomp$s2)

myfunc <- function(params){
  -1*pne_likelihood_by_hand(params, x, s, annot = A)
}
numDeriv::grad(myfunc, x = params)

ll_first_part <- pne_nllik(par = params, x, s, par_init = par_init,
                           fix_par = rep(FALSE, 4),
                           n0 = precomp$n0,
                           n1 = precomp$n1,
                           sum1 = precomp$sum1,
                           n2 = precomp$n2,
                           s2= precomp$s2,
                           z = precomp$z,
                           sum_z = precomp$sum_z,
                           annot = precomp$annot,
                           calc_grad = TRUE, calc_hess = FALSE)
