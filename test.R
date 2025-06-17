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
pne_likelihood_by_hand(params = c(params[1:2], c(0, 0)),
                      x = x, s = s, annot = A)

par_init <- list(bs = c(0, 0), beta = 0, mu = 0)
precomp <- pne_precomp(x, s, par_init = par_init, fix_par =  c(FALSE, TRUE, TRUE),
                       annot = matrix(A, ncol = 1))


ll_first_part <- pne_nllik(par = params[1:2], x, s, par_init = par_init,
                          fix_par = c(FALSE, TRUE, TRUE),
                          n0 = precomp$n0,
                          n1 = precomp$n1,
                          sum1 = precomp$sum1,
                          n2 = precomp$n2,
                          s2= precomp$s2,
                          z = precomp$z,
                          sum_z = precomp$sum_z,
                          annot = precomp$annot,
                          k = precomp$k,
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




## Test PNE functions

set.seed(0616)
n <- 10000
annot <- matrix(rbinom(n = n*1, size = 1, prob = 0.3), nrow = n)
annot_with_intercept <- cbind(rep(1, n), annot)
betas <- c(1.5, -1)
pi0 <- plogis(rowSums(t(t(annot_with_intercept)*betas) ))
sd <- 10
Z <- rbinom(n = n, size = 1, prob = 1-pi0)
theta <- numeric(length = n)
theta[Z == 1] <- rnorm(n = sum(Z), mean = 0, sd = sd)
s <- rep(1, n)
x <- rnorm(n = n, mean = theta, sd = s)

###
par_init <- pne_initpar(NULL, 0, "estimate", TRUE, x, s, annot)
fix_par <- c(FALSE, FALSE, TRUE)
precomp <- pne_precomp(x, s, par_init = par_init, fix_par = fix_par,
                       annot = annot)
params <- unlist(par_init[!fix_par])

### MLE function

optres <- mle_parametric(x = x,
                         s = s,
                         annot = annot,
                         par_init = par_init,
                         fix_par = fix_par,
                         scalepar_fn = pn_scalepar,
                         precomp_fn = pne_precomp,
                         nllik_fn = pne_nllik,
                         postcomp_fn = pne_postcomp,
                         optmethod = "nlm",
                         control = list(),
                         use_grad = FALSE,
                         use_hess = FALSE)

true_params <- c(betas, log(sd), 0)
pne_nllik(par = unlist(optres$par[!fix_par]), x, s, par_init = par_init,
                           fix_par = fix_par,
                           n0 = precomp$n0,
                           n1 = precomp$n1,
                           sum1 = precomp$sum1,
                           n2 = precomp$n2,
                           s2= precomp$s2,
                           z = precomp$z,
                           sum_z = precomp$sum_z,
                           annot = precomp$annot,
                           k = precomp$k,
                           calc_grad = FALSE, calc_hess = FALSE) %>%
        pn_llik_from_optval(., precomp$n1, precomp$n2, precomp$s2)

pne_nllik(par = true_params[1:3], x, s, par_init = par_init,
          fix_par = fix_par,
          n0 = precomp$n0,
          n1 = precomp$n1,
          sum1 = precomp$sum1,
          n2 = precomp$n2,
          s2= precomp$s2,
          z = precomp$z,
          sum_z = precomp$sum_z,
          annot = precomp$annot,
          k = precomp$k,
          calc_grad = FALSE, calc_hess = FALSE) %>%
  pn_llik_from_optval(., precomp$n1, precomp$n2, precomp$s2)


ebnm_result <- ebnm(x = x, s = s, annot = annot,
                   #output = ebnm_output_all(),
                   prior_family = "point_normal_enrichment",
                   optmethod = "nograd_nlm")

ebnm_result$fitted_g

plot(theta, ebnm_result$posterior$mean)
