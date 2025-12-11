#' Add sampler to an ebnm_object
#'
#' Adds a posterior sampler to a fitted \code{\link{ebnm}} object.
#'
#' @param ebnm_res The fitted \code{ebnm} object.
#'
#' @return The \code{ebnm} object with an additional field
#'   \code{posterior_sampler}.
#'
#' @export
#'
ebnm_add_sampler <- function(ebnm_res) {
  if (!inherits(ebnm_res, "ebnm")) {
    stop("Input argument must be an instance of class \"ebnm\".")
  }
  if (is.null(ebnm_res[[data_ret_str()]])) {
    stop("Data not found in ebnm object. Sampler cannot be added.")
  }
  if (is.null(ebnm_res[[g_ret_str()]])) {
    warning("Fitted prior not found in ebnm object. Sampler cannot be added.")
    incl_cdf <- FALSE
  }

  x <- ebnm_res[[data_ret_str()]][[obs_ret_str()]]
  names(x) <- rownames(ebnm_res[[data_ret_str()]])
  sampler <- ebnm(
    x = x,
    s = ebnm_res[[data_ret_str()]][[se_ret_str()]],
    g_init = ebnm_res[[g_ret_str()]],
    fix_g = TRUE,
    output = samp_arg_str()
  )
  ebnm_res[[samp_ret_str()]] <- sampler[[samp_ret_str()]]
  return(ebnm_res)
}


#'@export
ebnm_add_sampler_ep <- function(ebnm_res){
  if (!inherits(ebnm_res, "ebnm")) {
    stop("Input argument must be an instance of class \"ebnm\".")
  }
  if (is.null(ebnm_res[[data_ret_str()]])) {
    stop("Data not found in ebnm object. Sampler cannot be added.")
  }
  if (is.null(ebnm_res[["nllik_hess"]])) {
    stop("Hessian not found in ebnm object. Sampler cannot be added.")
  }
  if (is.null(ebnm_res[[g_ret_str()]])) {
    warning("Fitted prior not found in ebnm object. Sampler cannot be added.")
    incl_cdf <- FALSE
  }


  x <- ebnm_res[[data_ret_str()]][[obs_ret_str()]]
  names(x) <- rownames(ebnm_res[[data_ret_str()]])
  par <- ebnm_res[["par"]]
  npar <- length(par)
  fix_par <- ebnm_res[["fix_par"]]
  H <- ebnm_res[["nllik_hess"]]
  postsamp_fn <- ebnm_res[["postsamp_fn"]]

  nfix_par_i <- which(!fix_par)
  stopifnot(length(nfix_par_i) > 0)
  npar_nfix <- length(nfix_par_i)

  eig_I <- eigen(solve(H))
  I_sqrt_t <- with(eig_I, t(vectors) * sqrt(values))
  samplerep <- function(nsamp) {
    par_samp <- matrix(0, nrow = nsamp, ncol = npar)
    par_nfix_samp <- matrix(rnorm(n = npar_nfix*nsamp), nrow = nsamp)
    par_nfix_samp <- par_nfix_samp %*% I_sqrt_t
    par_samp[, nfix_par_i] <- par_nfix_samp
    par_samp <- t(t(par_samp) + unlist(par))
    samp <- sapply(seq(nsamp), function(i){
      my_par <- as.list(par_samp[i,])
      names(my_par) <- names(par)
      postsamp_fn(x, ebnm_res[[data_ret_str()]][[se_ret_str()]], my_par, 1)
    }) |> t()
    colnames(samp) <- names(x)
    return(samp)
    #return(list(samples = samp, params = par_samp))
  }

  ebnm_res[[sampep_ret_str()]] <- samplerep
  return(ebnm_res)
}
