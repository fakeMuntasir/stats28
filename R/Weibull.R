#' Title Weibull distibution parameters estimate.
#'
#' @param data data set / sample ( as a list ) to estimate the parameters and confidence interval.
#' @param cL confidence level.
#'
#' @return Point estimates the shape and scale parameter, Confidence interval for the shape and scale parameter.
#' @export
#'
#' @examples
#' (rweibull(0.5,.5),.95)
#' (rweibull(10,2),.90)
Weibull <- function(data,cL){


  weibull_log_likelihood <- function(params, x) {
    shape <- params[1]
    scale <- params[2]

    n <- length(x)
    log_likelihood <- -sum(log(dweibull(x, shape = shape, scale = scale)))

    return(log_likelihood)  # Return log-likelihood (not negative)
  }

  initial_params <- c(shape = 1, scale = exponential(data,cL)$Lamda)


  mle_result <- optim(par = initial_params, fn = weibull_log_likelihood,
                      x = data, method = "L-BFGS-B", hessian = T)

  shape_estimate <- mle_result$par[1]
  scale_estimate <- mle_result$par[2]

  hessian_matrix <- mle_result$hessian #numDeriv::hessian(weibull_log_likelihood, params = mle_result$par, x = data)
  cov_matrix <- solve(hessian_matrix)
  se_shape <- sqrt(cov_matrix[1, 1])
  se_scale <- sqrt(cov_matrix[2, 2])

  z_value <-qnorm((1-cL)/2,lower.tail = F)

  ci_shape <- c(shape_estimate - z_value * se_shape, shape_estimate + z_value * se_shape)
  ci_scale <- c(scale_estimate - z_value * se_scale, scale_estimate + z_value * se_scale)

  result <- c( shape_estimate,
               scale_estimate,
               CI_shape_lower = ci_shape[1],
               CI_shape_upper = ci_shape[2],
               CI_scale_lower = ci_scale[1],
               CI_scale_upper = ci_scale[2]
  )

  return(result)


}
