#' @title Power function
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param c                     A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-c,c)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to a probability.
power_TOST = function(alpha, theta, sigma_nu, nu, c, ...){
  tval = qt(1 - alpha, df = nu)
  mu1 = (theta + c)/sigma_nu
  mu2 = (theta - c)/sigma_nu
  R = (c*sqrt(nu))/(tval*sigma_nu)
  p1 = OwensQ(nu, tval, mu1, 0, R)
  p2 = OwensQ(nu, -tval, mu2, 0, R)
  pw = p2-p1
  pw[pw < 0] = 0
  pw
}

#' @title The size
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param c                     A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-c,c)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to a probability.
size_TOST = function(alpha, sigma_nu, nu, c, ...){
  power_TOST(alpha = alpha, theta = c, sigma_nu = sigma_nu,
             nu = nu, c = c)
}

#' @title Objective function to optimise
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param test                  A \code{numeric} value specifying the significance level to optimise.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param c                     A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-c,c)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value for the objective function.
obj_fun_alpha_star = function(test, alpha, sigma_nu, nu, c, ...){
  size = size_TOST(alpha = test, sigma_nu = sigma_nu, nu = nu, c = c)
  (size - alpha)^2
}

#' @title Get alpha star
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param c                     A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-c,c)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to the solution of the optimisation.
get_alpha_star = function(alpha, sigma_nu, nu, c, ...){
  out = optimize(obj_fun_alpha_star, c(alpha, 0.5),
                 alpha = alpha, sigma_nu = sigma_nu,
                 nu = nu, c = c)

  # size_out = size_TOST(alpha = out$minimum, sigma_nu = sigma_nu, nu = nu, c = c)

  out$minimum
}

#' @title Confidence Intervals
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param ...                   Additional parameters.
#' @return The function returns a numerical \code{vector} with the lower and upper bound of the confidence intervals.
ci = function(alpha, theta, sigma_nu, nu, ...){
  tval=qt(p=alpha,df=nu)
  lower <- theta+tval*sigma_nu
  upper <- theta-tval*sigma_nu
  cbind(lower,upper)
}

#' @title Assess Equivalence
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param ci                    A \code{vector} with the lower and upper bounds of the confidence interval.
#' @param c                     A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-c,c)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to 0 is there is no equivalence, or 1 if there is equivalence.
BE = function(ci, c, ...){
  return(prod((ci[1]>-c)*(ci[2]<c)))
}


#' @title Assess Equivalence in the univariate framework for the TOST and the alpha-TOST
#' @description The function returns many useful components in the assessment of bioequivalence
#' in the univariate frameowrk for each of the TOST and alphaTOST, namely, the estimated parameter,
#' the estimated standard error, confidence intervals for both methods, the equivalence limit,
#' the original significance level and the corrected one.
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param c                     A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-c,c)
#' @param ...                   Additional parameters.
#' @return The function returns an object of class `matrix`.
#' @export
cTOST = function(alpha, theta, sigma_nu, nu, c, ...){
  alpha_star = get_alpha_star(alpha=alpha, sigma_nu=sigma_nu, nu=nu, c=c)
  original_ci = ci(alpha=alpha, theta=theta, sigma_nu=sigma_nu, nu=nu)
  corrected_ci = ci(alpha=alpha_star, theta=theta, sigma_nu=sigma_nu, nu=nu)
  BE_TOST = BE(ci=original_ci, c=c)
  BE_cTOST = BE(ci=corrected_ci, c=c)
  res = matrix(NA,nrow=2,ncol=7)
  colnames(res) = c("theta","sigma_nu","Level","CI - low.","CI - up.","c","Decision")
  rownames(res) = c("TOST","cTOST")
  res[1,] = c(theta,sigma_nu,alpha,original_ci[1],original_ci[2],c,BE_TOST)
  res[2,] = c(theta,sigma_nu,alpha_star,corrected_ci[1],corrected_ci[2],c,BE_cTOST)
  return(round(res,2))
}



