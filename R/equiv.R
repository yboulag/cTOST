#' @title Power function
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to a probability.
power_TOST = function(alpha, theta, sigma_nu, nu, delta, ...){
  tval = qt(1 - alpha, df = nu)
  mu1 = (theta + delta)/sigma_nu
  mu2 = (theta - delta)/sigma_nu
  R = (delta*sqrt(nu))/(tval*sigma_nu)
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
#' @param delta                 A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to a probability.
size_TOST = function(alpha, sigma_nu, nu, delta, ...){
  power_TOST(alpha = alpha, theta = delta, sigma_nu = sigma_nu,
             nu = nu, delta = delta)
}

#' @title Objective function to optimise
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param test                  A \code{numeric} value specifying the significance level to optimise.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value for the objective function.
obj_fun_alpha_star = function(test, alpha, sigma_nu, nu, delta, ...){
  size = size_TOST(alpha = test, sigma_nu = sigma_nu, nu = nu, delta = delta)
  (size - alpha)^2
}

#' @title Get alpha star
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to the solution of the optimisation.
get_alpha_star = function(alpha, sigma_nu, nu, delta, ...){
  out = optimize(obj_fun_alpha_star, c(alpha, 0.5),
                 alpha = alpha, sigma_nu = sigma_nu,
                 nu = nu, delta = delta)

  # size_out = size_TOST(alpha = out$minimum, sigma_nu = sigma_nu, nu = nu, delta = delta)

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
#' @param delta                   A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @return The function returns a \code{numeric} value that corresponds to 0 is there is no equivalence, or 1 if there is equivalence.
BE = function(ci, delta, ...){
  return(prod((ci[1]>-delta)*(ci[2]<delta)))
}


#' @title Equivalence Assessment in the univariate framework for the TOST and the alpha-TOST
#' @description The function returns many useful components in the assessment of bioequivalence
#' in the univariate frameowrk for each of the TOST and alphaTOST, namely, the estimated parameter,
#' the estimated standard error, confidence intervals for both methods, the equivalence limit,
#' the original significance level and the corrected one.
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to equivalence limit. We assume symmetry, i.e, the equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @return The function returns an object of class `matrix`.
#' @export
#' @examples
#' theta <- diff(apply(skin,2,mean))
#' nu <- nrow(skin)-1
#' sigma_nu <- sd(apply(skin,1,diff))/sqrt(nu)
#' alpha <- 0.05
#' delta <- log(1.25)
#' cTOST(alpha, theta, sigma_nu, nu, delta)
cTOST = function(alpha, theta, sigma_nu, nu, delta, ...){
  alpha_star = get_alpha_star(alpha=alpha, sigma_nu=sigma_nu, nu=nu, delta=delta)
  original_ci = ci(alpha=alpha, theta=theta, sigma_nu=sigma_nu, nu=nu)
  corrected_ci = ci(alpha=alpha_star, theta=theta, sigma_nu=sigma_nu, nu=nu)
  BE_TOST = BE(ci=original_ci, delta=delta)
  BE_cTOST = BE(ci=corrected_ci, delta=delta)
  # res = matrix(NA,nrow=2,ncol=7)
  # colnames(res) = c("theta","sigma_nu","Level","CI - low.","CI - up.","delta","Decision")
  # rownames(res) = c("TOST","cTOST")
  # res[1,] = c(theta,sigma_nu,alpha,original_ci[1],original_ci[2],delta,BE_TOST)
  # res[2,] = c(theta,sigma_nu,alpha_star,corrected_ci[1],corrected_ci[2],delta,BE_cTOST)
  # out = res
  out = list(theta=theta,
             sigma_nu=sigma_nu,
             nu=nu,
             delta=delta,
             alpha=alpha,
             alpha_star=alpha_star,
             TOST_ci=original_ci,
             TOST_decision=BE_TOST,
             aTOST_ci=corrected_ci,
             aTOST_decision=BE_cTOST)
  class(out) = "cTOST"
  invisible(out)
}

#' @title Extract relevant info to format
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
coef.cTOST = function(x, ...){
  method = c("TOST","aTOST")
  decision = c(x$TOST_decision,x$aTOST_decision)
  ci = rbind(x$TOST_ci,x$aTOST_ci)
  df = cbind(ci, rep(x$delta, length(method)), decision)
  rownames(df) = method
  colnames(df) = c("CI - low.", "CI - up.", "Equiv. lim.", "Equiv.")
  df
}

#' @title Print results
#' @author Younes Boulaguiem, Stéphane Guerrier, etc.
print.cTOST = function(x, ...){
  cTOST_coef = coef(x)[,1:3]
  nom = rownames(cTOST_coef)
  var_name = "Procedure"
  BE = TRUE
  m = max(c(nchar(nom), nchar(var_name)))
  n = length(nom)

  space = "   "
  cTOST_coef_to_print = format(cTOST_coef, digits = 2, justify = "right", width = 11)
  col_nom = colnames(cTOST_coef_to_print)
  cat(sprintf(paste("%", m, "s", sep = ""), var_name))

  for (i in 1:(length(col_nom))){
    cat(space)
    cat(sprintf(paste("%", max(nchar(cTOST_coef_to_print[,i])), "s", sep = ""), col_nom[i]))
  }
  cat(space)
  cat(space)
  cat(" Equiv.")
  cat("\n")

  decision = c(x$TOST_decision,x$aTOST_decision)
  decision_label = sapply(decision, function(y) if(y) "Yes" else "No")

  for (i in 1:n){
    cat(sprintf(paste("%", m, "s", sep = ""), nom[i]))
    for (j in 1:(length(col_nom))){
      cat(space)
      cat(cTOST_coef_to_print[i,j])
    }
    cat(space)
    cat(space)

    abs_max_ci = max(abs(cTOST_coef[i,1:2]))

    if (abs_max_ci > cTOST_coef[i,3]){
      cat(" ")
      BE = FALSE
      cli_text(col_red("{symbol$cross}"))
    }else{
      cat(" ")
      cli_text(col_green(" {symbol$tick} "))
    }
  }
  cli_text("{.emph Estimate : }", round(x$theta,3), "; {.emph Std. error : }", round(x$sigma_nu,3),
           "; {.emph Signif. level : }", alpha, "; {.emph Corrected level : }", round(alpha_star,3))
  cli_text("{.emph Signif. codes:} {.red {symbol$cross}} Not Equivalent; {.green {symbol$tick}} Equivalent")
}

print(x)










