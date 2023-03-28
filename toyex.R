library(PowerTOST)

power_TOST = function(alpha, theta, sigma_nu, nu, c){
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

size_TOST = function(alpha, sigma_nu, nu, c){
  power_TOST(alpha = alpha, theta = c, sigma_nu = sigma_nu,
             nu = nu, c = c)
}

obj_fun_alpha_star = function(test, alpha, sigma_nu, nu, c){
  size = size_TOST(alpha = test, sigma_nu = sigma_nu, nu = nu, c = c)
  (size - alpha)^2
}

get_alpha_star = function(alpha, sigma_nu, nu, c){
  out = optimize(obj_fun_alpha_star, c(alpha, 0.5),
                 alpha = alpha, sigma_nu = sigma_nu,
                 nu = nu, c = c)

  # size_out = size_TOST(alpha = out$minimum, sigma_nu = sigma_nu, nu = nu, c = c)

  out$minimum
}

ci = function(alpha, theta, sigma_nu, nu){
  tval=qt(p=alpha,df=nu)
  lower <- theta+tval*sigma_nu
  upper <- theta-tval*sigma_nu
  cbind(lower,upper)
}

BE = function(ci, c){
  return(prod((ci[1]>-c)*(ci[2]<c)))
}

cTOST = function(alpha, theta, sigma_nu, nu, c){
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

alpha = 0.05
nu = 10
c = log(1.25)
sigma_nu = 2
theta=0

cTOST(alpha, theta, sigma_nu, nu, c)


