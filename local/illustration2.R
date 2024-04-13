build_Sigma = function(sigmas,rho,type=NULL){
  if (is.null(type)){
    R = diag(length(sigmas))
    R[R==0] = rho
    outer(sigmas,sigmas) * R
  } else if (type=="ar1") {
    ar1_cov(sigmas, rho)
  } else {
    error("Sigma type not available")
  }
}


find_sup_xtost = function(Sigma,delta,c_of_0){
  p=ncol(Sigma)
  if(p<=2){
    method="Brent"
  }else{
    method="Nelder-Mead"
    # these are the only possible values for lower and upper
    lower=-Inf;upper=Inf
  }
  # dimss2free yields all the combinations for the dimensions that we would let free.
  dims2free_ = combn(1:p, p-1)

  ## The solution stays sometimes too long at the same place unnecessarily, this is because during optimization
  # the difference between numerical fluctuation of the objective function at the same point
  # are not smaller than tol. --> solution= reltol.

  # we optimize for each combination of free dimensions. tmp_ outputs the parameter values
  # of the argsup's free dimension, and the value of the objective function (negative size)
  tmp_=t(apply(dims2free_,2, function(dims2free){
    if(p<=2){
      lower=-delta[dims2free];upper=delta[dims2free]
    }
    tmp = optim(par=rep(0,p-1),fn=objfun4sup_xtost,
                dims2free=dims2free,Sigma=Sigma,delta=delta, c_of_0=c_of_0,
                method=method, lower=lower,upper=upper,
                control=list(reltol=.Machine$double.eps^0.5))
    c(unlist(tmp[1]),unlist(tmp[2])) # the argsups fluctuate a bit but the objective function value is the very close each time
  }))
  # we look from among all the combinations, the optimal yielding the min negative size
  # tmp_[,ncol(tmp_)]=round(tmp_[,ncol(tmp_)],4)
  ind_max_size=which.min(tmp_[,ncol(tmp_)])
  dims2free=dims2free_[,ind_max_size]
  thetas=delta
  thetas[dims2free]=as.vector(tmp_[ind_max_size,-ncol(tmp_)])
  return(thetas)
}

objfun4sup_xtost = function(theta,dims2free,Sigma,delta,c_of_0){
  p=ncol(Sigma)
  # we define the potential argsup as a vector where we fix one dimension at
  # its corresponding delta, and let the other dimensions be fed by the values
  # considered by the optimisation
  thetas=delta # this is only a useful trick to get directly the dimension to fix at its corresponding delta
  thetas[dims2free]=theta # then we replace the dimensions we let free (to optimize) by the values considered by the optimization
  # cat(thetas,"\n")
  10^3*(-power_xTOST_mv(theta = thetas, Sigma = Sigma,
                        delta = c_of_0)[1])
}

power_xTOST_mv = function(theta, delta, Sigma, alpha=1/2, seed=10^5){
  #delta is what we optimise over
  set.seed(seed)
  mu = rep(0, ncol(Sigma))
  Sig_diag = sqrt(diag(Sigma))
  lb = -delta/Sig_diag-theta/Sig_diag + qnorm(1-alpha)
  ub = delta/Sig_diag-theta/Sig_diag - qnorm(1-alpha)
  if (sum(lb>ub)>0){
    return(0)
  } else {
    tmp_cdf_mvtnorm = mvtnorm::pmvnorm(lower=lb,
                                       upper=ub,
                                       mean=mu,corr=cov2cor(Sigma))[1]
    # tmp_cdf_TN = TruncatedNormal::pmvnorm(mu=mu,sigma=cov2cor(Sigma),lb=lb,
    #                                       ub=ub)[1]
    return(tmp_cdf_mvtnorm)
    # return(list(mvtnorm=tmp_cdf_mvtnorm, TruncatedN=tmp_cdf_TN) )
  }
}

###########################################

sig1 = 0.03
sig2 = 0.06
# Sigma = matrix(c(sig1^2, sig1*sig2*0.5, sig1*sig2*0.5, sig2^2), 2, 2)


sigmas = c(sig1, sig2)
rho = 0.8
Sigma = build_Sigma(sigmas,rho)

cte = rep(log(1.25), 2)

theta_sup_TOST_known_sigma = find_sup_xtost(Sigma,cte,cte)
lambda = theta_sup_TOST_known_sigma

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

  ## Get RGB values for named color
  rgb.val <- col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

  ## Save the color
  invisible(t.col)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(3)

reds = t_col(cols[1], c(0, 20, 40, 85))
greens = t_col(cols[2], c(0, 20, 40, 85))
blues = t_col(cols[3], c(0, 20, 40, 80))
col1 = "grey"
grey = t_col(col1, 80)
grey2 = t_col(col1, 30)
sig1 = 0.04
sig2 = 0.05
delta1 = sig1*qnorm(1-0.05)
delta2 = sig2*qnorm(1-0.05)
c_val = log(1.25)
alphas = c(0.2, 0.05, 0.01)
K = length(alphas)
eps = 0.0075

tikz("local/tex/fig2_illustration.tex", width = 4, height = 4, standAlone = TRUE,
     packages = c("\\usepackage{tikz}",
                  "\\usepackage[active,tightpage,psfixbb]{preview}",
                  "\\PreviewEnvironment{pgfpicture}",
                  "\\setlength\\PreviewBorder{0pt}",
                  "\\usepackage{amssymb}",
                  "\\usepackage{bm}","\\usepackage{amsthm}","\\usepackage{amsbsy}"
                  ,"\\usepackage{amsbsy}"
                  ,"\\usepackage{amsbsy}"
                  ,"\\usepackage{amsfonts}"))
par(mai = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))

plot(NA, xlim = c(-0.25, 0.25), ylim = c(-0.25, 0.25),
     axes = FALSE, xlab = "", ylab = "")

rect(-c_val, -c_val, c_val, c_val, col = grey, border = grey2)
rect(-c_val + delta1, -c_val + delta2, c_val - delta1, c_val - delta2, col = "white", border = grey2)


lines(c(0, 0), c(-1, 1))
lines(c(-1, 1), c(0, 0))
lines(c(c_val, c_val), c(-eps, eps))
lines(c(-c_val, -c_val), c(-eps, eps))
lines(c(c_val - delta1, c_val - delta1), c(-eps, eps))
lines(c(-c_val + delta1, -c_val + delta1), c(-eps, eps))

lines(c(-eps, eps), c(c_val, c_val))
lines(c(-eps, eps), c(-c_val, -c_val))
lines(c(-eps, eps), c(c_val - delta2, c_val - delta2))
lines(c(-eps, eps), c(-c_val + delta2, -c_val + delta2))



Sigma = matrix(c(sig1^2, sig1*sig2*0.5, sig1*sig2*0.5, sig2^2), 2, 2)

points(-c_val, 0.15, pch = 16, col = greens[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(-c_val, 0.15),
                    sigma = Sigma,
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = greens[i])
}

a = mixtools::ellipse(mu=c(-c_val, 0.15),
                      sigma = Sigma,
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = (a[,1] > (-c_val + delta1)) & (a[,2] < (c_val - delta2))
polygon(c(a[target,1], rep(-c_val + delta1, sum(target))),
        c(a[target,2], rep(c_val - delta2, sum(target))), col = greens[4],
        border = NA)

points(c_val, c_val, pch = 16, col = blues[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(c_val, c_val),
                    sigma = Sigma,
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = blues[i])
}

a = mixtools::ellipse(mu=c(c_val, c_val),
                      sigma = Sigma,
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
a = a[c(500:1000, 1:499),]
target = (a[,1] < (c_val - delta1)) & (a[,2] < (c_val - delta2))
polygon(c(a[target,1], rep(c_val - delta1, sum(target))),
        c(a[target,2], rep(c_val - delta2, sum(target))), col = blues[4],
        border = NA)

points(lambda[1], lambda[2], pch = 16, col = reds[1])
for (i in 1:K){
  mixtools::ellipse(mu=lambda,
                    sigma = Sigma,
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = reds[i])
}

a = mixtools::ellipse(mu=lambda,
                      sigma = Sigma,
                      alpha = alphas[K], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = (a[,2] < (c_val - delta2))
polygon(c(a[target,1], rev(a[target,1])),
        c(a[target,2], rep(-c_val + delta2, sum(target))), col = reds[4],
        border = NA)

points(-lambda[1], -lambda[2], pch = 16, col = reds[1])
for (i in 1:K){
  mixtools::ellipse(mu=-lambda,
                    sigma = Sigma,
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = reds[i])
}

a = mixtools::ellipse(mu=-lambda,
                      sigma = Sigma,
                      alpha = alphas[K], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = (a[,2] > (-c_val + delta2))
polygon(c(a[target,1], rep(c_val - delta1, sum(target))),
        c(a[target,2], rev(a[target,2])), col = reds[4],
        border = NA)

text(c_val, -3*eps, "$c$")
text(-c_val, -3*eps, "$-c$")
text(c_val - delta, -3*eps, "$c - \\delta_1$")
text(-c_val + delta, -3*eps, "$-c + \\delta_1$")

text(-2.5*eps, c_val, "$c$")
text(-2.5*eps, -c_val,  "$-c$")
text(-5*eps, c_val - delta,  "$c - \\delta_2$")
text(-6*eps, -c_val + delta, "$-c + \\delta_2$")
dev.off()

