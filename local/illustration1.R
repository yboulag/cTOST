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
cols = gg_color_hue(2)

reds = t_col(cols[1], c(0, 20, 40, 85))
blues = t_col(cols[2], c(0, 20, 40, 80))
col1 = "grey"
grey = t_col(col1, 80)
grey2 = t_col(col1, 30)
sig = 0.05
delta = sig*qnorm(1-0.05)
c_val = log(1.25)
alphas = c(0.2, 0.05, 0.01)
K = length(alphas)
eps = 0.0075


tikz("local/tex/fig1_illustration.tex", width = 4, height = 4, standAlone = TRUE,
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
rect(-c_val + delta, -c_val + delta, c_val - delta, c_val - delta, col = "white", border = grey2)

lines(c(-1, 1), c(0, 0))
lines(c(0, 0), c(-1, 1))
lines(c(c_val, c_val), c(-eps, eps))
lines(c(-c_val, -c_val), c(-eps, eps))
lines(c(c_val - delta, c_val - delta), c(-eps, eps))
lines(c(-c_val + delta, -c_val + delta), c(-eps, eps))

lines(c(-eps, eps), c(c_val, c_val))
lines(c(-eps, eps), c(-c_val, -c_val))
lines(c(-eps, eps), c(c_val - delta, c_val - delta))
lines(c(-eps, eps), c(-c_val + delta, -c_val + delta))




points(c_val,0, pch = 16, col = reds[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(c_val,0),
                    sigma = diag(c(sig,sig)^2),
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = reds[i])
}

a = mixtools::ellipse(mu=c(c_val,0),
                      sigma = diag(c(sig,sig)^2),
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = a[,1] < c_val - delta
polygon(c(a[target,1], rev(a[target,1])),
        c(a[target,2], rep(c_val - delta, sum(target))), col = reds[4],
        border = NA)

points(-c_val,0, pch = 16, col = reds[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(-c_val,0),
                    sigma = diag(c(sig,sig)^2),
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = reds[i])
}

a = mixtools::ellipse(mu=c(-c_val,0),
                      sigma = diag(c(sig,sig)^2),
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = a[,1] > -c_val + delta
polygon(c(a[target,1], rev(a[target,1])),
        c(a[target,2], rep(-c_val + delta, sum(target))), col = reds[4],
        border = NA)

points(0, c_val, pch = 16, col = reds[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(0, c_val),
                    sigma = diag(c(sig,sig)^2),
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = reds[i])
}

a = mixtools::ellipse(mu=c(0, c_val),
                      sigma = diag(c(sig,sig)^2),
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = a[,2] < c_val - delta
polygon(c(a[target,1], rev(a[target,1])),
        c(a[target,2], rep(-c_val + delta, sum(target))), col = reds[4],
        border = NA)

points(0, -c_val, pch = 16, col = reds[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(0, -c_val),
                    sigma = diag(c(sig,sig)^2),
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = reds[i])
}

a = mixtools::ellipse(mu=c(0, -c_val),
                      sigma = diag(c(sig,sig)^2),
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = a[,2] > -c_val + delta
polygon(c(a[target,1], rev(a[target,1])),
        c(a[target,2], rep(-c_val + delta, sum(target))), col = reds[4],
        border = NA)

points(c_val, c_val, pch = 16, col = blues[1])
for (i in 1:K){
  mixtools::ellipse(mu=c(c_val, c_val),
                    sigma = diag(c(sig,sig)^2),
                    alpha = alphas[i], npoints = 1000, newplot = FALSE,
                    draw = TRUE, col = blues[i])
}

a = mixtools::ellipse(mu=c(c_val, c_val),
                      sigma = diag(c(sig,sig)^2),
                      alpha = alphas[i], npoints = 1000, newplot = FALSE,
                      draw = FALSE, col = reds[i])
target = (a[,1] < (c_val - delta)) & (a[,2] < (c_val - delta))
polygon(c(a[target,1], rep(c_val - delta, sum(target))),
        c(a[target,2], rep(c_val - delta, sum(target))), col = blues[4],
        border = NA)

text(c_val, -2*eps, "$c$")
text(-c_val, -2*eps, "$-c$")
text(c_val - delta, -2*eps, "$c - \\delta$")
text(-c_val + delta, -2*eps, "$-c + \\delta$")

text(-2.5*eps, c_val, "$c$")
text(-2.5*eps, -c_val,  "$-c$")
text(-4*eps, c_val - delta,  "$c - \\delta$")
text(-5*eps, -c_val + delta, "$-c + \\delta$")


dev.off()
