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

library(tikzDevice)
tikz("local/tex/simu_table.tex", width = 9, height = 6, standAlone = TRUE,
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
# Bias

plot(NA, xlim = c(-0.32,4), ylim = c(-4.9,8),
     xlab = "", ylab = "", axes = FALSE)
#abline(h = 0:9)
#abline(v = 0:3)
text(2, 7.8, "\\textbf{Simulation}", cex = 1.2)
lines(c(0.1, 3.9), c(7.55, 7.55), col = "grey40", lwd = 1.5)
text(0.5, 7.25, "\\textbf{1}", cex = 1.15)
text(1.5, 7.25, "\\textbf{2} (Appendix E)", cex = 1.15)
text(2.5, 7.25, "\\textbf{3} (Appendix E)", cex = 1.15)
text(3.5, 7.25, "\\textbf{4} (Appendix E)", cex = 1.15)


col_grid = "grey80"
lines(c(0, 4), c(7, 7), col = col_grid)
lines(c(0, 4), c(6, 6), col = col_grid)
lines(c(0, 4), c(5, 5), col = col_grid)
lines(c(0, 4), c(4, 4), col = col_grid)
lines(c(0, 4), c(3, 3), col = col_grid)
lines(c(0, 4), c(2, 2), col = col_grid)
lines(c(0, 4), c(1, 1), col = col_grid)
lines(c(0, 4), c(0, 0), col = col_grid)
lines(c(0, 4), c(-1, -1), col = col_grid)
lines(c(0, 4), c(-2, -2), col = col_grid)
lines(c(0, 4), c(-3, -3), col = col_grid)
lines(c(0, 4), c(-4, -4), col = col_grid)
lines(c(0, 4), c(-5, -5), col = col_grid)

lines(c(0, 0), c(-5, 7), col = col_grid)
#lines(c(1, 1), c(1, 7), col = col_grid)
#lines(c(2, 2), c(1, 7), col = col_grid)
lines(c(4, 4), c(-5, 7), col = col_grid)

lines(c(0.1, 0.9), c(7, 7), col = "grey40", lwd = 1.5)
lines(c(1.1, 1.9), c(7, 7), col = "grey40", lwd = 1.5)
lines(c(2.1, 2.9), c(7, 7), col = "grey40", lwd = 1.5)
lines(c(3.1, 3.9), c(7, 7), col = "grey40", lwd = 1.5)

deltax = 0.03
deltay = 0.09
col1 = "#2e86c1"
col2 = "#f39c12"

col1b = t_col(col1, 60)
col2b = t_col(col2, 60)

col1 = t_col(col1, 80)
col2 = t_col(col2, 80)

grey_transp = t_col("grey40", 95)

rect(0 + deltax, 6 + deltay, 4-deltax, 7 - deltay, col = col1, border = col1b)

rect(0 + deltax, 5 + deltay, 2-deltax, 6 - deltay, col = col2, border = col2b)
rect(2 + deltax, 5 + deltay, 4-deltax, 6 - deltay, col = col2, border = col2b)

rect(0 + deltax, 4 + deltay, 1-deltax, 5 - deltay, col = col1, border = col1b)
rect(1 + deltax, 4 + deltay, 2-deltax, 5 - deltay, col = col1, border = col1b)
rect(2 + deltax, 4 + deltay, 3-deltax, 5 - deltay, col = col1, border = col1b)
rect(3 + deltax, 4 + deltay, 4-deltax, 5 - deltay, col = col1, border = col1b)

rect(0 + deltax, 3 + deltay, 2-deltax, 4 - deltay, col = col2, border = col2b)
rect(2 + deltax, 3 + deltay, 4-deltax, 4 - deltay, col = col2, border = col2b)

rect(0 + deltax, 2 + deltay, 2-deltax, 3 - deltay, col = col1, border = col1b)
rect(2 + deltax, 2 + deltay, 4-deltax, 3 - deltay, col = col1, border = col1b)

rect(0 + deltax, 1 + deltay, 2-deltax, 2 - deltay, col = col2, border = col2b)
rect(2 + deltax, 1 + deltay, 4-deltax, 2 - deltay, col = col2, border = col2b)

rect(0 + deltax, 0 + deltay, 2-deltax, 1 - deltay, col = col1, border = col1b)
rect(2 + deltax, 0 + deltay, 4-deltax, 1 - deltay, col = col1, border = col1b)

rect(0 + deltax, -1 + deltay, 2-deltax, 0 - deltay, col = col2, border = col2b)
rect(2 + deltax, -1 + deltay, 4-deltax, 0 - deltay, col = col2, border = col2b)

rect(0 + deltax, -2 + deltay, 2-deltax, -1 - deltay, col = col1, border = col1b)
rect(2 + deltax, -2 + deltay, 4-deltax, -1 - deltay, col = col1, border = col1b)

rect(0 + deltax, -3 + deltay, 4-deltax, -2 - deltay, col = col2, border = col2b)

rect(0 + deltax, -4 + deltay, 4-deltax, -3 - deltay, col = col1, border = col1b)

rect(0 + deltax, -5 + deltay, 1-deltax, -4 - deltay, col = col2, border = col2b)
rect(1 + deltax, -5 + deltay, 2-deltax, -4 - deltay, col = col2, border = col2b)
rect(2 + deltax, -5 + deltay, 3-deltax, -4 - deltay, col = col2, border = col2b)
rect(3 + deltax, -5 + deltay, 4-deltax, -4 - deltay, col = col2, border = col2b)

text(0, 6.5, "$c$", cex = 1.25, pos = 2)
text(0, 5.5, "$m$", cex = 1.25, pos = 2)
text(0, 4.5, "$\\nu$", cex = 1.25, pos = 2)
text(0, 3.5, "diag$(\\Sigma_1)$", cex = 1.25, pos = 2)
text(0, 2.5, "diag$(\\Sigma_2)$", cex = 1.25, pos = 2)
text(0, 1.5, "diag$(\\Sigma_3)$", cex = 1.25, pos = 2)
text(0, 0.5, "diag$(\\Sigma_4)$", cex = 1.25, pos = 2)
text(0, -0.5, "diag$(\\Sigma_5)$", cex = 1.25, pos = 2)
#text(-0.1, 2.5, "$\\text{diag}(\\boldsymbold{\\Sigma}_2)$", cex = 1.25)
#text(-0.1, 1.5, "$\\text{diag}(\\boldsymbold{\\Sigma}_3)$", cex = 1.25)
#text(-0.1, 0.5, "$\\text{diag}(\\boldsymbold{\\Sigma}_4)$", cex = 1.25)
#text(-0.1, -0.5, "$\\text{diag}(\\boldsymbold{\\Sigma}_5)$", cex = 1.25)
text(0, -1.5, "Covariance", cex = 1.25, pos = 2)
text(-0.1, -2.5, "$\\alpha$", cex = 1.25)
text(-0.1, -3.5, "$B$", cex = 1.25)
text(0, -4.5, "Results in", cex = 1.25, pos = 2)

text(2, 6.5, "$\\log(1.25) \\approx 0.2231$")
text(1, 5.5, "$2$")
text(3, 5.5, "$4$")
text(0.5, 4.5, "$20$")
text(1.5, 4.5, "$40$")
text(2.5, 4.5, "$20$")
text(3.5, 4.5, "$40$")
text(1, 3.5, "$[0.05^2,\\, 0.05^2]$")
text(3, 3.5, "$[0.05^2,\\, 0.05^2,\\, 0.05^2,\\, 0.05^2]$")

text(1, 2.5, "$[0.1^2,\\, 0.1^2]$")
text(3, 2.5, "$[0.1^2,\\, 0.1^2,\\, 0.1^2,\\, 0.1^2]$")

text(1, 1.5, "$[0.15^2,\\, 0.15^2]$")
text(3, 1.5, "$[0.15^2,\\, 0.15^2,\\, 0.15^2,\\, 0.15^2]$")

text(1, 0.5, "$[0.05^2,\\, 0.1^2]$")
text(3, 0.5, "$[0.05^2,\\, 0.05^2,\\, 0.1^2,\\, 0.1^2]$")

text(1, -0.5, "$[0.05^2,\\, 0.15^2]$")
text(3, -0.5, "$[0.05^2,\\, 0.05^2,\\, 0.15^2,\\, 0.15^2]$")

text(1, -1.5, "$\\sigma_{i,j} = \\rho \\sigma_i \\sigma_j$, for $\\rho \\in \\{0, \\, 0.5, \\, 0.9\\}$")
text(3, -1.5, "$\\sigma_{i,j} = \\rho^{|i-j|} \\sigma_i \\sigma_j$, for $\\rho \\in \\{0, \\, 0.5, \\, 0.9\\}$")

text(0.5, -4.5, "Figure 4")
text(1.5, -4.5, "Figure A.1")
text(2.5, -4.5, "Figure A.2")
text(3.5, -4.5, "Figure A.3")

text(2, -2.5, "0.05")

text(2, -3.5, "$5 \\times 10^4$")




#text(3.5, 5.5, "$45$", col = grey(0.4))
#text(0.5, 4.5, "$0.08,\\, 0.12,\\, 0.16$")
#text(3.5, 4.5, "$0.08,\\, 0.12,\\, 0.16$", col = grey(0.4))
#text(2, 4.5, "$100$ evenly spaced values between $0.01$ and $0.3$")

#text(0.5, 3.7, "$30$ evenly spaced values")
#text(0.5, 3.3, "between $0$ and $0.26$")
#text(1.5, 3.5, "$c$")
#text(2.5, 3.5, "$0$")

#text(3.5, 3.7, "$30$ evenly spaced values", col = grey(0.4))
#text(3.5, 3.3, "between $0$ and $0.26$", col = grey(0.4))

#text(2, 2.5, "$0.05$")
#text(2, 1.5, "$10^5$")
#text(1.5, 0.5, "General (canonical form)")
#text(3.5, 0.5, "Paired", col = grey(0.4))
#text(1.5, -0.5, "TOST, $\\alpha$-TOST and $\\delta$-TOST")
#text(3.5, -0.3, "TOST, $\\alpha$-TOST, $\\delta$-TOST,", col = grey(0.4))
#text(3.5, -0.7, "SABE and cSABE", col = grey(0.4))

dev.off()
