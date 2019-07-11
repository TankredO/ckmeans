## generate data
x1 = c(rnorm(10), rnorm(14, 5, 1) + 2, rnorm(30, -5, 1))
x2 = c(rnorm(10), rnorm(14, 5, 1) + 2, rnorm(30, -5, 1))
x3 = c(rnorm(10)-3, rnorm(14, 2, 1) + 2, rnorm(30, 1, 1))
x = matrix(c(x1, x2, x3), ncol = 3, dimnames = list(1:54, c('x1', 'x2', 'x3')))

pairs(x)

## run ckmeans for a single K
ckm = cKmeans(x, 3, n_rep = 100, p_samp = 0.5, p_pred = 0.5)

# plot consensus matrix with color coded clusters
plotDist(ckm)

# plot(x, col=ckm$cc, pch=c(rep(1, 10), rep(2, 14)))


## run ckmeans for multiple K
ckms = multickmeans(x, 1:7, n_rep = 100, p_samp = 0.8, p_pred = 0.5)
plot(ckms$bics, type='l')

ord = 1:54
for (i in 1:length(ckms$ckms)) {
  plotDist(ckms$ckms[[i]], ord=ord)
}
