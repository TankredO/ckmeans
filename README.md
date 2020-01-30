# ckmeans
A consensus K-Means clustering R package.

## Installation
```R
# install.packages('devtools')
library(devtools)

# On UNIX you might to run the following line if you are not a sudo user:
# Sys.setenv(TAR = "/bin/tar")
install_github('https://github.com/TankredO/ckmeans')
```

## Example
```R
library(ckmeans)

ckms <- multickmeans(iris[,1:4], ks = 1:10, p_pred = 0.8, p_samp = 0.8, n_rep = 500)

plotDist(ckms$ckms[[3]])
plot(ckms$bics, type='l')
```
