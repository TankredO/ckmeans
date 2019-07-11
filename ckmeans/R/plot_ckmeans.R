#' Generic plotting function for distances (distance matrices)
#' @title Plot distance matrix
#' @author Tankred Ott
#' @param x object to plot
#' @param ... further arguments passed to the class specific plotDist functions
#' @import graphics
#' @import grDevices
#' @export
plotDist <- function(x, ...) UseMethod('plotDist', x)

#' Function to plot a distance matrix
#' @title Plot distance matrix with clusters
#' @description Plots a distance matrix with color coded clusters
#' @author Tankred Ott
#' @param x n*n matrix
#' @param cl cluster membership
#' @param col vector of colors (optional)
#' @param col_cl vector of colors for the clusters (optional)
#' @param ord vectors of indices for ordering the matrix (optional).
#' @param returnOrd logical determining whether the indices used to order the plot should be returned (useful for plotting multiple distance plots with the same order)
#' @param main plot title
#' @param ... further arguments passed to the class specific plotDist functions
#' @export
plotDist.matrix <- function(x, cl = NULL , col = NULL, col_cl = NULL, ord = NULL, returnOrd = FALSE, main = "", ...) {
  if (is.null(ord)) ord <- seriation::seriate(as.dist(x), method = 'GW')[[1]]$order
  x_ord <- x[ord, ord]
  if (!is.null(cl)) cl_ord <- cl[ord]

  renameClusters <- function(cl) {
    1:length(unique(cl))
    tmp <- cl

    cls <- 1:length(unique(cl))
    i <- 1
    for (c in unique(cl)) {
      tmp[cl == c] <- cls[i]
      i <- i + 1
    }
    tmp
  }

  if (!is.null(cl)) cl_ord <- renameClusters(cl_ord)

  lerp <- function(a, b, .t) (1-.t) %*% t(a) + (.t) %*% t(b)

  if (is.null(col)) {
    .t <- seq(0, 1, length.out = 256)
    rgbs <- lerp(c(1,1,1), c(18/255, 136/255, 240/255),.t)
    col <- rev(rgb(rgbs[,1], rgbs[,2], rgbs[,3]))
  }

  getCol <- function(val) {
    col[ceiling(val * (length(col) - 1)) + 1]
  }

  nR <- nrow(x)
  nC <- ncol(x)

  w <- 1 / nC
  h <- 1 / nR
  m <- w / 20

  op <- par(xpd=TRUE, mar=c(1, 6, 6, 1), oma=c(0,0,0,0))
  on.exit(par(op))

  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = NA, ylab = NA, main = main)

  xs <- seq(0, 1, length.out = nC + 1)
  ys <- seq(0, 1, length.out = nR + 1)

  x12 <- cbind(xs[1:(length(xs)-1)], xs[2:length(xs)])
  y12 <- cbind(ys[1:(length(ys)-1)], ys[2:length(ys)])

  x_ax <- apply(x12, 1, mean)
  y_ax <- apply(y12, 1, mean)
  axis(3, x_ax, labels = colnames(x_ord), las = 2, lwd = 0.25, cex.axis = 0.75)
  axis(2, y_ax, labels = rev(colnames(x_ord)), las = 2, lwd = 0.25, cex.axis = 0.75)

  if (is.null(col_cl)) col_cl <- RColorBrewer::brewer.pal(7, 'Set1')

  for (.x in 1:nC) {

    for (.y in 1:nR) {
      .col <- getCol(x_ord[.y, .x])
      rect(xleft = xs[.y], xright = xs[.y] + w, ytop = 1 - ys[.x], ybottom = 1 - ys[.x] - h,
           border = NA, col = .col)
    }
  }

  for (.x in 1:nC) {
    lines(x = rep(xs[.x], 2), y = c(0,1), col = 'white', lwd = m)
    if (!is.null(cl)) {
      .col_cl <- col_cl[cl_ord[.x]]
      rect(xleft = xs[.x], xright = xs[.x] + w, ybottom = 1.01, ytop = 1.02, col = .col_cl, border = NA)
    }
  }
  for (.y in 1:nR) {
    lines(y = rep(ys[.y], 2), x = c(0,1), col = 'white', lwd = m)
    if (!is.null(cl)) {
      .col_cl <- col_cl[cl_ord[.y]]
      rect(ybottom = 1 - ys[.y], ytop = 1 - ys[.y] - h, xright = -0.01, xleft = -0.02, col = .col_cl, border = NA)
    }
  }

  if (returnOrd == TRUE) ord
}

#' Distance (and consensus cluster) plot for cKmeans
#' @title Plot ckmeans object as distance matrix
#' @description Plots the consensus matrix as distance matrix
#' @param x cKmeans object
#' @param col vector of colors (optional)
#' @param ord vectors of indices for ordering the matrix (optional).
#' @param returnOrd logical, determining whether the indices used to order the plot should be returned (useful for plotting multiple distance plots with the same order)
#' @param plotCC logical, determining whether the consensus cluster should be plotted
#' @param col_cl vector of colors for the clusters (optional)
#' @param main title of the plot
#' @param ... further arguments passed to the class specific plotDist functions
#' @export
plotDist.ckmeans <- function(x, col = NULL, ord = NULL, returnOrd = FALSE, plotCC = TRUE, col_cl = NULL, main = NULL, ...) {
  # standard plotDist, get order
  .ord <- plotDist(x = 1-x$pcc, cl = x$cc, col = col, ord = ord, returnOrd = returnOrd, col_cl = col_cl, main = main)
  #
  # if (plotCC == TRUE) {
  #   k <- x$k
  #   cc <- x$cc[.ord] # get (ordered) consensus clusters
  #
  #   if (is.null(colCC)) colCC <- .cols[5:(5+k)]
  #
  #   # prepare plotting
  #   nR <- nC <- ncol(x$pcc)
  #   h <- 1 / nR
  #
  #   xs <- seq(0, 1, length.out = nC + 1)
  #   x12 <- cbind(xs[1:(length(xs)-1)], xs[2:length(xs)])
  #
  #   # plot rectangles representing consensus cluster membership
  #   for (i in 1:nrow(x12)) {
  #     rect(xleft = x12[i,1], xright = x12[i,2], ybottom = 1.015, ytop = 1.015 + 2 * h, lty='blank',
  #          col = colCC[cc[i]])
  #   }
  # }
  #
  if (returnOrd == TRUE) .ord
}

.cols <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
           "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
           "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
           "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
           "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
           "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
           "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
           "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

           "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
           "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
           "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
           "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
           "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
           "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
           "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
           "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58")
