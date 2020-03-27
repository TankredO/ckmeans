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
#' @param x n*n (distance) matrix
#' @param cl vector determining cluster membership
#' @param value_range vector with two elements c(d_min, d_max) determining the possible value range within the matrix. By default this will be the range of values in x.
#' @param ord vectors of indices for ordering of the matrix or boolean determining whether the matrix should be ordered
#' @param col vector of colors for the distance matrix
#' @param col_cl vector of colors or color ramp function for the clusters
#' @param plot_colorbar logical determining whether a color bar should be plotted
#' @param ... further arguments
#' @export
plotDist.matrix <- function(x, cl=NULL, value_range=NULL, ord=TRUE, col=NULL, col_cl=NULL, plot_colorbar=FALSE,...) {
  # unpack ...
  dot <- list(...)

  # prepare par
  cb_mar_r <- if(is.null(dot$cb_mar_r)) 5 else dot$cb_mar_r
  mar <- if(is.null(dot$mar)) c(0,3,3,ifelse(plot_colorbar, cb_mar_r, 0)) else dot$mar

  old_par <- par(mar=mar, xpd=TRUE)
  on.exit(par(old_par))


  # prepare x
  n <- nrow(x)

  # if x has no row names set row indices as row names
  if (is.null(row.names(x))) row.names(x) <- 1:n

  # order
  if (length(ord) > 1) {
    x <- x[ord, ord]
  } else if (ord) {
    ord <- seriation::seriate(as.dist(x), method = 'GW')[[1]]$order
    x <- x[ord, ord]
  } else ord <- 1:n
  if (!is.null(cl)) cl <- cl[ord]
  names(ord) <- row.names(x)

  # color Ramp for distance matrix
  col <- if(is.null(col)) c('#24526E', 'white') else col
  if (length(col) < 2) stop('Passed a single color as col argument but at least two colors are required!')
  cRamp <- colorRamp(col)


  # create empty plot
  plot(NULL, xlim = c(0, n), ylim = c(0, n), frame.plot = FALSE, axes = FALSE, xlab = '', ylab = '')

  # figure
  # if value range is NULL calculate it from the input matrix
  if (is.null(value_range)) value_range <- c(min(x), max(x))
  x_scaled  <- if (value_range[1] == value_range[2]) x
               else (x - value_range[1]) / (value_range[2]  - value_range[1])
  dist_raster <- as.raster(
    matrix(
      # apply(cRamp(x_scaled), 1, function(col) rgb(col[1], col[2], col[3], maxColorValue = 255)),
      rgb(cRamp(x_scaled), maxColorValue = 255),
      ncol = n, nrow = n
    )
  )
  rasterImage(dist_raster, 0, 0, n, n, interpolate = FALSE)


  # plot cluster membership
  if(!is.null(cl)) {
    col_cl <- if(is.null(col_cl)) pals::kelly(22)[5:22] else col_cl
    if (is.function(col_cl)) {
      col_cl <- col_cl(length(unique(cl)))
    }

    if (length(unique(col_cl)) < length(unique(cl))) stop('Received less distinct colors than unique values in cl.')

    names(col_cl) <- unique(cl)

    cl_raster <- as.raster(col_cl[as.character(cl)])
    rasterImage(cl_raster, -1.5, 0, -0.5, n, interpolate = FALSE)
    rasterImage(t(cl_raster), 0, n+0.5, n, n+1.5, interpolate = FALSE)
  }


  # plot lines to create "boxes" instead of a simple raster image
  abline(v = 1:(n-1), col='white', lwd=1)
  abline(h = 1:(n-1), col='white', lwd=1)


  # plot color bar
  if(plot_colorbar) {
    cbar_raster <- as.raster(rev(rgb(t(sapply(seq(0, 1, length.out = 50), cRamp)), maxColorValue = 255)))

    cb_round_d <- if(is.null(dot$cb_round_d)) 3 else dot$cb_round_d

    cb_x1 <- n + 2.5
    cb_x2 <- n + 4
    cb_tick_at <- 0:10
    cb_n_ticks <- length(cb_tick_at)
    cb_tick_labels <- round(value_range[1] + cb_tick_at / (cb_n_ticks - 1) * value_range[2], cb_round_d)

    rasterImage(cbar_raster, cb_x1, 0, cb_x2, n, interpolate = TRUE)

    sapply(
      cb_tick_at,
      function(y) lines(c(cb_x1, cb_x2 + 0.5), rep(y * n / (cb_n_ticks - 1), 2), col='black')
    )
    text(cb_x2 + 1.0, cb_tick_at * n / (cb_n_ticks - 1), cb_tick_labels, adj=c(0, 0.55))
  }

  # axis
  at <- 1:n - 0.5
  labels <- rownames(x)

  axis(2, outer = F, at = rev(at), labels = labels, cex.axis=0.5, pos=(ifelse(is.null(cl), 0, -1.5)), las=1, lwd=0, lwd.ticks = 1)
  axis(3, outer = F, at = at, labels = labels, cex.axis=0.5, pos=ifelse(is.null(cl), n, n + 1.5), las=2, lwd=0, lwd.ticks = 1)

  # cl lines
  ## TODO: get start and end indices of ranges of the same value withing cl; plot v and h lines

  return(ord)
}

#' Distance (and consensus cluster) plot for cKmeans
#' @title Plot ckmeans object as distance matrix
#' @description Plots the consensus matrix as distance matrix
#' @param x cKmeans object
#' @param col vector of colors (optional)
#' @param ord vectors of indices for ordering the matrix (optional).
#' @param col_cl vector of colors for the clusters (optional)
#' @param plot_colorbar logical determining whether a color bar should be plotted
#' @param ... further arguments passed to the class specific plotDist functions
#' @export
plotDist.ckmeans <- function(x, col = NULL, ord = TRUE, col_cl = NULL, plot_colorbar=TRUE, ...) {
  # standard plotDist, get order
  .ord <- plotDist(x = 1-x$pcc, cl = x$cc, col = col, ord = ord, col_cl = col_cl, plot_colorbar=plot_colorbar, value_range=c(0,1))

  return(.ord)
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
