require(gplots)
require(RColorBrewer)
require(cba)



hcopt <- function(d, HC=NULL, method = "ward", members = NULL)
{
   if ( is.null(HC) ) {
      HC <- hclust(d,method=method,members=members)
   }
   ORD <- order.optimal(d,merge=HC$merge)
   HC$merge <- ORD$merge
   HC$order <- ORD$order
   HC
}

scale02<-function(x,cutoff=4){
   rm <- rowMeans(x, na.rm = F)
   x <- sweep(x, 1, rm)
   sx <- apply(x, 1, sd, na.rm = F)
   x <- sweep(x, 1, sx, "/")
   x[x<(-cutoff)]<--cutoff
   x[x>cutoff]<-cutoff
   
   return(x)
}

heatmap.2g<-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                      distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                         "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
                                                                                                                                       w), symm = FALSE, scale = c("none", "row", "column"), 
                      na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                      symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", 
                      col = "heat.colors", colsep, rowsep, sepcolor = "white", 
                      sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
                      na.color = par("bg"), trace = c("column", "row", "both", 
                                                      "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
                      linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
                      cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                      labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
                                                                              NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
                      key = TRUE, keysize = 1.5, density.info = c("histogram", 
                                                                  "density", "none"), denscol = tracecol, symkey = min(x < 
                                                                                                                          0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL, 
                      key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL, 
                      key.par = list(), main = NULL, xlab = NULL, ylab = NULL, 
                      lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, ...) 
{
   scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
   }
   retval <- list()
   scale <- if (symm && missing(scale)) 
      "none"
   else match.arg(scale)
   dendrogram <- match.arg(dendrogram)
   trace <- match.arg(trace)
   density.info <- match.arg(density.info)
   if (length(col) == 1 && is.character(col)) 
      col <- get(col, mode = "function")
   if (!missing(breaks) && (scale != "none")) 
      warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
   if (is.null(Rowv) || is.na(Rowv)) 
      Rowv <- FALSE
   if (is.null(Colv) || is.na(Colv)) 
      Colv <- FALSE
   else if (Colv == "Rowv" && !isTRUE(Rowv)) 
      Colv <- FALSE
   if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
      stop("`x' must be a numeric matrix")
   nr <- di[1]
   nc <- di[2]
   if (nr <= 1 || nc <= 1) 
      stop("`x' must have at least 2 rows and 2 columns")
   if (!is.numeric(margins) || length(margins) != 2) 
      stop("`margins' must be a numeric vector of length 2")
   if (missing(cellnote)) 
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
   if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                      c("both", "row"))) {
         if (is.logical(Colv) && (Colv)) 
            dendrogram <- "column"
         else dendrogram <- "none"
         warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                 dendrogram, "'. Omitting row dendogram.")
      }
   }
   if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                      c("both", "column"))) {
         if (is.logical(Rowv) && (Rowv)) 
            dendrogram <- "row"
         else dendrogram <- "none"
         warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                 dendrogram, "'. Omitting column dendogram.")
      }
   }
   if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
      if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
                                        nr)) 
         stop("Rowv dendrogram doesn't match size of x")
   }
   else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorderfun(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd)) 
         stop("row dendrogram ordering gave index of wrong length")
   }
   else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorderfun(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd)) 
         stop("row dendrogram ordering gave index of wrong length")
   }
   else {
      rowInd <- nr:1
   }
   if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
      if (length(colInd) > nc || any(colInd < 1 | colInd > 
                                        nc)) 
         stop("Colv dendrogram doesn't match size of x")
   }
   else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
         stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
         ddc <- ddr
         colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
   }
   else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm) 
         x
         else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorderfun(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd)) 
         stop("column dendrogram ordering gave index of wrong length")
   }
   else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm) 
         x
         else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorderfun(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd)) 
         stop("column dendrogram ordering gave index of wrong length")
   }
   else {
      colInd <- 1:nc
   }
   retval$rowInd <- rowInd
   retval$colInd <- colInd
   retval$call <- match.call()
   x <- x[rowInd, colInd]
   x.unscaled <- x
   cellnote <- cellnote[rowInd, colInd]
   if (is.null(labRow)) 
      labRow <- if (is.null(rownames(x))) 
         (1:nr)[rowInd]
   else rownames(x)
   else labRow <- labRow[rowInd]
   if (is.null(labCol)) 
      labCol <- if (is.null(colnames(x))) 
         (1:nc)[colInd]
   else colnames(x)
   else labCol <- labCol[colInd]
   if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
   }
   else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
   }
   if (missing(breaks) || is.null(breaks) || length(breaks) < 
          1) {
      if (missing(col) || is.function(col)) 
         breaks <- 16
      else breaks <- length(col) + 1
   }
   if (length(breaks) == 1) {
      if (!symbreaks) 
         breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                       length = breaks)
      else {
         extreme <- max(abs(x), na.rm = TRUE)
         breaks <- seq(-extreme, extreme, length = breaks)
      }
   }
   nbr <- length(breaks)
   ncol <- length(breaks) - 1
   if (class(col) == "function") 
      col <- col(ncol)
   min.breaks <- min(breaks)
   max.breaks <- max(breaks)
   x[x < min.breaks] <- min.breaks
   x[x > max.breaks] <- max.breaks
   if (missing(lhei) || is.null(lhei)) 
      lhei <- c(keysize, 4)
   if (missing(lwid) || is.null(lwid)) 
      lwid <- c(keysize, 4)
   if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      if (!missing(ColSideColors)) {
         if (!is.character(ColSideColors) || ncol(ColSideColors) != 
                nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
         lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                          1)
         lhei <- c(lhei[1], 0.2*nrow(ColSideColors), lhei[2])
      }
      if (!missing(RowSideColors)) {
         if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
         lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                               1), 1), lmat[, 2] + 1)
         lwid <- c(lwid[1], 0.2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
   }
   if (length(lhei) != nrow(lmat)) 
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
   if (length(lwid) != ncol(lmat)) 
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
   op <- par(no.readonly = TRUE)
   on.exit(par(op))
   layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
   if (!missing(RowSideColors)) {
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
   }
   if (!missing(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2]))
      image(matrix(1:length(ColSideColors),ncol=2),
            col = as.vector(t(ColSideColors[,colInd])), 
            axes = FALSE)
   }
   par(mar = c(margins[1], 0, 0, margins[2]))
   x <- t(x)
   cellnote <- t(cellnote)
   if (revC) {
      iy <- nr:1
      if (exists("ddr")) 
         ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
   }
   else iy <- 1:nr
   image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
            c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
         breaks = breaks, ...)
   retval$carpet <- x
   if (exists("ddr")) 
      retval$rowDendrogram <- ddr
   if (exists("ddc")) 
      retval$colDendrogram <- ddc
   retval$breaks <- breaks
   retval$col <- col
   if (any(is.na(x))) {
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
   }
   if (is.null(srtCol)) 
      axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
              offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
           padj = adjCol[2])
   else {
      if (is.numeric(srtCol)) {
         if (missing(adjCol) || is.null(adjCol)) 
            adjCol = c(1, NA)
         xpd.orig <- par("xpd")
         par(xpd = NA)
         xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                      tick = 0)
         text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                 strheight("M"), labels = labCol, adj = adjCol, 
              cex = cexCol, srt = srtCol)
         par(xpd = xpd.orig)
      }
      else warning("Invalid value for srtCol ignored.")
   }
   if (is.null(srtRow)) {
      axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
           tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
   }
   else {
      if (is.numeric(srtRow)) {
         xpd.orig <- par("xpd")
         par(xpd = NA)
         ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                      line = -0.5, tick = 0)
         text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
              y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
              srt = srtRow)
         par(xpd = xpd.orig)
      }
      else warning("Invalid value for srtRow ignored.")
   }
   if (!is.null(xlab)) 
      mtext(xlab, side = 1, line = margins[1] - 1.25)
   if (!is.null(ylab)) 
      mtext(ylab, side = 4, line = margins[2] - 1.25)
   if (!missing(add.expr)) 
      eval(substitute(add.expr))
   if (!missing(colsep)) 
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                                xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                   1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
   if (!missing(rowsep)) 
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                         1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                           1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                                col = sepcolor, border = sepcolor)
   min.scale <- min(breaks)
   max.scale <- max(breaks)
   x.scaled <- scale01(t(x), min.scale, max.scale)
   if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
         if (!is.null(vline)) {
            abline(v = i - 0.5 + vline.vals, col = linecol, 
                   lty = 2)
         }
         xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
         xv <- c(xv[1], xv)
         yv <- 1:length(xv) - 0.5
         lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
   }
   if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
         if (!is.null(hline)) {
            abline(h = i - 0.5 + hline.vals, col = linecol, 
                   lty = 2)
         }
         yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
         yv <- rev(c(yv[1], yv))
         xv <- length(yv):1 - 0.5
         lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
   }
   if (!missing(cellnote)) 
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
           col = notecol, cex = notecex)
   par(mar = c(margins[1], 0, 0, 0))
   if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
   }
   else plot.new()
   par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
   if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
   }
   else plot.new()
   if (!is.null(main)) 
      title(main, cex.main = 1.5 * op[["cex.main"]])
   if (key) {
      mar <- c(5, 4, 2, 1)
      if (!is.null(key.xlab) && is.na(key.xlab)) 
         mar[1] <- 2
      if (!is.null(key.ylab) && is.na(key.ylab)) 
         mar[2] <- 2
      if (!is.null(key.title) && is.na(key.title)) 
         mar[3] <- 1
      par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
      if (length(key.par) > 0) 
         do.call(par, key.par)
      tmpbreaks <- breaks
      if (symkey) {
         max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
         min.raw <- -max.raw
         tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
         tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
         min.raw <- min(x, na.rm = TRUE)
         max.raw <- max(x, na.rm = TRUE)
      }
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      if (is.null(key.xtickfun)) {
         lv <- pretty(breaks)
         xv <- scale01(as.numeric(lv), min.raw, max.raw)
         xargs <- list(at = xv, labels = lv)
      }
      else {
         xargs <- key.xtickfun()
      }
      xargs$side <- 1
      do.call(axis, xargs)
      if (is.null(key.xlab)) {
         if (scale == "row") 
            key.xlab <- "Row Z-Score"
         else if (scale == "column") 
            key.xlab <- "Column Z-Score"
         else key.xlab <- "Value"
      }
      if (!is.na(key.xlab)) {
         mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
      }
      if (density.info == "density") {
         dens <- density(x, adjust = densadj, na.rm = TRUE)
         omit <- dens$x < min(breaks) | dens$x > max(breaks)
         dens$x <- dens$x[-omit]
         dens$y <- dens$y[-omit]
         dens$x <- scale01(dens$x, min.raw, max.raw)
         lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
               lwd = 1)
         if (is.null(key.ytickfun)) {
            yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                             0.95, labels = pretty(dens$y))
         }
         else {
            yargs <- key.ytickfun()
         }
         yargs$side <- 2
         do.call(axis, yargs)
         if (is.null(key.title)) 
            key.title <- "Color Key\nand Density Plot"
         if (!is.na(key.title)) 
            title(key.title)
         par(cex = 0.5)
         if (is.null(key.ylab)) 
            key.ylab <- "Density"
         if (!is.na(key.ylab)) 
            mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5)
      }
      else if (density.info == "histogram") {
         h <- hist(x, plot = FALSE, breaks = breaks)
         hx <- scale01(breaks, min.raw, max.raw)
         hy <- c(h$counts, h$counts[length(h$counts)])
         lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
               col = denscol)
         if (is.null(key.ytickfun)) {
            yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                          labels = pretty(hy))
         }
         else {
            yargs <- key.ytickfun()
         }
         yargs$side <- 2
         do.call(axis, yargs)
         if (is.null(key.title)) 
            key.title <- "Color Key\nand Histogram"
         if (!is.na(key.title)) 
            title(key.title)
         par(cex = 0.5)
         if (is.null(key.ylab)) 
            key.ylab <- "Count"
         if (!is.na(key.ylab)) 
            mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5)
      }
      else title("Color Key")
      if (trace %in% c("both", "column")) {
         vline.vals <- scale01(vline, min.raw, max.raw)
         if (!is.null(vline)) {
            abline(v = vline.vals, col = linecol, lty = 2)
         }
      }
      if (trace %in% c("both", "row")) {
         hline.vals <- scale01(hline, min.raw, max.raw)
         if (!is.null(hline)) {
            abline(v = hline.vals, col = linecol, lty = 2)
         }
      }
   }
   else plot.new()
   retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                   high = retval$breaks[-1], color = retval$col)
   if (!is.null(extrafun)) 
      extrafun()
   invisible(retval)
}

barplot.g<-function (height, width = 1, space = NULL, names.arg = NULL, 
                     legend.text = NULL, beside = FALSE, horiz = FALSE, density = NULL, 
                     angle = 45, col = NULL, border = par("fg"), main = NULL, 
                     sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
                     xpd = TRUE, log = "", axes = TRUE, axisnames = TRUE, cex.axis = par("cex.axis"), 
                     cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
                     axis.lty = 0, offset = 0, add = F, args.legend = NULL, ...) 
{
  if (!missing(inside)) 
    .NotYetUsed("inside", error = FALSE)
  if (is.null(space)) 
    space <- if (is.matrix(height) && beside) 
      c(0, 1)
  else 0.2
  space <- space * mean(width)
  if (plot && axisnames && is.null(names.arg)) 
    names.arg <- if (is.matrix(height)) 
      colnames(height)
  else names(height)
  if (is.vector(height) || (is.array(height) && (length(dim(height)) == 
                                                 1))) {
    height <- cbind(height)
    beside <- TRUE
    if (is.null(col)) 
      col <- "grey"
  }
  else if (is.matrix(height)) {
    if (is.null(col)) 
      col <- gray.colors(nrow(height))
  }
  else stop("'height' must be a vector or a matrix")
  if (is.logical(legend.text)) 
    legend.text <- if (legend.text && is.matrix(height)) 
      rownames(height)
  stopifnot(is.character(log))
  logx <- logy <- FALSE
  if (log != "") {
    logx <- length(grep("x", log)) > 0L
    logy <- length(grep("y", log)) > 0L
  }
  if ((logx || logy) && !is.null(density)) 
    stop("Cannot use shading lines in bars when log scale is used")
  NR <- nrow(height)
  NC <- ncol(height)
  if (beside) {
    if (length(space) == 2) 
      space <- rep.int(c(space[2L], rep.int(space[1L], 
                                            NR - 1)), NC)
    width <- rep_len(width, NR)
  }
  else {
    width <- rep_len(width, NC)
  }
  offset <- rep_len(as.vector(offset), length(width))
  delta <- width/2
  w.r <- cumsum(space + width)
  w.m <- w.r - delta
  w.l <- w.m - delta
  log.dat <- (logx && horiz) || (logy && !horiz)
  if (log.dat) {
    if (min(height + offset, na.rm = TRUE) <= 0) 
      stop("log scale error: at least one 'height + offset' value <= 0")
    if (logx && !is.null(xlim) && min(xlim) <= 0) 
      stop("log scale error: 'xlim' <= 0")
    if (logy && !is.null(ylim) && min(ylim) <= 0) 
      stop("log scale error: 'ylim' <= 0")
    rectbase <- if (logy && !horiz && !is.null(ylim)) 
      ylim[1L]
    else if (logx && horiz && !is.null(xlim)) 
      xlim[1L]
    else 0.9 * min(height, na.rm = TRUE)
  }
  else rectbase <- 0
  if (!beside) 
    height <- rbind(rectbase, apply(height, 2L, cumsum))
  rAdj <- offset + (if (log.dat) 
    0.9 * height
    else -0.01 * height)
  delta <- width/2
  w.r <- cumsum(space + width)
  w.m <- w.r - delta
  w.l <- w.m - delta
  if (horiz) {
    if (is.null(xlim)) 
      xlim <- range(rAdj, height + offset, na.rm = TRUE)
    if (is.null(ylim)) 
      ylim <- c(min(w.l), max(w.r))
  }
  else {
    if (is.null(xlim)) 
      xlim <- c(min(w.l), max(w.r))
    if (is.null(ylim)) 
      ylim <- range(rAdj, height + offset, na.rm = TRUE)
  }
  if (beside) 
    w.m <- matrix(w.m, ncol = NC)
  if (plot) {
    dev.hold()
    opar <- if (horiz) 
      par(xaxs = "i", xpd = xpd)
    else par(yaxs = "i", xpd = xpd)
    on.exit({
      dev.flush()
      par(opar)
    })
    if (!add) {
      plot.new()
      par(xaxs='i')
      par(yaxs='i')
      plot.window(xlim, ylim, log = log, ...)
    }

        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
                       ...) {
      if (horizontal) 
        rect(x1, y1, x2, y2, ...)
      else rect(y1, x1, y2, x2, ...)
    }
    if (beside){
      xyrect(rectbase + offset, w.l, c(height) + offset, 
             w.r, horizontal = horiz, angle = angle, density = density, 
             col = col, border = border)
    }else {
      for (i in 1L:NC) {
        xyrect(height[1L:NR, i] + offset[i], w.l[i], 
               height[-1, i] + offset[i], w.r[i], horizontal = horiz, 
               angle = angle, density = density, col = col, 
               border = border)
      }
    }
    if (axisnames && !is.null(names.arg)) {
      at.l <- if (length(names.arg) != length(w.m)) {
        if (length(names.arg) == NC) 
          colMeans(w.m)
        else stop("incorrect number of names")
      }
      else w.m
      axis(if (horiz) 
        2
        else 1, at = at.l, labels = names.arg, lty = axis.lty, 
        cex.axis = cex.names, ...)
    }
    if (!is.null(legend.text)) {
      legend.col <- rep_len(col, length(legend.text))
      if ((horiz & beside) || (!horiz & !beside)) {
        legend.text <- rev(legend.text)
        legend.col <- rev(legend.col)
        density <- rev(density)
        angle <- rev(angle)
      }
      xy <- par("usr")
      if (is.null(args.legend)) {
        legend(xy[2L] - xinch(0.1), xy[4L] - yinch(0.1), 
               legend = legend.text, angle = angle, density = density, 
               fill = legend.col, xjust = 1, yjust = 1)
      }
      else {
        args.legend1 <- list(x = xy[2L] - xinch(0.1), 
                             y = xy[4L] - yinch(0.1), legend = legend.text, 
                             angle = angle, density = density, fill = legend.col, 
                             xjust = 1, yjust = 1)
        args.legend1[names(args.legend)] <- args.legend
        do.call("legend", args.legend1)
      }
    }
    title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
          ...)
    if (axes) 
      axis(if (horiz) 
        1
        else 2, cex.axis = cex.axis, ...)
    invisible(w.m)
  }
  else w.m
}

