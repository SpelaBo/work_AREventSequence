traitData.plot.spela <- function (y, phy, col.label = "red", col.tree = "black", 
          col.hist = "navy", cex.plot = 0.7, cex.tips = 0.7, 
          show.tips = FALSE, include.hist = FALSE, n.split = 5, lwd.traits = 1, 
          show.axis = TRUE, axis.text = NULL, transform.axis.label = NULL, 
          at = NULL, labels = NULL, axis.text.line = 1, offset.bars = 1, 
          ...) 
{

  tree.height <- nodeTimes(phy)[1, 1] * offset.bars
  times <- tree.height/(max(y))
  range.fun <- function(x) (x - min(x))/(max(x) - min(x))
  if (is.null(at)) {
    splits <- (seq(0, max(y), length.out = n.split + 1) * 
                 times * 0.9) + (tree.height)
  }
  else {
    splits <- (at * times * 0.9) + (tree.height)
  }
  if (is.null(labels)) {
    label.splits <- seq(min(y), max(y), length.out = n.split + 
                          1)
  }
  else {
    label.splits <- labels
  }
  plot(ladderize(phy), show.tip.label = FALSE, edge.col = col.tree, 
       x.lim = c(0, tree.height * 2))
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x_cord <- lastPP$xx
  y_cord <- lastPP$yy
  tip.x.orig <- tip.x <- x_cord[1:Ntip(phy)]
  tip.y.orig <- tip.y <- y_cord[1:Ntip(phy)]
  t.data.scale.orig <- t.data.scale <- max(x_cord) + (range.fun(y) * 
                                                        (tree.height * 0.9))
  recycle.color <- rep(col.label, Ntip(phy))
  ord.y <- order(tip.y)
  tip.y <- tip.y[ord.y]
  tip.x <- tip.x[ord.y]
  recycle.color <- recycle.color[ord.y]
  t.data.scale <- t.data.scale[ord.y]
  zero.height <- tip.y[-length(tip.y)] + diff(tip.y)/2
  all.y <- rep(c(tip.y, zero.height), each = 3)
  all.x <- rep(tip.x[1], length(all.y))
  all.x[seq(2, length(all.y), 6)] <- t.data.scale
  for (uu in seq(0, length(all.x) - 7, 6)) all.x[c(4:6) + uu] <- NA
  ord.y2 <- order(all.y)
  all.y <- all.y[ord.y2]
  recycle.color <- rep(recycle.color, 6)
  for (x in seq(1, length(splits), 2)) polygon(rep(splits[(x):(x + 
                                                                 1)], each = 2), c(Ntip(phy), 1, 1, Ntip(phy)), col = "#00000020", 
                                               border = FALSE, xpd = TRUE)
  polygon(all.x * offset.bars, all.y, xpd = TRUE, col = NA, 
          border = recycle.color, lwd = lwd.traits)
  if (!is.null(transform.axis.label)) {
    if (transform.axis.label == "exp") 
      label.splits <- exp(label.splits)
    if (transform.axis.label == "exp10") 
      label.splits <- 10^(label.splits)
  }
  if (show.axis) 
    axis(3, at = splits, labels = signif(label.splits, 3), 
         xpd = TRUE, line = -0.5, ...)
  if (show.tips) 
    text(t.data.scale.orig, tip.y.orig, phy$tip.label, pos = 4, 
         cex = cex.tips, xpd = TRUE)
  if (!is.null(axis.text)) 
    mtext(axis.text, 3, at = mean(splits), line = axis.text.line, 
          cex = cex.plot)
  if (include.hist) {
    par(mar = c(3, 3, 3, 3))
    if (is.null(colnames(y))) {
      name.trait <- "trait"
    }
    else {
      name.trait <- colnames(y)
    }
    hist(y, col = col.hist, xaxs = "i", yaxt = "n", 
         border = "white", main = "", xlab = "", 
         ylab = "", cex.axis = cex.plot, las = 2, axes = TRUE, 
         mgp = c(3, 1, 0))
    box()
    axis(2, cex.axis = cex.plot, mgp = c(3, 1, 0))
    mtext("frequency", 2, line = axis.text.line, cex = cex.plot)
    mtext(name.trait, 1, line = axis.text.line, cex = cex.plot)
  }
}