library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(TDbook)
## ref: http://www.phytools.org/eqg2015/asr.html
##
## load `tree_anole` and `df_svl` from 'TDbook'
svl <- as.matrix(df_svl)[,1]
fit <- phytools::fastAnc(tree_anole, svl, vars=TRUE, CI=TRUE)

td <- data.frame(node = nodeid(tree_anole, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

comb_b_gg <- as.data.frame(comb_b)
comb_b_gg$node <- c(1:nrow(comb_b_gg))

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree_b, comb_b_gg, by = 'node')

p1 <- ggtree(tree, aes(color=gpII6_size), layout = 'circular', 
             ladderize = FALSE, continuous = 'colour', size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(hjust = -.1) + 
  theme(legend.position = c(.05, .85)) 

p1

ggtree(tree, aes(color=gpII6_size), continuous = 'colour', yscale = "gpII6_size") + 
  scale_color_viridis_c() + theme_minimal()

# plot reconstructed phenograms ####
par(mfrow=c(4,4))
for (i in 1:ncol(comb_b)) {
  phenogram(tree=tree_b,x=comb_b[,i], ftype="off", ylim=c(min(comb_b[,i]),max(comb_b[,i])), xlab = colnames(comb_b)[i])
}

# plot reconstructed phenograms ####
par(mfrow=c(4,4))
for (i in 1:ncol(comb_b)) {
  ggtree(tree, aes(color=tree), continuous = 'colour', yscale = "trait") + 
    scale_color_viridis_c() + theme_minimal()
}




plot(tree_mean_new)