# testing BM, EB and OU model for each of the trait groups, clade B #########
# performance traits: 
## trophic niche: gpI6_size in gpII6_size; MORPHO_GLS[[i]][,20:21]
## body shape: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width; MORPHO_GLS[[i]][,c(7:11)]
## locomotion: body_length, pV_length, pVI_length, pVII_length; MORPHO_GLS[[i]][,c(1,4:6)]
## sensoric: antenna I and II; MORPHO_GLS[[i]][,c(2:3)]

# preparing the data
## tree
tree_a1 <- CLADES_MORPHO_GLS[[1]]
is.ultrametric(tree_a1)
## distance to the root (will need for reconstructions)
dist_to_root_a1 <- node.depth.edgelength(tree_a1)
## morpho data
data_a1 <- MORPHO_GLS[[1]]
## sort as node labels
data_a1<-data_a1[tree_a1$tip.label,]
## check
all(rownames(data_a1)==tree_a1$tip.label)
name.check(tree_a1, data_a1)

# models trophic niche ####
model_a1M_gnato_a1 <- mvBM(tree_a1,data_a1[,20:21])
model_EB_gnato_a1 <- mvEB(tree_a1,data_a1[,20:21])
model_OU_gnato_a1 <- mvOU(tree_a1,data_a1[,20:21])
# merge into list
models_gnato_a1 <- list(model_a1M_gnato_a1,model_EB_gnato_a1,model_OU_gnato_a1)
# choose the best model
gnato_a1_aic <- aicw(models_gnato_a1)
gnato_a1_aic

# models body shape ####
model_a1M_body_a1 <- mvBM(tree_a1,data_a1[,c(7:11)])
model_EB_body_a1 <- mvEB(tree_a1,data_a1[,c(7:11)])
model_OU_body_a1 <- mvOU(tree_a1,data_a1[,c(7:11)])
# merge into list
models_body_a1 <- list(model_a1M_body_a1,model_EB_body_a1,model_OU_body_a1)
# choose the best model
body_a1_aic <- aicw(models_body_a1)
body_a1_aic

# models locomotion ####
model_a1M_loco_a1 <- mvBM(tree_a1,data_a1[,c(1,4:6)])
model_EB_loco_a1 <- mvEB(tree_a1,data_a1[,c(1,4:6)])
model_OU_loco_a1 <- mvOU(tree_a1,data_a1[,c(1,4:6)])
# merge into list
models_loco_a1 <- list(model_a1M_loco_a1,model_EB_loco_a1,model_OU_loco_a1)
# choose the best model
loco_a1_aic <- aicw(models_loco_a1)
loco_a1_aic

# models sensoric ####
model_a1M_sens_a1 <- mvBM(tree_a1,data_a1[,c(2:3)])
model_EB_sens_a1 <- mvEB(tree_a1,data_a1[,c(2:3)])
model_OU_sens_a1 <- mvOU(tree_a1,data_a1[,c(2:3)])
# merge into list
models_sens_a1 <- list(model_a1M_sens_a1,model_EB_sens_a1,model_OU_sens_a1)
# choose the best model
sens_a1_aic <- aicw(models_sens_a1)
sens_a1_aic

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato_a1<-estim(tree_a1, data_a1[,20:21], model_OU_gnato_a1, asr=TRUE)
comb_gnato_a1 <- rbind(data_a1[,20:21],ASR_gnato_a1$estim)

ASR_body_a1<-estim(tree_a1, data_a1[,c(7:11)], model_OU_body_a1, asr=TRUE)
comb_body_a1 <- rbind(data_a1[,c(7:11)],ASR_body_a1$estim)

ASR_loco_a1<-estim(tree_a1, data_a1[,c(1,4:6)], model_OU_loco_a1, asr=TRUE)
comb_loco_a1 <- rbind(data_a1[,c(1,4:6)],ASR_loco_a1$estim)

ASR_sens_a1<-estim(tree_a1, data_a1[,c(2:3)], model_OU_sens_a1, asr=TRUE)
comb_sens_a1 <- rbind(data_a1[,c(2:3)],ASR_sens_a1$estim)

comb_a1 <- cbind(comb_body_a1,comb_loco_a1,comb_sens_a1,comb_gnato_a1)

# za vsak stolpec v comb b dobi koordinate 50
T50_a1lade_a1 <- c()
for (i in 1:ncol(comb_a1)) {
  test <- coordinatesTimeslices(tree_a1, comb_a1[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T50_a1lade_a1[i] <- c(coordinatesXXMorphospace(test_grouped,50))$sample
}

T90_a1lade_a1 <- c()
for (i in 1:ncol(comb_a1)) {
  test <- coordinatesTimeslices(tree_a1, comb_a1[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T90_a1lade_a1[i] <- c(coordinatesXXMorphospace(test_grouped,90))$sample
  
}

par(mfrow=c(4,4))
for (i in 1:ncol(comb_a1)) {
  phenogram(tree=tree_a1,x=comb_a1[,i], ftype="off", ylim=c(min(comb_a1[,i]),max(comb_a1[,i])), xlab = colnames(comb_a1)[i])
  abline(v=T50_a1lade_a1[i], col = "blue")
  abline(v=T90_a1lade_a1[i], col = "red")
}
