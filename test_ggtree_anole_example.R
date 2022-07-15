# testing BM, EB and OU model for each of the trait groups, clade B #########
# performance traits: 
## trophic niche: gpI6_size in gpII6_size; MORPHO_GLS[[i]][,20:21]
## body shape: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width; MORPHO_GLS[[i]][,c(7:11)]
## locomotion: body_length, pV_length, pVI_length, pVII_length; MORPHO_GLS[[i]][,c(1,4:6)]
## sensoric: antenna I and II; MORPHO_GLS[[i]][,c(2:3)]

# preparing the data
## tree
tree_a2 <- CLADES_MORPHO_GLS[[2]]
is.ultrametric(tree_a2)
## distance to the root (will need for reconstructions)
dist_to_root_a2 <- node.depth.edgelength(tree_a2)
## morpho data
data_a2 <- MORPHO_GLS[[2]]
## sort as node labels
data_a2<-data_a2[tree_a2$tip.label,]
## check
all(rownames(data_a2)==tree_a2$tip.label)
name.check(tree_a2, data_a2)

# models trophic niche ####
model_a2M_gnato_a2 <- mvBM(tree_a2,data_a2[,20:21])
model_EB_gnato_a2 <- mvEB(tree_a2,data_a2[,20:21])
model_OU_gnato_a2 <- mvOU(tree_a2,data_a2[,20:21])
# merge into list
models_gnato_a2 <- list(model_a2M_gnato_a2,model_EB_gnato_a2,model_OU_gnato_a2)
# choose the best model
gnato_a2_aic <- aicw(models_gnato_a2)
gnato_a2_aic

# models body shape ####
model_a2M_body_a2 <- mvBM(tree_a2,data_a2[,c(7:11)])
model_EB_body_a2 <- mvEB(tree_a2,data_a2[,c(7:11)])
model_OU_body_a2 <- mvOU(tree_a2,data_a2[,c(7:11)])
# merge into list
models_body_a2 <- list(model_a2M_body_a2,model_EB_body_a2,model_OU_body_a2)
# choose the best model
body_a2_aic <- aicw(models_body_a2)
body_a2_aic

# models locomotion ####
model_a2M_loco_a2 <- mvBM(tree_a2,data_a2[,c(1,4:6)])
model_EB_loco_a2 <- mvEB(tree_a2,data_a2[,c(1,4:6)])
model_OU_loco_a2 <- mvOU(tree_a2,data_a2[,c(1,4:6)])
# merge into list
models_loco_a2 <- list(model_a2M_loco_a2,model_EB_loco_a2,model_OU_loco_a2)
# choose the best model
loco_a2_aic <- aicw(models_loco_a2)
loco_a2_aic

# models sensoric ####
model_a2M_sens_a2 <- mvBM(tree_a2,data_a2[,c(2:3)])
model_EB_sens_a2 <- mvEB(tree_a2,data_a2[,c(2:3)])
model_OU_sens_a2 <- mvOU(tree_a2,data_a2[,c(2:3)])
# merge into list
models_sens_a2 <- list(model_a2M_sens_a2,model_EB_sens_a2,model_OU_sens_a2)
# choose the best model
sens_a2_aic <- aicw(models_sens_a2)
 

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato_a2<-estim(tree_a2, data_a2[,20:21], model_OU_gnato_a2, asr=TRUE)
comb_gnato_a2 <- rbind(data_a2[,20:21],ASR_gnato_a2$estim)

ASR_body_a2<-estim(tree_a2, data_a2[,c(7:11)], model_OU_body_a2, asr=TRUE)
comb_body_a2 <- rbind(data_a2[,c(7:11)],ASR_body_a2$estim)

ASR_loco_a2<-estim(tree_a2, data_a2[,c(1,4:6)], model_OU_loco_a2, asr=TRUE)
comb_loco_a2 <- rbind(data_a2[,c(1,4:6)],ASR_loco_a2$estim)

ASR_sens_a2<-estim(tree_a2, data_a2[,c(2:3)], model_OU_sens_a2, asr=TRUE)
comb_sens_a2 <- rbind(data_a2[,c(2:3)],ASR_sens_a2$estim)

comb_a2 <- cbind(comb_body_a2,comb_loco_a2,comb_sens_a2,comb_gnato_a2)

# za vsak stolpec v comb b dobi koordinate 50
T50_a2lade_a2 <- c()
for (i in 1:ncol(comb_a2)) {
  test <- coordinatesTimeslices(tree_a2, comb_a2[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T50_a2lade_a2[i] <- c(coordinatesXXMorphospace(test_grouped,50))$sample
}

T90_a2lade_a2 <- c()
for (i in 1:ncol(comb_a2)) {
  test <- coordinatesTimeslices(tree_a2, comb_a2[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T90_a2lade_a2[i] <- c(coordinatesXXMorphospace(test_grouped,90))$sample
  
}

par(mfrow=c(4,4))
for (i in 1:ncol(comb_a2)) {
  phenogram(tree=tree_a2,x=comb_a2[,i], ftype="off", ylim=c(min(comb_a2[,i]),max(comb_a2[,i])), xlab = colnames(comb_a2)[i])
  abline(v=T50_a2lade_a2[i], col = "blue")
  abline(v=T90_a2lade_a2[i], col = "red")
}
