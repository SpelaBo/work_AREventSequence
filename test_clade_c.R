# testing BM, EB and OU model for each of the trait groups, clade B #########
# performance traits: 
## trophic niche: gpI6_size in gpII6_size; MORPHO_GLS[[i]][,20:21]
## body shape: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width; MORPHO_GLS[[i]][,c(7:11)]
## locomotion: body_length, pV_length, pVI_length, pVII_length; MORPHO_GLS[[i]][,c(1,4:6)]
## sensoric: antenna I and II; MORPHO_GLS[[i]][,c(2:3)]

# preparing the data
## tree
tree_c <- CLADES_MORPHO_GLS[[4]]
is.ultrametric(tree_c)
## distance to the root (will need for reconstructions)
dist_to_root_c <- node.depth.edgelength(tree_c)
## morpho data
data_c <- MORPHO_GLS[[4]]
## sort as node labels
data_c<-data_c[tree_c$tip.label,]
## check
all(rownames(data_c)==tree_c$tip.label)
name.check(tree_c, data_c)

# models trophic niche ####
model_cM_gnato_c <- mvBM(tree_c,data_c[,20:21])
model_EB_gnato_c <- mvEB(tree_c,data_c[,20:21])
model_OU_gnato_c <- mvOU(tree_c,data_c[,20:21])
# merge into list
models_gnato_c <- list(model_cM_gnato_c,model_EB_gnato_c,model_OU_gnato_c)
# choose the best model
gnato_c_aic <- aicw(models_gnato_c)
gnato_c_aic

# models body shape ####
model_cM_cody_c <- mvBM(tree_c,data_c[,c(7:11)])
model_EB_cody_c <- mvEB(tree_c,data_c[,c(7:11)])
model_OU_cody_c <- mvOU(tree_c,data_c[,c(7:11)])
# merge into list
models_cody_c <- list(model_cM_cody_c,model_EB_cody_c,model_OU_cody_c)
# choose the best model
body_c_aic <- aicw(models_cody_c)
body_c_aic

# models locomotion ####
model_cM_loco_c <- mvBM(tree_c,data_c[,c(1,4:6)])
model_EB_loco_c <- mvEB(tree_c,data_c[,c(1,4:6)])
model_OU_loco_c <- mvOU(tree_c,data_c[,c(1,4:6)])
# merge into list
models_loco_c <- list(model_cM_loco_c,model_EB_loco_c,model_OU_loco_c)
# choose the best model
loco_c_aic <- aicw(models_loco_c)
loco_c_aic

# models sensoric ####
model_cM_sens_c <- mvBM(tree_c,data_c[,c(2:3)])
model_EB_sens_c <- mvEB(tree_c,data_c[,c(2:3)])
model_OU_sens_c <- mvOU(tree_c,data_c[,c(2:3)])
# merge into list
models_sens_c <- list(model_cM_sens_c,model_EB_sens_c,model_OU_sens_c)
# choose the best model
sens_c_aic <- aicw(models_sens_c)
sens_c_aic

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato_c<-estim(tree_c, data_c[,20:21], model_OU_gnato_c, asr=TRUE)
comb_gnato_c <- rbind(data_c[,20:21],ASR_gnato_c$estim)

ASR_cody_c<-estim(tree_c, data_c[,c(7:11)], model_OU_cody_c, asr=TRUE)
comb_cody_c <- rbind(data_c[,c(7:11)],ASR_cody_c$estim)

ASR_loco_c<-estim(tree_c, data_c[,c(1,4:6)], model_OU_loco_c, asr=TRUE)
comb_loco_c <- rbind(data_c[,c(1,4:6)],ASR_loco_c$estim)

ASR_sens_c<-estim(tree_c, data_c[,c(2:3)], model_OU_sens_c, asr=TRUE)
comb_sens_c <- rbind(data_c[,c(2:3)],ASR_sens_c$estim)

comb_c <- cbind(comb_cody_c,comb_loco_c,comb_sens_c,comb_gnato_c)

# za vsak stolpec v comb b dobi koordinate 50
T50_clade_c <- c()
for (i in 1:ncol(comb_c)) {
  test <- coordinatesTimeslices(tree_c, comb_c[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T50_clade_c[i] <- c(coordinatesXXMorphospace(test_grouped,50))$sample
}

T90_clade_c <- c()
for (i in 1:ncol(comb_c)) {
  test <- coordinatesTimeslices(tree_c, comb_c[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T90_clade_c[i] <- c(coordinatesXXMorphospace(test_grouped,90))$sample
  
}

par(mfrow=c(4,4))
for (i in 1:ncol(comb_c)) {
  phenogram(tree=tree_c,x=comb_c[,i], ftype="off", ylim=c(min(comb_c[,i]),max(comb_c[,i])), xlab = colnames(comb_c)[i])
  abline(v=T50_clade_c[i], col = "blue")
  abline(v=T90_clade_c[i], col = "red")
}
