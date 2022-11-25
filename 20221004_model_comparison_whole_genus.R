### analiza diverzifikacije, rekonstrukcija na celotnem drevesu ####

#### Data - see 20220713_model_comparison for import####
# final tree, max credibility tree, mean heights
tree_mean
# Morphology, continuous
morphology # whole, including NAs
morphology_noNA # Remove NA-s
tree_morpho # Adjusted tree to data
is.ultrametric(tree_morpho)
tree_morpho_l # ladderized
morpho_all_gls # Phylogenetically corrected GLS
morpho_all_gls<-morpho_all_gls[order(match(rownames(morpho_all_gls), tree_morpho_clades$tip.label)), , drop = FALSE]
all(rownames(morpho_all_gls)==tree_morpho$tip.label)

tree_morpho_clades <-tree_morpho 
is.ultrametric(tree_morpho_clades)
tree_morpho_clades<-paintSubTree(tree_morpho_clades,node=mrcaA1,state="A1",stem=T)
tree_morpho_clades<-paintSubTree(tree_morpho_clades,node=mrcaA2,state="A2",stem=T)
tree_morpho_clades<-paintSubTree(tree_morpho_clades,node=mrcaB,state="B",stem=T)
tree_morpho_clades<-paintSubTree(tree_morpho_clades,node=mrcaC,state="C",stem=T)
cols<-c("black","deepskyblue3","chartreuse3","gold", "pink"); names(cols)<-c(1,"A1","A2","B","C")
plotSimmap(tree_morpho_clades,cols,lwd=2,pts=F, ftype="off")

#tree_morpho_nv=drop.tip(tree_morpho, "N_vadimi")
#morpho_all_gls_nv <- subset(morpho_all_gls, rownames(morpho_all_gls) %in% tree_morpho_nv$tip.label)

# testing BM, EB and OU model for each of the trait groups, whole niphargus #########
# performance traits: 
## trophic niche: gpI6_size in gpII6_size; morpho_all_gls[,20:21]
## body shape: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width; morpho_all_gls[,c(7:11)]
## locomotion: body_length, pV_length, pVI_length, pVII_length; morpho_all_gls[,c(1,4:6)]
## sensoric: antenna I and II; morpho_all_gls[,c(2:3)]
# models trophic niche ####
tic()
model_BM_gnato <- mvBM(tree_morpho,morpho_all_gls[,20:21])
toc() # 3.5s, for pseudoinverse 121.41
model_BMM_gnato <- mvBM(tree_morpho_clades,morpho_all_gls[,20:21],model = "BMM", param = list(root = F))
tic()
model_EB_gnato <- mvEB(tree_morpho,morpho_all_gls[,20:21])
toc() # 4.09s
tic()
model_OU_gnato <- mvOU(tree_morpho,morpho_all_gls[,20:21], param = list(root = F))
toc() # 25.88 s, for pseudoinverse 545.64s
tic()
model_OUM_gnato <- mvOU(tree_morpho_clades,morpho_all_gls[,20:21],model = "OUM", param = list(root = F))
toc()
# merge into list
models_gnato <- list(model_BM_gnato,model_BMM_gnato,model_EB_gnato,model_OU_gnato, model_OUM_gnato)
# choose the best model
gnato_aic <- aicw(models_gnato)
gnato_aic

# models body shape ####
tic()
model_BM_body <- mvBM(tree_morpho,morpho_all_gls[,c(7:11)])
toc() # 1218.84 s (20 min)
model_BMM_body <- mvBM(tree_morpho_clades,morpho_all_gls[,c(7:11)],model = "BMM", param = list(root = F))
tic()
model_EB_body <- mvEB(tree_morpho,morpho_all_gls[,c(7:11)], method = "sparse")
toc() # 1647.94 sec elapsed
tic()
model_OU_body <- mvOU(tree_morpho,morpho_all_gls[,c(7:11)], param = list(root = F))
toc() # 7947.5s
tic()
model_OUM_body <- mvOU(tree_morpho_clades,morpho_all_gls[,c(7:11)],model = "OUM", param = list(root = F))
toc()
# merge into list
models_body <- list(model_BM_body,model_BMM_body,model_EB_body,model_OU_body,model_OUM_body)
# choose the best model
body_aic <- aicw(models_body)
body_aic
save.image()

# models locomotion ####
model_BM_loco <- mvBM(tree_morpho,morpho_all_gls[,c(1,4:6)])
model_BMM_loco <- mvBM(tree_morpho_clades,morpho_all_gls[,c(1,4:6)],model = "BMM", param = list(root = F))
model_EB_loco <- mvEB(tree_morpho,morpho_all_gls[,c(1,4:6)])
model_OU_loco <- mvOU(tree_morpho,morpho_all_gls[,c(1,4:6)], param = list(root = F))
model_OUM_loco <- mvOU(tree_morpho_clades,morpho_all_gls[,c(1,4:6)],model = "OUM", param = list(root = F))

# merge into list
models_loco <- list(model_BM_loco,model_BMM_loco,model_EB_loco,model_OU_loco,model_OUM_loco)
# choose the best model
loco_aic <- aicw(models_loco)
loco_aic
save.image()

# models sensoric ####
model_BM_sens <- mvBM(tree_morpho,morpho_all_gls[,c(2:3)])
model_BMM_sens <- mvBM(tree_morpho_clades,morpho_all_gls[,c(2:3)],model = "BMM")
model_EB_sens <- mvEB(tree_morpho,morpho_all_gls[,c(2:3)])
model_OU_sens <- mvOU(tree_morpho,morpho_all_gls[,c(2:3)], param = list(root = F))
model_OUM_sens <- mvOU(tree_morpho_clades,morpho_all_gls[,c(2:3)],model = "OUM", param = list(root = F))

# merge into list
models_sens <- list(model_BM_sens,model_BMM_sens,model_EB_sens,model_OU_sens,model_OUM_sens)
# choose the best model
sens_aic <- aicw(models_sens)
sens_aic
save.image()

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato<-estim(tree_morpho_clades,morpho_all_gls[,20:21], model_OUM_gnato, asr=TRUE)
comb_gnato <- rbind(morpho_all_gls[,20:21],ASR_gnato$estim)

ASR_body<-estim(tree_morpho_clades,morpho_all_gls[,c(7:11)], model_BMM_body, asr=TRUE)
comb_body <- rbind(morpho_all_gls[,c(7:11)],ASR_body$estim)

ASR_loco<-estim(tree_morpho_clades,morpho_all_gls[,c(1,4:6)], model_OUM_loco, asr=TRUE)
comb_loco <- rbind(morpho_all_gls[,c(1,4:6)],ASR_loco$estim)

ASR_sens<-estim(tree_morpho_clades,morpho_all_gls[,c(2:3)], model_BMM_sens, asr=TRUE)
comb_sens <- rbind(morpho_all_gls[,c(2:3)],ASR_sens$estim)

comb <- cbind(comb_body,comb_loco,comb_sens,comb_gnato)
col.order <- colnames(morpho_all_gls[,c(1:11,20:21)])
comb <- comb[ , col.order]

# plot reconstructed phenograms ####
par(mfrow=c(4,4))
for (i in 1:ncol(comb)) {
  phenogram(tree=tree_morpho_clades,x=comb[,i], ftype="off", ylim=c(min(comb[,i]),max(comb[,i])), xlab = colnames(comb)[i])
}

par(mfrow=c(4,4))
for (i in 1:ncol(comb)) {
  phenogram(tree=tree_morpho_clades,x=morpho_all_gls[,i], ftype="off", ylim=c(min(comb[,i]),max(comb[,i])), xlab = colnames(comb)[i])
}


#subset of species
speciesA1_m <- subset(speciesA1, speciesA1 %in% tree_morpho$tip.label)
speciesA2_m <- subset(speciesA2, speciesA2 %in% tree_morpho$tip.label)
speciesC_m <- subset(speciesC, speciesC %in% tree_morpho$tip.label)
speciesB_m <- subset(speciesB, speciesB %in% tree_morpho$tip.label)

#mrca
mrcaA1 <- findMRCA(tree_morpho, speciesA1_m)
mrcaA2 <- findMRCA(tree_morpho, speciesA2_m)
mrcaC <- findMRCA(tree_morpho, speciesC_m)
mrcaB <- findMRCA(tree_morpho, speciesB_m)


tree_tibble <- as_tibble(tree_morpho_clades)
nodeHeights(tree_morpho_clades)
nodeheight(tree_morpho_clades,330)
tree_morpho_clades$edge
branching.times(tree_morpho_clades)

ASR_mrca_clades <- subset(comb,rownames(comb)%in% c(mrcaA1,mrcaA2,mrcaB,mrcaC))