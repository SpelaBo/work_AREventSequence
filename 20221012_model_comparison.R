### analiza diverzifikacije, rekonstrukcija na celotnem drevesu ####

# functions for modeling ####
allModels <-  function(tree, morpho_df, traits) {
  model_BM <- mvBM(tree,morpho_df[,traits],model = "BM1")
  model_BMM <- mvBM(tree,morpho_df[,traits],model = "BMM", param = list(root = F))
  model_EB <- mvEB(tree,morpho_df[,traits])
  model_OU <- mvOU(tree,morpho_df[,traits],model = "OU1", param = list(root = F))
  model_OUM <- mvOU(tree,morpho_df[,traits],model = "OUM", param = list(root = F))
  models <- list(model_BM,model_BMM,model_EB,model_OU, model_OUM)
  return(models)
}
testBestModel <- function(models){
  aic <- aicw(models)
  best <- models[[which.max(aic$aicweights)]]
  return(best)
}

# Modelling, genus level ####
models_list <- list()
for (i in 1:length(traits)) {
  trait <- traits[[i]]
  models <- allModels(tree_morpho_painted,morpho_all_gls,trait)
  models_list[[i]] <- models
}
names(models_list) <- names(traits)

best_model_list <- list()
for (i in 1:length(models_list)) {
  models <- models_list[[i]]
  best_model <- testBestModel(models)
  best_model_list[[i]] <- best_model
}
names(best_model_list) <- names(models_list)

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_list <-  list()
for (i in 1:length(traits)) {
  model <- best_model_list[[i]]
  trait <- traits[[i]]
  ASR<-estim(tree_morpho_painted,morpho_all_gls[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(morpho_all_gls[,trait],ASR$estim))
  comb$species <- rownames(comb)
  
  ASR_list[[i]] <- comb
}
ASR_df <-  Reduce(inner_join,ASR_list)
rownames(ASR_df) <- ASR_df$species
ASR_df[c("species")] <- NULL


# clade subsetting ####
## Node translation and subsets of reconstructed morphologies for clades
selectTipNode <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}
nodeTranslation <- function(trees){
  node_translation_df <- list()
  for (i in 1:length(trees)) {
    tree <- trees[[i]]
    t <- data.frame(
      "node" = tree$edge[,1],
      "label" = sapply(tree$edge[,1],
                       selectTipNode,
                       tree = tree))
    t <- t[!duplicated(t), ]
    t <- t[order(t$node),]
    node_translation_df[[i]] <- t
  }
  return(node_translation_df)
}
#trees=tree_morpho_clades,translation = node_translation,traits = comb, descendants = desc_clades, nameclades = nameclades
reconstructionsCladeSubsets <- function(trees,translation,traits, descendants,nameclades) {
  comb_clades <- list()
  for (i in 1:length(trees)) {
    tree <- trees[[i]]
    nnode <- tree$Nnode
    ntip <- nnode+1
    node_seq <- c((ntip+1):(ntip+nnode))
    desc <- descendants[[i]]
    comb_sub <- as.data.frame(traits[desc,])
    comb_sub$node <- rownames(comb_sub)
    nodes <- translation[[i]]
    comb_sub[node_seq,]$node <- nodes$node
    rownames(comb_sub) <- comb_sub$node
    comb_sub <- subset(comb_sub, select = -node)
    comb_clades[[i]] <- comb_sub
  }
  names(comb_clades) <- nameclades
  return(comb_clades)
}


node_translation <- nodeTranslation(tree_morpho_clades)

# ASR for each clade ####
ASR_clades <- reconstructionsCladeSubsets(trees = tree_morpho_clades,
                                                 translation = node_translation,
                                                 traits = ASR_df, 
                                                 descendants = desc_clades,
                                                 nameclades = nameclades)

design <- matrix(c(1,2,3,4,5,6,7,8,0,9,10,0,11,12,0), nrow=3, byrow = F)
for (j in 1:length(nameclades)) {
  c=as.matrix(ASR_clades[[j]])
  #c <- c[,sort]
  t=tree_morpho_clades[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno_new.pdf"))
  layout(design)
  par(mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}


design <- matrix(c(1,2,3,4,5,6,7,8,0,9,10,0,11,12,0), nrow=3, byrow = F)
design3 <- matrix(c(1,2,3,4,5,6:9,0,10,11,0,0,0,12,13,0,0,0,14,15,0,0,0), nrow=5, byrow = F)

for (j in 1:length(nameclades)) {
  c=as.matrix(ASR_clades[[j]])
  #c <- c[,sort]
  t=tree_morpho_clades[[j]]
  tips=t$tip.label
  c=morpho_all_gls[rownames(morpho_all_gls)%in%tips,]
  title <- nameclades[j]
  layout(design3)
  par(mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
}




















ASR_body<-estim(tree_morpho_painted,morpho_all_gls[,body_traits], model_BMM_body, asr=TRUE)
comb_body <- rbind(morpho_all_gls[,body_traits],ASR_body$estim)

ASR_loco<-estim(tree_morpho_painted,morpho_all_gls[,loco_traits], model_OUM_loco, asr=TRUE)
comb_loco <- rbind(morpho_all_gls[,loco_traits],ASR_loco$estim)

ASR_sens<-estim(tree_morpho_painted,morpho_all_gls[,sens_traits], model_BMM_sens, asr=TRUE)
comb_sens <- rbind(morpho_all_gls[,sens_traits],ASR_sens$estim)

comb <- cbind(comb_body,comb_loco,comb_sens,comb_gnato)
comb <- comb[ ,colnames(morpho_all_gls)]












# Modelling, genus level ####
# testing BM, EB and OU model for each of the trait groups, whole niphargus
## trophic_traits: gpI6_size in gpII6_size
## body_traits: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width
## loco_traits: body_length, pV_length, pVI_length, pVII_length
## sens_traits: antenna I and II
# .trophic niche ####
model_BM_gnato <- mvBM(tree_morpho_painted,morpho_all_gls[,trophic_traits],model = "BM1")
model_BMM_gnato <- mvBM(tree_morpho_painted,morpho_all_gls[,trophic_traits],model = "BMM", param = list(root = F))
model_EB_gnato <- mvEB(tree_morpho,morpho_all_gls[,trophic_traits])
model_OU_gnato <- mvOU(tree_morpho,morpho_all_gls[,trophic_traits], param = list(root = F))
model_OUM_gnato <- mvOU(tree_morpho_painted,morpho_all_gls[,trophic_traits],model = "OUM", param = list(root = F))
# merge into list
models_gnato <- list(model_BM_gnato,model_BMM_gnato,model_EB_gnato,model_OU_gnato, model_OUM_gnato)
# choose the best model
gnato_aic <- mvMORPH::aicw(models_gnato)
gnato_aic

# .body shape ####
tic()
model_BM_body <- mvBM(tree_morpho,morpho_all_gls[,body_traits])
toc() # 1218.84 s (20 min)
model_BMM_body <- mvBM(tree_morpho_painted,morpho_all_gls[,body_traits],model = "BMM", param = list(root = F))
tic()
model_EB_body <- mvEB(tree_morpho,morpho_all_gls[,body_traits], method = "sparse")
toc() # 1647.94 sec elapsed
tic()
model_OU_body <- mvOU(tree_morpho,morpho_all_gls[,body_traits], param = list(root = F))
toc() # 7947.5s
tic()
model_OUM_body <- mvOU(tree_morpho_painted,morpho_all_gls[,body_traits],model = "OUM", param = list(root = F))
toc()
# merge into list
models_body <- list(model_BM_body,model_BMM_body,model_EB_body,model_OU_body,model_OUM_body)
# choose the best model
body_aic <- aicw(models_body)
body_aic
save.image()

# .locomotion ####
model_BM_loco <- mvBM(tree_morpho,morpho_all_gls[,loco_traits])
model_BMM_loco <- mvBM(tree_morpho_painted,morpho_all_gls[,loco_traits],model = "BMM", param = list(root = F))
model_EB_loco <- mvEB(tree_morpho,morpho_all_gls[,loco_traits])
model_OU_loco <- mvOU(tree_morpho,morpho_all_gls[,loco_traits], param = list(root = F))
model_OUM_loco <- mvOU(tree_morpho_painted,morpho_all_gls[,loco_traits],model = "OUM", param = list(root = F))
save.image()
# merge into list & choose the best model
models_loco <- list(model_BM_loco,model_BMM_loco,model_EB_loco,model_OU_loco,model_OUM_loco)
loco_aic <- aicw(models_loco)
best_loco <- models_loco[[which.max(loco_aic$aicweights)]]

# .sensoric ####
model_BM_sens <- mvBM(tree_morpho,morpho_all_gls[,sens_traits])
model_BMM_sens <- mvBM(tree_morpho_painted,morpho_all_gls[,sens_traits],model = "BMM")
model_EB_sens <- mvEB(tree_morpho,morpho_all_gls[,sens_traits])
model_OU_sens <- mvOU(tree_morpho,morpho_all_gls[,sens_traits], param = list(root = F))
model_OUM_sens <- mvOU(tree_morpho_painted,morpho_all_gls[,sens_traits],model = "OUM", param = list(root = F))

# merge into list
models_sens <- list(model_BM_sens,model_BMM_sens,model_EB_sens,model_OU_sens,model_OUM_sens)
# choose the best model
sens_aic <- aicw(models_sens)
sens_aic
save.image()

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato<-estim(tree_morpho_painted,morpho_all_gls[,trophic_traits], model_OUM_gnato, asr=TRUE)
comb_gnato <- rbind(morpho_all_gls[,trophic_traits],ASR_gnato$estim)

ASR_body<-estim(tree_morpho_painted,morpho_all_gls[,body_traits], model_BMM_body, asr=TRUE)
comb_body <- rbind(morpho_all_gls[,body_traits],ASR_body$estim)

ASR_loco<-estim(tree_morpho_painted,morpho_all_gls[,loco_traits], model_OUM_loco, asr=TRUE)
comb_loco <- rbind(morpho_all_gls[,loco_traits],ASR_loco$estim)

ASR_sens<-estim(tree_morpho_painted,morpho_all_gls[,sens_traits], model_BMM_sens, asr=TRUE)
comb_sens <- rbind(morpho_all_gls[,sens_traits],ASR_sens$estim)

comb <- cbind(comb_body,comb_loco,comb_sens,comb_gnato)
comb <- comb[ ,colnames(morpho_all_gls)]

# plot reconstructed phenograms ####
par(mfrow=c(4,4))
for (i in 1:ncol(comb)) {
  phenogram(tree=tree_morpho_painted,x=comb[,i], ftype="off", ylim=c(min(comb[,i]),max(comb[,i])), xlab = colnames(comb)[i])
}

par(mfrow=c(4,4))
for (i in 1:ncol(comb)) {
  phenogram(tree=tree_morpho_painted,x=morpho_all_gls[,i], ftype="off", ylim=c(min(comb[,i]),max(comb[,i])), xlab = colnames(comb)[i])
}


for (j in 1:4) {
  c=as.matrix(comb_clades[[j]])
  c <- c[,sort]
  t=tree_morpho_clades[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno.pdf"))
  layout(design)
  par(mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}


