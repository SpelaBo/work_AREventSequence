## PCAs ####
#run PCAs
trophic_pca <- phyl.pca(tree = tree_morpho, morpho_all_gls[,trophic_traits], method = "BM", mode = "cov")
body_pca <- phyl.pca(tree = tree_morpho, morpho_all_gls[,body_traits], method = "BM", mode = "cov")
loco_pca <- phyl.pca(tree = tree_morpho, morpho_all_gls[,loco_traits], method = "BM", mode = "cov")
sens_pca <- phyl.pca(tree = tree_morpho, morpho_all_gls[,sens_traits], method = "BM", mode = "cov")

pc1 <- as.data.frame(cbind(body_pca$S[,1], loco_pca$S[,1], sens_pca$S[,1],trophic_pca$S[,1]))
colnames(pc1) <- c("PC1_body", "PC1_loco","PC1_sens","PC1_trophic")

# trophic_pca <- prcomp(morpho_all_gls[,trophic_traits])
# body_pca <- prcomp(morpho_all_gls[,body_traits])
# loco_pca <- prcomp(morpho_all_gls[,loco_traits])
# sens_pca <- prcomp(morpho_all_gls[,sens_traits])
# 
# pc1 <- as.data.frame(cbind(body_pca$x[,1], loco_pca$x[,1], sens_pca$x[,1],trophic_pca$x[,1]))
# colnames(pc1) <- c("PC1_body", "PC1_loco","PC1_sens","PC1_trophic")


## MODELS on PC1 ####
# .trophic niche ####
model_BM_gnato_PC1 <- mvBM(tree_morpho,pc1[,"PC1_trophic", drop = FALSE])
model_BMM_gnato_PC1 <- mvBM(tree_morpho_painted,pc1[,"PC1_trophic", drop = FALSE],model = "BMM", param = list(root = F))
model_EB_gnato_PC1 <- mvEB(tree_morpho,pc1[,"PC1_trophic", drop = FALSE])
model_OU_gnato_PC1 <- mvOU(tree_morpho,pc1[,"PC1_trophic", drop = FALSE], param = list(root = F))
model_OUM_gnato_PC1 <- mvOU(tree_morpho_painted,pc1[,"PC1_trophic", drop = FALSE],model = "OUM", param = list(root = F))
# merge into list
models_gnato_PC1 <- list(model_BM_gnato_PC1,model_BMM_gnato_PC1,model_EB_gnato_PC1,model_OU_gnato_PC1, model_OUM_gnato_PC1)
gnato_aic_PC1 <- aicw(models_gnato_PC1)
best_gnato_PC1 <- models_gnato_PC1[[which.max(gnato_aic_PC1$aicweights)]]


# .body shape ####
model_BM_body_PC1 <- mvBM(tree_morpho,pc1[,"PC1_body", drop = FALSE])
model_BMM_body_PC1 <- mvBM(tree_morpho_painted,pc1[,"PC1_body", drop = FALSE],model = "BMM", param = list(root = F))
model_EB_body_PC1 <- mvEB(tree_morpho,pc1[,"PC1_body", drop = FALSE], method = "sparse")
model_OU_body_PC1 <- mvOU(tree_morpho,pc1[,"PC1_body", drop = FALSE], param = list(root = F))
model_OUM_body_PC1 <- mvOU(tree_morpho_painted,pc1[,"PC1_body", drop = FALSE],model = "OUM", param = list(root = F))
# merge into list
models_body_PC1 <- list(model_BM_body_PC1,model_BMM_body_PC1,model_EB_body_PC1,model_OU_body_PC1,model_OUM_body_PC1)
body_aic_PC1 <- aicw(models_body_PC1)
best_body_PC1 <- models_body_PC1[[which.max(body_aic_PC1$aicweights)]]


# .locomotion ####
model_BM_loco_PC1 <- mvBM(tree_morpho,pc1[,"PC1_loco", drop = FALSE])
model_BMM_loco_PC1 <- mvBM(tree_morpho_painted,pc1[,"PC1_loco", drop = FALSE],model = "BMM", param = list(root = F))
model_EB_loco_PC1 <- mvEB(tree_morpho,pc1[,"PC1_loco", drop = FALSE])
model_OU_loco_PC1 <- mvOU(tree_morpho,pc1[,"PC1_loco", drop = FALSE], param = list(root = F))
model_OUM_loco_PC1 <- mvOU(tree_morpho_painted,pc1[,"PC1_loco", drop = FALSE],model = "OUM", param = list(root = F))

# merge into list
models_loco_PC1 <- list(model_BM_loco_PC1,model_BMM_loco_PC1,model_EB_loco_PC1,model_OU_loco_PC1,model_OUM_loco_PC1)
loco_aic_PC1 <- aicw(models_loco_PC1)
best_loco_PC1 <- models_loco_PC1[[which.max(loco_aic_PC1$aicweights)]]


# .sensoric ####
model_BM_sens_PC1 <- mvBM(tree_morpho,pc1[,"PC1_sens", drop = FALSE])
model_BMM_sens_PC1 <- mvBM(tree_morpho_painted,pc1[,"PC1_sens", drop = FALSE],model = "BMM")
model_EB_sens_PC1 <- mvEB(tree_morpho,pc1[,"PC1_sens", drop = FALSE])
model_OU_sens_PC1 <- mvOU(tree_morpho,pc1[,"PC1_sens", drop = FALSE], param = list(root = F))
model_OUM_sens_PC1 <- mvOU(tree_morpho_painted,pc1[,"PC1_sens", drop = FALSE],model = "OUM", param = list(root = F))

# merge into list
models_sens_PC1 <- list(model_BM_sens_PC1,model_BMM_sens_PC1,model_EB_sens_PC1,model_OU_sens_PC1,model_OUM_sens_PC1)
sens_aic_PC1 <- aicw(models_sens_PC1)
best_sens_PC1 <- models_sens_PC1[[which.max(sens_aic_PC1$aicweights)]]


# Ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato_PC1<-estim(tree_morpho_painted,pc1[,"PC1_trophic", drop = FALSE], best_gnato_PC1, asr=TRUE)
comb_gnato_PC1 <- rbind(pc1[,"PC1_trophic", drop = FALSE],ASR_gnato_PC1$estim)

ASR_body_PC1<-estim(tree_morpho_painted,pc1[,"PC1_body", drop = FALSE], best_body_PC1, asr=TRUE)
comb_body_PC1 <- rbind(pc1[,"PC1_body", drop = FALSE],ASR_body_PC1$estim)

ASR_loco_PC1<-estim(tree_morpho_painted,pc1[,"PC1_loco", drop = FALSE], best_loco_PC1, asr=TRUE)
comb_loco_PC1 <- rbind(pc1[,"PC1_loco", drop = FALSE],ASR_loco_PC1$estim)

ASR_sens_PC1<-estim(tree_morpho_painted,pc1[,"PC1_sens", drop = FALSE], best_sens_PC1, asr=TRUE)
comb_sens_PC1 <- rbind(pc1[,"PC1_sens", drop = FALSE],ASR_sens_PC1$estim)

comb_PC1 <- cbind(comb_body_PC1,comb_loco_PC1,comb_sens_PC1,comb_gnato_PC1)

comb_clades_PC1 <- reconstructionsCladeSubsets(trees = tree_morpho_clades,translation = node_translation,traits = comb_PC1, descendants = desc_clades, nameclades = nameclades)

## Morphospace through time ####
morphospaceTime_df_PC1 <- morphospaceThroughTime(tree_morpho,0.5,comb_PC1)
morphospaceTime_df_means_PC1 <- morphospaceTime_df_PC1 %>% 
  group_by(trait) %>%
  summarise(mean=max(range)/2)

ggplot(morphospaceTime_df_PC1, aes(time, range))+
  geom_point()+
  geom_hline(data = morphospaceTime_df_means_PC1, mapping = aes(yintercept = mean)) +
  facet_wrap(~trait, scales = "free")


morphospaceTime_clades_PC1 <- list()
for (z in 1:length(comb_clades_PC1)) {
  comb_clade <- comb_clades_PC1[[z]]
  tree_morpho_clade <- tree_morpho_clades[[z]]
  morphospaceTime_df <- morphospaceThroughTime(tree_morpho_clade,0.5,comb_clade)
  morphospaceTime_clades_PC1[[z]] <- morphospaceTime_df
}

morphospaceTime_means_PC1 <- list()
for (i in 1:length(morphospaceTime_clades_PC1)) {
  morphospaceTime_df_means_PC1 <- morphospaceTime_clades_PC1[[i]] %>% 
    group_by(trait) %>%
    summarise(mean=max(range)/2)
  morphospaceTime_means_PC1[[i]] <- morphospaceTime_df_means_PC1
}

names(morphospaceTime_clades_PC1) <- nameclades
names(morphospaceTime_means_PC1) <- nameclades

plots_PC1 <- list()
for (i in 1:length(morphospaceTime_clades_PC1)) {
  title <- names(morphospaceTime_clades_PC1[i])
  df <- morphospaceTime_clades_PC1[[i]]
  means <- morphospaceTime_means_PC1[[i]]
  p <- ggplot(df, aes(time, range))+
    geom_point()+
    geom_hline(data = means, mapping = aes(yintercept = mean)) +
    labs(title=title)+
    facet_wrap(trait~., scales = "free")
  plots_PC1[[i]] <- p
}
plots_PC1

for (i in 1:length(plots_PC1)) {
  ggsave(plots_PC1[[i]], filename = paste0(nameclades[i],"_morphospace_PC1_phylpca.pdf"), scale = 5)
}


for (j in 1:length(comb_clades_PC1)) {
  c=as.matrix(comb_clades_PC1[[j]])
  t=tree_morpho_clades[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno_PC1_phylpca.pdf"))
  par(mfrow=c(2,2), mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}


