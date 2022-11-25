vadimi <- c("N_vadimi")
morpho_all_gls_nv <- morpho_all_gls[!rownames(morpho_all_gls) %in% vadimi, ]

tree_morpho_nv=drop.tip(tree_morpho, "N_vadimi")

## ...mrca ####
mrcaA1_nv <- findMRCA(tree_morpho_nv, speciesA1_m[!speciesA1_m %in% c("N_vadimi")])
mrcaA2_nv <- findMRCA(tree_morpho_nv, speciesA2_m)
mrcaB_nv <- findMRCA(tree_morpho_nv, speciesB_m)
mrcaC_nv <- findMRCA(tree_morpho_nv, speciesC_m)

tree_morpho_painted_nv <-tree_morpho_nv
is.ultrametric(tree_morpho_painted_nv)
tree_morpho_painted_nv<-paintSubTree(tree_morpho_painted_nv,node=mrcaA1_nv,state="A1",stem=T)
tree_morpho_painted_nv<-paintSubTree(tree_morpho_painted_nv,node=mrcaA2_nv,state="A2",stem=T)
tree_morpho_painted_nv<-paintSubTree(tree_morpho_painted_nv,node=mrcaB_nv,state="B",stem=T)
tree_morpho_painted_nv<-paintSubTree(tree_morpho_painted_nv,node=mrcaC_nv,state="C",stem=T)
plotSimmap(tree_morpho_painted_nv,cols,lwd=2,pts=F, ftype="off")
axisPhylo()

tree_morpho_n_nv <- tree_morpho_nv
tree_morpho_n_nv$node.label <- c((tree_morpho_n_nv$Nnode+2):((tree_morpho_n_nv$Nnode*2)+1))
#use labeled tree for subtrees
tree_morpho_A1_nv=keep.tip(tree_morpho_n_nv, speciesA1_m[!speciesA1_m %in% c("N_vadimi")])
tree_morpho_A2_nv=keep.tip(tree_morpho_n_nv, speciesA2_m)
tree_morpho_B_nv=keep.tip(tree_morpho_n_nv, speciesB_m)
tree_morpho_C_nv=keep.tip(tree_morpho_n_nv, speciesC_m)
tree_morpho_clades_nv <- list(tree_morpho_A1_nv,tree_morpho_A2_nv,tree_morpho_B_nv,tree_morpho_C_nv)
names(tree_morpho_clades_nv) <- nameclades

node_translation_nv <- nodeTranslation(tree_morpho_clades_nv)

## ...descendants ####
descA1_nv <- sort(c(getDescendants(tree_morpho_painted_nv,mrcaA1_nv),mrcaA1_nv))
descA2_nv <- sort(c(getDescendants(tree_morpho_painted_nv,mrcaA2_nv),mrcaA2_nv))
descB_nv <- sort(c(getDescendants(tree_morpho_painted_nv,mrcaB_nv),mrcaB_nv))
descC_nv <- sort(c(getDescendants(tree_morpho_painted_nv,mrcaC_nv),mrcaC_nv))
desc_clades_nv <- list(descA1_nv,descA2_nv,descB_nv,descC_nv)
names(desc_clades_nv) <- nameclades


## PCAs ####
#run PCAs
trophic_pca_nv <- prcomp(morpho_all_gls_nv[,trophic_traits], scale. = TRUE)
trophic2_pca_nv <- prcomp(morpho_all_gls_nv[,trophic_traits2], scale. = TRUE)
body_pca_nv <- prcomp(morpho_all_gls_nv[,body_traits], scale. = TRUE)
loco_pca_nv <- prcomp(morpho_all_gls_nv[,loco_traits], scale. = TRUE)
sens_pca_nv <- prcomp(morpho_all_gls_nv[,sens_traits], scale. = TRUE)

 
pc1_nv <- as.matrix(cbind(body_pca_nv$x[,1], loco_pca_nv$x[,1], sens_pca_nv$x[,1],
                          trophic_pca_nv$x[,1],trophic2_pca_nv$x[,1]))
colnames(pc1_nv) <- PC_traits

all(rownames(pc1_nv)==tree_morpho_painted_nv$tip.label)
## MODELS on PC1 ####
models_PC1_list_nv <- list()
for (i in 1:length(PC_traits)) {
  trait <- PC_traits[i]
  models_PC1 <- allModels(tree_morpho_painted_nv,pc1_nv,trait)
  models_PC1_list_nv[[i]] <- models_PC1
}
names(models_PC1_list_nv) <- names(traits)

best_model_PC1_list_nv <- list()
for (i in 1:length(models_PC1_list_nv)) {
  models <- models_PC1_list_nv[[i]]
  best_model_PC1 <- testBestModel(models)
  best_model_PC1_list_nv[[i]] <- best_model_PC1
}
names(best_model_PC1_list_nv) <- names(models_PC1_list_nv)
sink("modelling_PCA_novadimi.txt")
print(best_model_PC1_list_nv)
sink()
# Ancestral states ####
ASR_PC1_list_nv <-  list()
for (i in 1:length(PC_traits)) {
  model <- best_model_PC1_list_nv[[i]][[2]]
  trait <- PC_traits[i]
  ASR<-estim(tree_morpho_painted_nv,pc1_nv[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(as.matrix(pc1_nv[,trait]),ASR$estim))
  colnames(comb) <- trait
  comb$species <- rownames(comb)
  ASR_PC1_list_nv[[i]] <- comb
}
ASR_PC1_df_nv <-  Reduce(inner_join,ASR_PC1_list_nv)
rownames(ASR_PC1_df_nv) <- ASR_PC1_df_nv$species
ASR_PC1_df_nv[c("species")] <- NULL

# plot phenogram ####
dev.off()
cols<-c("#00000000","deepskyblue3","purple","red", "orange"); names(cols)<-c(1,"A1","A2","B","C")
par(mfrow=c(3,2),mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
for (i in 1:ncol(ASR_PC1_df_nv)) {
  phenogram(tree=tree_morpho_painted_nv,colors=cols, x=as.matrix(ASR_PC1_df_nv)[,i],ftype="off",
            ylim=c(min(as.matrix(ASR_PC1_df_nv)[,i]),max(as.matrix(ASR_PC1_df_nv)[,i])),
            xlab = colnames(ASR_PC1_df_nv)[i], lwd=0.6)
}

# clades ####

comb_clades_PC1_nv <- reconstructionsCladeSubsets(trees = tree_morpho_clades_nv,
                                               translation = node_translation_nv,
                                               traits = ASR_PC1_df_nv,
                                               descendants = desc_clades_nv,
                                               nameclades = nameclades)

## Morphospace through time ####
morphospaceTime_df_PC1_nv <- morphospaceThroughTime(tree_morpho_painted_nv,0.5,ASR_PC1_df_nv)
morphospaceTime_df_means_PC1_nv <- morphospaceTime_df_PC1_nv %>% 
  group_by(trait) %>%
  summarise(mean=max(range)/2)

ggplot(morphospaceTime_df_PC1_nv, aes(time, range))+
  geom_point()+
  geom_hline(data = morphospaceTime_df_means_PC1_nv, mapping = aes(yintercept = mean)) +
  facet_wrap(~trait, scales = "free")


morphospaceTime_clades_PC1_nv <- list()
for (z in 1:length(comb_clades_PC1_nv)) {
  comb_clade <- comb_clades_PC1_nv[[z]]
  tree_morpho_clade <- tree_morpho_clades_nv[[z]]
  morphospaceTime_df <- morphospaceThroughTime(tree_morpho_clade,0.5,comb_clade)
  morphospaceTime_clades_PC1_nv[[z]] <- morphospaceTime_df
}

morphospaceTime_means_PC1_nv <- list()
for (i in 1:length(morphospaceTime_clades_PC1_nv)) {
  morphospaceTime_df_means_PC1 <- morphospaceTime_clades_PC1_nv[[i]] %>% 
    group_by(trait) %>%
    summarise(mean=max(range)/2)
  morphospaceTime_means_PC1_nv[[i]] <- morphospaceTime_df_means_PC1
}

names(morphospaceTime_clades_PC1_nv) <- nameclades
names(morphospaceTime_means_PC1_nv) <- nameclades


plots_scaled_PC1_nv <- list()
for (i in 1:length(morphospaceTime_clades_PC1_nv)) {
  title <- names(morphospaceTime_clades_PC1_nv[i])
  df <- morphospaceTime_clades_PC1_nv[[i]]
  df <- df %>% group_by(trait) %>% mutate(range = range01(range))
  p <- ggplot(df, aes(time, range, group = trait, colour = trait))+
    geom_line()+
    scale_color_manual(values=colors3)+
    theme_bw()+
    labs(title=title)
  plots_scaled_PC1_nv[[i]] <- p
}

multiplot_scaled_PC1_nv <- ggarrange(plotlist=plots_scaled_PC1_nv,ncol=2,nrow=2,common.legend = TRUE)
multiplot_scaled_PC1_nv
ggexport(multiplot_scaled_PC1_nv, filename = "clades_morphospace_scaledPCA_nv.pdf")


for (i in 1:length(plots_PC1)) {
  ggsave(plots_PC1[[i]], filename = paste0(nameclades[i],"_morphospace_PC1_phylpca.pdf"), scale = 5)
}


for (j in 1:length(comb_clades_PC1_nv)) {
  c=as.matrix(comb_clades_PC1_nv[[j]])
  t=tree_morpho_clades_nv[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno_scaledPCA_nv.pdf"))
  par(mfrow=c(2,3), mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}


