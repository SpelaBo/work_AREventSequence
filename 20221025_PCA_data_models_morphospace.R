## PCAs ####
#run PCAs
trophic_pca <- prcomp(morpho_all_gls[,trophic_traits], scale. = TRUE)
trophic2_pca <- prcomp(morpho_all_gls[,trophic_traits2], scale. = TRUE)
body_pca <- prcomp(morpho_all_gls[,body_traits], scale. = TRUE)
loco_pca <- prcomp(morpho_all_gls[,loco_traits], scale. = TRUE)
sens_pca <- prcomp(morpho_all_gls[,sens_traits], scale. = TRUE)

summary(trophic_pca)
summary(trophic2_pca)
summary(body_pca)
summary(loco_pca)
summary(sens_pca)
 
pc1 <- as.matrix(cbind(body_pca$x[,1], loco_pca$x[,1], sens_pca$x[,1],trophic_pca$x[,1],trophic2_pca$x[,1]))
PC_traits <- names(traits)
 # c("PC1_body", "PC1_loco","PC1_sens","PC1_trophic","PC1_trophic2")
colnames(pc1) <- names(traits)

all(rownames(pc1)==tree_morpho_painted$tip.label)
## MODELS on PC1 ####
models_PC1_list <- list()
for (i in 1:length(PC_traits)) {
  trait <- PC_traits[i]
  models_PC1 <- allModels(tree_morpho_painted,pc1,trait)
  models_PC1_list[[i]] <- models_PC1
}
names(models_PC1_list) <- names(traits)

best_model_PC1_list <- list()
for (i in 1:length(models_PC1_list)) {
  models <- models_PC1_list[[i]]
  best_model_PC1 <- testBestModel(models)
  best_model_PC1_list[[i]] <- best_model_PC1
}
names(best_model_PC1_list) <- names(models_PC1_list)
sink("modelling_PCA.txt")
print(best_model_PC1_list)
sink()
# Ancestral states ####
ASR_PC1_list <-  list()
for (i in 1:length(PC_traits)) {
  model <- best_model_PC1_list[[i]][[2]]
  trait <- PC_traits[i]
  ASR<-estim(tree_morpho_painted,pc1[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(as.matrix(pc1[,trait]),ASR$estim))
  colnames(comb) <- trait
  comb$species <- rownames(comb)
  ASR_PC1_list[[i]] <- comb
}
ASR_PC1_df <-  Reduce(inner_join,ASR_PC1_list)
rownames(ASR_PC1_df) <- ASR_PC1_df$species
ASR_PC1_df[c("species")] <- NULL

# plot phenogram ####
dev.off()
cols<-c("#00000000","deepskyblue3","purple","red", "orange"); names(cols)<-c(1,"A1","A2","B","C")
par(mfrow=c(3,2),mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
for (i in 1:ncol(ASR_PC1_df)) {
  phenogram(tree=tree_morpho_painted,colors=cols, x=as.matrix(ASR_PC1_df)[,i],ftype="off",
            ylim=c(min(as.matrix(ASR_PC1_df)[,i]),max(as.matrix(ASR_PC1_df)[,i])),
            xlab = colnames(ASR_PC1_df)[i], lwd=0.6)
}

#clades ####
comb_clades_PC1 <- reconstructionsCladeSubsets(trees = tree_morpho_clades,
                                               translation = node_translation,
                                               traits = ASR_PC1_df,
                                               descendants = desc_clades,
                                               nameclades = nameclades)

## Morphospace through time ####
morphospaceTime_df_PC1 <- morphospaceThroughTime(tree_morpho,0.5,ASR_PC1_df)
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

colors3=c("brown","red","orange","darkblue","#1B96CC")

plots_scaled_PC1 <- list()
for (i in 1:length(morphospaceTime_clades_PC1)) {
  title <- names(morphospaceTime_clades_PC1[i])
  df <- morphospaceTime_clades_PC1[[i]]
  df <- df %>% group_by(trait) %>% mutate(range = range01(range))
  #df[c(2)] <- scale(df[c(2)])
  p <- ggplot(df, aes(time, range, group = trait, colour = trait))+
    geom_line()+
    scale_color_manual(values=colors3)+
    theme_bw()+
    labs(title=title)
  plots_scaled_PC1[[i]] <- p
}

multiplot_scaled_PC1 <- ggarrange(plotlist=plots_scaled_PC1,ncol=2,nrow=2,common.legend = TRUE)
ggexport(multiplot_scaled_PC1, filename = "clades_morphospace_scaledPCA.pdf")


for (j in 1:length(comb_clades_PC1)) {
  c=as.matrix(comb_clades_PC1[[j]])
  t=tree_morpho_clades[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno_scaledPCA.pdf"))
  par(mfrow=c(2,3), mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}


