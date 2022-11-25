### Analysis of diversification, reconstruction on whole tree ####
# Modelling, genus level ####
models_list_nv <- list()
for (i in 1:length(traits)) {
  trait <- traits[[i]]
  models <- allModels(tree_morpho_painted_nv,morpho_all_gls_nv,trait)
  models_list_nv[[i]] <- models
}
names(models_list_nv) <- names(traits)

best_model_list_nv <- list()
for (i in 1:length(models_list_nv)) {
  models <- models_list_nv[[i]]
  best_model <- testBestModel(models)
  best_model_list_nv[[i]] <- best_model
}
names(best_model_list_nv) <- names(models_list_nv)
sink("modelling_GLS_novadimi.txt")
print(best_model_list_nv)
sink()

# ancestral states ####
ASR_list_nv <-  list()
for (i in 1:length(traits)) {
  model <- best_model_list_nv[[i]][[2]]
  trait <- traits[[i]]
  ASR<-estim(tree_morpho_painted_nv,morpho_all_gls_nv[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(morpho_all_gls_nv[,trait],ASR$estim))
  comb$species <- rownames(comb)
  ASR_list_nv[[i]] <- comb
}
ASR_df_nv <-  Reduce(inner_join,ASR_list_nv)
rownames(ASR_df_nv) <- ASR_df_nv$species
ASR_df_nv[c("species")] <- NULL

# plot phenogram ####
dev.off()
par(mfrow=c(3,4),mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
for (i in 1:ncol(ASR_df_nv)) {
  phenogram(tree=tree_morpho_painted_nv,colors=cols, x=as.matrix(ASR_df_nv)[,i], ftype="off", ylim=c(min(as.matrix(ASR_df_nv)[,i]),max(as.matrix(ASR_df_nv)[,i])), xlab = colnames(ASR_df_nv)[i], lwd=0.6)
}
# clade subsetting ####
# node_translation_nv and desc_clades_nv defined in PCA no vadimi file

# ASR for each clade ####
ASR_clades_nv <- reconstructionsCladeSubsets(trees = tree_morpho_clades_nv,
                                                 translation = node_translation_nv,
                                                 traits = ASR_df_nv, 
                                                 descendants = desc_clades_nv,
                                                 nameclades = nameclades)

# plot phenograms ####
design <- matrix(c(1,2,3,4,5,6,7,8,0,9,10,0,11,12,0), nrow=3, byrow = F)
for (j in 1:length(nameclades)) {
  c=as.matrix(ASR_clades_nv[[j]])
  #c <- c[,sort]
  t=tree_morpho_clades_nv[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno_nv.pdf"))
  layout(design)
  par(mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}

# morphospace through time: genus ####
morphospaceTime_df_nv <- morphospaceThroughTime(tree_morpho_nv,0.5,ASR_df_nv)
morphospaceTime_df_means_nv <- morphospaceTime_df_nv %>% 
  group_by(trait) %>%
  summarise(mean=max(range)/2)

ggplot(morphospaceTime_df_nv, aes(time, range))+
  geom_point()+
  geom_hline(data = morphospaceTime_df_means_nv, mapping = aes(yintercept = mean)) +
  facet_wrap(~trait, scales = "free")


# morphospace through time: clades ####
morphospaceTime_clades_nv <- list()
for (z in 1:length(ASR_clades_nv)) {
  ASR_clade <- ASR_clades_nv[[z]]
  tree_morpho_clade <- tree_morpho_clades_nv[[z]]
  morphospaceTime_df <- morphospaceThroughTime(tree_morpho_clade,0.5,ASR_clade)
  morphospaceTime_clades_nv[[z]] <- morphospaceTime_df
}

morphospaceTime_means_nv <- list()
for (i in 1:length(morphospaceTime_clades_nv)) {
  morphospaceTime_df_means <- morphospaceTime_clades_nv[[i]] %>% 
    group_by(trait) %>%
    summarise(mean=max(range)/2)
  morphospaceTime_means_nv[[i]] <- morphospaceTime_df_means
}

names(morphospaceTime_clades_nv) <- nameclades
names(morphospaceTime_means_nv) <- nameclades

# plot morphospace ####
# facet manual plots are in old script (just morphospace)
sort <- c(body_traits,loco_traits,sens_traits,trophic_traits,trophic_traits2)
design2 <- matrix(c(1:8,NA,9,10,NA,11,12,NA), nrow=3, byrow = F)

colors=c("#234d20","#77ab59","#c9df8a","#ffff66","#ffcc66","#ff6600","#ff00ff","#cc66ff","#663300","#996600","#0099ff","#0000ff")
colors2=c("brown","brown","brown","red","red","red","orange","orange","darkblue","darkblue","#1B96CC","#1B96CC")

plots_scaled_nv <- list()
for (i in 1:length(morphospaceTime_clades_nv)) {
  title <- names(morphospaceTime_clades_nv[i])
  df <- morphospaceTime_clades_nv[[i]]
  df <- df %>% group_by(trait) %>% mutate(range = range01(range))
  means <- morphospaceTime_means_nv[[i]]
  p <- ggplot(transform(df, trait = factor(trait, levels = sort)), aes(time, range, group = trait, colour = trait))+
    geom_line()+
    scale_color_manual(values=colors2)+
    theme_bw()+
    labs(title=title)
  plots_scaled_nv[[i]] <- p
}
multiplot_scaled_nv <- ggarrange(plotlist=plots_scaled_nv,ncol=2,nrow=2,common.legend = TRUE)
ggexport(multiplot_scaled_nv, filename = "clades_morphospace_scaled_nv.pdf")

