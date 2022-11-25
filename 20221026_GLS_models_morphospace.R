### Analysis of diversification, reconstruction on whole tree ####
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
sink("modelling_GLS.txt")
print(best_model_list)
sink()

# ancestral states ####
ASR_list <-  list()
for (i in 1:length(traits)) {
  model <- best_model_list[[i]][[2]]
  trait <- traits[[i]]
  ASR<-estim(tree_morpho_painted,morpho_all_gls[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(morpho_all_gls[,trait],ASR$estim))
  comb$species <- rownames(comb)
  ASR_list[[i]] <- comb
}
ASR_df <-  Reduce(inner_join,ASR_list)
rownames(ASR_df) <- ASR_df$species
ASR_df[c("species")] <- NULL

# plot phenogram ####
dev.off()
par(mfrow=c(3,4),mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
for (i in 1:ncol(ASR_df)) {
  phenogram(tree=tree_morpho_painted,colors=cols, x=as.matrix(ASR_df)[,i], ftype="off", ylim=c(min(as.matrix(ASR_df)[,i]),max(as.matrix(ASR_df)[,i])), xlab = colnames(ASR_df)[i], lwd=0.6)
}

# clade subsetting ####
node_translation <- nodeTranslation(tree_morpho_clades)

# ASR for each clade ####
ASR_clades <- reconstructionsCladeSubsets(trees = tree_morpho_clades,
                                                 translation = node_translation,
                                                 traits = ASR_df, 
                                                 descendants = desc_clades,
                                                 nameclades = nameclades)

# plot phenograms ####
design <- matrix(c(1,2,3,4,5,6,7,8,0,9,10,0,11,12,0), nrow=3, byrow = F)
for (j in 1:length(nameclades)) {
  c=as.matrix(ASR_clades[[j]])
  #c <- c[,sort]
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

# morphospace through time: genus ####
morphospaceTime_df1 <- morphospaceThroughTime(tree_morpho,0.5,comb)
morphospaceTime_df_means <- morphospaceTime_df %>% 
  group_by(trait) %>%
  summarise(mean=max(range)/2)

ggplot(morphospaceTime_df, aes(time, range))+
  geom_point()+
  geom_hline(data = morphospaceTime_df_means, mapping = aes(yintercept = mean)) +
  facet_wrap(~trait, scales = "free")


# morphospace through time: clades ####
morphospaceTime_clades <- list()
for (z in 1:length(ASR_clades)) {
  ASR_clade <- ASR_clades[[z]]
  tree_morpho_clade <- tree_morpho_clades[[z]]
  morphospaceTime_df <- morphospaceThroughTime(tree_morpho_clade,0.5,ASR_clade)
  morphospaceTime_clades[[z]] <- morphospaceTime_df
}

morphospaceTime_means <- list()
for (i in 1:length(morphospaceTime_clades)) {
  morphospaceTime_df_means <- morphospaceTime_clades[[i]] %>% 
    group_by(trait) %>%
    summarise(mean=max(range)/2)
  morphospaceTime_means[[i]] <- morphospaceTime_df_means
}

names(morphospaceTime_clades) <- nameclades
names(morphospaceTime_means) <- nameclades

# plot morphospace ####
# facet manual plots are in old script (just morphospace)
sort <- c(body_traits,loco_traits,sens_traits,trophic_traits,trophic_traits2)
design2 <- matrix(c(1:8,NA,9,10,NA,11,12,NA), nrow=3, byrow = F)

colors=c("#234d20","#77ab59","#c9df8a","#ffff66","#ffcc66","#ff6600","#ff00ff","#cc66ff","#663300","#996600","#0099ff","#0000ff")
colors2=c("brown","brown","brown","red","red","red","orange","orange","darkblue","darkblue","#1B96CC","#1B96CC")

plots_scaled <- list()
for (i in 1:length(morphospaceTime_clades)) {
  title <- names(morphospaceTime_clades[i])
  df <- morphospaceTime_clades[[i]]
  df <- df %>% group_by(trait) %>% mutate(range = range01(range))
  means <- morphospaceTime_means[[i]]
  p <- ggplot(transform(df, trait = factor(trait, levels = sort)), aes(time, range, group = trait, colour = trait))+
    geom_line()+
    scale_color_manual(values=colors2)+
    theme_bw()+
    labs(title=title)
  plots_scaled[[i]] <- p
}
multiplot_scaled <- ggarrange(plotlist=plots_scaled,ncol=2,nrow=2,common.legend = TRUE)
ggexport(multiplot_scaled, filename = "clades_morphospace_scaled.pdf")

