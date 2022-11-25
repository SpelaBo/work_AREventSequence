# .morphology #### 
morphology_individuals <- read_excel("./data/20221018_morphology_imput_complete.xlsx", na = "NA")
morphology_individuals_top3 <- morphology_individuals %>%# Top N highest values by group
  arrange(desc(body_length)) %>% 
  group_by(species) %>%
  slice(1:3)

morpho_AT3 <- morphology_individuals_top3 %>%# average
  arrange(desc(body_length)) %>% 
  group_by(species) %>%
  summarise(across(-MOTU_ind, mean, na.rm = TRUE))


morpho_AT3 <- as.data.frame(morpho_AT3)
row.names(morpho_AT3) <- morpho_AT3$species
morpho_AT3[1] <- NULL
## calculate gnathopod size and remove unnecessary columns
morpho_AT3$gpI6_size <- morpho_AT3$gI6_length+morpho_AT3$gI6_width+morpho_AT3$gI6_diag
morpho_AT3$gpII6_size <- morpho_AT3$gII6_length+morpho_AT3$gII6_width+morpho_AT3$gII6_diag

morpho_AT3$Icosinus <- ((morpho_AT3$gI6_length)^2+(morpho_AT3$gI6_width)^2-(morpho_AT3$gI6_diag)^2)/(2*morpho_AT3$gI6_length*morpho_AT3$gI6_width)
morpho_AT3$IIcosinus <- ((morpho_AT3$gII6_length)^2+(morpho_AT3$gII6_width)^2-(morpho_AT3$gII6_diag)^2)/(2*morpho_AT3$gII6_length*morpho_AT3$gII6_width)
morpho_AT3$Ialpha <- 180*(acos(morpho_AT3$Icosinus))/pi
morpho_AT3$IIalpha <- 180*(acos(morpho_AT3$IIcosinus))/pi
morpho_AT3 <- morpho_AT3[,c(1:11,20:21,24:25)]
morpho_AT3 <- as.matrix(morpho_AT3[ order(row.names(morpho_AT3)), ])
## Remove NA-s
morpho_AT3_noNA <- morpho_AT3[complete.cases(morpho_AT3), ]
morpho_AT3_noNA<-morpho_AT3_noNA[order(match(rownames(morpho_AT3_noNA), tree_morpho$tip.label)), , drop = FALSE]


# Phylogenetically corrected GLS
phyl_gls_morpho_avg_top3 <- phyl.resid(tree_morpho, morpho_AT3_noNA[,1], morpho_AT3_noNA[,c(2:13)])
residuals_gls_avg_top3 <- as.matrix(phyl_gls_morpho_avg_top3$resid)
body_length_avg_top3 <- as.matrix(morpho_AT3_noNA[,1], )
gnato_angles_avg_top3 <- as.matrix(morpho_AT3_noNA[,14:15], )
colnames(body_length_avg_top3)<-c("body_length")
morpho_all_gls_AT3 <- cbind(body_length_avg_top3, residuals_gls_avg_top3, gnato_angles_avg_top3)

all(rownames(morpho_all_gls_AT3)==tree_morpho$tip.label)
name.check(tree_morpho,morpho_all_gls_AT3)

### Analysis of diversification, reconstruction on whole tree ####
# Modelling, genus level ####
models_list_AT3 <- list()
for (i in 1:length(traits)) {
  trait <- traits[[i]]
  models <- allModels(tree_morpho_painted,morpho_all_gls_AT3,trait)
  models_list_AT3[[i]] <- models
}
names(models_list_AT3) <- names(traits)

best_model_list_AT3 <- list()
for (i in 1:length(models_list)) {
  models <- models_list_AT3[[i]]
  best_model <- testBestModel(models)
  best_model_list_AT3[[i]] <- best_model
}
names(best_model_list_AT3) <- names(models_list_AT3)

# ancestral states ####
ASR_list_AT3 <-  list()
for (i in 1:length(traits)) {
  model <- best_model_list_AT3[[i]]
  trait <- traits[[i]]
  ASR<-estim(tree_morpho_painted,morpho_all_gls_AT3[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(morpho_all_gls_AT3[,trait],ASR$estim))
  comb$species <- rownames(comb)
  ASR_list_AT3[[i]] <- comb
}
ASR_df_AT3 <-  Reduce(inner_join,ASR_list_AT3)
rownames(ASR_df_AT3) <- ASR_df_AT3$species
ASR_df_AT3[c("species")] <- NULL

# ASR for each clade ####
ASR_clades_AT3 <- reconstructionsCladeSubsets(trees = tree_morpho_clades,
                                          translation = node_translation,
                                          traits = ASR_df_AT3, 
                                          descendants = desc_clades,
                                          nameclades = nameclades)

# morpho through time clade ####
morphospaceTime_clades_AT3 <- list()
for (z in 1:length(ASR_clades_AT3)) {
  ASR_clade <- ASR_clades_AT3[[z]]
  tree_morpho_clade <- tree_morpho_clades[[z]]
  morphospaceTime_df <- morphospaceThroughTime(tree_morpho_clade,0.1,ASR_clade)
  morphospaceTime_clades_AT3[[z]] <- morphospaceTime_df
}

morphospaceTime_means_AT3 <- list()
for (i in 1:length(morphospaceTime_clades_AT3)) {
  morphospaceTime_df_means <- morphospaceTime_clades_AT3[[i]] %>% 
    group_by(trait) %>%
    summarise(mean=max(range)/2)
  morphospaceTime_means_AT3[[i]] <- morphospaceTime_df_means
}

names(morphospaceTime_clades_AT3) <- nameclades
names(morphospaceTime_means_AT3) <- nameclades

# plot ####
design <- matrix(c(1,2,3,4,5,6,7,8,0,9,10,0,11,12,0), nrow=3, byrow = F)
for (j in 1:length(nameclades)) {
  c=as.matrix(ASR_clades_AT3[[j]])
  #c <- c[,sort]
  t=tree_morpho_clades[[j]]
  title <- nameclades[j]
  #pdf(file = paste0(title,"_pheno_newAT3.pdf"))
  layout(design)
  par(mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  #title(title)
  dev.off()
}





plots_scaled_AT3 <- list()
for (i in 1:length(morphospaceTime_clades_AT3)) {
  title <- names(morphospaceTime_clades_AT3[i])
  df <- morphospaceTime_clades_AT3[[i]]
  df <- df %>% group_by(trait) %>% mutate(range = range01(range))
  #df[c(2)] <- scale(df[c(2)])
  means <- morphospaceTime_means_AT3[[i]]
  p <- ggplot(transform(df, trait = factor(trait, levels = sort)), aes(time, range, group = trait, colour = trait))+
    geom_line()+
    scale_color_manual(values=colors2)+
    labs(title=title)
  plots_scaled_AT3[[i]] <- p
}
plots_scaled_AT3
for (i in 1:length(plots_scaled_AT3)) {
  ggsave(plots_scaled_AT3[[i]], filename = paste0(nameclades[i],"_morphospace_scaled_AT3.pdf"), scale = 5)
}


