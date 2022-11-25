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

# ancestral states ####
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
node_translation <- nodeTranslation(tree_morpho_clades)

# ASR for each clade ####
ASR_clades <- reconstructionsCladeSubsets(trees = tree_morpho_clades,
                                                 translation = node_translation,
                                                 traits = ASR_df, 
                                                 descendants = desc_clades,
                                                 nameclades = nameclades)

# plot ####
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

