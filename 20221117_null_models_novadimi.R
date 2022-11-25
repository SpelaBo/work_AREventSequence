# NULL no vadimi ####
#multirate BM model, each clade its own rate, models_PC1_list[[i]][[2]]

# 1 simulate
sim_traits_nv <- list()
for (i in 1:length(models_PC1_list_nv)) {
  obj <- models_PC1_list_nv[[i]][[2]]
  sim <- simulate(obj, nsim=100, tree=tree_morpho_painted_nv)
  sim_traits_nv[[i]] <- sim
}
names(sim_traits_nv) <- names(models_PC1_list_nv)

# 2 reconstruct inner nodes
sim_ASR_traits_nv <- list()
for (i in 1:length(sim_traits_nv)) {
  sim_trait_nv <- sim_traits_nv[[i]]
  trait_nv <- names(sim_traits_nv[i])
  model_nv <- models_PC1_list_nv[[i]][[3]]
  sim_ASRs_nv <- list()
  for (j in 1:ncol(sim_trait_nv)) {
    sim_nodes_nv <- sim_trait_nv[,j]
    ASR_nv<-estim(tree_morpho_painted_nv,sim_nodes_nv, model_nv, asr=TRUE)
    comb_nv <- as.data.frame(rbind(as.matrix(sim_nodes_nv),ASR_nv$estim))
    colnames(comb_nv) <- trait_nv
    comb_nv$species <- rownames(comb_nv)
    sim_ASRs_nv[[j]] <- comb_nv
  }
  sim_ASR_traits_nv[[i]] <- sim_ASRs_nv
}

#vzamem iz vsakega lista 1 element in jih joinam
sim_ASR_traits_df_nv <- list()
for (i in 1:100) {
  x1_nv = lapply(sim_ASR_traits_nv, function(l) l[[i]])
  df_nv <-  Reduce(inner_join,x1_nv)
  rownames(df_nv) <- df_nv$species
  df_nv[c("species")] <- NULL
  sim_ASR_traits_df_nv[[i]] <- df_nv
}



# 3 calculate morphospace filling for clades
sim_comb_PC1_nv <- list()
for (i in 1:100) {
  traits_nv <- sim_ASR_traits_df_nv[[i]]
  sim_comb_PC1_nv[[i]]<- reconstructionsCladeSubsets(trees = tree_morpho_clades_nv,
                              translation = node_translation_nv,
                              traits = traits_nv,
                              descendants = desc_clades_nv,
                              nameclades = nameclades)
}

## Morphospace through time for clades and for sims
sim_morphospaceTime_nv <- list()
for (i in 1:100) {
  sim_comb_clade_nv <- sim_comb_PC1_nv[[i]]
  sim_morphospaceTime_clades_PC1_nv <- list()
  for (z in 1:length(sim_comb_clade_nv)) {
    comb_clade_nv <- sim_comb_clade_nv[[z]]
    tree_morpho_clade_nv <- tree_morpho_clades_nv[[z]]
    morphospaceTime_df_nv <- morphospaceThroughTime(tree_morpho_clade_nv,0.5,comb_clade_nv)
    sim_morphospaceTime_clades_PC1_nv[[z]] <- morphospaceTime_df_nv
  }
  names(sim_morphospaceTime_clades_PC1_nv) <- nameclades
  sim_morphospaceTime_nv[[i]] <- sim_morphospaceTime_clades_PC1_nv
  
}

# 4 compare to real data
# sim_morphospaceTime has 100 lists with 4 dfs with morphospace through time  for each clade, for all traits.
# merge lists of clades in df and add column clade
sim_morphospaceTime_df_nv <- list()
for (i in 1:100) {
  x1_nv = sim_morphospaceTime_nv[[i]]
  df_nv <-  ldply(x1_nv)
  names(df_nv)[names(df_nv) == '.id'] <- 'clade'
  sim_morphospaceTime_df_nv[[i]] <- df_nv
}
names(sim_morphospaceTime_df_nv) <- sprintf("seq%s",seq(1:100))

# merge lists of simulations in df and add column .id of simulation
sims_df_nv <- ldply(sim_morphospaceTime_df_nv)

#original data
real_nv <- ldply(morphospaceTime_clades_PC1_nv)
names(real_nv)[names(real_nv) == '.id'] <- 'clade'
real_nv$.id <- c("real")

colors3=c("brown","red","orange","darkblue","#1B96CC")
#plots of morphospace filling per clade and per trait null and real
plots_scaled_PC1_sims_nv <- list()
for (i in 1:length(nameclades)) {
  real_nv2 <- real_nv%>% filter(clade==nameclades[i]) %>% group_by(trait) %>% mutate(range = range01(range))
  sims_df_nv2 <-  sims_df_nv %>% filter(clade==nameclades[i])%>% group_by(.id,trait) %>% mutate(range = range01(range))
  p <- ggplot(sims_df_nv2, aes(time, range, group = .id, colour = trait))+
    geom_line()+
    geom_line(data = real_nv2, aes(time, range), colour = "black")+
    scale_color_manual(values=alpha(colors3,0.1))+
    theme_bw()+
    facet_wrap(~trait, scales = "fixed", nrow = 3, ncol = 2)+
    labs(title=nameclades[i])
  plots_scaled_PC1_sims_nv[[i]] <- p
}

multiplot_null_PC1 <- ggarrange(plotlist=plots_scaled_PC1_sims_nv,ncol=2,nrow=2,common.legend = TRUE)
ggexport(multiplot_null_PC1, filename = "clades_morphospace_scaledPCA.pdf")

for (j in 1:length(nameclades)) {
  ggsave(filename = paste0(nameclades[j],"_null_PCA_novadimi.pdf"),
         plot=plots_scaled_PC1_sims_nv[[j]],
         width = 8,
         height= 8)
}

plots_scaled_PC1_sims_test

