# NULL with vadimi ####
#multirate BM model, each clade its own rate, models_PC1_list[[i]][[2]]

# 1 simulate
sim_traits <- list()
for (i in 1:length(models_PC1_list)) {
  obj <- models_PC1_list[[i]][[2]]
  sim <- simulate(obj, nsim=100, tree=tree_morpho_painted)
  sim_traits[[i]] <- sim
}
names(sim_traits) <- names(models_PC1_list)

# 2 reconstruct inner nodes
sim_ASR_traits <- list()
for (i in 1:length(sim_traits)) {
  sim_trait <- sim_traits[[i]]
  trait <- names(sim_traits[i])
  model <- models_PC1_list[[i]][[3]]
  sim_ASRs <- list()
  for (j in 1:ncol(sim_trait)) {
    sim_nodes <- sim_trait[,j]
    ASR<-estim(tree_morpho_painted,sim_nodes, model, asr=TRUE)
    comb <- as.data.frame(rbind(as.matrix(sim_nodes),ASR$estim))
    colnames(comb) <- trait
    comb$species <- rownames(comb)
    sim_ASRs[[j]] <- comb
  }
  sim_ASR_traits[[i]] <- sim_ASRs
}

#vzamem iz vsakega lista 1 element in jih joinam
sim_ASR_traits_df <- list()
for (i in 1:100) {
  x1 = lapply(sim_ASR_traits, function(l) l[[i]])
  df <-  Reduce(inner_join,x1)
  rownames(df) <- df$species
  df[c("species")] <- NULL
  sim_ASR_traits_df[[i]] <- df
}



# 3 calculate morphospace filling for clades
sim_comb_PC1 <- list()
for (i in 1:100) {
  traits <- sim_ASR_traits_df[[i]]
  sim_comb_PC1[[i]]<- reconstructionsCladeSubsets(trees = tree_morpho_clades,
                              translation = node_translation,
                              traits = traits,
                              descendants = desc_clades,
                              nameclades = nameclades)
}

## Morphospace through time for clades and for sims
sim_morphospaceTime <- list()
for (i in 1:100) {
  sim_comb_clade <- sim_comb_PC1[[i]]
  sim_morphospaceTime_clades_PC1 <- list()
  for (z in 1:length(sim_comb_clade)) {
    comb_clade <- sim_comb_clade[[z]]
    tree_morpho_clade <- tree_morpho_clades[[z]]
    morphospaceTime_df <- morphospaceThroughTime(tree_morpho_clade,0.5,comb_clade)
    sim_morphospaceTime_clades_PC1[[z]] <- morphospaceTime_df
  }
  names(sim_morphospaceTime_clades_PC1) <- nameclades
  sim_morphospaceTime[[i]] <- sim_morphospaceTime_clades_PC1
  
}

# 4 compare to real data
# sim_morphospaceTime has 100 lists with 4 dfs with morphospace through time  for each clade, for all traits.
# merge lists of clades in df and add column clade
sim_morphospaceTime_df <- list()
for (i in 1:100) {
  x1 = sim_morphospaceTime[[i]]
  df <-  ldply(x1)
  names(df)[names(df) == '.id'] <- 'clade'
  sim_morphospaceTime_df[[i]] <- df
}
names(sim_morphospaceTime_df) <- sprintf("seq%s",seq(1:100))

# merge lists of simulations in df and add column .id of simulation
sims_df <- ldply(sim_morphospaceTime_df)

colors3=c("brown","red","orange","darkblue","#1B96CC")

#original data
real <- ldply(morphospaceTime_clades_PC1)
names(real)[names(real) == '.id'] <- 'clade'
real$.id <- c("real")


#plots of morphospace filling per clade and per trait null and real
plots_scaled_PC1_sims_test <- list()
for (i in 1:length(nameclades)) {
  real2 <- real%>% filter(clade==nameclades[i]) %>% group_by(trait) %>% mutate(range = range01(range))
  sims_df2 = sims_df %>% filter(clade==nameclades[i])%>% group_by(.id,trait) %>% mutate(range = range01(range))
  p <- ggplot(sims_df2, aes(time, range, group = .id, colour = trait))+
    geom_line()+
    geom_line(data = real2, aes(time, range), colour = "black")+
    scale_color_manual(values=alpha(colors3,0.1))+
    theme_bw()+
    facet_wrap(~trait, scales = "fixed", nrow = 3, ncol = 2)+
    labs(title=nameclades[i])
  plots_scaled_PC1_sims_test[[i]] <- p
}
for (j in 1:length(nameclades)) {
  ggsave(filename = paste0(nameclades[j],"_null_PCA.pdf"),
         plot=plots_scaled_PC1_sims_test[[j]],
         width = 8,
         height= 8)
}

# Slopes
#k = [y(t+1)-y(t)] / 0.5 
#k = df$range[i]-df$range[i-1]/(df$time[i]-df$time[i-1])
#i = 2:nrow(df)

# real slope
real_slope <- real %>% 
    group_by(clade, trait) %>% 
    mutate(range = range01(range))%>%
    mutate(slope = map_dbl(row_number(), function(i) {
      ifelse(i==1, 0, (range[i]-range[i-1])/(time[i]-time[i-1]))
    }))

# simulated slopes
sims_slope <- sims_df %>% 
  group_by(.id,clade, trait) %>% 
  mutate(range = range01(range))%>%
  mutate(slope = map_dbl(row_number(), function(i) {
    ifelse(i==1, 0, (range[i]-range[i-1])/(time[i]-time[i-1]))
  }))

slopes <- left_join(sims_slope, real_slope %>% select(clade, trait, time, slope), by = c("clade", "trait", "time"), suffix = c("_sim","_real"))
slopes <- slopes %>% mutate(slope_dif = slope_real- slope_sim)

plots_slopes <- list()
for (i in 1:length(nameclades)) {
  slopes2 <- slopes%>% filter(clade==nameclades[i]) %>% group_by(trait)
  p <- ggplot(slopes2, aes(time, slope_dif, group = .id, colour = trait))+
    geom_line()+
    geom_hline(yintercept = 0, linetype='dotted')+
    scale_color_manual(values=alpha(colors3,0.1))+
    theme_bw()+
    facet_wrap(~trait, scales = "fixed", nrow = 3, ncol = 2)+
    labs(title=nameclades[i])
  plots_slopes[[i]] <- p
}

for (j in 1:length(nameclades)) {
  ggsave(filename = paste0(nameclades[j],"_slopes.pdf"),
         plot=plots_slopes[[j]],
         width = 8,
         height= 8)
}
