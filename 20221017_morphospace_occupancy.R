# whole genus ####
morphospaceTime_df1 <- morphospaceThroughTime(tree_morpho,0.5,comb)
morphospaceTime_df_means <- morphospaceTime_df %>% 
  group_by(trait) %>%
  summarise(mean=max(range)/2)

ggplot(morphospaceTime_df, aes(time, range))+
  geom_point()+
  geom_hline(data = morphospaceTime_df_means, mapping = aes(yintercept = mean)) +
  facet_wrap(~trait, scales = "free")


# clades ####
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

sort <- c(body_traits,loco_traits,sens_traits,trophic_traits,trophic_traits2)

design2 <- matrix(c(1:8,NA,9,10,NA,11,12,NA), nrow=3, byrow = F)

plots <- list()
for (i in 1:length(morphospaceTime_clades)) {
  title <- names(morphospaceTime_clades[i])
  df <- morphospaceTime_clades[[i]]
  means <- morphospaceTime_means[[i]]
  p <- ggplot(transform(df, trait = factor(trait, levels = sort)), aes(time, range))+
    geom_point()+
    geom_hline(data = transform(means, trait = factor(trait, levels = sort)), mapping = aes(yintercept = mean)) +
      labs(title=title)+
    facet_manual(trait~., scales = "free", design2, trim_blank = FALSE)
  plots[[i]] <- p
}
plots
for (i in 1:length(plots)) {
  ggsave(plots[[i]], filename = paste0(nameclades[i],"_morphospace_new.pdf"), scale = 5)
}

ggplot(transform(df, trait = factor(trait, levels = sort)), aes(time, range))+
  +
  geom_hline(data = transform(means, trait = factor(trait, levels = sort)), mapping = aes(yintercept = mean)) +
  labs(title=title)+
  facet_manual(trait~., scales = "free", design2, trim_blank = FALSE)


colors=c("#234d20","#77ab59","#c9df8a","#ffff66","#ffcc66","#ff6600","#ff00ff","#cc66ff","#663300","#996600","#0099ff","#0000ff")
colors2=c("brown","brown","brown","red","red","red","orange","orange","darkblue","darkblue","#1B96CC","#1B96CC")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
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



for (i in 1:length(plots_scaled)) {
  ggsave(plots_scaled[[i]], filename = paste0(nameclades[i],"_morphospace_scaled.pdf"), scale = 3)
}

