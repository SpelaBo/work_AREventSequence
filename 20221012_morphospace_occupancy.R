H = nodeHeights(tree_morpho, root.edge = F)  ### get age of each node
root=max(H[,2]) ##   get age of the root

##  define time slices: start at the root and increase in 0.5 Million years + additionally the tips
n= 0.5
v= c(seq(from=0, to=root, by= n) , root)
v2= v[-(length(v)-1)]
xx = node.depth.edgelength(tree_morpho) # get age of each node

###  combine states of the nodes with node ages
ASR_age =  cbind(comb, xx)

###  seting reconstructed phenotype with x coolrdinate (age) of corresponding node
morphospaceTime_list <- list()
for (j in 1:(ncol(ASR_age)-1)) {
  yy.bayPC <- ASR_age[,c(j,ncol(ASR_age))]#get ith trait
  Yt= list()
  for ( i in c(1: (length(v)-1)) ) {      
    tmp.bayPC3      = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
    tmp.bayPC.node3 = tree_morpho $edge[tmp.bayPC3,]                   ### get nodes of all branches between root and timepoint 
    start           = yy.bayPC[tmp.bayPC.node3[,1],]
    end             = yy.bayPC[tmp.bayPC.node3[,2],]
    t               = rep(v[i],times= length(end[,1]))
    y3              =(((end[,1]-start[,1])/(end[,2]-start[,2])) * (v[i] - start[,2])) + start[,1]
    Yt[[i]]         =cbind(y3, t)
  }
  out3=NULL
  for ( k in 1: length(v2)){
    tmp3= unlist(Yt[[k]])
    tmp3.2= data.frame(tmp3, "sample" = v2[k])
    out3= rbind(out3, tmp3.2)
  }
  d=out3[,c(1,3)]
  d$sample2 = as.factor(d$sample)
  ######### prepare the output for the different stats per timeslice over the trait axes
  statsOtime= as.data.frame(matrix(NA, length(levels(d$sample2)) , 1 ))
  statsOtime[,1]= as.numeric(as.character(levels(d$sample2)))
  names(statsOtime)="time"
  ######### get range from each timeslice:
  minv= aggregate( d[,1], list(d$sample2), min)
  maxv= aggregate( d[,1], list(d$sample2), max)
  ### range of the axes
  statsOtime$range= abs(minv[,-1]-maxv[,-1])
  ###################### get number of linnegaes per time slice
  count_per_bin= aggregate( d[,1], list(d$sample2), length )
  ## combine
  statsOtime$numLineages = count_per_bin[,2]
  statsOtime$trait <- colnames(yy.bayPC)[1]
  
  morphospaceTime_list[[j]] <- statsOtime
}


morphospaceTime_df <- bind_rows(morphospaceTime_list)
morphospaceTime_df_means <- morphospaceTime_df %>% 
  group_by(trait) %>%
  summarise(mean=max(range)/2)

ggplot(morphospaceTime_df, aes(time, range))+
  geom_point()+
  geom_hline(data = morphospaceTime_df_means, mapping = aes(yintercept = mean)) +
  facet_wrap(~trait, scales = "free")




morphospaceTime_clades <- list()
for (z in 1:length(comb_clades)) {
  comb_clade <- comb_clades[[z]]
  tree_morpho_clade <- tree_morpho_clades[[z]]
  H = nodeHeights(tree_morpho_clade, root.edge = F)
  root=max(H[,2])
  n= 0.5
  v= c(seq(from=0, to=root, by= n) , root)
  v2= v[-(length(v)-1)]
  xx = node.depth.edgelength(tree_morpho_clade)
  ASR_clade =  cbind(comb_clade, xx)
  morphospaceTime_list <- list()
  for (j in 1:(ncol(ASR_clade)-1)) {
    yy.bayPC <- ASR_clade[,c(j,ncol(ASR_clade))]#get ith trait
    Yt= list()
    for ( i in c(1: (length(v)-1)) ) {      
      tmp.bayPC3      = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
      tmp.bayPC.node3 = tree_morpho_clade $edge[tmp.bayPC3,]                   ### get nodes of all branches between root and timepoint 
      start           = yy.bayPC[tmp.bayPC.node3[,1],]
      end             = yy.bayPC[tmp.bayPC.node3[,2],]
      t               = rep(v[i],times= length(end[,1]))
      y3              =(((end[,1]-start[,1])/(end[,2]-start[,2])) * (v[i] - start[,2])) + start[,1]
      Yt[[i]]         =cbind(y3, t)
    }
    out3=NULL
    for ( k in 1: length(v2)){
      tmp3= unlist(Yt[[k]])
      tmp3.2= data.frame(tmp3, "sample" = v2[k])
      out3= rbind(out3, tmp3.2)
    }
    d=out3[,c(1,3)]
    d$sample2 = as.factor(d$sample)
    ######### prepare the output for the different stats per timeslice over the trait axes
    statsOtime= as.data.frame(matrix(NA, length(levels(d$sample2)) , 1 ))
    statsOtime[,1]= as.numeric(as.character(levels(d$sample2)))
    names(statsOtime)="time"
    ######### get range from each timeslice:
    minv= aggregate( d[,1], list(d$sample2), min)
    maxv= aggregate( d[,1], list(d$sample2), max)
    ### range of the axes
    statsOtime$range= abs(minv[,-1]-maxv[,-1])
    ###################### get number of linnegaes per time slice
    count_per_bin= aggregate( d[,1], list(d$sample2), length )
    ## combine
    statsOtime$numLineages = count_per_bin[,2]
    statsOtime$trait <- colnames(yy.bayPC)[1]
    
    morphospaceTime_list[[j]] <- statsOtime
  }
  morphospaceTime_df <- bind_rows(morphospaceTime_list)
  morphospaceTime_clades[[z]] <- morphospaceTime_df
}

morphospaceTime_clades

morphospaceTime_means <- list()
for (i in 1:length(morphospaceTime_clades)) {
  morphospaceTime_df_means <- morphospaceTime_clades[[i]] %>% 
    group_by(trait) %>%
    summarise(mean=max(range)/2)
  morphospaceTime_means[[i]] <- morphospaceTime_df_means
}
names(morphospaceTime_clades) <- nameclades
names(morphospaceTime_means) <- nameclades

sort <- c(body_traits,loco_traits,sens_traits,trophic_traits)

design <- matrix(c(1,2,3,4,5,6:9,0,10,11,0,0,0,12,13,0,0,0), nrow=5, byrow = F)
design2 <- matrix(c(1,2,3,4,5,6:9,NA,10,11,NA,NA,NA,12,13,NA,NA,NA), nrow=5, byrow = F)

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
for (i in 1:length(plots)) {
  ggsave(plots[[i]], filename = paste0(nameclades[i],".pdf"), scale = 5)
}




layout(design)
layout.show(n=13)
