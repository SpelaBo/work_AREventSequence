# functions for modeling ####
allModels <-  function(tree, morpho_df, traits) {
  model_BM <- mvBM(tree,morpho_df[,traits],model = "BM1")
  model_BMM <- mvBM(tree,morpho_df[,traits],model = "BMM", param = list(root = F))
  model_EB <- mvEB(tree,morpho_df[,traits])
  model_OU <- mvOU(tree,morpho_df[,traits],model = "OU1", param = list(root = F))
  model_OUM <- mvOU(tree,morpho_df[,traits],model = "OUM", param = list(root = F))
  models <- list(model_BM,model_BMM,model_EB,model_OU, model_OUM)
  return(models)
}
testBestModel <- function(models){
  aic <- aicw(models)
  best <- models[[which.max(aic$aicweights)]]
  return(list(aic,best))
}

## Node translation and subsets of reconstructed morphologies for clades
selectTipNode <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}
nodeTranslation <- function(trees){
  node_translation_df <- list()
  for (i in 1:length(trees)) {
    tree <- trees[[i]]
    t <- data.frame(
      "node" = tree$edge[,1],
      "label" = sapply(tree$edge[,1],
                       selectTipNode,
                       tree = tree))
    t <- t[!duplicated(t), ]
    t <- t[order(t$node),]
    node_translation_df[[i]] <- t
  }
  return(node_translation_df)
}
#trees=tree_morpho_clades,translation = node_translation,traits = comb, descendants = desc_clades, nameclades = nameclades
reconstructionsCladeSubsets <- function(trees,translation,traits, descendants,nameclades) {
  comb_clades <- list()
  for (i in 1:length(trees)) {
    tree <- trees[[i]]
    nnode <- tree$Nnode
    ntip <- nnode+1
    node_seq <- c((ntip+1):(ntip+nnode))
    desc <- descendants[[i]]
    comb_sub <- as.data.frame(traits[desc,])
    comb_sub$node <- rownames(comb_sub)
    nodes <- translation[[i]]
    comb_sub[node_seq,]$node <- nodes$node
    rownames(comb_sub) <- comb_sub$node
    comb_sub <- subset(comb_sub, select = -node)
    comb_clades[[i]] <- comb_sub
  }
  names(comb_clades) <- nameclades
  return(comb_clades)
}



#function that returns list of matrices per each time slice with columns: trait value and time slice 
traitPerTimeSlice <- function(tree,yy.bayPC,v,H){
  Yt= list()
  for ( i in c(1: (length(v)-1)) ) {      
    tmp.bayPC3 = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))# get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
    tmp.bayPC.node3 = tree$edge[tmp.bayPC3,]# get nodes of all branches between root and timepoint 
    start = yy.bayPC[tmp.bayPC.node3[,1],]
    end = yy.bayPC[tmp.bayPC.node3[,2],]
    t = rep(v[i],times= length(end[,1]))
    y3 =(((end[,1]-start[,1])/(end[,2]-start[,2])) * (v[i] - start[,2])) + start[,1]
    Yt[[i]] = cbind(y3, t)
  }
  return(Yt)
}



# v funkcijo hočeš podat: drevo, time slice, states on the nodes
## n: define time slices: start at the root and increase in 0.5 Million years + additionally the tips

###  seting reconstructed phenotype with x coolrdinate (age) of corresponding node
morphospaceThroughTime <- function(tree, n, traits){
  if(is.ultrametric(tree)==FALSE) stop('tree is not ultrametric')
  H = nodeHeights(tree, root.edge = F)  ### get age of each node
  root=max(H[,2]) ##   get age of the root
  v= c(seq(from=0, to=root, by= n) , root)
  v2= v[-(length(v)-1)]
  xx = node.depth.edgelength(tree) # get age of each node
  ASR_age =  cbind(traits, xx)# combine states of the nodes with node ages
morphospaceTime_list <- list()
for (j in 1:(ncol(ASR_age)-1)) {
  yy.bayPC <- ASR_age[,c(j,ncol(ASR_age))]#get ith trait
  Yt <- traitPerTimeSlice(tree,yy.bayPC,v,H)
  out3=NULL
  for ( k in 1: length(v2)){
    tmp3= unlist(Yt[[k]])
    tmp3.2= data.frame(tmp3, "sample" = v2[k])
    out3= rbind(out3, tmp3.2)
  }
  d=out3[,c(1,3)]
  d$sample2 = as.factor(d$sample)
  statsOtime= as.data.frame(matrix(NA, length(levels(d$sample2)) , 1 )) # prepare the output with stats per time slice
  statsOtime[,1]= as.numeric(as.character(levels(d$sample2)))
  names(statsOtime)="time"
  minv= aggregate( d[,1], list(d$sample2), min)# get range from each timeslice: min
  maxv= aggregate( d[,1], list(d$sample2), max)# get range from each timeslice: max
  statsOtime$range= abs(minv[,-1]-maxv[,-1])# range of the axes
  count_per_bin= aggregate( d[,1], list(d$sample2), length )# get number of linnegaes per time slice
  statsOtime$numLineages = count_per_bin[,2]# combine nr. of linneages and trait name
  statsOtime$trait <- colnames(yy.bayPC)[1]
  morphospaceTime_list[[j]] <- statsOtime
}
morphospaceTime_df <- bind_rows(morphospaceTime_list)
}


range01 <- function(x){(x-min(x))/(max(x)-min(x))}