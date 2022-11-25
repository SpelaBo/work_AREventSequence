## ladderize the tree and make sure is ultrametric
tree_morpho_l <- ladderize(tree_morpho, right = F)
is.ultrametric(tree_morpho_l)

H = nodeHeights(tree_morpho_l, root.edge = F)  ### get age of each node
root=max(H[,2]) ##   get age of the root  ]]

##  define time slices: start at the root and increase in 0.5 Million years + additionally the tips
n= 0.5
v= c(seq(from=0, to=root, by= n) , root)
v2= v[-(length(v)-1)]

comb2 <- comb[order(match(rownames(comb), tree_morpho_l$tip.label)), , drop = FALSE] # order reconstructions to match tree tips
xx = node.depth.edgelength(tree_morpho_l) # get age of each node



body_ASR <- comb2[,2]


###  combine states of the nodes with node ages

yy.bayPC =  cbind(body_ASR, xx)            ###  seting reconstructed phenotype with x coolrdinate (age) of corresponding node

##### for each timepoint reconstruct the state linearly to the distance between the two node stages
Yt= list()
for ( i in c(1: (length(v)-1)) ) {      
  tmp.bayPC3      = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
  tmp.bayPC.node3 = tree_morpho_l $edge[tmp.bayPC3,]                   ### get nodes of all branches between root and timepoint 
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
ggplot(statsOtime, aes(time, range))+geom_line()
