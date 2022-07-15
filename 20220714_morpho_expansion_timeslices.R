################################ ################################ ################################ 
################################   sample Ancetral states in time slices  
################################ ################################ ################################


tree_b
comb_body_b[,1]


##### for each timepoint reconstruct the state linearly to the distance between the two node stages
coordinatesTimeslices <- function(tree, data) { #data is named vector of reconstructed morpho states, sorted by nodes
  H = nodeHeights(tree)  ### get age of each node
  root=max(H[,2]) ##   get age of the root
  root2= floor(root*1000)/1000 ## round down the root age to avoid problems with rounding artefacts
  v= c(seq(from=0, to=root2, by= 0.5) , root2)##  define timeslices: start at the root and increase in 1.5 Million years + additionally the tips
  xx = node.depth.edgelength(tree)  ### get age of each node
  Nnodes <- length(tree$tip.label) + tree$Nnode
  yy.bayPC =  as.data.frame(cbind(c(1:Nnodes),data, xx)) ## seting reconstructed phenotype with x coordinate (age) of corresponding node
  Yt= list()
  for ( i in c(1: (length(v)-1)) ) {
    tmp.bayPC3      = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
    tmp.bayPC.node3 = tree $edge[tmp.bayPC3,]                   ### get nodes of all branches between root and timepoint 
    start           = yy.bayPC[tmp.bayPC.node3[,1],]
    end             = yy.bayPC[tmp.bayPC.node3[,2],]
    t               = rep(v[i],times= length(end[,1]))
    k = ((end[,2]-start[,2])/(end[,3]-start[,3]))
    x = (v[i] - start[,3])
    n = start[,2]
    y3              =( k * x ) + n
    Yt[[i]]         =cbind(y3, t)
  }
  v2= v[-(length(v)-1)]
  out3=NULL
  for ( k in 1: length(v2)){
    tmp3= unlist(Yt[[k]])
    tmp3.2= data.frame(tmp3, "sample" = v2[k])
    out3= rbind(out3, tmp3.2)
  }
  out4=out3[,c(1,3)]
  return(out4)
}

coordinatesGroupedBySlice <- function(data){# input from function coordinatesTimeslices
  range_max <-  diff(range(data$y3))
  out_grouped <- data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(diff(range(y3))/range_max*100)
  return(out_grouped)
}

coordinatesXXMorphospace <- function(data, perc){#input from function coordinatesGroupedBySlice
  t <- data[(min(which(data[,2] >= perc))),1]
  return(t)
}

# za vsak stolpec v comb b dobi koordinate 50
T50_clade_b <- c()
for (i in 1:ncol(comb_b)) {
  test <- coordinatesTimeslices(tree_b, comb_b[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T50_clade_b[i] <- c(coordinatesXXMorphospace(test_grouped,50))$sample
}

T90_clade_b <- c()
for (i in 1:ncol(comb_b)) {
  test <- coordinatesTimeslices(tree_b, comb_b[,i])
  test_grouped <- coordinatesGroupedBySlice(test)
  T90_clade_b[i] <- c(coordinatesXXMorphospace(test_grouped,90))$sample
  
}

par(mfrow=c(4,4))
for (i in 1:ncol(comb_b)) {
  phenogram(tree=tree_b,x=comb_b[,i], ftype="off", ylim=c(min(comb_b[,i]),max(comb_b[,i])), xlab = colnames(comb_b)[i])
  abline(v=T50_clade_b[i], col = "blue")
  abline(v=T90_clade_b[i], col = "red")
}