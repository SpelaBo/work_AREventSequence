#### Packages ####
library(readxl)
library(vegan)
library(phytools)
library(geiger)
library(mvMORPH)
library(RColorBrewer)
library(ggplot2)
library(plyr)
library(tidyverse)
library(ggh4x)
#### Data ####
# .tree #### 
tree_mean <- read.nexus("./data/20220721_niphargus.tree") # max credibility tree, mean heights 
tree_mean <- drop.tip(tree_mean, "Nl_nolli")
is.ultrametric(tree_mean)
#tree_mean <- force.ultrametric(tree_mean,method="extend")
plot(tree_mean, show.tip.label=F)
axisPhylo()

# .morphology #### 
morphology <- read_excel("./data/20221018_morphology_imput.xlsx", na = "NA")
morphology <- as.data.frame(morphology)
row.names(morphology) <- morphology$species
morphology[1] <- NULL
## calculate gnathopod size and remove unnecessary columns
morphology$gpI6_size <- morphology$gI6_length+morphology$gI6_width+morphology$gI6_diag
morphology$gpII6_size <- morphology$gII6_length+morphology$gII6_width+morphology$gII6_diag
morphology$Icosinus <- ((morphology$gI6_length)^2+(morphology$gI6_width)^2-(morphology$gI6_diag)^2)/(2*morphology$gI6_length*morphology$gI6_width)
morphology$IIcosinus <- ((morphology$gII6_length)^2+(morphology$gII6_width)^2-(morphology$gII6_diag)^2)/(2*morphology$gII6_length*morphology$gII6_width)
morphology$Ialpha <- 180*(acos(morphology$Icosinus))/pi
morphology$IIalpha <- 180*(acos(morphology$IIcosinus))/pi
morphology <- morphology[,c(1:11,20:21,24:25)]
morphology <- as.matrix(morphology[ order(row.names(morphology)), ])
## Remove NA-s
morphology_noNA <- morphology[complete.cases(morphology), ]

## .morpho_tree and gls ####
to_keep_morpho=rownames(morphology_noNA)
tree_morpho=keep.tip(tree_mean, to_keep_morpho)
name.check(tree_morpho,morphology_noNA) 
is.ultrametric(tree_morpho)
## ladderized morpho tree
tree_morpho_l <- ladderize(tree_morpho, right = F)
is.ultrametric(tree_morpho_l)
morphology_noNA<-morphology_noNA[order(match(rownames(morphology_noNA), tree_morpho$tip.label)), , drop = FALSE]
all(rownames(morphology_noNA)==tree_morpho$tip.label)

# Phylogenetically corrected GLS
phyl_gls_morpho <- phyl.resid(tree_morpho, morphology_noNA[,1], morphology_noNA[,c(2:13)])
residuals_gls <- as.matrix(phyl_gls_morpho$resid)
body_length <- as.matrix(morphology_noNA[,1], )
gnato_angles <- as.matrix(morphology_noNA[,14:15], )
colnames(body_length)<-c("body_length")
morpho_all_gls <- cbind(body_length, residuals_gls, gnato_angles)
all(rownames(morpho_all_gls)==tree_morpho$tip.label)

## .trait groups ####
trophic_traits2 <- c("Ialpha","IIalpha")
trophic_traits <- c("gpI6_size","gpII6_size")
body_traits <- c("cxII_length","pV2_width","pVI2_width") # zmanjšam nabor traitov za hitrejše računanje (cov2cor(model_OUM_body$sigma)), odvzamem "cxIII_length","pVII2_width"
loco_traits <-  c("body_length", "pV_length", "pVII_length") # enako kot zgoraj, odvzamem "pVI_length"
sens_traits <- c("antennaI","antennaII")
traits <- list(body_traits,loco_traits,sens_traits,trophic_traits,trophic_traits2)
names(traits) <- c("body_traits","loco_traits","sens_traits","trophic_traits","trophic_traits2")

## .clades ####
#pannonian
speciesA1=c("N_spn_NC576","N_spn_NC586","N_spn_NC580","N_spn_NC660","N_spn_Sukhaja_balka","N_inermis","N_spn_NC575","N_spn_NC561","N_spn_NC613","N_alasonius","N_spn_NC724","N_spn_NC707","N_spn_NC711","N_spn_NC716","N_vadimi","N_spn_Gumbrini","N_spn_NC702","N_spn_NC630","N_ablaskiri","N_dimorphus","N_gebhardti","N_spn_ND429","N_spn_ND565","N_spn_Huda_Luknja","N_spn_ND209","N_aquilex_T81","N_spn_Borgnone","N_ambulator","N_tauri","N_fongi","N_sp_A_iz_wolfi","N_spn_ND113","N_cf_tauri_3","N_carniolicus","N_wolfi_A","N_aquilex_T94","N_dobati","N_spn_NC850","N_spn_ND770","N_spn_Podutik","N_spn_tauri_Brezno3src","N_spn_Sitarjevec","N_bihorensis_A","N_bihorensis_B","N_cvetkovi","N_spn_ND053","N_pontoruffoi","N_racovitzai","N_daniali","N_spn_Sbirkovska_cave")
#pontic
speciesA2=c("N_minor","N_carpathicus","N_alpinus","N_spn_Jelovica","N_spn_NC786","N_wolfi_B","N_spn_NC217","N_spn_Vodni_kevder","N_tamaninii","N_spn_ND439","N_galvagnii","N_spn_ND301","N_spn_ND540","N_spn_ND188","N_spn_ND924","N_spn_ND426","N_spn_ND705","N_spn_ND182","N_labacensis_A","N_lattingerae","N_longidactylus","N_spn_NC508","N_spn_ND670","N_aquilex_B","N_spn_ND216","N_danielopoli","N_strouhali","N_armatus","N_labacensis_B","N_spn_ND669","N_spn_ND211","N_barbatus","N_similis","N_aberrans","N_multipennatus","N_serbicus","N_spn_ND605","N_grandii","N_italicus","N_spn_ND213","N_bajuvaricus","N_microcerberus_A","N_microcerberus_B","N_dissonus","N_transitivus","N_parapupetta","N_pupetta","C_paradoxa","N_pectinicauda","N_karkabounasi","N_spn_Arkadi","N_spn_NC624","N_transsylvanicus_C","N_moldavicus_cf","N_petrosani","N_sp_T110","N_decui","N_andropus","N_transsylvanicus_A","N_transsylvanicus_B")
#south dinaric
speciesB=c("N_bilecanus","N_vjetrenicensis_A","N_spn_NB613","N_trullipes","N_spn_ND463","N_vjetrenicensis_B","N_kusceri","N_hercegovinensis","N_balcanicus","N_polymorphus","N_dabarensis","N_navotinus","N_spn_ND374","N_spn_ND270","N_spn_Vilina_pecina","N_podgoricensis","N_buturovici","N_spn_ND690","N_spn_ND752","N_spn_ND833","N_spn_ND193","N_spn_ND201","N_spn_ND198","N_brevicuspis_B","N_spn_ND204","N_sketi","N_spn_ND195","N_brevicuspis_A","N_spn_NF987","N_factor","N_spn_ND749","N_spn_NE178","N_spn_NF592","N_spn_NF550","N_spn_NE188","N_spn_NF575","N_spn_NC759","N_spn_ND096","N_aulicus","N_lunaris_B","N_hvarensis_B","N_spn_ND501","N_spn_ND558","N_hvarensis_C","N_miljeticus","N_sp_A","N_salernianus_A","N_salernianus_B","N_salernianus_C","N_salernianus_D","N_boskovici","N_spn_NC812","N_sp_C","N_spn_ND739","N_sp_E","N_zagorae","N_spn_ND398","N_sp_D")
#west balkan
speciesC=c("N_doli","N_brevirostris","N_spn_ND818","N_cornicolanus","N_hebereri_1","N_parenzani","N_stefanellii","N_messanai","N_stefanellii_C","N_ictus","N_stefanellii_A","N_spn_Zenadija","N_hebereri_2","N_hebereri_3","N_rejici","N_antipodes","N_arbiter","N_spn_ND469","N_alpheus","N_arethusa","N_anchialinus","N_salonitanus","N_fjakae","N_pincinovae","N_spn_NB849","N_spn_ND338","N_spn_ND139","N_longiflagellum","N_pachytelson","N_spn_ND119","N_spn_ND160","N_spn_Prodisce_Sava","N_lunaris_A","N_spn_Lusci_Palanka","N_stenopus_A","N_spn_ND130","N_spn_tertius","N_mirocensis","N_stenopus_B","N_spn_ND711","N_pectencoronatae","N_kolombatovici","N_spn_ND256","N_steueri","N_liburnicus_A","N_liburnicus_B","N_croaticus","N_patrizii","N_subtypicus","Nb_orophobata","Cn_lubuskensis","N_spn_NF561","N_orcinus")
nameclades <- c("A1_pannonian", "A2_pontic", "B_south_dinaric","C_west_balkan")

#subset of species in clades with morpho data
speciesA1_m <- subset(speciesA1, speciesA1 %in% tree_morpho$tip.label)
speciesA2_m <- subset(speciesA2, speciesA2 %in% tree_morpho$tip.label)
speciesB_m <- subset(speciesB, speciesB %in% tree_morpho$tip.label)
speciesC_m <- subset(speciesC, speciesC %in% tree_morpho$tip.label)

## ...mrca ####
mrcaA1 <- findMRCA(tree_morpho_painted, speciesA1_m)
mrcaA2 <- findMRCA(tree_morpho, speciesA2_m)
mrcaB <- findMRCA(tree_morpho, speciesB_m)
mrcaC <- findMRCA(tree_morpho_painted, speciesC_m)

## ...descendants ####
descA1 <- sort(c(getDescendants(tree_morpho,mrcaA1),mrcaA1))
descA2 <- sort(c(getDescendants(tree_morpho,mrcaA2),mrcaA2))
descB <- sort(c(getDescendants(tree_morpho,mrcaB),mrcaB))
descC <- sort(c(getDescendants(tree_morpho,mrcaC),mrcaC))
desc_clades <- list(descA1,descA2,descB,descC)
names(desc_clades) <- nameclades

# ...tree w. painted clades ####
tree_morpho_painted <-tree_morpho_l
is.ultrametric(tree_morpho_painted)
tree_morpho_painted<-paintSubTree(tree_morpho_painted,node=mrcaA1,state="A1",stem=T)
tree_morpho_painted<-paintSubTree(tree_morpho_painted,node=mrcaA2,state="A2",stem=T)
tree_morpho_painted<-paintSubTree(tree_morpho_painted,node=mrcaB,state="B",stem=T)
tree_morpho_painted<-paintSubTree(tree_morpho_painted,node=mrcaC,state="C",stem=T)
cols<-c("gray","deepskyblue3","purple","red", "orange"); names(cols)<-c(1,"A1","A2","B","C")
plotSimmap(tree_morpho_painted,cols,lwd=2,pts=F, ftype="i", fsize = 0.5)
axisPhylo()

## ...clade trees ####
# label nodes for future reference of original node numbers
tree_morpho_n <- tree_morpho
tree_morpho_n$node.label <- c((tree_morpho_n$Nnode+2):((tree_morpho_n$Nnode*2)+1))
#use labeled tree for subtrees
tree_morpho_A1=keep.tip(tree_morpho_n, speciesA1_m)
tree_morpho_A2=keep.tip(tree_morpho_n, speciesA2_m)
tree_morpho_B=keep.tip(tree_morpho_n, speciesB_m)
tree_morpho_C=keep.tip(tree_morpho_n, speciesC_m)
tree_morpho_clades <- list(tree_morpho_A1,tree_morpho_A2,tree_morpho_B,tree_morpho_C)
names(tree_morpho_clades) <- nameclades

species_nr <- c()
for (i in 1:length(tree_morpho_clades)) {
  x=length(tree_morpho_clades[[i]]$tip.label)
  name=names(tree_morpho_clades[i])
  species_nr[i] <- (paste(name,": ",x))
}
species_nr <- append(species_nr,paste("whole genus:",length(tree_morpho$tip.label)))
write(species_nr, "number_of_species.txt")