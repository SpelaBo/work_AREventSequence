#### Packages ####
library(readxl)
library(vegan)
library(phytools)
library(geiger)
library(mvMORPH)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggh4x)
#### Data ####
# .tree #### 
tree_mean_old <- read.nexus("./data/20220325_niphargus_calibrated_edit.tree") # max credibility tree, mean heights 
tree_mean_old <- drop.tip(tree_mean_old, "Nl_nolli")
is.ultrametric(tree_mean_old)
tree_mean <- force.ultrametric(tree_mean_old,method="extend")
plot(tree_mean_old, show.tip.label=F)
axisPhylo()

# .morphology #### 
morphology_old <- read_excel("./data/morphology.xlsx", na = "NA")
morphology_old <- as.data.frame(morphology_old)
row.names(morphology_old) <- morphology_old$species
morphology_old[1] <- NULL
## calculate gnathopod size and remove unnecessary columns
morphology_old$gpI6_size <- morphology_old$gI6_length+morphology_old$gI6_width+morphology_old$gI6_diag
morphology_old$gpII6_size <- morphology_old$gII6_length+morphology_old$gII6_width+morphology_old$gII6_diag
morphology_old$Icosinus <- ((morphology_old$gI6_length)^2+(morphology_old$gI6_width)^2-(morphology_old$gI6_diag)^2)/(2*morphology_old$gI6_length*morphology_old$gI6_width)
morphology_old$IIcosinus <- ((morphology_old$gII6_length)^2+(morphology_old$gII6_width)^2-(morphology_old$gII6_diag)^2)/(2*morphology_old$gII6_length*morphology_old$gII6_width)
morphology_old$Ialpha <- 180*(acos(morphology_old$Icosinus))/pi
morphology_old$IIalpha <- 180*(acos(morphology_old$IIcosinus))/pi
morphology_old <- morphology_old[,c(1:11,20:21,24:25)]
morphology_old <- as.matrix(morphology_old[ order(row.names(morphology_old)), ])
## Remove NA-s
morphology_old_noNA <- morphology_old[complete.cases(morphology_old), ]

## .morpho_tree and gls ####
to_keep_morpho_old=rownames(morphology_old_noNA)
tree_morpho_old=keep.tip(tree_mean_old, to_keep_morpho_old)
name.check(tree_morpho_old,morphology_old_noNA) 
is.ultrametric(tree_morpho_old)
## ladderized morpho tree
tree_morpho_old_l <- ladderize(tree_morpho_old, right = F)
is.ultrametric(tree_morpho_old_l)
morphology_old_noNA<-morphology_old_noNA[order(match(rownames(morphology_old_noNA), tree_morpho_old$tip.label)), , drop = FALSE]
all(rownames(morphology_old_noNA)==tree_morpho_old$tip.label)

# Phylogenetically corrected GLS
phyl_gls_morpho_old <- phyl.resid(tree_morpho_old, morphology_old_noNA[,1], morphology_old_noNA[,c(2:13)])
residuals_gls_old <- as.matrix(phyl_gls_morpho_old$resid)
body_length_old <- as.matrix(morphology_old_noNA[,1], )
gnato_angles_old <- as.matrix(morphology_old_noNA[,14:15], )
colnames(body_length_old)<-c("body_length")
morpho_all_gls_old <- cbind(body_length_old, residuals_gls_old, gnato_angles_old)
all(rownames(morpho_all_gls_old)==tree_morpho_old$tip.label)

## .trait groups ####
trophic_traits2 <- c("Ialpha","IIalpha")
trophic_traits <- c("gpI6_size","gpII6_size")
body_traits <- c("cxII_length","pV2_width","pVI2_width") # zmanjšam nabor traitov za hitrejše računanje (cov2cor(model_OUM_body$sigma)), odvzamem "cxIII_length","pVII2_width"
loco_traits <-  c("body_length", "pV_length", "pVII_length") # enako kot zgoraj, odvzamem "pVI_length"
sens_traits <- c("antennaI","antennaII")
traits <- list(body_traits,loco_traits,sens_traits,trophic_traits,trophic_traits2)
names(traits) <- c("body_traits","loco_traits","sens_traits","trophic_traits","trophic_traits2")

## .clades ####
# define subclades
#pannonian
speciesA1_old=c("N_ablaskiri","N_ambulator","N_aquilex_T81","N_aquilex_T94","N_bihorensis_A","N_bihorensis_B","N_carniolicus","N_cf_tauri_3","N_cvetkovi","N_daniali","N_dimorphus","N_dobati","N_fongi","N_gebhardti","N_inermis","N_pontoruffoi","N_racovitzai","N_sp_A_iz_wolfi","N_spn_Gumbrini","N_spn_Huda_Luknja","N_spn_NC561","N_spn_NC575","N_spn_NC576","N_spn_NC580","N_spn_NC586","N_spn_NC613","N_spn_NC630","N_spn_NC660","N_spn_NC702","N_spn_NC707","N_spn_NC711","N_spn_NC716","N_spn_NC718","N_spn_NC724","N_spn_NC850","N_spn_ND053","N_spn_ND113","N_spn_ND209","N_spn_ND429","N_spn_ND565","N_spn_ND770","N_spn_Podutik","N_spn_Sbirkovska_cave","N_spn_Sitarjevec","N_spn_Sukhaja_balka","N_spn_tauri_Brezno3src","N_tauri","N_vadimi","N_wolfi_A")
#pontic
speciesA2_old=c("C_paradoxa","N_aberrans","N_alpinus","N_andropus","N_aquilex_B","N_armatus","N_bajuvaricus","N_barbatus","N_carpathicus","N_danielopoli","N_decui","N_dissonus","N_galvagnii","N_grandii","N_italicus","N_karkabounasi","N_labacensis_A","N_labacensis_B","N_lattingerae","N_longidactylus","N_microcerberus_A","N_microcerberus_B","N_minor","N_moldavicus_cf","N_multipennatus","N_parapupetta","N_pectinicauda","N_petrosani","N_pupetta","N_serbicus","N_similis","N_sp_T110","N_spn_Arkadi","N_spn_Jelovica","N_spn_NC217","N_spn_NC624","N_spn_NC786","N_spn_ND182","N_spn_ND188","N_spn_ND211","N_spn_ND213","N_spn_ND216","N_spn_ND301","N_spn_ND426","N_spn_ND439","N_spn_ND540","N_spn_ND605","N_spn_ND669","N_spn_ND670","N_spn_ND705","N_spn_ND924","N_spn_Vodni_kevder","N_strouhali","N_tamaninii","N_transitivus","N_transsylvanicus_A","N_transsylvanicus_B","N_transsylvanicus_C","N_wolfi_B")
#south dinaric
speciesB_old=c("N_aulicus","N_balcanicus","N_bilecanus","N_boskovici","N_brevicuspis_A","N_brevicuspis_B","N_buturovici","N_dabarensis","N_factor","N_hercegovinensis","N_hvarensis_B","N_hvarensis_C","N_kusceri","N_lunaris_B","N_miljeticus","N_navotinus","N_podgoricensis","N_polymorphus","N_salernianus_A","N_salernianus_B","N_salernianus_C","N_salernianus_D","N_sketi","N_sp_A","N_sp_C","N_sp_D","N_sp_E","N_spn_NB613","N_spn_NC759","N_spn_NC812","N_spn_ND096","N_spn_ND193","N_spn_ND195","N_spn_ND198","N_spn_ND201","N_spn_ND204","N_spn_ND270","N_spn_ND374","N_spn_ND398","N_spn_ND463","N_spn_ND501","N_spn_ND558","N_spn_ND690","N_spn_ND739","N_spn_ND749","N_spn_ND752","N_spn_ND833","N_spn_Vilina_pecina","N_trullipes","N_vjetrenicensis_A","N_vjetrenicensis_B","N_zagorae")
#west balkan
speciesC_old=c("Cn_lubuskensis","N_alpheus","N_anchialinus","N_antipodes","N_arbiter","N_arethusa","N_brevirostris","N_cornicolanus","N_croaticus","N_doli","N_fjakae","N_hebereri_1","N_hebereri_2","N_hebereri_3","N_ictus","N_kolombatovici","N_liburnicus_A","N_liburnicus_B","N_longiflagellum","N_lunaris_A","N_messanai","N_mirocensis","N_orcinus","N_pachytelson","N_parenzani","N_patrizii","N_pectencoronatae","N_pincinovae","N_rejici","N_rhenorhodanensis_1","N_rhenorhodanensis_2","N_rhenorhodanensis_3","N_rhenorhodanensis_4","N_rhenorhodanensis_6","N_rhenorhodanensis_9","N_salonitanus","N_spn_Lusci_Palanka","N_spn_NB849","N_spn_ND119","N_spn_ND130","N_spn_ND139","N_spn_ND160","N_spn_ND256","N_spn_ND338","N_spn_ND469","N_spn_ND711","N_spn_ND818","N_spn_Prodisce_Sava","N_spn_Zenadija","N_spn_tertius","N_stefanellii","N_stefanellii_A","N_stefanellii_C","N_stenopus_A","N_stenopus_B","N_steueri","N_subtypicus","Nb_orophobata")
nameclades <- c("A1_pannonian", "A2_pontic", "B_south_dinaric","C_west_balkan")

#subset of species in clades with morpho data
speciesA1_m_old <- subset(speciesA1_old, speciesA1_old %in% tree_morpho_old$tip.label)
speciesA2_m_old <- subset(speciesA2_old, speciesA2_old %in% tree_morpho_old$tip.label)
speciesB_m_old <- subset(speciesB_old, speciesB_old %in% tree_morpho_old$tip.label)
speciesC_m_old <- subset(speciesC_old, speciesC_old %in% tree_morpho_old$tip.label)
tree_morpho_old_painted <-tree_morpho_old_l
is.ultrametric(tree_morpho_old_painted)
## ...mrca ####
mrcaA1_old <- findMRCA(tree_morpho_old, speciesA1_m_old)
mrcaA2_old <- findMRCA(tree_morpho_old, speciesA2_m_old)
mrcaB_old <- findMRCA(tree_morpho_old, speciesB_m_old)
mrcaC_old <- findMRCA(tree_morpho_old, speciesC_m_old)

## ...descendants ####
descA1_old <- sort(c(getDescendants(tree_morpho_old,mrcaA1_old),mrcaA1_old))
descA2_old <- sort(c(getDescendants(tree_morpho_old,mrcaA2_old),mrcaA2_old))
descB_old <- sort(c(getDescendants(tree_morpho_old,mrcaB_old),mrcaB_old))
descC_old <- sort(c(getDescendants(tree_morpho_old,mrcaC_old),mrcaC_old))
desc_clades_old <- list(descA1_old,descA2_old,descB_old,descC_old)
names(desc_clades_old) <- nameclades

# ...tree w. painted clades ####

tree_morpho_old_painted<-paintSubTree(tree_morpho_old_painted,node=mrcaA1_old,state="A1",stem=T)
tree_morpho_old_painted<-paintSubTree(tree_morpho_old_painted,node=mrcaA2_old,state="A2",stem=T)
tree_morpho_old_painted<-paintSubTree(tree_morpho_old_painted,node=mrcaB_old,state="B",stem=T)
tree_morpho_old_painted<-paintSubTree(tree_morpho_old_painted,node=mrcaC_old,state="C",stem=T)
cols<-c("black","deepskyblue3","chartreuse3","gold", "pink"); names(cols)<-c(1,"A1","A2","B","C")
plotSimmap(tree_morpho_old_painted,cols,lwd=2,pts=F, ftype="off")
axisPhylo()

## ...clade trees ####
# label nodes for future reference of original node numbers
tree_morpho_old_n <- tree_morpho_old
tree_morpho_old_n$node.label <- c((tree_morpho_old_n$Nnode+2):((tree_morpho_old_n$Nnode*2)+1))
#use labeled tree for subtrees
tree_morpho_old_A1=keep.tip(tree_morpho_old_n, speciesA1_m_old)
tree_morpho_old_A2=keep.tip(tree_morpho_old_n, speciesA2_m_old)
tree_morpho_old_B=keep.tip(tree_morpho_old_n, speciesB_m_old)
tree_morpho_old_C=keep.tip(tree_morpho_old_n, speciesC_m_old)
tree_morpho_old_clades <- list(tree_morpho_old_A1,tree_morpho_old_A2,tree_morpho_old_B,tree_morpho_old_C)
names(tree_morpho_old_clades) <- nameclades

### analiza diverzifikacije, rekonstrukcija na celotnem drevesu ####
# Modelling, genus level ####
models_list_old <- list()
for (i in 1:length(traits)) {
  trait <- traits[[i]]
  models <- allModels(tree_morpho_old_painted,morpho_all_gls_old,trait)
  models_list_old[[i]] <- models
}
names(models_list_old) <- names(traits)

best_model_list_old <- list()
for (i in 1:length(models_list_old)) {
  models <- models_list_old[[i]]
  best_model <- testBestModel(models)
  best_model_list_old[[i]] <- best_model
}
names(best_model_list_old) <- names(models_list_old)

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_list_old <-  list()
for (i in 1:length(traits)) {
  model <- best_model_list_old[[i]]
  trait <- traits[[i]]
  ASR<-estim(tree_morpho_old_painted,morpho_all_gls_old[,trait], model, asr=TRUE)
  comb <- as.data.frame(rbind(morpho_all_gls_old[,trait],ASR$estim))
  comb$species <- rownames(comb)
  
  ASR_list_old[[i]] <- comb
}
ASR_df_old <-  Reduce(inner_join,ASR_list_old)
rownames(ASR_df_old) <- ASR_df_old$species
ASR_df_old[c("species")] <- NULL


# clade subsetting ####
node_translation_old <- nodeTranslation(tree_morpho_old_clades)

# ASR for each clade ####
ASR_clades_old <- reconstructionsCladeSubsets(trees = tree_morpho_old_clades,
                                          translation = node_translation_old,
                                          traits = ASR_df_old, 
                                          descendants = desc_clades_old,
                                          nameclades = nameclades)

# plot ####
design <- matrix(c(1,2,3,4,5,6,7,8,0,9,10,0,11,12,0), nrow=3, byrow = F)
for (j in 1:length(nameclades)) {
  c=as.matrix(ASR_clades_old[[j]])
  #c <- c[,sort]
  t=tree_morpho_old_clades[[j]]
  title <- nameclades[j]
  pdf(file = paste0(title,"_pheno_old.pdf"))
  layout(design)
  par(mar = c(4, 2, 1, 1), oma = c(1,0, 0, 0))
  for (i in 1:ncol(c)) {
    phenogram(tree=t,x=c[,i], ftype="off", ylim=c(min(c[,i]),max(c[,i])), xlab = colnames(c)[i])
  }
  title(title)
  dev.off()
}

phenogram(tree=tree_morpho_C,x=as.matrix(ASR_clades[[j]])[,6], xlab = colnames(ASR_clades[[j]][6]))
phenogram(tree=tree_morpho_old_C,x=as.matrix(ASR_clades_old[[j]])[,6], xlab = colnames(ASR_clades_old[[j]][6]))
  
