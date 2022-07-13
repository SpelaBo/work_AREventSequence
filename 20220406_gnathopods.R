#### Packages ####
library(readxl)
library(vegan)
library(phytools)
library(geiger)
library(mvMORPH)

library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(plyr)

library(motmot)

#needed for DTT envelope test after Murrell (2018)
library(devtools)
install_github('myllym/spptest', ref = 'no_fastdepth')
library(spptest)
#corrected DTT, source code from https://github.com/mwpennell/geiger-v2/blob/master/R/disparity.R
source("https://raw.githubusercontent.com/mwpennell/geiger-v2/master/R/disparity.R")
#Modified dtt function to call modified MDI function and allow user to change y-axis
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/dtt1.R")
#Modified code for two sided MDI test. Note this is the CORRECTED version of the MDI test in geiger
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/getMDI1.R")

#This function takes a dtt produced object and runs it through the rank envelope test
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/rank_dtt.R")

#### Data ####
# final tree, max credibility tree, mean heights
tree_mean <- read.nexus("./data/20220325_niphargus_calibrated_edit.tree")
tree_mean <- drop.tip(tree_mean, "Nl_nolli")

# Morphology, continuous
morphology <- read_excel("./data/morphology.xlsx", na = "NA")
morphology <- as.data.frame(morphology)
row.names(morphology) <- morphology$species
morphology[1] <- NULL
## calculate gnathopod size and angles
morphology$gpI6_size <- morphology$gI6_length+morphology$gI6_width+morphology$gI6_diag
morphology$gpII6_size <- morphology$gII6_length+morphology$gII6_width+morphology$gII6_diag
morphology$Icosinus <- ((morphology$gI6_length)^2+(morphology$gI6_width)^2-(morphology$gI6_diag)^2)/(2*morphology$gI6_length*morphology$gI6_width)
morphology$IIcosinus <- ((morphology$gII6_length)^2+(morphology$gII6_width)^2-(morphology$gII6_diag)^2)/(2*morphology$gII6_length*morphology$gII6_width)
morphology$Ialpha <- 180*(acos(morphology$Icosinus))/pi
morphology$IIalpha <- 180*(acos(morphology$IIcosinus))/pi

morphology <- as.matrix(morphology[ order(row.names(morphology)), ])


## Remove NA-s
morphology_noNA <- morphology[complete.cases(morphology), ]
morphology_gnato <- morphology[,c(1,12:25)]
morphology_gnato_noNA <- morphology_gnato[complete.cases(morphology_gnato),]

morphology_body <- morphology[,c(1:11)]
morphology_body_noNA <- morphology_body[complete.cases(morphology_body),]


# Adjust tree to data
to_keep_morpho=rownames(morphology_noNA)
tree_morpho=keep.tip(tree_mean, to_keep_morpho)
name.check(tree_morpho,morphology_noNA) 

tree_morpho_l <- ladderize(tree_morpho, right = F)

# Phylogenetically corrected GLS
phyl_gls_morpho=phyl.resid(tree_morpho, morphology_noNA[,1], morphology_noNA[,c(2:21)])
residuals_gls=as.matrix(phyl_gls_morpho$resid)
residuals_gls <- residuals_gls[ order(row.names(residuals_gls)), ]

body_length=as.matrix(morphology_noNA[,1], )
colnames(body_length)<-c("body_length")

morpho_all_gls=cbind(body_length, residuals_gls,morphology_noNA[,c(22:25)])

morpho_gnato_gls <- morpho_all_gls[,c(1,12:25)]
morpho_body_gls <- morpho_all_gls[,c(1:11)]

morpho_gnato_gls_scaled <- scale(morpho_gnato_gls)


## Sort trait data and plot it on phylogeny ####
#loop sorting
sortedData <- vector("list",ncol(morpho_all_gls))
for (x in 1:ncol(morpho_all_gls)) {
  x_mat <- as.matrix(morpho_all_gls[,x])
  names_x <- paste(colnames(morpho_all_gls)[x])
  colnames(x_mat) <- names_x
  sortedData[[x]] <- sortTraitData(phy = tree_morpho_l, y = x_mat, 
                                   data.name = names_x, pass.ultrametric = TRUE, log.trait = F)
  names(sortedData)[x] <- names_x
}
# plot all sorted traits on phylo
par(mfrow=c(5,5), mar = c(1, 0, 2, 0), oma = c(0,0, 2, 0))
for (x in 1:ncol(morpho_all_gls)) {
  sortedData_x <- sortedData[[x]]
  traitData.plot.spela(sortedData_x$trait,sortedData_x$phy,show.tips = T,cex.tips = 0.3,
                       axis.text = paste(names(sortedData)[x]),cex.plot = 1)
}

## DTT plots - all parameters, whole phylogeny ####
par(mfrow=c(5,5), mar = c(1, 1, 2, 2), oma = c(0,1, 0, 0))
DDT0=vector("list", 19)
RANK1=vector("list", 19)
for (i in 1:21) {
  DDT0[[i]]<-dtt1(tree_morpho, morpho_all_gls[,i], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1[[i]]<-rank_env_dtt(DDT0[[i]], Plot=F)
  plot(c(0,1), c(0,2), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main=colnames(morpho_all_gls)[i])
  polygon(c(RANK1[[i]]$r, rev(RANK1[[i]]$r)), c(RANK1[[i]]$upper, rev(RANK1[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1[[i]]$r, RANK1[[i]]$data_curve, lwd=1)
  lines(RANK1[[i]]$r, RANK1[[i]]$central_curve, lty=2)
  abline(v=c(0.2,0.4,0.6,0.8), col = brewer.pal(5, "Blues"))
}



## some random stuff
pheno_bodylength <- phenogram(tree_morpho, morpho_all_gls[,1], colors =c("blue"))
cont_bodylength <- contMap(tree_morpho, morpho_all_gls[,1], plot = F)
plotTree.wBars(tree_morpho_l, morpho_all_gls[,1],type="phylogram")
plot(cont_bodylength)
phenogram(tree_morpho, morpho_all_gls[,20], add=T)



par(mfrow=c(5,1))
## DTT locomotion whole genus----
d1_loco<-dtt1(tree_morpho, morpho_all_gls[,c(1,4:6)], plot=F, nsim=1000, calculateMDIp=T)	
ylim_loco<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1_loco<-rank_env_dtt(d1_loco, Plot=F)
# Plot DTT --> body whole genus
plot(c(0,1), ylim_loco[1:2], yaxp=c(ylim_loco[1],ylim_loco[2],ylim_loco[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="Locomotion - whole genus")  
x_loco<-r1_loco$r
y1_loco<-r1_loco$upper
y2_loco<-r1_loco$lower
polygon(c(r1_loco$r, rev(r1_loco$r)), c(r1_loco$upper, rev(r1_loco$lower)), col="grey60", border=NA)
lines(r1_loco$r, r1_loco$data_curve, lwd=1)
lines(r1_loco$r, r1_loco$central_curve, lty=2)

## DTT shape whole genus----
d1_shape<-dtt1(tree_morpho, morpho_all_gls[,c(7:11)], plot=F, nsim=1000, calculateMDIp=T)	
ylim_shape<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1_shape<-rank_env_dtt(d1_shape, Plot=F)
# Plot DTT --> body whole genus
plot(c(0,1), ylim_shape[1:2], yaxp=c(ylim_shape[1],ylim_shape[2],ylim_shape[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="Shape - whole genus")  
x_shape<-r1_shape$r
y1_shape<-r1_shape$upper
y2_shape<-r1_shape$lower
polygon(c(r1_shape$r, rev(r1_shape$r)), c(r1_shape$upper, rev(r1_shape$lower)), col="grey60", border=NA)
lines(r1_shape$r, r1_shape$data_curve, lwd=1)
lines(r1_shape$r, r1_shape$central_curve, lty=2)

## DTT senzoric whole genus----
d1_sensoric<-dtt1(tree_morpho, morpho_all_gls[,c(2,3)], plot=F, nsim=1000, calculateMDIp=T)	
ylim_sensoric<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1_sensoric<-rank_env_dtt(d1_sensoric, Plot=F)
# Plot DTT --> body whole genus
plot(c(0,1), ylim_sensoric[1:2], yaxp=c(ylim_sensoric[1],ylim_sensoric[2],ylim_sensoric[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="senzoric - whole genus")  
x_sensoric<-r1_sensoric$r
y1_sensoric<-r1_sensoric$upper
y2_sensoric<-r1_sensoric$lower
polygon(c(r1_sensoric$r, rev(r1_sensoric$r)), c(r1_sensoric$upper, rev(r1_sensoric$lower)), col="grey60", border=NA)
lines(r1_sensoric$r, r1_sensoric$data_curve, lwd=1)
lines(r1_sensoric$r, r1_sensoric$central_curve, lty=2)

## DTT trophic whole genus----
d1_gnatho<-dtt1(tree_morpho, morpho_all_gls[,20:21], plot=F, nsim=1000, calculateMDIp=T)	
ylim_gnatho<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1_gnatho<-rank_env_dtt(d1_gnatho, Plot=F)
# Plot DTT --> body whole genus
plot(c(0,1), ylim_gnatho[1:2], yaxp=c(ylim_gnatho[1],ylim_gnatho[2],ylim_gnatho[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="trophic - whole genus")  
x_gnatho<-r1_gnatho$r
y1_gnatho<-r1_gnatho$upper
y2_gnatho<-r1_gnatho$lower
polygon(c(r1_gnatho$r, rev(r1_gnatho$r)), c(r1_gnatho$upper, rev(r1_gnatho$lower)), col="grey60", border=NA)
lines(r1_gnatho$r, r1_gnatho$data_curve, lwd=1)
lines(r1_gnatho$r, r1_gnatho$central_curve, lty=2)

#########Analysis on subclades#########

# define subclades
#pannonian
speciesA1=c("N_ablaskiri","N_ambulator","N_aquilex_T81","N_aquilex_T94","N_bihorensis_A","N_bihorensis_B","N_carniolicus","N_cf_tauri_3","N_cvetkovi","N_daniali","N_dimorphus","N_dobati","N_fongi","N_gebhardti","N_inermis","N_pontoruffoi","N_racovitzai","N_sp_A_iz_wolfi","N_spn_Gumbrini","N_spn_Huda_Luknja","N_spn_NC561","N_spn_NC575","N_spn_NC576","N_spn_NC580","N_spn_NC586","N_spn_NC613","N_spn_NC630","N_spn_NC660","N_spn_NC702","N_spn_NC707","N_spn_NC711","N_spn_NC716","N_spn_NC718","N_spn_NC724","N_spn_NC850","N_spn_ND053","N_spn_ND113","N_spn_ND209","N_spn_ND429","N_spn_ND565","N_spn_ND770","N_spn_Podutik","N_spn_Sbirkovska_cave","N_spn_Sitarjevec","N_spn_Sukhaja_balka","N_spn_tauri_Brezno3src","N_tauri","N_vadimi","N_wolfi_A")
#pontic
speciesA2=c("C_paradoxa","N_aberrans","N_alpinus","N_andropus","N_aquilex_B","N_armatus","N_bajuvaricus","N_barbatus","N_carpathicus","N_danielopoli","N_decui","N_dissonus","N_galvagnii","N_grandii","N_italicus","N_karkabounasi","N_labacensis_A","N_labacensis_B","N_lattingerae","N_longidactylus","N_microcerberus_A","N_microcerberus_B","N_minor","N_moldavicus_cf","N_multipennatus","N_parapupetta","N_pectinicauda","N_petrosani","N_pupetta","N_serbicus","N_similis","N_sp_T110","N_spn_Arkadi","N_spn_Jelovica","N_spn_NC217","N_spn_NC624","N_spn_NC786","N_spn_ND182","N_spn_ND188","N_spn_ND211","N_spn_ND213","N_spn_ND216","N_spn_ND301","N_spn_ND426","N_spn_ND439","N_spn_ND540","N_spn_ND605","N_spn_ND669","N_spn_ND670","N_spn_ND705","N_spn_ND924","N_spn_Vodni_kevder","N_strouhali","N_tamaninii","N_transitivus","N_transsylvanicus_A","N_transsylvanicus_B","N_transsylvanicus_C","N_wolfi_B")
#south dinaric
speciesB=c("N_aulicus","N_balcanicus","N_bilecanus","N_boskovici","N_brevicuspis_A","N_brevicuspis_B","N_buturovici","N_dabarensis","N_factor","N_hercegovinensis","N_hvarensis_B","N_hvarensis_C","N_kusceri","N_lunaris_B","N_miljeticus","N_navotinus","N_podgoricensis","N_polymorphus","N_salernianus_A","N_salernianus_B","N_salernianus_C","N_salernianus_D","N_sketi","N_sp_A","N_sp_C","N_sp_D","N_sp_E","N_spn_NB613","N_spn_NC759","N_spn_NC812","N_spn_ND096","N_spn_ND193","N_spn_ND195","N_spn_ND198","N_spn_ND201","N_spn_ND204","N_spn_ND270","N_spn_ND374","N_spn_ND398","N_spn_ND463","N_spn_ND501","N_spn_ND558","N_spn_ND690","N_spn_ND739","N_spn_ND749","N_spn_ND752","N_spn_ND833","N_spn_Vilina_pecina","N_trullipes","N_vjetrenicensis_A","N_vjetrenicensis_B","N_zagorae")
#west balkan
speciesC=c("Cn_lubuskensis","N_alpheus","N_anchialinus","N_antipodes","N_arbiter","N_arethusa","N_brevirostris","N_cornicolanus","N_croaticus","N_doli","N_fjakae","N_hebereri_1","N_hebereri_2","N_hebereri_3","N_ictus","N_kolombatovici","N_liburnicus_A","N_liburnicus_B","N_longiflagellum","N_lunaris_A","N_messanai","N_mirocensis","N_orcinus","N_pachytelson","N_parenzani","N_patrizii","N_pectencoronatae","N_pincinovae","N_rejici","N_rhenorhodanensis_1","N_rhenorhodanensis_2","N_rhenorhodanensis_3","N_rhenorhodanensis_4","N_rhenorhodanensis_6","N_rhenorhodanensis_9","N_salonitanus","N_spn_Lusci_Palanka","N_spn_NB849","N_spn_ND119","N_spn_ND130","N_spn_ND139","N_spn_ND160","N_spn_ND256","N_spn_ND338","N_spn_ND469","N_spn_ND711","N_spn_ND818","N_spn_Prodisce_Sava","N_spn_Zenadija","N_spn_tertius","N_stefanellii","N_stefanellii_A","N_stefanellii_C","N_stenopus_A","N_stenopus_B","N_steueri","N_subtypicus","Nb_orophobata")
#north dinaric
speciesE=c("N_malagorae","N_spn_ND348","N_spn_ND441","N_brachytelson","N_spn_NC767","N_spn_ND302","N_illidzensis","N_sp_T262","N_spn_ND837","N_spn_ND377","N_dalmatinus","N_kapelanus","N_elegans","N_vinodolensis","N_spn_ND937","N_debilis","N_pretneri_A","N_pretneri_B","N_spn_NB236","N_spn_NC839","N_spn_NC768","N_spn_ND822","N_kordunensis","N_spn_NA057","N_spn_Vlasic","N_zagrebensis","N_chagankae","N_gottscheeanensis","N_spn_NC281","N_novomestanus","N_likanus","N_podpecanus","N_slovenicus","N_sp_B","N_spn_ND273","N_spn_ND750","N_spn_ND790","N_pedemontanus","N_rhenorhodanensis_8","N_rhenorhodanensis_7","N_spn_NC838","N_sp_T261","N_spn_NB659","N_iskae","N_spoeckeri","N_hadzii","N_cvajcki","N_spn_Iska_vas","N_spn_ND934")

# prune trees
clade_A1=keep.tip(tree_mean, speciesA1)
clade_A2=keep.tip(tree_mean, speciesA2)
clade_B=keep.tip(tree_mean, speciesB)
clade_C=keep.tip(tree_mean, speciesC)
clade_E=keep.tip(tree_mean, speciesE)

CLADES = as.list(c(clade_A1, clade_A2, clade_B, clade_C, clade_E))
names_clades = c("A1_pannonian", "A2_pontic", "B_south_dinaric","C_west_balkan", "E_north_dinaric")

# subset of morphology
MORPHO <- vector("list", length(CLADES)) 
for (i in 1:length(CLADES)) {
  MORPHO[[i]]<- subset(morphology_noNA, rownames(morphology_noNA) %in% CLADES[[i]]$tip.label)
}
MORPHO_GLS <- vector("list", length(CLADES)) 
for (i in 1:length(CLADES)) {
  MORPHO_GLS[[i]]<- subset(morpho_all_gls, rownames(morpho_all_gls) %in% CLADES[[i]]$tip.label)
}
CLADES_MORPHO=vector("list", length(CLADES)) 
for (i in 1:length(CLADES)) {
  CLADES_MORPHO[[i]]=keep.tip(CLADES[[i]], rownames(MORPHO[[i]]))
}
CLADES_MORPHO_GLS=vector("list", length(CLADES)) 
for (i in 1:length(CLADES)) {
  CLADES_MORPHO_GLS[[i]]=keep.tip(CLADES[[i]], rownames(MORPHO_GLS[[i]]))
}

# DTT for each clade
for (i in 1:length(CLADES_MORPHO)) {
  print(name.check(CLADES_MORPHO[[i]], MORPHO_GLS[[i]]))
}
## DTT plots --> shape parameters
par(mfrow=c(4,5), mar = c(1, 1, 2, 2), oma = c(3,3, 1, 1))
DDT0_body=vector("list", length(CLADES))
RANK1_body=vector("list", length(CLADES))
for (i in 1:length(CLADES)) {
  DDT0_body[[i]]<-dtt1(CLADES_MORPHO[[i]], MORPHO_GLS[[i]][,c(7:11)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_body[[i]]<-rank_env_dtt(DDT0_body[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main=names_clades[i])
  polygon(c(RANK1_body[[i]]$r, rev(RANK1_body[[i]]$r)), c(RANK1_body[[i]]$upper, rev(RANK1_body[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_body[[i]]$r, RANK1_body[[i]]$data_curve, lwd=1)
  lines(RANK1_body[[i]]$r, RANK1_body[[i]]$central_curve, lty=2)
}
mtext(text="Gnathopods                          Locomotion                          Sensoric                          Shape",side=2,line=1,outer=TRUE)

## DTT plots --> sensoric parameters
DDT0_sensoric=vector("list", length(CLADES))
RANK1_sensoric=vector("list", length(CLADES))
for (i in 1:length(CLADES)) {
  DDT0_sensoric[[i]]<-dtt1(CLADES_MORPHO[[i]], MORPHO_GLS[[i]][,c(2,3)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_sensoric[[i]]<-rank_env_dtt(DDT0_sensoric[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_sensoric[[i]]$r, rev(RANK1_sensoric[[i]]$r)), c(RANK1_sensoric[[i]]$upper, rev(RANK1_sensoric[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_sensoric[[i]]$r, RANK1_sensoric[[i]]$data_curve, lwd=1)
  lines(RANK1_sensoric[[i]]$r, RANK1_sensoric[[i]]$central_curve, lty=2)
}


## DTT plots --> locomotion parameters
DDT0_loco=vector("list", length(CLADES))
RANK1_loco=vector("list", length(CLADES))
for (i in 1:length(CLADES)) {
  DDT0_loco[[i]]<-dtt1(CLADES_MORPHO[[i]], MORPHO_GLS[[i]][,c(1,4:6)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_loco[[i]]<-rank_env_dtt(DDT0_loco[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_loco[[i]]$r, rev(RANK1_loco[[i]]$r)), c(RANK1_loco[[i]]$upper, rev(RANK1_loco[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_loco[[i]]$r, RANK1_loco[[i]]$data_curve, lwd=1)
  lines(RANK1_loco[[i]]$r, RANK1_loco[[i]]$central_curve, lty=2)
}

## DTT plots --> gnathopod parameters
DDT0_gnatho=vector("list", length(CLADES))
RANK1_gnatho=vector("list", length(CLADES))
for (i in 1:length(CLADES)) {
  DDT0_gnatho[[i]]<-dtt1(CLADES_MORPHO[[i]], MORPHO_GLS[[i]][,20:21], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim_gnatho<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_gnatho[[i]]<-rank_env_dtt(DDT0_gnatho[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim_gnatho[1],ylim_gnatho[2],ylim_gnatho[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_gnatho[[i]]$r, rev(RANK1_gnatho[[i]]$r)), c(RANK1_gnatho[[i]]$upper, rev(RANK1_gnatho[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_gnatho[[i]]$r, RANK1_gnatho[[i]]$data_curve, lwd=1)
  lines(RANK1_gnatho[[i]]$r, RANK1_gnatho[[i]]$central_curve, lty=2)
}


## Phenograms ####

PHENOGRSMS_all_traits=vector("list", length(CLADES))
fancyTree(tree_morpho_l,type="scattergram",X=morpho_all_gls[,c(1,2,7,20)], label="off", ftype="off")

par(mfrow=c(5,5), mar = c(1, 4, 1, 1), oma = c(1,0, 1, 1))
PHENOGRSMS_all_traits=vector("list", 21)
for (i in 1:21) {
  PHENOGRSMS_all_traits[[i]] <- phenogram(tree_morpho, morpho_all_gls[,i], ftype="off",
                                          ylab=colnames(morpho_all_gls)[i])
  abline(v=c(17,27,37), col = brewer.pal(5, "Blues"))
}


dev.off()


par(mfrow=c(5,5), mar = c(1, 4, 1, 1), oma = c(1,0, 1, 1))
PHENOGRSMS_all_traits=vector("list", 21)
for (i in 1:21) {
  PHENOGRSMS_all_traits[[i]] <- phenogram(CLADES_MORPHO[[5]], MORPHO_GLS[[5]][,i], ftype="off",
                                          ylab=colnames(morpho_all_gls)[i])
  abline(v=c(17,27,37), col = brewer.pal(5, "Blues"))
}
tree_morpho_fix<-di2multi(tree_morpho)
nodeA1 <- findMRCA(tree_morpho, c("N_malagorae","N_spn_ND348","N_spn_ND441","N_brachytelson","N_spn_NC767","N_spn_ND302",
                                  "N_illidzensis","N_sp_T262")
)

nodeA1=c(626)
#pontic
nodeA2=c(600)
#south dinaric
nodeB=c(467)
#west balkan
nodeC=c(425)
#north dinaric
nodeE=c(342)

tree_morpho_paint<-paintSubTree(tree_morpho,node=nodeA1,state="2",anc.state="1")
tree_morpho_paint<-paintSubTree(tree_morpho,node=nodeA2,state="3", stem = FALSE)
tree_morpho_paint<-paintSubTree(tree_morpho,node=nodeB,state="4")
tree_morpho_paint<-paintSubTree(tree_morpho,node=nodeC,state="5")
tree_morpho_paint<-paintSubTree(tree_morpho,node=nodeE,state="6")
cols<-c("black","blue","red", "yellow", "green","gray")
names(cols)<-1:6
plotSimmap(tree_morpho_paint,cols,pts=F,lwd=3,node.numbers=T)
par(mfrow=c(5,5), mar = c(1, 4, 1, 1), oma = c(1,0, 1, 1))
PHENOGRSMS_all_traits=vector("list", 21)
for (i in 1:21) {
  PHENOGRSMS_all_traits[[i]] <- phenogram(tree_morpho_paint, morpho_all_gls[,i], ftype="off",
                                          ylab=colnames(morpho_all_gls)[i], colors = cols)
}
