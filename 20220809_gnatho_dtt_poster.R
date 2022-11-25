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
library(parallel)


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
tree_mean_nv <- drop.tip(tree_mean, "N_vadimi")

# Morphology, continuous
morphology <- read_excel("./data/morphology.xlsx", na = "NA")
morphology <- as.data.frame(morphology)
row.names(morphology) <- morphology$species
morphology[1] <- NULL
## calculate gnathopod size and angles
morphology$gpI6_size <- morphology$gI6_length+morphology$gI6_width+morphology$gI6_diag
morphology$gpII6_size <- morphology$gII6_length+morphology$gII6_width+morphology$gII6_diag
morphology <- as.matrix(morphology[ order(row.names(morphology)), ])

## Remove NA-s and N_vadimi
morphology_noNA <- morphology[complete.cases(morphology), ]
morphology_noNA_nv <- morphology_noNA[row.names(morphology_noNA) != "N_vadimi",]
# Adjust tree to data
to_keep_morpho_nv=rownames(morphology_noNA_nv)
tree_morpho_nv=keep.tip(tree_mean_nv, to_keep_morpho_nv)
name.check(tree_morpho_nv,morphology_noNA_nv) 

# Phylogenetically corrected GLS
phyl_gls_morpho_nv=phyl.resid(tree_morpho_nv, morphology_noNA_nv[,1], morphology_noNA_nv[,c(2:21)])
residuals_gls_nv=as.matrix(phyl_gls_morpho_nv$resid)
residuals_gls_nv <- residuals_gls_nv[ order(row.names(residuals_gls_nv)), ]

body_length_nv=as.matrix(morphology_noNA_nv[,1], )
colnames(body_length_nv)<-c("body_length")

morpho_all_gls_nv=cbind(body_length_nv, residuals_gls_nv)

par(mfrow=c(3,1))
## DTT body whole genus----
d1_loco<-dtt1(tree_morpho_nv, morpho_all_gls_nv[,c(1,4:11)], plot=F, nsim=1000, calculateMDIp=T)	
ylim_loco<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1_loco<-rank_env_dtt(d1_loco, Plot=F)
# Plot DTT --> body whole genus
plot(c(0,1), ylim_loco[1:2], yaxp=c(ylim_loco[1],ylim_loco[2],ylim_loco[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="body - whole genus")  
x_loco<-r1_loco$r
y1_loco<-r1_loco$upper
y2_loco<-r1_loco$lower
polygon(c(r1_loco$r, rev(r1_loco$r)), c(r1_loco$upper, rev(r1_loco$lower)), col="grey60", border=NA)
lines(r1_loco$r, r1_loco$data_curve, lwd=1, col = "#009E73")
lines(r1_loco$r, r1_loco$central_curve, lty=2, col = "#009E73")


## DTT senzoric whole genus----
d1_sensoric<-dtt1(tree_morpho_nv, morpho_all_gls_nv[,c(2,3)], plot=F, nsim=1000, calculateMDIp=T)	
ylim_sensoric<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1_sensoric<-rank_env_dtt(d1_sensoric, Plot=F)
# Plot DTT --> body whole genus
plot(c(0,1), ylim_sensoric[1:2], yaxp=c(ylim_sensoric[1],ylim_sensoric[2],ylim_sensoric[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="senzoric - whole genus")  
x_sensoric<-r1_sensoric$r
y1_sensoric<-r1_sensoric$upper
y2_sensoric<-r1_sensoric$lower
polygon(c(r1_sensoric$r, rev(r1_sensoric$r)), c(r1_sensoric$upper, rev(r1_sensoric$lower)), col = "grey", border=NA)
lines(r1_sensoric$r, r1_sensoric$data_curve, lwd=1, col = "#0072B2")
lines(r1_sensoric$r, r1_sensoric$central_curve, lty=2, col = "#0072B2")

## DTT trophic whole genus----
d1_gnatho<-dtt1(tree_morpho_nv, morpho_all_gls_nv[,20:21], plot=F, nsim=1000, calculateMDIp=T)	
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
lines(r1_gnatho$r, r1_gnatho$data_curve, lwd=1, col = "#D55E00")
lines(r1_gnatho$r, r1_gnatho$central_curve, lty=2, col = "#D55E00")

#########Analysis on subclades#########

# define subclades
#pannonian
speciesA1_nv=c("N_ablaskiri","N_ambulator","N_aquilex_T81","N_aquilex_T94","N_bihorensis_A","N_bihorensis_B","N_carniolicus","N_cf_tauri_3","N_cvetkovi","N_daniali","N_dimorphus","N_dobati","N_fongi","N_gebhardti","N_inermis","N_pontoruffoi","N_racovitzai","N_sp_A_iz_wolfi","N_spn_Gumbrini","N_spn_Huda_Luknja","N_spn_NC561","N_spn_NC575","N_spn_NC576","N_spn_NC580","N_spn_NC586","N_spn_NC613","N_spn_NC630","N_spn_NC660","N_spn_NC702","N_spn_NC707","N_spn_NC711","N_spn_NC716","N_spn_NC718","N_spn_NC724","N_spn_NC850","N_spn_ND053","N_spn_ND113","N_spn_ND209","N_spn_ND429","N_spn_ND565","N_spn_ND770","N_spn_Podutik","N_spn_Sbirkovska_cave","N_spn_Sitarjevec","N_spn_Sukhaja_balka","N_spn_tauri_Brezno3src","N_tauri","N_wolfi_A")
#pontic
speciesA2=c("C_paradoxa","N_aberrans","N_alpinus","N_andropus","N_aquilex_B","N_armatus","N_bajuvaricus","N_barbatus","N_carpathicus","N_danielopoli","N_decui","N_dissonus","N_galvagnii","N_grandii","N_italicus","N_karkabounasi","N_labacensis_A","N_labacensis_B","N_lattingerae","N_longidactylus","N_microcerberus_A","N_microcerberus_B","N_minor","N_moldavicus_cf","N_multipennatus","N_parapupetta","N_pectinicauda","N_petrosani","N_pupetta","N_serbicus","N_similis","N_sp_T110","N_spn_Arkadi","N_spn_Jelovica","N_spn_NC217","N_spn_NC624","N_spn_NC786","N_spn_ND182","N_spn_ND188","N_spn_ND211","N_spn_ND213","N_spn_ND216","N_spn_ND301","N_spn_ND426","N_spn_ND439","N_spn_ND540","N_spn_ND605","N_spn_ND669","N_spn_ND670","N_spn_ND705","N_spn_ND924","N_spn_Vodni_kevder","N_strouhali","N_tamaninii","N_transitivus","N_transsylvanicus_A","N_transsylvanicus_B","N_transsylvanicus_C","N_wolfi_B")
#south dinaric
speciesB=c("N_aulicus","N_balcanicus","N_bilecanus","N_boskovici","N_brevicuspis_A","N_brevicuspis_B","N_buturovici","N_dabarensis","N_factor","N_hercegovinensis","N_hvarensis_B","N_hvarensis_C","N_kusceri","N_lunaris_B","N_miljeticus","N_navotinus","N_podgoricensis","N_polymorphus","N_salernianus_A","N_salernianus_B","N_salernianus_C","N_salernianus_D","N_sketi","N_sp_A","N_sp_C","N_sp_D","N_sp_E","N_spn_NB613","N_spn_NC759","N_spn_NC812","N_spn_ND096","N_spn_ND193","N_spn_ND195","N_spn_ND198","N_spn_ND201","N_spn_ND204","N_spn_ND270","N_spn_ND374","N_spn_ND398","N_spn_ND463","N_spn_ND501","N_spn_ND558","N_spn_ND690","N_spn_ND739","N_spn_ND749","N_spn_ND752","N_spn_ND833","N_spn_Vilina_pecina","N_trullipes","N_vjetrenicensis_A","N_vjetrenicensis_B","N_zagorae")
#west balkan
speciesC=c("Cn_lubuskensis","N_alpheus","N_anchialinus","N_antipodes","N_arbiter","N_arethusa","N_brevirostris","N_cornicolanus","N_croaticus","N_doli","N_fjakae","N_hebereri_1","N_hebereri_2","N_hebereri_3","N_ictus","N_kolombatovici","N_liburnicus_A","N_liburnicus_B","N_longiflagellum","N_lunaris_A","N_messanai","N_mirocensis","N_orcinus","N_pachytelson","N_parenzani","N_patrizii","N_pectencoronatae","N_pincinovae","N_rejici","N_rhenorhodanensis_1","N_rhenorhodanensis_2","N_rhenorhodanensis_3","N_rhenorhodanensis_4","N_rhenorhodanensis_6","N_rhenorhodanensis_9","N_salonitanus","N_spn_Lusci_Palanka","N_spn_NB849","N_spn_ND119","N_spn_ND130","N_spn_ND139","N_spn_ND160","N_spn_ND256","N_spn_ND338","N_spn_ND469","N_spn_ND711","N_spn_ND818","N_spn_Prodisce_Sava","N_spn_Zenadija","N_spn_tertius","N_stefanellii","N_stefanellii_A","N_stefanellii_C","N_stenopus_A","N_stenopus_B","N_steueri","N_subtypicus","Nb_orophobata")

# prune trees
clade_A1_nv=keep.tip(tree_mean_nv, speciesA1_nv)
clade_A2_nv=keep.tip(tree_mean_nv, speciesA2)
clade_B_nv=keep.tip(tree_mean_nv, speciesB)
clade_C_nv=keep.tip(tree_mean_nv, speciesC)

CLADES_nv = as.list(c(clade_A1_nv, clade_A2_nv, clade_B_nv, clade_C_nv))
names_clades = c("A1_pannonian", "A2_pontic", "B_south_dinaric","C_west_balkan")
names(CLADES_nv) <- names_clades

# subset of morphology
MORPHO_GLS_nv <- vector("list", length(CLADES_nv)) 
for (i in 1:length(CLADES)) {
  MORPHO_GLS_nv[[i]]<- subset(morpho_all_gls_nv, rownames(morpho_all_gls_nv) %in% CLADES_nv[[i]]$tip.label)
}
CLADES_MORPHO_nv=vector("list", length(CLADES_nv)) 
for (i in 1:length(CLADES_nv)) {
  CLADES_MORPHO_nv[[i]]=keep.tip(CLADES_nv[[i]], rownames(MORPHO_GLS_nv[[i]]))
}
names(CLADES_MORPHO_nv) <- names_clades

# DTT for each clade
for (i in 1:length(CLADES_MORPHO_nv)) {
  print(name.check(CLADES_MORPHO_nv[[i]], MORPHO_GLS_nv[[i]]))
}
## DTT plots --> 3 trait groups
dev.off()
dev.new()
par(mfrow=c(3,4), mar = c(1, 1, 2, 2), oma = c(3,3, 1, 1))

## DTT plots --> body parameters
DDT0_loco=vector("list", length(CLADES_nv))
RANK1_loco=vector("list", length(CLADES_nv))
for (i in 1:length(CLADES_nv)) {
  DDT0_loco[[i]]<-dtt1(CLADES_MORPHO_nv[[i]], MORPHO_GLS_nv[[i]][,c(1,4:11)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_loco[[i]]<-rank_env_dtt(DDT0_loco[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main=names(CLADES_MORPHO_nv)[i])
  polygon(c(RANK1_loco[[i]]$r, rev(RANK1_loco[[i]]$r)), c(RANK1_loco[[i]]$upper, rev(RANK1_loco[[i]]$lower)), col = "#00d99e", border=NA)
  lines(RANK1_loco[[i]]$r, RANK1_loco[[i]]$data_curve, lwd=1, col = "#009E73")
  lines(RANK1_loco[[i]]$r, RANK1_loco[[i]]$central_curve, lty=2, col = "#009E73")
}
mtext(text="Gnathopods                      Sensoric                  Body",side=2,line=1,outer=TRUE)

## DTT plots --> sensoric parameters
DDT0_sensoric=vector("list", length(CLADES_nv))
RANK1_sensoric=vector("list", length(CLADES_nv))
YLIM_sensoric=vector("list", length(CLADES_nv))
for (i in 1:length(CLADES_nv)) {
  DDT0_sensoric[[i]]<-dtt1(CLADES_MORPHO_nv[[i]], MORPHO_GLS_nv[[i]][,c(2,3)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  YLIM_sensoric[[i]]<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_sensoric[[i]]<-rank_env_dtt(DDT0_sensoric[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(YLIM_sensoric[[i]][1],YLIM_sensoric[[i]][2],YLIM_sensoric[[i]][3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_sensoric[[i]]$r, rev(RANK1_sensoric[[i]]$r)), c(RANK1_sensoric[[i]]$upper, rev(RANK1_sensoric[[i]]$lower)), border=NA, col = rgb(0/255, 114/255, 178/255, alpha = 0.1))
  lines(RANK1_sensoric[[i]]$r, RANK1_sensoric[[i]]$data_curve, lwd=1, col="#0072B2")
  lines(RANK1_sensoric[[i]]$r, RANK1_sensoric[[i]]$central_curve, lty=2, col = "#0072B2")
}

## DTT plots --> trophic parameters
DDT0_gnatho=vector("list", length(CLADES_nv))
RANK1_gnatho=vector("list", length(CLADES_nv))
for (i in 1:length(CLADES_nv)) {
  DDT0_gnatho[[i]]<-dtt1(CLADES_MORPHO_nv[[i]], MORPHO_GLS_nv[[i]][,20:21], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim_gnatho<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_gnatho[[i]]<-rank_env_dtt(DDT0_gnatho[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim_gnatho[1],ylim_gnatho[2],ylim_gnatho[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_gnatho[[i]]$r, rev(RANK1_gnatho[[i]]$r)), c(RANK1_gnatho[[i]]$upper, rev(RANK1_gnatho[[i]]$lower)), col=rgb(213/255, 94/255, 0/255, alpha = 0.1), border=NA)
  lines(RANK1_gnatho[[i]]$r, RANK1_gnatho[[i]]$data_curve, lwd=1, col = "#D55E00")
  lines(RANK1_gnatho[[i]]$r, RANK1_gnatho[[i]]$central_curve, lty=2, col = "#D55E00")
}


