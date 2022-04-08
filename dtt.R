######Packages####
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


# final tree, max credibility tree, mean heights
tree_mean <- read.nexus("./data/beast_meanH.tre")
tree_mean <- drop.tip(tree_mean, "Nl_nolli")

#Morphology, continuous
morphology <- read_excel("./data/morphology.xlsx", na = "NA")
morphology <- as.data.frame(morphology)
row.names(morphology) <- morphology$species
morphology[1] <- NULL
morphology$gpI6_size <- morphology$gI6_length+morphology$gI6_width+morphology$gI6_diag
morphology$gpII6_size <- morphology$gII6_length+morphology$gII6_width+morphology$gII6_diag
morphology$Icosinus <- ((morphology$gI6_length)^2+(morphology$gI6_width)^2-(morphology$gI6_diag)^2)/(2*morphology$gI6_length*morphology$gI6_width)
morphology$IIcosinus <- ((morphology$gII6_length)^2+(morphology$gII6_width)^2-(morphology$gII6_diag)^2)/(2*morphology$gII6_length*morphology$gII6_width)
morphology$Ialpha <- 180*(acos(morphology$Icosinus))/pi
morphology$IIalpha <- 180*(acos(morphology$IIcosinus))/pi
#morphology$IIa_deg <- (180 * morphology$a_rad) / pi
#morphology[c(22,23)] <- NULL
morphology <- as.matrix(morphology)
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

# Phylogenetically corrected GLS
phyl_gls_morpho=phyl.resid(tree_morpho, morphology_noNA[,1], morphology_noNA[,c(2:21)])
rstandard(phyl_gls_morpho$resid)
residuals_gls=as.matrix(phyl_gls_morpho$resid)
residuals_gls <- residuals_gls[ order(row.names(residuals_gls)), ]
body_length=as.matrix(morphology_noNA[,1], )
body_length <- as.matrix(body_length[ order(row.names(body_length)), ])
colnames(body_length)<-c("body_length")
Icosinus <- as.matrix(morphology_noNA[,22], )
Icosinus <- as.matrix(Icosinus[ order(row.names(Icosinus)), ])
colnames(Icosinus)<-c("Icosinus")
IIcosinus <- as.matrix(morphology_noNA[,23], )
IIcosinus <- as.matrix(IIcosinus[ order(row.names(IIcosinus)), ])
colnames(IIcosinus)<-c("IIcosinus")
Ialpha <- as.matrix(morphology_noNA[,24], )
Ialpha <- as.matrix(Ialpha[ order(row.names(Ialpha)), ])
colnames(Ialpha)<-c("Ialpha")
IIalpha <- as.matrix(morphology_noNA[,25], )
IIalpha <- as.matrix(IIalpha[ order(row.names(IIalpha)), ])
colnames(IIalpha)<-c("IIalpha")
morpho_all_gls=cbind(body_length,residuals_gls,Icosinus,IIcosinus,Ialpha,IIalpha)
gpI6_size <- as.matrix(morpho_gnato_gls[,10], )
colnames(gpI6_size)<-c("gpI6_size")
gpII6_size <- as.matrix(morpho_gnato_gls[,11], )
colnames(gpII6_size)<-c("gpII6_size")
body_length <- as.matrix(morpho_gnato_gls[,1], )
colnames(body_length)<-c("body_length")

morpho_gnato_gls <- morpho_all_gls[,c(1,12:25)]
morpho_body_gls <- morpho_all_gls[,c(1:11)]


morpho_gnato_gls_scaled <- scale(morpho_gnato_gls)





## DTT GNATO SIZE1 ----
#velikost gp I (=gpI dolžina + gp I palm + gp I diagonal)
#velikost 2 (=gpII dolžina + gpII palm + gpII diagonal); meri gpI/5 in gpII/5 zavržemo.
d1gnato<-dtt1(tree_morpho, morpho_gnato_gls[,10], plot=F, nsim=1000, calculateMDIp=T)	
cat("MDI-pvalue= ", d1gnato$MDIp, "\n")
d1gnato$MDI
d1gnato$MDIpVal
ylim<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1gnato<-rank_env_dtt(d1gnato, Plot=F)
# Note that the p-value is a range because the ranks will almost always lead to some ties
r1gnato$p_interval
# Plot DTT --> gnatopod 1 size whole genus
par(mfcol=c(5,1))
plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="Gnathopod I size phyl.resid")  
xgnato<-r1gnato$r
y1gnato<-r1gnato$upper
y2gnato<-r1gnato$lower
polygon(c(xgnato, rev(xgnato)), c(y1gnato, rev(y2gnato)), col="grey60", border=NA)
lines(xgnato, r1gnato$data_curve, lwd=1)
lines(xgnato, r1gnato$central_curve, lty=2)


## DTT GNATO SIZE2 ----
#velikost 2 (=gpII dolžina + gpII palm + gpII diagonal); meri gpI/5 in gpII/5 zavržemo.
d2gnato<-dtt1(tree_morpho, morpho_gnato_gls[,11], plot=F, nsim=1000, calculateMDIp=T)	
cat("MDI-pvalue= ", d2gnato$MDIp, "\n")
d2gnato$MDI
d2gnato$MDIpVal
ylim2<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r2gnato<-rank_env_dtt(d2gnato, Plot=F)
# Note that the p-value is a range because the ranks will almost always lead to some ties
r2gnato$p_interval
# Plot DTT --> gnatopod 2 size whole genus
plot(c(0,1), ylim2[1:2], yaxp=c(ylim2[1],ylim2[2],ylim2[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="Gnathopod II size phyl.resid")  
xgnato2<-r2gnato$r
y1gnato2<-r2gnato$upper
y2gnato2<-r2gnato$lower
polygon(c(xgnato2, rev(xgnato2)), c(y1gnato2, rev(y2gnato2)), col="grey60", border=NA)
lines(xgnato2, r2gnato$data_curve, lwd=1)
lines(xgnato2, r2gnato$central_curve, lty=2)

## DTT GNATO cosinus ----
# morpho_gnato_gls_df <- as.data.frame(morpho_gnato_gls)
# morphology_noNA_df <- as.data.frame(morphology_noNA)
# ggplot(morphology_noNA_df, aes(x = body_length, y = Icosinus))+
#   geom_point()
# 
tree_morpho_L=ladderize(tree_morpho, right = FALSE)
# is_tip <- tree_morpho_L$edge[,2] <= length(tree_morpho_L$tip.label)
# ordered_tips <- tree_morpho_L$edge[is_tip, 2]
# ordered_labels <- tree_morpho_L$tip.label[ordered_tips]
# morpho_gnato_gls_ord <- morpho_gnato_gls[order(match(rownames(morpho_gnato_gls), ordered_labels)), , drop = FALSE]
# plotTree.wBars(tree_morpho_L,(morpho_gnato_gls[,12]),
#                method="plotTree", lwd=1, col="black", width = 0.7, scale = 0.4)
# plotTree.barplot(tree_morpho_L,(morpho_gnato_gls_ord[,12]),args.plotTree=list(ftype="off"))
# 
# 
# 
Ialpha
IIalpha
gpI6_size
gpII6_size
body_length
sortedData <- sortTraitData(phy = tree_morpho, y = IIalpha, 
                             data.name = "alpha", pass.ultrametric = TRUE, log.trait = F)
phy <- sortedData$phy
alpha <- sortedData$trait
traitData.plot(alpha,phy,show.tips = T,cex.tips = 0.3,
                axis.text = "gpII alpha",cex.plot = 1)

sortedDatagpI6_size <- sortTraitData(phy = tree_morpho, y = gpI6_size, 
                            data.name = "gpI6_size", pass.ultrametric = TRUE, log.trait = F)
gpI6_size <- sortedDatagpI6_size$trait
traitData.plot(gpI6_size,phy,show.tips = T,cex.tips = 0.3,
               axis.text = "gpI6_size",cex.plot = 1)

sortedDatagpII6_size <- sortTraitData(phy = tree_morpho, y = gpII6_size, 
                                     data.name = "gpII6_size", pass.ultrametric = TRUE, log.trait = F)
gpII6_size <- sortedDatagpII6_size$trait
traitData.plot(gpII6_size,phy,show.tips = T,cex.tips = 0.3,
               axis.text = "gpII6_size",cex.plot = 1)

sortedDatabody_length <- sortTraitData(phy = tree_morpho, y = body_length, 
                                      data.name = "body_length", pass.ultrametric = TRUE, log.trait = F)
body_length_s <- sortedDatabody_length$trait
traitData.plot(body_length_s,phy,show.tips = T,cex.tips = 0.3,
               axis.text = "body_length",cex.plot = 1)
# gpI6_size,gpII6_size,alpha,gpI5_length,gpII5length
d3gnato<-dtt1(tree_morpho, morpho_gnato_gls[,c(2,6,10:13)], plot=F, nsim=1000, calculateMDIp=T)	
cat("MDI-pvalue= ", d3gnato$MDIp, "\n")
d3gnato$MDI
d3gnato$MDIpVal
ylim3<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r3gnato<-rank_env_dtt(d3gnato, Plot=F)
# Note that the p-value is a range because the ranks will almost always lead to some ties
r3gnato$p_interval
# Plot DTT --> gnatopod 2 size whole genus
plot(c(0,1), ylim3[1:2], yaxp=c(ylim3[1],ylim3[2],ylim3[3]), type="n",
     xlab="relative time", frame.plot=F, ylab="disparity",
     main="gpI,gpII size, gpI5,gpII5 length: phyl.resid; gpI, gpII cosinus")  
xgnato3<-r3gnato$r
y1gnato3<-r3gnato$upper
y2gnato3<-r3gnato$lower
polygon(c(xgnato3, rev(xgnato3)), c(y1gnato3, rev(y2gnato3)), col="grey60", border=NA)
lines(xgnato3, r3gnato$data_curve, lwd=1)
lines(xgnato3, r3gnato$central_curve, lty=2)

## DTT gnathopod alpha only ----
d4gnato<-dtt1(tree_morpho, morpho_gnato_gls[,c(12,13)], plot=F, nsim=1000, calculateMDIp=T)	
cat("MDI-pvalue= ", d4gnato$MDIp, "\n")
d4gnato$MDI
d4gnato$MDIpVal
ylim4<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r4gnato<-rank_env_dtt(d4gnato, Plot=F)
# Note that the p-value is a range because the ranks will almost always lead to some ties
r4gnato$p_interval
# Plot DTT --> gnatopod 2 size whole genus
plot(c(0,1), c(0,7), yaxp=c(ylim4[1],ylim4[2],ylim4[3]), type="n",
     xlab="relative time", frame.plot=F, ylab="disparity",
     main="gpI cosinus, gpII cosinus")  
xgnato4<-r4gnato$r
y1gnato4<-r4gnato$upper
y2gnato4<-r4gnato$lower
polygon(c(xgnato4, rev(xgnato4)), c(y1gnato4, rev(y2gnato4)), col="grey60", border=NA)
lines(xgnato4, r4gnato$data_curve, lwd=1)
lines(xgnato4, r4gnato$central_curve, lty=2)



## DTT BODY ----
d1body<-dtt1(tree_morpho, morpho_body_gls, plot=F, nsim=1000, calculateMDIp=T)	
cat("MDI-pvalue= ", d1body$MDIp, "\n")
d1body$MDI
d1body$MDIpVal
ylim<-par("yaxp")
# Compute the rank envelope after Murrell (2018)
r1body<-rank_env_dtt(d1body, Plot=F)
# Note that the p-value is a range because the ranks will almost always lead to some ties
r1body$p_interval
# Plot DTT --> body whole genus
plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="Body all traits, phyl.resid")  
xbody<-r1body$r
y1body<-r1body$upper
y2body<-r1body$lower
polygon(c(xbody, rev(xbody)), c(y1body, rev(y2body)), col="grey60", border=NA)
lines(xbody, r1body$data_curve, lwd=1)
lines(xbody, r1body$central_curve, lty=2)

## DTT plots (only Murell (2018) variant shown) --> gnathopod parameters
par(mfcol=c(5,4))
DDT0=vector("list", 19)
RANK1=vector("list", 19)
for (i in 1:19) {
  DDT0[[i]]<-dtt1(tree_morpho, morpho_all_gls[,i], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1[[i]]<-rank_env_dtt(DDT0[[i]], Plot=F)
  plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main=colnames(morpho_all_gls)[i])
  polygon(c(RANK1[[i]]$r, rev(RANK1[[i]]$r)), c(RANK1[[i]]$upper, rev(RANK1[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1[[i]]$r, RANK1[[i]]$data_curve, lwd=1)
  lines(RANK1[[i]]$r, RANK1[[i]]$central_curve, lty=2)
}

#########Analysis on subclades#########

# define subclades
speciesA1=c("N_ablaskiri","N_ambulator","N_aquilex_T81","N_aquilex_T94","N_bihorensis_A_H","N_bihorensis_B_H","N_carniolicus","N_cf_tauri_3","N_cvetkovi","N_daniali","N_dimorphus","N_dobati","N_fongi","N_gebhardti_H","N_inermis","N_sp_A__iz_wolfi__H","N_spn_Brezno3src","N_spn_Gumbrini","N_spn_Huda_Luknja","N_spn_Podutik","N_spn_Sbirkovska_cave","N_spn_Sitarjevec","N_spn_Sukhaja_balka_H","N_tauri","N_vadimi","N_wolfi_A","PN_racovitzai","PN_ruffoi")
speciesA2=c("C_paradoxa","N_aberrans","N_alpinus_H","N_andropus","N_aquilex_B","N_armatus","N_bajuvaricus","N_barbatus","N_carpathicus","N_danielopoli","N_decui","N_galvagnii","N_grandii","N_italicus","N_karkabounasi_H","N_labacensis_A","N_labacensis_B","N_lattingerae","N_longidactylus","N_microcerberus_A","N_microcerberus_B","N_minor","N_moldavicus_cf","N_multipennatus","N_parapupetta","N_pectinicauda_H","N_petrosani","N_pupetta","N_serbicus","N_similis","N_sp_T110","N_spn_Arkadi","N_spn_Jelovica","N_spn_Vodni_kevder","N_strouhali","N_tamaninii","N_transitivus","N_transsylvanicus_A","N_transsylvanicus_B","N_transsylvanicus_C","N_wolfi_B")
speciesB=c("N_bilecanus","N_vjetrenicensis_A","N_trullipes","N_kusceri_H","N_hercegovinensis","N_vjetrenicensis_B","N_balcanicus_H","N_polymorphus","N_dabarensis","N_spn_Vilina_pecina","N_podgoricensis","N_brevicuspis","N_factor","N_buturovici","N_aulicus","N_lunaris","N_hvarensis_C","N_miljeticus","N_hvarensis_B","N_boskovici","N_salernianus_A","N_salernianus_B","N_salernianus_C","N_salernianus_D","N_sp_A","N_sp_C","N_sp_E","N_zagorae","N_sp_D")
speciesC=c("N_hebereri_1","N_hebereri_2","N_hebereri_3","N_spn_Zenadija","N_cornicolanus","N_messanai","N_stefanellii_A","N_parenzani","N_stefanellii_C","N_ictus","N_stefanellii","N_rejici","N_antipodes","N_arbiter_H","N_doli_H","N_fjakae_H","N_pincinovae","N_alpheus_H","N_arethusa","N_anhihalinus_H","N_salonitanus_H","N_longiflagellum","N_pachytelson_H","N_pectencoronatae","N_spn_Lusci_Palanka","N_spn_Suhaca_jama","N_stenopus_A","N_spn_tertius","N_brevirostris","N_stenopus_B_H","N_mirocensis","N_spn_Prodisce_Sava","N_kolombatovici_H","N_steueri_H","N_liburnicus_B","N_patrizii","N_croaticus_H","N_subtypicus_H","Nb_orophobata_H","N_orcinus","Nb_sp")
speciesD=c("N_longicaudatus_T39","N_sodalis","N_longicaudatus_T40","N_longicaudatus_T42","N_longicaudatus_T43","N_longicaudatus_T38","N_longicaudatus_T36","N_longicaudatus_T37","N_longicaudatus_T45","N_longicaudatus_T47","N_longicaudatus_T46","N_pasquinii_H","N_sibillinianus","N_longicaudatus_T51","N_longicaudatus_B","N_longicaudatus_C","N_corsicanus","N_longicaudatus_T58","N_longicaudatus_T54","N_longicaudatus_T55","N_longicaudatus_T61","N_longicaudatus_T62","N_frasassianus","N_longicaudatus_T65","N_longicaudatus_T50","N_longicaudatus_T64","N_cvijici","N_aitolosi_H","N_longicaudatus_T32","N_versluysi","N_longicaudatus_T27","N_longicaudatus_A","N_longicaudatus_E","N_longicaudatus_F","N_longicaudatus_D","N_radzai","N_sp_T163","N_sp_T164")
speciesE=c("N_brachytelson","N_chagankae","N_cvajcki_H","N_dalmatinus","N_elegans","N_gottscheeanensis_H","N_hadzii","N_illidzensis","N_iskae","N_kapelanus","N_kordunensis","N_likanus_H","N_malagorae_H","N_novomestanus","N_podpecanus_H","N_pretneri_A","N_pretneri_B","N_slovenicus","N_sp_B","N_sp_T262","N_spn_Iska_vas","N_spn_Vlasic","N_spoeckeri","N_vinodolensis","N_zagrebensis")

# prune trees
clade_A1=keep.tip(tree_mean, speciesA1)
clade_A2=keep.tip(tree_mean, speciesA2)
clade_B=keep.tip(tree_mean, speciesB)
clade_C=keep.tip(tree_mean, speciesC)
clade_D=keep.tip(tree_mean, speciesD)
clade_E=keep.tip(tree_mean, speciesE)

CLADES = as.list(c(clade_A1, clade_A2, clade_B, clade_C, clade_D, clade_E))
names_clades = c("clade_A1", "clade_A2", "clade_B","clade_C", "clade_D", "clade_E")

# subset of morphology
MORPHO <- vector("list", 6) 
for (i in 1:length(CLADES)) {
  MORPHO[[i]]<- subset(morphology_noNA, rownames(morphology_noNA) %in% CLADES[[i]]$tip.label)
}
MORPHO_GLS <- vector("list", 6) 
for (i in 1:length(CLADES)) {
  MORPHO_GLS[[i]]<- subset(morpho_all_gls, rownames(morpho_all_gls) %in% CLADES[[i]]$tip.label)
}
CLADES_MORPHO=vector("list", 6) 
for (i in 1:length(CLADES)) {
  CLADES_MORPHO[[i]]=keep.tip(CLADES[[i]], rownames(MORPHO[[i]]))
}
CLADES_MORPHO_GLS=vector("list", 6) 
for (i in 1:length(CLADES)) {
  CLADES_MORPHO_GLS[[i]]=keep.tip(CLADES[[i]], rownames(MORPHO_GLS[[i]]))
}

# DTT for each clade
for (i in 1:length(CLADES_MORPHO)) {
  print(name.check(CLADES_MORPHO[[i]], MORPHO[[i]]))
}
## DTT plots (only Murell (2018) variant shown) --> body parameters
par(mfcol=c(2,6))
DDT0=vector("list", 6)
RANK1=vector("list", 6)
for (i in 1:6) {
  DDT0[[i]]<-dtt1(CLADES_MORPHO[[i]], MORPHO[[i]][,1:11], plot=T, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1[[i]]<-rank_env_dtt(DDT0[[i]], Plot=F)
  plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main=names_clades[i])
  polygon(c(RANK1[[i]]$r, rev(RANK1[[i]]$r)), c(RANK1[[i]]$upper, rev(RANK1[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1[[i]]$r, RANK1[[i]]$data_curve, lwd=1)
  lines(RANK1[[i]]$r, RANK1[[i]]$central_curve, lty=2)
}

## DTT plots (only Murell (2018) variant shown) --> gnathopod parameters
par(mfcol=c(1,6))
DDT0=vector("list", 6)
RANK1=vector("list", 6)
for (i in 1:6) {
  DDT0[[i]]<-dtt1(CLADES_MORPHO[[i]], MORPHO[[i]][,12:21], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
  ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1[[i]]<-rank_env_dtt(DDT0[[i]], Plot=F)
  plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1[[i]]$r, rev(RANK1[[i]]$r)), c(RANK1[[i]]$upper, rev(RANK1[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1[[i]]$r, RANK1[[i]]$data_curve, lwd=1)
  lines(RANK1[[i]]$r, RANK1[[i]]$central_curve, lty=2)
}



## DTT statistics --> Table 4
for (i in 1:6) {
  print(DDT0[[i]]$MDI)
  print(DDT0[[i]]$MDIpVal)
  print(RANK1[[i]]$p_interval)
}

