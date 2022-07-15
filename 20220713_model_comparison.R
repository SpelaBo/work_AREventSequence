### analiza diverzifikacije, test 1 ####

#### Data ####
# final tree, max credibility tree, mean heights
tree_mean <- read.nexus("./data/20220325_niphargus_calibrated_edit.tree")
tree_mean <- drop.tip(tree_mean, "Nl_nolli")
is.ultrametric(tree_mean)
tree_mean <- force.ultrametric(tree_mean)

# Morphology, continuous
morphology <- read_excel("./data/morphology.xlsx", na = "NA")
morphology <- as.data.frame(morphology)
row.names(morphology) <- morphology$species
morphology[1] <- NULL
## calculate gnathopod size and angles
morphology$gpI6_size <- morphology$gI6_length+morphology$gI6_width+morphology$gI6_diag
morphology$gpII6_size <- morphology$gII6_length+morphology$gII6_width+morphology$gII6_diag
morphology <- as.matrix(morphology[ order(row.names(morphology)), ])

## Remove NA-s
morphology_noNA <- morphology[complete.cases(morphology), ]

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

morpho_all_gls=cbind(body_length, residuals_gls)

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
MORPHO_GLS <- vector("list", length(CLADES)) 
for (i in 1:length(CLADES)) {
  MORPHO_GLS[[i]]<- subset(morpho_all_gls, rownames(morpho_all_gls) %in% CLADES[[i]]$tip.label)
}
CLADES_MORPHO_GLS=vector("list", length(CLADES)) 
for (i in 1:length(CLADES)) {
  CLADES_MORPHO_GLS[[i]]=keep.tip(CLADES[[i]], rownames(MORPHO_GLS[[i]]))
}

# testing BM, EB and OU model for each of the trait groups, clade B #########
# performance traits: 
## trophic niche: gpI6_size in gpII6_size; MORPHO_GLS[[i]][,20:21]
## body shape: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width; MORPHO_GLS[[i]][,c(7:11)]
## locomotion: body_length, pV_length, pVI_length, pVII_length; MORPHO_GLS[[i]][,c(1,4:6)]
## sensoric: antenna I and II; MORPHO_GLS[[i]][,c(2:3)]

# preparing the data
## tree
tree_b <- CLADES_MORPHO_GLS[[3]]
is.ultrametric(tree_b)
## distance to the root (will need for reconstructions)
dist_to_root_b <- node.depth.edgelength(tree_b)
## morpho data
data_b <- MORPHO_GLS[[3]]
## sort as node labels
data_b<-data_b[tree_b$tip.label,]
## check
all(rownames(data_b)==tree_b$tip.label)
name.check(tree_b, data_b)

# models trophic niche ####
model_BM_gnato_b <- mvBM(tree_b,data_b[,20:21])
model_EB_gnato_b <- mvEB(tree_b,data_b[,20:21])
model_OU_gnato_b <- mvOU(tree_b,data_b[,20:21])
# merge into list
models_gnato_b <- list(model_BM_gnato_b,model_EB_gnato_b,model_OU_gnato_b)
# choose the best model
gnato_b_aic <- aicw(models_gnato_b)
gnato_b_aic

# models body shape ####
model_BM_body_b <- mvBM(tree_b,data_b[,c(7:11)])
model_EB_body_b <- mvEB(tree_b,data_b[,c(7:11)])
model_OU_body_b <- mvOU(tree_b,data_b[,c(7:11)])
# merge into list
models_body_b <- list(model_BM_body_b,model_EB_body_b,model_OU_body_b)
# choose the best model
body_b_aic <- aicw(models_body_b)
body_b_aic

# models locomotion ####
model_BM_loco_b <- mvBM(tree_b,data_b[,c(1,4:6)])
model_EB_loco_b <- mvEB(tree_b,data_b[,c(1,4:6)])
model_OU_loco_b <- mvOU(tree_b,data_b[,c(1,4:6)])
# merge into list
models_loco_b <- list(model_BM_loco_b,model_EB_loco_b,model_OU_loco_b)
# choose the best model
loco_b_aic <- aicw(models_loco_b)
loco_b_aic

# models sensoric ####
model_BM_sens_b <- mvBM(tree_b,data_b[,c(2:3)])
model_EB_sens_b <- mvEB(tree_b,data_b[,c(2:3)])
model_OU_sens_b <- mvOU(tree_b,data_b[,c(2:3)])
# merge into list
models_sens_b <- list(model_BM_sens_b,model_EB_sens_b,model_OU_sens_b)
# choose the best model
sens_b_aic <- aicw(models_sens_b)
sens_b_aic

# ancestral states ####
## We want the ancestral states values at each nodes:
ASR_gnato_b<-estim(tree_b, data_b[,20:21], model_OU_gnato_b, asr=TRUE)
comb_gnato_b <- rbind(data_b[,20:21],ASR_gnato_b$estim)

ASR_body_b<-estim(tree_b, data_b[,c(7:11)], model_OU_body_b, asr=TRUE)
comb_body_b <- rbind(data_b[,c(7:11)],ASR_body_b$estim)

ASR_loco_b<-estim(tree_b, data_b[,c(1,4:6)], model_OU_loco_b, asr=TRUE)
comb_loco_b <- rbind(data_b[,c(1,4:6)],ASR_loco_b$estim)

ASR_sens_b<-estim(tree_b, data_b[,c(2:3)], model_OU_sens_b, asr=TRUE)
comb_sens_b <- rbind(data_b[,c(2:3)],ASR_sens_b$estim)

comb_b <- cbind(comb_body_b,comb_loco_b,comb_sens_b,comb_gnato_b)

# plot reconstructed phenograms ####
par(mfrow=c(4,4))
for (i in 2:ncol(comb_b)) {
  phenogram(tree=tree_b,x=comb_b[,i], ftype="off", ylim=c(min(comb_b[,i]),max(comb_b[,i])), xlab = colnames(comb_b)[i])
}

# calculate occupied morphospace ####






## ggtree ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
library(ggtree)

ggtree(tree_b, aes(color=trait), continuous = 'colour', yscale = "trait") + 
  scale_color_viridis_c() + theme_minimal()



GNATO_GLS <- vector("list", length(MORPHO_GLS)) 
for (i in 1:length(MORPHO_GLS)) {
  GNATO_GLS[[i]]<- MORPHO_GLS[[i]][,20:21]
}
gnato_gls_B <- GNATO_GLS[[3]]


## drevesa: zbrana v listu CLADES, posamično pa: clade_A1, clade_A2, clade_B, clade_C, clade_E


## test modelov za shape/locomotion/senzoric/trophic traits posebej/ali za vsak trait posebej?
# klad B




## We want the ancestral states values at each nodes:
ASR_gnato_b<-estim(tree_b, data_gnato_b, model_OU_gnato_b, asr=TRUE)
ASR_gnato_b$estim
comb_data_gnato_b <- rbind(data_gnato_b,ASR_gnato_b$estim)
# Check the ancestral states
## default gpI6_size
phenogram(tree=tree_b,x=data_gnato_b[,1], ftype="off", ylim=c(-0.2,2))
## best model gpI6_size  
pheno_gpI6_size <- phenogram(tree=tree_b,x=comb_data_gnato_b[,1], ftype="off", ylim=c(-0.2,2))
## default gpII6_size
phenogram(tree=tree_b,x=data_gnato_b[,2], ftype="off", ylim=c(-0.4,2))
## best model gpII6_size
phenogram(tree=tree_b,x=comb_data_gnato_b[,2], ftype="off", ylim=c(-0.4,2))


diff(range(comb_data_gnato_b[,1]))

# Fitting the models
tree <- CLADES_MORPHO_GLS[[3]]
data <- MORPHO_GLS[[3]][,c(1:11,20:21)]
name.check(tree, data)
# BM1 - (Equal rate matrix)
model_1<-mvBM(tree, data, model="BM1", diagnostic=FALSE, echo=FALSE)
model_1$theta

# BMM - (Proportional rate matrices)
model_2<-mvBM(tree, data, param=list(constraint="diagonal"), diagnostic=FALSE, echo=FALSE)
# BMM - (Shared eigenvectors between rate matrices)
model_3<-mvBM(tree, data, param=list(constraint="equaldiagonal"), diagnostic=FALSE, echo=FALSE)
# BMM - (Similar correlations between rate matrices)
model_4<-mvBM(tree, data, param=list(constraint="equal"), diagnostic=FALSE, echo=FALSE)

# Compare the models with AIC
AIC(model_1)
AIC(model_2)
AIC(model_3)
AIC(model_4)
AIC(model_5)
AIC(model_6)
# Test significance with LRT
LRT(model_6,model_5)
LRT(model_6,model_4)
LRT(model_6,model_3)
LRT(model_6,model_2)
LRT(model_6,model_1)