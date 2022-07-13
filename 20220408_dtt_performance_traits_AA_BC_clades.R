#pannonian
speciesA=c("N_ablaskiri","N_ambulator","N_aquilex_T81","N_aquilex_T94","N_bihorensis_A","N_bihorensis_B","N_carniolicus","N_cf_tauri_3","N_cvetkovi","N_daniali","N_dimorphus","N_dobati","N_fongi","N_gebhardti","N_inermis","N_pontoruffoi","N_racovitzai","N_sp_A_iz_wolfi","N_spn_Gumbrini","N_spn_Huda_Luknja","N_spn_NC561","N_spn_NC575","N_spn_NC576","N_spn_NC580","N_spn_NC586","N_spn_NC613","N_spn_NC630","N_spn_NC660","N_spn_NC702","N_spn_NC707","N_spn_NC711","N_spn_NC716","N_spn_NC718","N_spn_NC724","N_spn_NC850","N_spn_ND053","N_spn_ND113","N_spn_ND209","N_spn_ND429","N_spn_ND565","N_spn_ND770","N_spn_Podutik","N_spn_Sbirkovska_cave","N_spn_Sitarjevec","N_spn_Sukhaja_balka","N_spn_tauri_Brezno3src","N_tauri","N_vadimi","N_wolfi_A","C_paradoxa","N_aberrans","N_alpinus","N_andropus","N_aquilex_B","N_armatus","N_bajuvaricus","N_barbatus","N_carpathicus","N_danielopoli","N_decui","N_dissonus","N_galvagnii","N_grandii","N_italicus","N_karkabounasi","N_labacensis_A","N_labacensis_B","N_lattingerae","N_longidactylus","N_microcerberus_A","N_microcerberus_B","N_minor","N_moldavicus_cf","N_multipennatus","N_parapupetta","N_pectinicauda","N_petrosani","N_pupetta","N_serbicus","N_similis","N_sp_T110","N_spn_Arkadi","N_spn_Jelovica","N_spn_NC217","N_spn_NC624","N_spn_NC786","N_spn_ND182","N_spn_ND188","N_spn_ND211","N_spn_ND213","N_spn_ND216","N_spn_ND301","N_spn_ND426","N_spn_ND439","N_spn_ND540","N_spn_ND605","N_spn_ND669","N_spn_ND670","N_spn_ND705","N_spn_ND924","N_spn_Vodni_kevder","N_strouhali","N_tamaninii","N_transitivus","N_transsylvanicus_A","N_transsylvanicus_B","N_transsylvanicus_C","N_wolfi_B")
#south dinaric
speciesBC=c("N_aulicus","N_balcanicus","N_bilecanus","N_boskovici","N_brevicuspis_A","N_brevicuspis_B","N_buturovici","N_dabarensis","N_factor","N_hercegovinensis","N_hvarensis_B","N_hvarensis_C","N_kusceri","N_lunaris_B","N_miljeticus","N_navotinus","N_podgoricensis","N_polymorphus","N_salernianus_A","N_salernianus_B","N_salernianus_C","N_salernianus_D","N_sketi","N_sp_A","N_sp_C","N_sp_D","N_sp_E","N_spn_NB613","N_spn_NC759","N_spn_NC812","N_spn_ND096","N_spn_ND193","N_spn_ND195","N_spn_ND198","N_spn_ND201","N_spn_ND204","N_spn_ND270","N_spn_ND374","N_spn_ND398","N_spn_ND463","N_spn_ND501","N_spn_ND558","N_spn_ND690","N_spn_ND739","N_spn_ND749","N_spn_ND752","N_spn_ND833","N_spn_Vilina_pecina","N_trullipes","N_vjetrenicensis_A","N_vjetrenicensis_B","N_zagorae","Cn_lubuskensis","N_alpheus","N_anchialinus","N_antipodes","N_arbiter","N_arethusa","N_brevirostris","N_cornicolanus","N_croaticus","N_doli","N_fjakae","N_hebereri_1","N_hebereri_2","N_hebereri_3","N_ictus","N_kolombatovici","N_liburnicus_A","N_liburnicus_B","N_longiflagellum","N_lunaris_A","N_messanai","N_mirocensis","N_orcinus","N_pachytelson","N_parenzani","N_patrizii","N_pectencoronatae","N_pincinovae","N_rejici","N_rhenorhodanensis_1","N_rhenorhodanensis_2","N_rhenorhodanensis_3","N_rhenorhodanensis_4","N_rhenorhodanensis_6","N_rhenorhodanensis_9","N_salonitanus","N_spn_Lusci_Palanka","N_spn_NB849","N_spn_ND119","N_spn_ND130","N_spn_ND139","N_spn_ND160","N_spn_ND256","N_spn_ND338","N_spn_ND469","N_spn_ND711","N_spn_ND818","N_spn_Prodisce_Sava","N_spn_Zenadija","N_spn_tertius","N_stefanellii","N_stefanellii_A","N_stefanellii_C","N_stenopus_A","N_stenopus_B","N_steueri","N_subtypicus","Nb_orophobata")

# prune trees
clade_A=keep.tip(tree_mean, speciesA)
clade_BC=keep.tip(tree_mean, speciesBC)

CLADES_double = as.list(c(clade_A, clade_BC))
names_clades_double = c("A_pannonian__pontic", "BC_south_dinaric_west_balkan")

# subset of morphology
MORPHO_GLS_double <- vector("list", length(CLADES_double)) 
for (i in 1:length(CLADES_double)) {
  MORPHO_GLS_double[[i]]<- subset(morpho_all_gls, rownames(morpho_all_gls) %in% CLADES_double[[i]]$tip.label)
}
CLADES_MORPHO_double=vector("list", length(CLADES_double)) 
for (i in 1:length(CLADES_double)) {
  CLADES_MORPHO_double[[i]]=keep.tip(CLADES_double[[i]], rownames(MORPHO_GLS_double[[i]]))
}

# DTT for each clade
for (i in 1:length(CLADES_MORPHO_double)) {
  print(name.check(CLADES_MORPHO_double[[i]], MORPHO_GLS_double[[i]]))
}
## DTT plots --> shape parameters
par(mfrow=c(4,2), mar = c(1, 1, 2, 2), oma = c(3,3, 1, 1))
DDT0_body_double=vector("list", length(CLADES_double))
RANK1_body_double=vector("list", length(CLADES_double))
for (i in 1:length(CLADES_double)) {
   DDT0_body_double[[i]]<-dtt1(CLADES_MORPHO_double[[i]], MORPHO_GLS_double[[i]][,c(7:11)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
   ylim<-par("yaxp")
  #Compute the rank envelope after Murell
   RANK1_body_double[[i]]<-rank_env_dtt(DDT0_body_double[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main=names_clades_double[i])
  polygon(c(RANK1_body_double[[i]]$r, rev(RANK1_body_double[[i]]$r)), c(RANK1_body_double[[i]]$upper, rev(RANK1_body_double[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_body_double[[i]]$r, RANK1_body_double[[i]]$data_curve, lwd=1)
  lines(RANK1_body_double[[i]]$r, RANK1_body_double[[i]]$central_curve, lty=2)
}
mtext(text="Gnathopods                          Locomotion                          Sensoric                          Shape",side=2,line=1,outer=TRUE)

## DTT plots --> sensoric parameters
DDT0_loco_double=vector("list", length(CLADES_double))
RANK1_loco_double=vector("list", length(CLADES_double))
for (i in 1:length(CLADES_double)) {
   DDT0_loco_double[[i]]<-dtt1(CLADES_MORPHO_double[[i]], MORPHO_GLS_double[[i]][,c(2,3)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
   ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_loco_double[[i]]<-rank_env_dtt(DDT0_loco_double[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_loco_double[[i]]$r, rev(RANK1_loco_double[[i]]$r)), c(RANK1_loco_double[[i]]$upper, rev(RANK1_loco_double[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_loco_double[[i]]$r, RANK1_loco_double[[i]]$data_curve, lwd=1)
  lines(RANK1_loco_double[[i]]$r, RANK1_loco_double[[i]]$central_curve, lty=2)
}


## DTT plots --> locomotion parameters
DDT0_loco_double=vector("list", length(CLADES_double))
RANK1_loco_double=vector("list", length(CLADES_double))
for (i in 1:length(CLADES_double)) {
   DDT0_loco_double[[i]]<-dtt1(CLADES_MORPHO_double[[i]], MORPHO_GLS_double[[i]][,c(1,4:6)], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
   ylim<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_loco_double[[i]]<-rank_env_dtt(DDT0_loco_double[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_loco_double[[i]]$r, rev(RANK1_loco_double[[i]]$r)), c(RANK1_loco_double[[i]]$upper, rev(RANK1_loco_double[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_loco_double[[i]]$r, RANK1_loco_double[[i]]$data_curve, lwd=1)
  lines(RANK1_loco_double[[i]]$r, RANK1_loco_double[[i]]$central_curve, lty=2)
}

## DTT plots --> gnathopod parameters
DDT0_gnatho_double=vector("list", length(CLADES_double))
RANK1_gnatho_double=vector("list", length(CLADES_double))
for (i in 1:length(CLADES_double)) {
   DDT0_gnatho_double[[i]]<-dtt1(CLADES_MORPHO_double[[i]], MORPHO_GLS_double[[i]][,20:21], plot=F, nsim=1000, calculateMDIp=T, Ylim = c(0,2.5))
   ylim_gnatho<-par("yaxp")
  #Compute the rank envelope after Murell
  RANK1_gnatho_double[[i]]<-rank_env_dtt(DDT0_gnatho_double[[i]], Plot=F)
  plot(c(0,1), c(0,3), yaxp=c(ylim_gnatho[1],ylim_gnatho[2],ylim_gnatho[3]), type="n", xlab="relative time", frame.plot=F, ylab="disparity", main="")
  polygon(c(RANK1_gnatho_double[[i]]$r, rev(RANK1_gnatho_double[[i]]$r)), c(RANK1_gnatho_double[[i]]$upper, rev(RANK1_gnatho_double[[i]]$lower)), col="grey60", border=NA)
  lines(RANK1_gnatho_double[[i]]$r, RANK1_gnatho_double[[i]]$data_curve, lwd=1)
  lines(RANK1_gnatho_double[[i]]$r, RANK1_gnatho_double[[i]]$central_curve, lty=2)
}


par(mfrow=c(5,5), mar = c(1, 4, 1, 1), oma = c(1,0, 1, 1))
PHENOGRSMS_all_traits=vector("list", 21)
for (i in 1:21) {
  PHENOGRSMS_all_traits[[i]] <- phenogram(CLADES_MORPHO_double[[2]], MORPHO_GLS_double[[2]][,i], ftype="off",
                                          ylab=colnames(morpho_all_gls)[i])
}