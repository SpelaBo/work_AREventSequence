library(BAMMtools)
library(coda)
# tree
to_keep_morpho=rownames(morpho_all_gls)
tree_bamm=keep.tip(tree_mean, to_keep_morpho)
name.check(tree_morpho,morpho_all_gls) 

is.ultrametric(tree_bamm)
# if not ultrmetric: tree_bamm <- force.ultrametric(tree_bamm)
is.binary.tree(tree_bamm)
sum(tree_bamm$edge.length < 0)
sum(tree_bamm$edge.length == 0)

# for plotting for poster - exclude vadimi
tree_bamm_novadimi=drop.tip(tree_bamm, "N_vadimi")
write.tree(tree_bamm,"E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/tree.tre")
write.tree(tree_bamm_novadimi,"E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/tree.tre")

#subset of species with and without vadimi in A1 clade
speciesA1_bamm <- subset(speciesA1, speciesA1 %in% tree_bamm$tip.label)
speciesA2_bamm <- subset(speciesA2, speciesA2 %in% tree_bamm$tip.label)
speciesC_bamm <- subset(speciesC, speciesC %in% tree_bamm$tip.label)
speciesB_bamm <- subset(speciesB, speciesB %in% tree_bamm$tip.label)
speciesA1_novadimi_bamm <- subset(speciesA1_novadimi, speciesA1_novadimi %in% tree_bamm$tip.label)

#mrca without vadimi - note that node numbers are not the same as with it!!
mrcaA1_novadimi <- findMRCA(tree_bamm_novadimi, speciesA1_novadimi_bamm)
mrcaA2_novadimi <- findMRCA(tree_bamm_novadimi, speciesA2_bamm)
mrcaC_novadimi <- findMRCA(tree_bamm_novadimi, speciesC_bamm)
mrcaB_novadimi <- findMRCA(tree_bamm_novadimi, speciesB_bamm)

#mrca with vadimi - note that node numbers are not the same as without it!!
mrcaA1 <- findMRCA(tree_bamm, speciesA1_bamm)
mrcaA2 <- findMRCA(tree_bamm, speciesA2_bamm)
mrcaC <- findMRCA(tree_bamm, speciesC_bamm)
mrcaB <- findMRCA(tree_bamm, speciesB_bamm)

#check tree
plotTree(tree_bamm_novadimi, type = "fan",show.tip.label = F, ftype="off")
nodes<-labelnodes(text=c("A1","A2","B","C"),node=c(mrcaA1_novadimi,mrcaA2_novadimi,mrcaB_novadimi,mrcaC_novadimi), shape="ellipse",cex=0.8,interactive=FALSE)
# morpho PCA ordered as tree labels
morpho_all_gls<-morpho_all_gls[tree_bamm$tip.label,]
## check
all(rownames(morpho_all_gls)==tree_bamm$tip.label)
name.check(tree_bamm, morpho_all_gls)

#create data matrices for each set of traits (traits in columns)
int_body <- as.matrix(morpho_all_gls[,c(1,4:11)]) #body
int_sens <- as.matrix(morpho_all_gls[,c(2:3)]) #sensoric
int_gnatho <- as.matrix(morpho_all_gls[,c(20:21)]) #trophic

#run PCAs
int_body_pca <- phyl.pca(tree = tree_mean, int_body, method = "BM", mode = "cov")
int_sens_pca <- phyl.pca(tree = tree_mean, int_sens, method = "BM", mode = "cov")
int_gnatho_pca <- phyl.pca(tree = tree_mean, int_gnatho, method = "BM", mode = "cov")

int_pc1 <- as.data.frame(cbind(int_sens_pca$S[,1], int_body_pca$S[,1], int_gnatho_pca$S[,1]))
colnames(int_pc1) <- c("PC1_sens", "PC1_body", "PC1_gnatho")

# note that I deleted vadimi in txt files and not here - see separate folder in bamm folder
write.table(int_sens_pca$S[,1], file="PC1_sens.txt", quote = F, append = FALSE,
            sep = "\t", dec = ".", eol = "\n",row.names = T, col.names = F)
write.table(int_body_pca$S[,1], file="PC1_body.txt", quote = F, append = FALSE,
            sep = "\t", dec = ".", eol = "\n",row.names = T, col.names = F)
write.table(int_gnatho_pca$S[,1], file="PC1_gnatho.txt", quote = F, append = FALSE,
            sep = "\t", dec = ".", eol = "\n",row.names = T, col.names = F)


setBAMMpriors(phy = tree_morpho, traits = "PC1_gnatho.txt", outfile = "gnatho_priors.txt")
setBAMMpriors(phy = tree_morpho, traits = "PC1_body.txt", outfile = "body_priors.txt")
setBAMMpriors(phy = tree_morpho, traits = "PC1_sens.txt", outfile = "sens_priors.txt")

###
tree <- read.tree("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/tree.tre")
gnatho_edata <- getEventData(tree, eventdata = "E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/gnatho_event_data.txt", burnin=0.25, type = "trait")
gnatho_shift_probs <- summary(gnatho_edata)

gnatho_mcmcout <- read.csv("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/gnatho_mcmc_out.txt", header=T)
plot(gnatho_mcmcout$logLik ~ gnatho_mcmcout$generation)

gnatho_burnstart <- floor(0.25 * nrow(gnatho_mcmcout))
gnatho_postburn <- gnatho_mcmcout[gnatho_burnstart:nrow(gnatho_mcmcout), ]

effectiveSize(gnatho_postburn$N_shifts)
effectiveSize(gnatho_postburn$logLik)


gnatho_best <- getBestShiftConfiguration(gnatho_edata, expectedNumberOfShifts=1)
plot.bammdata(gnatho_best, lwd = 2, legend = T)
addBAMMshifts(gnatho_best, cex=1)
nodes<-labelnodes(text=c("A1","A2","B","C"),node=c(mrcaA1,mrcaA2,mrcaB,mrcaC),
                  shape="ellipse",cex=0.8,interactive=FALSE)

plot.new()
par(mfrow=c(2,3))
st <- max(branching.times(tree))
plotRateThroughTime(gnatho_edata, intervalCol="#D55E00", avgCol="#D55E00", start.time=st, ylim=c(0,1), cex.axis=1, cex.lab = 1)
text(x=40, y= 0.8, label="Niphargus", font=1, cex=1.5, pos=4)
plotRateThroughTime(gnatho_edata, intervalCol="#009E73", avgCol="#009E73", start.time=st, node=mrcaA1, ylim=c(0,1),cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="A1", font=1, cex=1.5, pos=4)
plotRateThroughTime(gnatho_edata, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, node=mrcaA2, ylim=c(0,1), cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="A2", font=1, cex=1.5, pos=4)
plotRateThroughTime(gnatho_edata, intervalCol="purple", avgCol="purple", start.time=st, node=mrcaB, ylim=c(0,1),cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="B", font=1, cex=1.5, pos=4)
plotRateThroughTime(gnatho_edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st, node=mrcaC, ylim=c(0,1), cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="C", font=1, cex=1.5, pos=4)

## BAMM results - body
body_edata <- getEventData(tree, eventdata = "E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/body_event_data.txt", burnin=0.25, type = "trait")
body_shift_probs <- summary(body_edata)

body_mcmcout <- read.csv("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/body_mcmc_out.txt", header=T)
plot(body_mcmcout$logLik ~ body_mcmcout$generation)

body_burnstart <- floor(0.25 * nrow(body_mcmcout))
body_postburn <- body_mcmcout[body_burnstart:nrow(body_mcmcout), ]

effectiveSize(body_postburn$N_shifts)
effectiveSize(body_postburn$logLik)

body_best <- getBestShiftConfiguration(body_edata, expectedNumberOfShifts=1)
plot.bammdata(body_best, lwd = 2, legend = T)
addBAMMshifts(body_best, cex=1)
nodes<-labelnodes(text=c("A1","A2","B","C"),node=c(mrcaA1,mrcaA2,mrcaB,mrcaC),
                  shape="ellipse",cex=0.8,interactive=FALSE)

plot.new()
par(mfrow=c(2,3))
st <- max(branching.times(tree))
plotRateThroughTime(body_edata, intervalCol="#D55E00", avgCol="#D55E00", start.time=st, ylim=c(0,20), cex.axis=1, cex.lab = 1)
text(x=40, y= 0.8, label="Niphargus", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#009E73", avgCol="#009E73", start.time=st, node=mrcaA1, ylim=c(0,20),cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="A1", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, node=mrcaA2, ylim=c(0,20), cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="A2", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="purple", avgCol="purple", start.time=st, node=mrcaB, ylim=c(0,20),cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="B", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st, node=mrcaC, ylim=c(0,20), cex.axis=1, cex.lab = 1)
text(x=10, y= 0.8, label="C", font=1, cex=1.5, pos=4)


## BAMM results - sens
sens_edata <- getEventData(tree, eventdata = "E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/sens_event_data.txt", burnin=0.25, type = "trait")
sens_shift_probs <- summary(sens_edata)

sens_mcmcout <- read.csv("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/sens_mcmc_out.txt", header=T)
plot(sens_mcmcout$logLik ~ sens_mcmcout$generation)

sens_burnstart <- floor(0.25 * nrow(sens_mcmcout))
sens_postburn <- sens_mcmcout[sens_burnstart:nrow(sens_mcmcout), ]

effectiveSize(sens_postburn$N_shifts)
effectiveSize(sens_postburn$logLik)

sens_best <- getBestShiftConfiguration(sens_edata, expectedNumberOfShifts=1)
plot.bammdata(sens_best, lwd = 2, legend = T)
addBAMMshifts(sens_best, cex=1)
nodes<-labelnodes(text=c("A1","A2","B","C"),node=c(mrcaA1,mrcaA2,mrcaB,mrcaC),
                  shape="ellipse",cex=0.8,interactive=FALSE)

plot.new()
par(mfrow=c(2,3))
st <- max(branching.times(tree))
plotRateThroughTime(sens_edata, intervalCol="#D55E00", avgCol="#D55E00", start.time=st, ylim=c(0,5), cex.axis=1, cex.lab = 1)
text(x=40, y= 5, label="Niphargus", font=1, cex=1.5, pos=4)
plotRateThroughTime(sens_edata, intervalCol="#009E73", avgCol="#009E73", start.time=st, node=mrcaA1, ylim=c(0,5),cex.axis=1, cex.lab = 1)
text(x=10, y= 5, label="A1", font=1, cex=1.5, pos=4)
plotRateThroughTime(sens_edata, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, node=mrcaA2, ylim=c(0,5), cex.axis=1, cex.lab = 1)
text(x=10, y= 5, label="A2", font=1, cex=1.5, pos=4)
plotRateThroughTime(sens_edata, intervalCol="purple", avgCol="purple", start.time=st, node=mrcaB, ylim=c(0,5),cex.axis=1, cex.lab = 1)
text(x=10, y= 5, label="B", font=1, cex=1.5, pos=4)
plotRateThroughTime(sens_edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st, node=mrcaC, ylim=c(0,5), cex.axis=1, cex.lab = 1)
text(x=10, y= 5, label="C", font=1, cex=1.5, pos=4)





## BAMM results, without vadimi
tree_novadimi <- read.tree("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/tree.tre")

body_novadimi_edata <- getEventData(tree_novadimi, eventdata = "E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/body_novadimi_event_data.txt", burnin=0.25, type = "trait")
gnatho_novadimi_edata <- getEventData(tree_novadimi, eventdata = "E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/gnatho_novadimi_event_data.txt", burnin=0.25, type = "trait")
sens_novadimi_edata <- getEventData(tree_novadimi, eventdata = "E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/sens_novadimi_event_data.txt", burnin=0.25, type = "trait")

# check convergence and effective size
body_shift_probs_nv <- summary(body_novadimi_edata)
gnatho_shift_probs_nv <- summary(gnatho_novadimi_edata)
sens_shift_probs_nv <- summary(sens_novadimi_edata)

body_mcmcout_nv <- read.csv("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/body_novadimi_mcmc_out.txt", header=T)
plot(body_mcmcout_nv$logLik ~ body_mcmcout_nv$generation)
gnatho_mcmcout_nv <- read.csv("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/gnatho_novadimi_mcmc_out.txt", header=T)
plot(gnatho_mcmcout_nv$logLik ~ gnatho_mcmcout_nv$generation)
sens_mcmcout_nv <- read.csv("E:/sluzba/phylo_programi/bamm-2.5.0-Windows/2022_niphargus/no_vadimi/sens_novadimi_mcmc_out.txt", header=T)
plot(sens_mcmcout_nv$logLik ~ sens_mcmcout_nv$generation)

body_burnstart_nv <- floor(0.25 * nrow(body_mcmcout_nv))
body_postburn_nv <- body_mcmcout_nv[body_burnstart_nv:nrow(body_mcmcout_nv), ]
effectiveSize(body_postburn_nv$N_shifts)
effectiveSize(body_postburn_nv$logLik)

gnatho_burnstart_nv <- floor(0.25 * nrow(gnatho_mcmcout_nv))
gnatho_postburn_nv <- gnatho_mcmcout_nv[gnatho_burnstart_nv:nrow(gnatho_mcmcout_nv), ]
effectiveSize(gnatho_postburn_nv$N_shifts)
effectiveSize(gnatho_postburn_nv$logLik)

sens_burnstart_nv <- floor(0.25 * nrow(sens_mcmcout_nv))
sens_postburn_nv <- sens_mcmcout_nv[sens_burnstart_nv:nrow(sens_mcmcout_nv), ]
effectiveSize(sens_postburn_nv$N_shifts)
effectiveSize(sens_postburn_nv$logLik)


## combined rates per clade, without vadimi
dev.off()
par(mfrow=c(3,4))
plotRateThroughTime(body_novadimi_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA1_novadimi, ylim=c(0,10), cex.axis=1, cex.lab = 1)
text(x=15, y= 8, label="A1) Pannonian", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_novadimi_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA2_novadimi, ylim=c(0,10), cex.axis=1, cex.lab = 1)
text(x=15, y= 8, label="A2) Pontic", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_novadimi_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaB_novadimi, ylim=c(0,10), cex.axis=1, cex.lab = 1)
text(x=12, y= 8, label="B) S Dinaric", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_novadimi_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaC_novadimi, ylim=c(0,10), cex.axis=1, cex.lab = 1)
text(x=14, y= 8, label="C) W Balkan", font=1, cex=1.5, pos=4)

plotRateThroughTime(sens_novadimi_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaA1_novadimi, ylim=c(0,4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_novadimi_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaA2_novadimi, ylim=c(0,4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_novadimi_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaB_novadimi, ylim=c(0,4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_novadimi_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaC_novadimi, ylim=c(0,4), cex.axis=1, cex.lab = 1)

plotRateThroughTime(gnatho_novadimi_edata,  node=mrcaA1_novadimi, intervalCol="#D55E00", avgCol="#D55E00", ylim=c(0,0.4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(gnatho_novadimi_edata,  node=mrcaA2_novadimi, intervalCol="#D55E00", avgCol="#D55E00", ylim=c(0,0.4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(gnatho_novadimi_edata,  node=mrcaB_novadimi, intervalCol="#D55E00", avgCol="#D55E00", ylim=c(0,0.4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(gnatho_novadimi_edata,  node=mrcaC_novadimi, intervalCol="#D55E00", avgCol="#D55E00", ylim=c(0,0.4), cex.axis=1, cex.lab = 1)

legend('topleft', legend = c("gnatho", "body", "sens"), col = c("#D55E00", "#009E73","#0072B2"),
       fill = c("#D55E00", "#009E73", "#0072B2"), border = FALSE, lty = 1, lwd = 2, merge = TRUE,
       seg.len=0.6)








dev.off()
dev.new()
par(mfrow=c(3,1))
plotRateThroughTime(body_novadimi_edata, intervalCol="#009E73", avgCol="#009E73", ylim=c(0,10), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_novadimi_edata, intervalCol="#0072B2", avgCol="#0072B2", ylim=c(0,4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(gnatho_novadimi_edata,intervalCol="#D55E00", avgCol="#D55E00",ylim=c(0,0.4), cex.axis=1, cex.lab = 1)
text(x=15, y= 0.36, label="Niphargus", font=1, cex=1.5, pos=4)
legend('topleft', legend = c("gnatho", "body", "sens"), col = c("#D55E00", "#009E73","#0072B2"),
       fill = c("#D55E00", "#009E73", "#0072B2"), border = FALSE, lty = 1, lwd = 2, merge = TRUE,
       seg.len=0.6)
































## combined rates
par(mfrow=c(3,1))
plotRateThroughTime(gnatho_edata, intervalCol="#D55E00", avgCol="#D55E00", start.time=30, ylim=c(0,0.6), cex.axis=1, cex.lab = 1)
text(x=25, y= 0.5, label="PC1 gnatho", font=1, cex=1.5, pos=4)
plotRateThroughTime(gnatho_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA1,add = T)
plotRateThroughTime(gnatho_edata, intervalCol="darkgreen", avgCol="darkgreen", node=mrcaA2, add = T)
plotRateThroughTime(gnatho_edata, intervalCol="purple", avgCol="purple",  node=mrcaB, add = T)
plotRateThroughTime(gnatho_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaC,add = T)
legend('topleft', legend = c("Niphargus","A1) Pannonian", "A2) Pontic", "B) S Dinaric","C) W Balkan"), col = c("#D55E00", "#009E73", "darkgreen","purple","#0072B2"),
       fill = c("#D55E00", "#009E73", "darkgreen","purple","#0072B2"), border = FALSE, lty = 1, lwd = 2, merge = TRUE,
       seg.len=0.6)

plotRateThroughTime(body_edata, intervalCol="#D55E00", avgCol="#D55E00", start.time=30, ylim=c(0,10), cex.axis=1, cex.lab = 1)
text(x=25, y= 8, label="PC1 body", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA1, add = T)
plotRateThroughTime(body_edata, intervalCol="darkgreen", avgCol="darkgreen", node=mrcaA2, add = T)
plotRateThroughTime(body_edata, intervalCol="purple", avgCol="purple", node=mrcaB, add = T)
plotRateThroughTime(body_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaC, add = T)

plotRateThroughTime(sens_edata, intervalCol="#D55E00", avgCol="#D55E00", start.time=30, ylim=c(0,4), cex.axis=1, cex.lab = 1)
text(x=25, y= 3, label="PC1 sens", font=1, cex=1.5, pos=4)
plotRateThroughTime(sens_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA1, add = T)
plotRateThroughTime(sens_edata, intervalCol="darkgreen", avgCol="darkgreen", node=mrcaA2, add = T)
plotRateThroughTime(sens_edata, intervalCol="purple", avgCol="purple", node=mrcaB, add = T)
plotRateThroughTime(sens_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaC, add = T)

## combined rates per clade
dev.off()
par(mfrow=c(4,3))
plotRateThroughTime(gnatho_edata,  node=mrcaA1, intervalCol="#D55E00", avgCol="#D55E00", start.time=15, ylim=c(0,0.2), cex.axis=1, cex.lab = 1)
text(x=15, y= 0.18, label="A1) Pannonian", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA1, start.time=15, ylim=c(0,4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaA1, start.time=15, ylim=c(0,0.6), cex.axis=1, cex.lab = 1)
legend('topleft', legend = c("gnatho", "body", "sens"), col = c("#D55E00", "#009E73","#0072B2"),
       fill = c("#D55E00", "#009E73", "#0072B2"), border = FALSE, lty = 1, lwd = 2, merge = TRUE,
       seg.len=0.6)

plotRateThroughTime(gnatho_edata,  node=mrcaA2, intervalCol="#D55E00", avgCol="#D55E00", start.time=15, ylim=c(0,0.2), cex.axis=1, cex.lab = 1)
text(x=15, y= 0.18, label="A2) Pontic", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaA2, start.time=15, ylim=c(0,4), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaA2, start.time=15, ylim=c(0,0.6), cex.axis=1, cex.lab = 1)

plotRateThroughTime(gnatho_edata,  node=mrcaB, intervalCol="#D55E00", avgCol="#D55E00", start.time=15, ylim=c(0,0.2), cex.axis=1, cex.lab = 1)
text(x=15, y= 0.18, label="B) S Dinaric", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaB,start.time=15, ylim=c(0,10), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaB, start.time=15, ylim=c(0,4), cex.axis=1, cex.lab = 1)

plotRateThroughTime(gnatho_edata,  node=mrcaC, intervalCol="#D55E00", avgCol="#D55E00", start.time=15, ylim=c(0,0.4), cex.axis=1, cex.lab = 1)
text(x=15, y= 0.36, label="C) W Balkan", font=1, cex=1.5, pos=4)
plotRateThroughTime(body_edata, intervalCol="#009E73", avgCol="#009E73", node=mrcaC, start.time=15, ylim=c(0,10), cex.axis=1, cex.lab = 1)
plotRateThroughTime(sens_edata, intervalCol="#0072B2", avgCol="#0072B2", node=mrcaC, start.time=15, ylim=c(0,4), cex.axis=1, cex.lab = 1)