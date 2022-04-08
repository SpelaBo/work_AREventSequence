
par(mfrow=c(3,1), mar = c(1, 1, 2, 2), oma = c(0,1, 0, 0))
#body shape: cxII_length,cxIII_length,pV2_width,pVI2_width,pVII2_width
#morpho_all_gls[,c(7:11)]
dtt_body_shape<-dtt1(tree_morpho, morpho_all_gls[,c(7:11)], plot=F, nsim=1000, calculateMDIp=T)	
ylim<-par("yaxp")
r1_body_shape<-rank_env_dtt(dtt_body_shape, Plot=F)
plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="body shape")  
x_body_shape<-r1_body_shape$r
y1_body_shape<-r1_body_shape$upper
y2_body_shape<-r1_body_shape$lower
polygon(c(x_body_shape, rev(x_body_shape)), c(y1_body_shape, rev(y2_body_shape)), col="grey60", border=NA)
lines(x_body_shape, r1_body_shape$data_curve, lwd=1)
lines(x_body_shape, r1_body_shape$central_curve, lty=2)

#locomotion: body_length, pV_length, pVI_length, pVII_length
#morpho_all_gls[,c(1,4:6)]
dtt_locomotion<-dtt1(tree_morpho, morpho_all_gls[,c(1,4:6)], plot=F, nsim=1000, calculateMDIp=T)	
ylim<-par("yaxp")
r1_locomotion<-rank_env_dtt(dtt_locomotion, Plot=F)
plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="locomotion")  
x_locomotion<-r1_locomotion$r
y1_locomotion<-r1_locomotion$upper
y2_locomotion<-r1_locomotion$lower
polygon(c(x_locomotion, rev(x_locomotion)), c(y1_locomotion, rev(y2_locomotion)), col="grey60", border=NA)
lines(x_locomotion, r1_locomotion$data_curve, lwd=1)
lines(x_locomotion, r1_locomotion$central_curve, lty=2)

#trophic: gpI6_size,gpII6_size
#morpho_all_gls[,c(20:21)]
dtt_trophic<-dtt1(tree_morpho, morpho_all_gls[,c(20:21)], plot=F, nsim=1000, calculateMDIp=T)	
ylim<-par("yaxp")
r1_trophic<-rank_env_dtt(dtt_trophic, Plot=F)
plot(c(0,1), ylim[1:2], yaxp=c(ylim[1],ylim[2],ylim[3]), type="n", xlab="relative time",
     frame.plot=F, ylab="disparity", main="trophic")  
x_trophic<-r1_trophic$r
y1_trophic<-r1_trophic$upper
y2_trophic<-r1_trophic$lower
polygon(c(x_trophic, rev(x_trophic)), c(y1_trophic, rev(y2_trophic)), col="grey60", border=NA)
lines(x_trophic, r1_trophic$data_curve, lwd=1)
lines(x_trophic, r1_trophic$central_curve, lty=2)
