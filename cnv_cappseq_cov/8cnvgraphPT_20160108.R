################################################################################
#Setting my working directory
################################################################################
wd <- '/Volumes/morinlab/home/sepideh/exome_coverage/exome_coverage/modifiedcoding'
if(!is.null(wd)) {setwd(wd) }
################################################################################
# two-variable pointplot with fine control for cappseq and exome coverage ratios across different patients colored based on Titan calls
################################################################################
pointplot2pdf <- function(cappseq.df, exome_2.df, colvec_titan,colvec_titan_short, colvec_seq){
xmax=nrow(cappseq.df)
rangeDf <- range(cappseq.df[(cappseq.df != Inf) & !is.na(cappseq.df)])
print(cappseq.df)
ymin=min(rangeDf)
ymax=max(rangeDf)
par(mar=c(4.3,5,1.5,8), font.axis=1, xpd=TRUE)
plot(cappseq.df[,1], col=colvec_titan, main=colnames(cappseq.df), alpha=80, pch=4, cex=1, bty="n", xaxt="n", xlab="", ylab="", las=1, mgp=c(0,0.5,0), tck=-0.02, ylim=c(-6,6), xlim=c(1.5,39))#, panel.first=grid(lty=1, nx=NA, ny=NULL, lwd=0.1)) #control the size of graph by ylim and xlim
points(cappseq.df[,1], col=colvec_titan_short, alpha=100, pch=20, cex= 2, add=TRUE)
points(exome_2.df[,1], col=colvec_titan, alpha=80, pch=4, cex=1, add=TRUE)
points(exome_2.df[,1], col=colvec_seq, alpha=80, pch=21, cex= 1.3, add=TRUE)
axis(side=1, line=0, labels=rownames(cappseq.df), at=1:xmax, las=3, cex.axis=0.7, mgp=c(3,0.6,0), tck=-0.02)#mgp the second one change the placement of the samples codes
for(i in seq(from=1, to=xmax, by=1)){ lines(c(i,i),c(-6,6), lwd=0.1, col="grey") }
for(i in seq(from=-6, to=6, by=2)){ lines(c(0.35,xmax+0.65), c(i,i), lwd=0.1, col="grey") }
points(cappseq.df[,1], col=colvec_titan, alpha=80, pch=4, cex=2, add=TRUE, na.col=NA)
points(cappseq.df[,1], col=colvec_titan_short, alpha=100, pch=20, cex= 2, add=TRUE)
points(exome_2.df[,1], col=colvec_titan, alpha=80, pch=4, cex=1, add=TRUE)
points(exome_2.df[,1], col=colvec_seq, alpha=80, pch=21, cex=1.3, add=TRUE, na.col=NA)
mtext("ZCoverage_Ratio", side=2, line=1.8, adj=0.5, cex=1, col="black")
legend(0.5,-7.9, inset=.02, title= "", c("CappseqCov.TitanFocal   X CappseqCov.300Kb-ExpandedTitanSegments     O Log2ExomeCov.Sequenza   x Log2ExomeCov.Titan"), pch=20, col="black", cex=0.8, pt.cex=2, box.lwd=0, bty="n", box.col="gray", horiz = F)
legend(42,6.3,inset=.02, title="CNS 2=HT",c('0','1','2','2NLOH','3','4','5','NA'),pch=15, col=c("cornflowerblue","deepskyblue3","gray0","yellowgreen","mediumorchid","orangered2","orangered4","plum1"), cex=0.7, pt.cex=1.5, box.lwd=0, bty="n", box.col="gray", horiz = F)
recordPlot()
}
cns = c("plum1", "cornflowerblue", "deepskyblue3", "gray0", "yellowgreen","mediumorchid", "orangered2","orangered4")
names(cns) = c("100", "0", "1", "2", "50", "3", "4", "5")
################################################################################
# two-variable pointplot with fine control to summarize titan and sequenza cnv variant calls
################################################################################
titanexome.df <- ptbafPT_AR_Ex.df
pointplot3pdf <- function(titanexome.df, colvec,titanexomeE.df, colvecE){
xmax=nrow(titanexome.df)
rangeDf <- range(titanexome.df[(titanexome.df != Inf) & !is.na(titanexome.df)])
ymin=min(rangeDf)
ymax=max(rangeDf)
par(mar=c(8,5,2,8), font.axis=1, xpd=TRUE)
plot(titanexome.df[,1], col= colvec, main="", pch=15, bty="n", cex= 1.8, xaxt="n", xlab="", ylab="", las=1, mgp=c(0,0.5,0), tck=-0.02, ylim=c(0,1), xlim=c(1.5,39))#, panel.first=grid(lty=1, nx=NA, ny=NULL, lwd=0.1))mgp=c(0,0.5,0)control the placement of Ytick lables
points(titanexomeE.df[,1], col=colvecE, alpha=80, pch=21, bg="white", add=TRUE)
axis(side=1, line=0, labels=rownames(titanexome.df), at=1:xmax, las=3, cex.axis=0.7, mgp=c(3,0.6,0), tck=-0.02)
for(i in seq(from=1, to=xmax, by=1)){ lines(c(i,i),c(0,1), lwd=0.1, col="grey") }
for(i in seq(from=0, to=1, by=0.1)){ lines(c(0.35,xmax+0.65), c(i,i), lwd=0.1, col="grey") }
points(titanexome.df[,1], col=colvec, alpha=100, pch=15, cex= 1.8, add=TRUE)
points(titanexomeE.df[,1], col=colvecE, alpha=80, pch=21, cex= 0.8, bg="white", add=TRUE)
mtext("AllelicRatio", side=2, line=1.9, adj=0.5, cex=1.2, col="black")
legend(0.6,-0.2, inset=.02, title= "", c("GeneCoordinates     o 300Kb-ExpandedGeneCoordinates"),pch=15, col="black", cex=0.8, pt.cex=1.4, box.lwd=0, bty="n", box.col="gray", horiz = F)
legend(42,1.02, inset=.02, title="TitanCall",c("DLOH", "HOMD", "NLOH", "HET", "GAIN", "ASCNA", "BCNA", "UBCNA","NA"),pch=15, col=c("blue2","deepskyblue3","yellowgreen", "gray0", "red", "darkred","darkred","darkred","white"), cex=0.7, pt.cex=1.5, box.lwd=0, bty="n", box.col="gray", horiz = F)
}
cns_1 = c("cornflowerblue","deepskyblue3","yellowgreen", "gray0", "orangered2","orangered4","orangered4","orangered4","white")
names(cns_1) = c("DLOH", "HOMD", "NLOH", "HET", "GAIN", "ASCNA", "BCNA", "UBCNA","NA")
cns_1[c("NLOH")]
as.character(cns_1["NA"])
################################################################################
# Combining graphs based on different patients.
################################################################################
for (patient in c(colnames(patientcovz.df))){
  ztemp <-as.data.frame(patientcovz.df[,patient,drop=FALSE])
  logtemp <- as.data.frame(zLog_tumNor.df[,patient,drop=FALSE])
  logtempcappcov <- as.data.frame(zLogcappcov_tumNor.df[,patient,drop=FALSE])
  ctemp_titan <- as.character(c(t(colvec_titan[,patient])))
  ctemp_titan_short <- as.character(c(t(colvec_titan_short[,patient])))
  ctemp_seq <- as.character(c(t(colvec_seq[,patient])))
  baf <- as.data.frame(ptbafPT_AR_Ex.df[,patient,drop=FALSE])
  is.na(baf)="NA" 
  cCNS <- as.character(c(t(ptbafPT_CNS_Ex.df[,patient])))
  is.na(cCNS)="NA"
  bafE <- as.data.frame(ptbafPT_AR_EEX.df[,patient,drop=FALSE])
  is.na(bafE)="NA" 
  cCNSE <- as.character(c(t(ptbafPT_CNS_EEx.df[,patient])))
  is.na(cCNSE)="NA"
#   Sequenza <- t(Tpurity_sequenza_EEx_clean.df[,patient])
#   Titan <- t(Tpurity_titan_EEx_clean.df[,patient])
  pdf(paste(patient,"_cnvpt.pdf",sep = ""), height=10, width=18)
  par(mfcol=c(2,1))
  pointplot2pdf(ztemp,logtempcappcov, cns[ctemp_titan], cns[ctemp_titan_short], cns[ctemp_seq])
  pointplot3pdf(baf,cns_1[cCNS], bafE,cns_1[cCNSE])
  # pointplot4pdf(Sequenza,Titan)
  dev.off()
}
