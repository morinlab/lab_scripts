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
  plot(cappseq.df[,1], col=colvec_titan, main=colnames(cappseq.df), alpha=80, pch=4, cex=1, bty="n", xaxt="n", xlab="", ylab="", las=1, mgp=c(0,0.5,0), tck=-0.02, ylim=c(-6,6), xlim=c(2.4,91))#, panel.first=grid(lty=1, nx=NA, ny=NULL, lwd=0.1)) #control the size of graph by ylim and xlim
  points(cappseq.df[,1], col=colvec_titan_short, alpha=100, pch=20, cex= 2, add=TRUE)
  points(exome_2.df[,1], col=colvec_titan, alpha=80, pch=4, cex=1, add=TRUE)
  points(exome_2.df[,1], col=colvec_seq, alpha=80, pch=21, cex= 1.3, add=TRUE)
  axis(side=1, line=0, labels=rownames(cappseq.df), at=1:xmax, las=3, cex.axis=0.6, mgp=c(3,0.6,0), tck=-0.02)#mgp the second one change the placement of the samples codes
  for(i in seq(from=1, to=xmax, by=1)){ lines(c(i,i),c(-6,6), lwd=0.1, col="grey") }
  for(i in seq(from=-6, to=6, by=2)){ lines(c(0.35,xmax+0.65), c(i,i), lwd=0.1, col="grey") }
  points(cappseq.df[,1], col=colvec_titan, alpha=80, pch=4, cex=2, add=TRUE, na.col=NA)
  points(cappseq.df[,1], col=colvec_titan_short, alpha=100, pch=20, cex= 2, add=TRUE)
  points(exome_2.df[,1], col=colvec_titan, alpha=80, pch=4, cex=1, add=TRUE)
  points(exome_2.df[,1], col=colvec_seq, alpha=80, pch=21, cex=1.3, add=TRUE, na.col=NA)
  mtext("ZCoverage_Ratio", side=2, line=1.8, adj=0.5, cex=1, col="black")
  legend(0.5,-7.9, inset=.02, title= "", c("CappseqCov.TitanFocal   X CappseqCov.300Kb-ExpandedTitanSegments     O Log2ExomeCov.Sequenza   x Log2ExomeCov.Titan"), pch=20, col="black", cex=0.8, pt.cex=2, box.lwd=0, bty="n", box.col="gray", horiz = F)
  legend(99,6.3,inset=.02, title="CNS 2=HT",c('0','1','2','2NLOH','3','4','5','NA'),pch=15, col=c("cornflowerblue","deepskyblue3","gray0","yellowgreen","mediumorchid","orangered2","orangered4","plum1"), cex=0.7, pt.cex=1.5, box.lwd=0, bty="n", box.col="gray", horiz = F)
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
  plot(titanexome.df[,1], col= colvec, main="", pch=15, bty="n", cex= 1.8, xaxt="n", xlab="", ylab="", las=1, mgp=c(0,0.5,0), tck=-0.02, ylim=c(0,1), xlim=c(2.4,91))#, panel.first=grid(lty=1, nx=NA, ny=NULL, lwd=0.1))mgp=c(0,0.5,0)control the placement of Ytick lables
  points(titanexomeE.df[,1], col=colvecE, alpha=80, pch=21, bg="white", add=TRUE)
  axis(side=1, line=0, labels=rownames(titanexome.df), at=1:xmax, las=3, cex.axis=0.6, mgp=c(3,0.6,0), tck=-0.02)
  for(i in seq(from=1, to=xmax, by=1)){ lines(c(i,i),c(0,1), lwd=0.1, col="grey") }
  for(i in seq(from=0, to=1, by=0.1)){ lines(c(0.35,xmax+0.65), c(i,i), lwd=0.1, col="grey") }
  points(titanexome.df[,1], col=colvec, alpha=100, pch=15, cex= 1.8, add=TRUE)
  points(titanexomeE.df[,1], col=colvecE, alpha=80, pch=21, cex= 0.8, bg="white", add=TRUE)
  mtext("AllelicRatio", side=2, line=1.9, adj=0.5, cex=1.2, col="black")
  legend(99,1.02, inset=.02, title="TitanCall",c("DLOH", "HOMD", "NLOH", "HET", "GAIN", "ASCNA", "BCNA", "UBCNA","NA"),pch=15, col=c("blue2","deepskyblue3","yellowgreen", "gray0", "red", "darkred","darkred","darkred","white"), cex=0.7, pt.cex=1.5, box.lwd=0, bty="n", box.col="gray", horiz = F)
  legend(0.6,-0.2, inset=.02, title= "", c("GeneCoordinates     o 300Kb-ExpandedGeneCoordinates"),pch=15, col="black", cex=0.8, pt.cex=1.4, box.lwd=0, bty="n", box.col="gray", horiz = F)
}
cns_1 = c("cornflowerblue","deepskyblue3","yellowgreen", "gray0", "orangered2","orangered4","orangered4","orangered4","white")
names(cns_1) = c("DLOH", "HOMD", "NLOH", "HET", "GAIN", "ASCNA", "BCNA", "UBCNA","NA")
cns_1[c("NLOH")]
as.character(cns_1["NA"])
# ################################################################################
# #pointplots with error-bars based on Sequeza cellurarity estimates and Titan Tumour purity
# ################################################################################
# pointplot4pdf <- function(Tpurity_sequenza.df, Tpurity_titan.df){
#   fix.df <- Tpurity_sequenza.df[grepl("_F", rownames(Tpurity_sequenza.df)),1]    
#   xmax <- length(fix.df)
#   rangeDf <- range(Tpurity_sequenza.df[(Tpurity_sequenza.df != Inf) & !is.na(Tpurity_sequenza.df)])
#   ymin=min(rangeDf)
#   ymax=max(rangeDf)
#   names(fix.df)=sub('_F', '', names(fix.df))
#   pt <- names(fix.df)
#   par(mar=c(5,5,5,8), mgp=c(3,0.4,0), font.axis=1, xpd=TRUE)
#   plot(Tpurity_sequenza.df[grepl("_F", rownames(Tpurity_sequenza.df)),1], col='mediumpurple1' , main='', pch=18, bty="n", cex= 1.8, xaxt="n", xlab="", ylab="", las=1, mgp=c(3,0.5,0), tck=-0.02, ylim=c(0,1), xlim=c(2.4,83))#, panel.first=grid(lty=1, nx=NA, ny=NULL, lwd=0.1))
#   points(Tpurity_titan.df[grepl("_T", rownames(Tpurity_titan.df)),1], col="gray72", alpha=80, pch=18, bg="white", cex=1.5, add=TRUE)
#   axis(side=1, line=0, labels=pt, at=1:xmax, las=3, cex.axis=0.7, mgp=c(3,0.6,0), tck=-0.02)
#   for(i in seq(from=1, to=xmax, by=1)){ lines(c(i,i),c(0,1), lwd=0.1, col="grey") }
#   for(i in seq(from=0, to=1, by=0.1)){ lines(c(0.35,xmax+0.65), c(i,i), lwd=0.1, col="grey") }
#   points(as.character(Tpurity_sequenza.df[grepl("_U", rownames(Tpurity_sequenza.df)),1]), col='gray28', pch="-", cex= 1)
#   points(as.character(Tpurity_sequenza.df[grepl("_F", rownames(Tpurity_sequenza.df)),1]), col='mediumpurple1', pch=18, cex=1.8)
#   points(as.character(Tpurity_sequenza.df[grepl("_L", rownames(Tpurity_sequenza.df)),1]), col='gray28', pch="-", cex=1)
#   points(as.character(Tpurity_titan.df[grepl("_T", rownames(Tpurity_titan.df)),1]),  col='gray72', pch=18, cex=1.5)
#   points(as.character(Tpurity_titan.df[grepl("_Z1", rownames(Tpurity_titan.df)),1]), col='gray72', pch=4, bg="white", cex=0.8)
#   points(as.character(Tpurity_titan.df[grepl("_Z2", rownames(Tpurity_titan.df)),1]), col='gray72', pch=4, bg="white", cex=0.8)
#   points(as.character(Tpurity_titan.df[grepl("_Z3", rownames(Tpurity_titan.df)),1]), col='gray72', pch=4, bg="white", cex=0.8)
#   points(as.character(Tpurity_titan.df[grepl("_Z4", rownames(Tpurity_titan.df)),1]), col='gray72', pch=4, bg="white", cex=0.8)
#   points(as.character(Tpurity_titan.df[grepl("_Z5", rownames(Tpurity_titan.df)),1]), col='gray72', pch=4, bg="white", cex=0.8)
#   mtext("TumourCellFractions", side=2, line=1.9, adj=0.5, cex=1.2, col="black")
#   legend(80,1.08, inset=.02, title= "", c("Sequenza-CI","TitanxSubclones"),pch=c(18,18), col=c("mediumpurple1","gray72"), cex=0.7, pt.cex=c(1.3,1.3), box.lwd=0, bty="n", box.col="gray", horiz = F)
# }
################################################################################
# Combining graphs based on different genes.
################################################################################
for (gene in c(rownames(patientz.df))){
  ztemp <- t(patientcovz.df[gene,])
  ztemp [c(is.na(ztemp))] =  "NA"
  logtemp <- t(zLog_tumNor.df[gene,])
  logtemp [c(is.na(logtemp))] =  "NA"
  logtempcappcov <- t(zLogcappcov_tumNor.df[gene,])
  logtempcappcov [c(is.na(logtempcappcov))] =  "NA"
  ctemp_titan <- t(colvec_titan[gene,])
  ctemp_titan_short <- t(colvec_titan_short[gene,])
  ctemp_seq <- t(colvec_seq[gene,])
  baf <- t(ptbafPT_AR_Ex.df[gene,])
  baf [c(is.na(baf))] =  "NA" 
  cCNS <- as.character(c(t(ptbafPT_CNS_Ex.df[gene,])))
  cCNS[c(is.na(cCNS))] =  "NA"
  bafE <- t(ptbafPT_AR_EEX.df[gene,])
  bafE [c(is.na(bafE))] =  "NA" 
  cCNSE <- as.character(c(t(ptbafPT_CNS_EEx.df[gene,])))
  cCNSE[c(is.na(cCNSE))] =  "NA"
#   Sequenza <- t(Tpurity_sequenza_EEx_clean.df[gene,])
#   Titan <- t(Tpurity_titan_EEx_clean.df[gene,])
  pdf(paste(gene,"_cnv.pdf",sep = ""), height=10, width=18)
  par(mfcol=c(2,1))
  pointplot2pdf(ztemp,logtempcappcov, cns[ctemp_titan], cns[ctemp_titan_short], cns[ctemp_seq])
  pointplot3pdf(baf,cns_1[cCNS], bafE,cns_1[cCNSE])
  # pointplot4pdf(Sequenza,Titan)
  dev.off()
}
