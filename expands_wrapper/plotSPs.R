

newPlotSPs<-function(dm, sampleID=NA,cex=0.5, legend="CN_Estimate", orderBy="chr", rawAF=F){
  
  keep=which(!is.na(dm[,"SP"]),); dm=dm[keep,];
  ia=order(dm[,"startpos"]);dm=dm[ia,];
  ia=order(dm[,orderBy]);dm=dm[ia,];
  ia=order(dm[,"SP"],decreasing = TRUE);dm=dm[ia,];
  
  maxPloidy=max(dm[,"PM"],na.rm=TRUE);
  yticklab=c(1:maxPloidy,seq(0,1,by=0.1));
  at=c(sort(c(1:maxPloidy)*-0.1),seq(0,1,by=0.1));
  
  if (!any(colnames(dm)=="AF_Tumor_Adjusted")){
  dm=.addColumn(dm,"AF_Tumor_Adjusted",NA);
  }
  if(!rawAF){
    dm[,"AF_Tumor_Adjusted"]=(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]-dm[,"PN_B"])/(dm[,"PM_B"]-dm[,"PN_B"])
    adjusted="Adjusted"
  }else{
    dm[,"AF_Tumor_Adjusted"]=dm[,"AF_Tumor"]
    adjusted=""
  }
  
  
  par(xpd=T, cex=cex, cex.axis=1/cex,cex.lab=1/cex, cex.main=1/cex,mar=par()$mar+c(0,0.5,0,4.2));
  #   plot.new(); 
  par(xpd=FALSE)
  
  plot(c(1:length(at)),at, col="white",pch=8,xlim=c(0,nrow(dm)),
       yaxt="n", bty="L", main=sampleID, xlab="Mutation", 
       ylab=paste("Copy-number <--->",adjusted,"Allele-frequency and SP size"));
  axis(2, at, yticklab) 
  
  legend1=.plotSPPerVar(dm,17,1,0,legend);
  
  x=gray.colors(100)
  norm=1/length(x);
  for (k in 1:nrow(dm)){
    ci=max(1,ceil(dm[k,"%maxP"]/norm));
    matpoints(k,dm[k,"SP"],pch=15,col=x[ci]);
    if (k==1){
      legend1$text[length(legend1$text)+1]="SP";
      legend1$col[length(legend1$text)]=x[ci];
      par(xpd=TRUE)
      legend("topright",legend1$text,fill=legend1$col,inset=c(-0.1,-0.02),cex=0.75/cex,bty = "n")
    }
  }
  .plotSPPerVar(dm,8,0,0,legend);
  .plotSPPerVar(dm,8,0,1,legend);
  .plotSPPerVar(dm,8,1,1,legend);
  lines(c(0,nrow(dm)),c(0,0),col="black");
  return(dm)
}


.plotSPPerVar<-function(dm,lineType,lohFlag,cnFlag,var){
  if(var=="CN_Estimate"){
    .plotSPPerCopyNumber(dm,lineType,lohFlag,cnFlag)
  }else if(var=="chr"){
    .plotSPPerChr(dm,lineType,lohFlag,cnFlag)
  }
}

.plotSPPerChr<-function(dm,lineType,lohFlag,cnFlag){
  x=rainbow(40);
  legend1=list("text"=c(),"col"=c());
  maxploidy=max(dm[,"PM"],na.rm=TRUE)+1;
  for (i in 1:22){
    idx=which(dm[,"chr"]==i & ((lohFlag & dm[,"PN_B"]==1) |(!lohFlag & dm[,"PN_B"]==0)) );
    if (!isempty(idx)){
      if (!cnFlag){
        matpoints(idx,dm[idx,"AF_Tumor_Adjusted"],pch=lineType,col=x[i]);
      }else{
        if (any(!is.na(dm[idx,"PM"]))){
          matpoints(idx,0.1*(dm[idx,"PM"]-maxploidy),pch=20,col=x[i]);
        }
      }
    }
    legend1$text[i]=paste("chr",i);
    legend1$col[i]=x[i];
  }
  return(legend1);
}

.addColumn<-function(M,newCol,initVal){
  if (!any(colnames(M)==newCol)){
    M=matrix(cbind(M,matrix(initVal,nrow(M),1)),nrow=nrow(M),ncol=ncol(M)+1,
             dimnames = list(rownames(M), c(colnames(M),newCol)));
  }
  return(M);
}

.plotSPPerCopyNumber<-function(dm,lineType,lohFlag,cnFlag){
  #x=c(terrain.colors(4),rainbow(4));
  library(RColorBrewer)
  x=brewer.pal(8, "Paired")
  legend1=list("text"=c(),"col"=c());
  maxploidy=max(dm[,"PM"],na.rm=TRUE)+1;
  for (i in 1:length(x)){
    idx=which(round(dm[,"CN_Estimate"])==i & ((lohFlag & dm[,"PN_B"]==1) |(!lohFlag & dm[,"PN_B"]==0)) );
    if (!isempty(idx)){
      if (!cnFlag){
        matpoints(idx,dm[idx,"AF_Tumor_Adjusted"],pch=lineType,col=x[i]);
      }else{
        if (any(!is.na(dm[idx,"PM"]))){
          matpoints(idx,0.1*(dm[idx,"PM"]-maxploidy),pch=20,col=x[i]);
        }
      }
    }
    legend1$text[i]=paste("CN",i);
    legend1$col[i]=x[i];
  }
  return(legend1);
}
