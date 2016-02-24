  


#this doesn't actually seem to work any better
  newClusterCellFrequencies <- function(densities, precision, nrep=30, min_CellFreq=0.1, p_cutoff=0.8){ ##, plotF=0
#   if(plotF>0 && !require(rgl)){
#     plotF=0;
#     message("Plot supressed:  Package rgl required for 3D plot of subpopulation clusters. Load this package befor using this option.")
#   }

  library(flexmix)
  library(matlab)
  library(mclust)
  library(moments)
  library(permute)
  freq=as.numeric(colnames(densities));
  print(paste("Clustering ",nrow(densities),"probability distributions..."))
  cols=c("red","yellow","green","pink","magenta","cyan","lightblue","blue");
  
  count=0;##counting xtick for cluster plot
  clusterMethod="average";
  print(paste("Clustering agglomeration method:",clusterMethod))
  
  ##Cluster probabilities based on Kullback-Leibler divergence
  D = KLdiv(t(densities),eps=10^-18);
  idxx=which(apply(is.finite(D),1,all) );
  print(paste(nrow(densities)-length(idxx),"SNVs excluded due to non-finite pdfs") );
  D=D[idxx,idxx];densities=densities[idxx,];
  Z = hclust(as.dist(D),method = clusterMethod);
  TC = cutree(Z,k=round(sqrt(nrow(densities))));
  clIdx=unique(TC);
  print("Done");
  print("Filtering Clusters...");
  allSPs=list();tRep=nrep;
  while (nrep>0){
    if (mod(nrep,3)==0){
      print(paste(100*(tRep-nrep)/tRep, "% completed"))
    }
    ##test each cluster for significance --> printlay SPs
    spCols=c("Max Size","Mean Size","Mean Weighted","wilcTest","wilcTest_Mean","wilcTest_Sum","wilcTest_Kurtosis","kurtosis","kurtosisNoise","kurtosisMean","kurtosisNoiseMean","nMutations","precision","score");
    SPs <- matrix(nrow = length(clIdx), ncol = length(spCols), byrow = TRUE, dimnames = list(paste(c(1:length(clIdx))), spCols))
    
    for (k in 1:length(clIdx)){
      meanClPeak<-peak<-wMean<-score<-NA; ##init
      
      ia=which(TC==clIdx[k]);
      if (length(ia)<2){
        next;
      }
      clusterM=densities[ia,];
      ##Sort
      #ia=apply(clusloerM,1,which.max);

      ia=apply(clusterM,1,localSum,threshold=p_cutoff);
      ia=order(freq[ia]);
      clusterM=clusterM[ia,];
      
      ##Extend around maxima
      
      #possibly replace weightedMean here too?/

      peak=weightedMean(clusterM,freq);
      #peak_old = .weightedMean(clusterM,freq)
      #print(paste("K=",k,"old_wm:",peak_old,"new_wm:",peak))
      peakMin=peak-0.05;peakMax=peak + 0.05;
      idx=find(freq>=peakMin & freq<=peakMax);
      stdOk=which(apply(densities[,idx],1,std)>10^-5);
      if (length(stdOk)<5){
        next;
      }
      densitiesOk=densities[stdOk,];##remove densities with small std
      
      Tx=kmeans(densitiesOk[,idx], centers=2, nstart =10);
      Tx=Tx$cluster;
      if (length(unique(Tx))!=2 || length(which(Tx==1))<=1 ||
            length(which(Tx==2))<=1){
        next; #second cluster step unsuccessfull for this range;
      }
      
      
      
      cl1=densitiesOk[which(Tx==1),]; 
      cl2=densitiesOk[which(Tx==2),]; 

      mean1 = apply(cl1,2,mean)
      mean2 = apply(cl2,2,mean)

      # Instead of global max, get weighted max (sum of P around peak)  
     
      ia_old=which.max(apply(cl1,2,mean));
      ia = localSum(mean1,threshold=p_cutoff)
      ib_old=which.max(apply(cl2,2,mean));
      ib = localSum(mean2,threshold=p_cutoff)
      
  

      meanCl=c(freq[ia],freq[ib]);

      #print(paste("1:",freq[ia],"2:",freq[ib],sep=" "))
      peakIdx=which.min(abs(meanCl-peak));
      peakCl=densitiesOk[which(Tx==peakIdx),];
      meanClPeak=meanCl[peakIdx]; 
      
      wMean=weightedMean(peakCl[,idx],freq[idx]);
      #print(peakCl[,idx])
      #print(freq[idx])
      tryCatch({   
        ##find peak range of cluster
        maxCl=apply(peakCl,2,mean,na.rm=T);
        x=apply(peakCl[,idx],1,max);
        y=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,max);
        zz=wilcox.test(x,y,conf.level=0.99,alternative="greater");
        
        x2=apply(peakCl[,idx],1,mean);
        y2=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,mean);  
        zz2=wilcox.test(x2,y2,conf.level=0.99,alternative="greater");
        
        x3=apply(peakCl[,idx],1,sum);
        y3=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,sum);              
        zz3=wilcox.test(x3,y3,conf.level=0.99,alternative="greater");
        kurt=apply(peakCl[,idx],1,kurtosis);
        kurtNoise=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,kurtosis);
        zzK=wilcox.test(kurt,kurtNoise,conf.level=0.99,alternative="greater");
        
        SPs[k,]=c(peak,meanClPeak,wMean, zz$p.value,zz2$p.value,zz3$p.value,zzK$p.value,max(kurt),max(kurtNoise),mean(kurt),mean(kurtNoise),nrow(peakCl),precision, score);
        ##calculate score
        score=SPs[k,"wilcTest"]+SPs[k,"wilcTest_Mean"]+SPs[k,"wilcTest_Sum"]+1/log(SPs[k,"nMutations"]);
        if(!is.na(SPs[k,"kurtosisNoiseMean"])){
          score=score+SPs[k,"kurtosisNoiseMean"]/500;
        }else{
          score=score+1;
        }
        SPs[k,"score"]=score;
        
#         #plot option
#         if (exists("plotF") && plotF>0 && nrep==tRep){
#            col=cols[mod(k,length(cols))+1];
#            count=.addTo3DPlot(count,clusterM,freq,col);
#         }
      },error = function(e) {
        print(e);
      })
    }    
    ##collapse similar
    ia=order(SPs[,"Mean Weighted"],decreasing=T);
    SPs=SPs[ia,];
    SPs=.collapseSimilar(SPs,precision);
    ##print(paste("Found ",size(SPs,1),"SPs."));
    if (size(SPs,1)>0){
      allSPs[[length(allSPs)+1]]=SPs;
    }
    nrep=nrep-1;
  }
if (length(allSPs)==0){
    return(NULL);
  }
  
#   ##plot option
#   if (plotF>0){
#     title3d("",label);
#   }
  
  robSPs=.chooseRobustSPs(allSPs,precision,min_CellFreq);
  SPs=.collapseSimilar(robSPs$SPs,precision);
  
  outcols=c("Mean Weighted","score","precision","nMutations");##printlay only these columns
  if (is.null(dim(SPs))){
    SPs=SPs[outcols];
  }else{
    SPs=SPs[,outcols];
  }
  print("Done.");
  ##print(SPs);
  ##out=list("SPs"=SPs,"spGrid"=robSPs$spGrid);
  return(SPs);
}

#appears in testing to work a lot better than the default function
newAssignMutations<-function( dm, finalSPs, max_PM=6, p_cutoff=0.8){
  library(moments)
  if (is.null(dim(finalSPs))) {
    spFreq = finalSPs[ "Mean Weighted"]
    precision=finalSPs["precision"]
  }  else {
    spFreq = finalSPs[, "Mean Weighted"]
    precision=finalSPs[1,"precision"]
  }
  spFreq=sort(spFreq);
  
  ##PM_B is the ploidy of the B-allele in SP; PM is the total ploidy in SP_cnv
  #addCols=c("%maxP","SP","PM_B","SP_cnv","PM","PM_cnv","scenario");
  #for (k in 1:length(addCols)){
  #  dm=.addColumn(dm,addCols[k],NA);
  #}
  added = matrix(nrow=length(dm[,1]),ncol=7)
  colnames(added) = c("%maxP","SP","PM_B","SP_cnv","PM","PM_cnv","scenario");
  dm=cbind(dm,added)

  dm[,"SP"]=NA;  dm[,"SP_cnv"]=NA; ##delete any potentially existing SP info
  freq=c()
  for (sp in spFreq){
     freq=c(freq,seq(sp-precision/2,sp+precision/2,by=precision/20))
  }
  success=0;
  densities=matrix(matrix(NA,nrow(dm),length(freq)),nrow=nrow(dm),ncol=length(freq),dimnames=list(1:nrow(dm),freq));
  for(k in 1:nrow(dm)){
    ##Joined fit
    snvJ=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=max_PM),silent=TRUE)
    ##Separate fit
    f_CNV=NA; pm=NA;
    cnv=try(cellfrequency_pdf(NA,dm[k,"CN_Estimate"],NA,freq, max_PM=max_PM, snv_cnv_flag=2),silent=TRUE)
    snvSbeforeC=NULL;
    if(class(cnv)!="try-error" && any(!is.na(cnv$p))){
      if (max(cnv$p, na.rm=T)>0){
        cnv_bestp_ind = localSum(cnv$p)[1]
        idx = which.min(abs(spFreq-freq[cnv_bestp_ind]))
        #idx=which.min(abs(spFreq-freq[which.max(cnv$p)]))
        f_CNV=spFreq[idx];
        idx=which.min(abs(cnv$fit[,"f"]-f_CNV));
        pm=cnv$fit[idx,"PM"]; ##pm=(dm[k,"CN_Estimate"]-(1-f_CNV)*2)/f_CNV;  pm=max(0,pm); pm=min(max_PM,pm,na.rm=T);
        ##Fit under the assumption that SNV happened before CNV
        snvSbeforeC=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=NA, snv_cnv_flag=4, SP_cnv=f_CNV, PM_cnv=pm),silent=TRUE);
      }
    }
    ##Max_PM is either 2 if SP with SNV is not a descendant of SP with CNV, and is PM of SP with CNV otherwise. Since we don't know which applies we choose the maximum value
    snvS=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=max(c(pm,2),na.rm=T), snv_cnv_flag=1),silent=TRUE)
    
    maxP_J=0; ##Maximum probability from joined fit
    maxP_S=0; ##Maximum probability from separate fit
    maxP_SbeforeC=0; ##Maximum probability from separate fit, under the assumption that CNV happened in descendant of SP with SNV
    
    if(class(snvJ)!="try-error" && any(!is.na(snvJ$p))){
      #maxP_J=max(snvJ$p,na.rm=T)
      
      maxP_J  = localSum(snvJ$p,threshold=p_cutoff)[2]
      #print(paste("LS on:",k))
      #print(snvJ$p)
      #print(snvJ$p)
      if(max(snvJ$p) == as.numeric(1)){
        print(paste(k,"IGNORING JOINT BECAUSE maxP==1"))
        maxP_J = 0
      } else{
        kurt_J = localSum(snvJ$p,threshold=p_cutoff,simple=FALSE)[4]
        if(is.na(kurt_J) | kurt_J < 3){
	        print(paste(k,"IGNORING JOINT BECAUSE",kurt_J))
          maxP_J = 0
        }
      }
    }
    if(class(snvS)!="try-error" && any(!is.na(snvS$p))){
      #maxP_S=max(snvS$p,na.rm=T)
      #print(paste("LS on:",k))
      #print(snvS$p)
      #print(snvS$p)
      
      maxP_S  = localSum(snvS$p,threshold=p_cutoff)[2]
      #kurt_S = localSum(snvS$p,simple=FALSE)[4]
    }
    if(!is.null(snvSbeforeC) && class(snvSbeforeC)!="try-error" && any(!is.na(snvSbeforeC$p))){
      #print(paste("LS on:",k))
      #print(snvSbeforeC$p)
      #print(snvSbeforeC$p)
      if(max(snvSbeforeC$p) == as.numeric(1)){
        print(paste(k,"IGNORING SBC BECAUSE maxP==1"))
        maxP_SbeforeC = 0
      } else{
        maxP_SbeforeC  = localSum(snvSbeforeC$p,threshold=p_cutoff)[2]
        kurt_SC = localSum(snvSbeforeC$p,threshold=p_cutoff,simple=FALSE)[4]
      }
	   if(is.na(kurt_SC) | kurt_SC < 3){
	     print(paste(k,"IGNORING SbeforeC BECAUSE",kurt_SC))
        maxP_SbeforeC = 0
      }
      #maxP_SbeforeC=max(snvSbeforeC$p,na.rm=T)
    }
    
    ##Skip if no solution found
    if (maxP_J==0 && maxP_S==0 && maxP_SbeforeC==0){
      dm[k,"SP"]=NA;
      next;
    }
    
    joinedFit=FALSE;
    
    if (class(snvJ)!="try-error" && (dm[k,"PN_B"]==1 || maxP_J >=max(maxP_S,maxP_SbeforeC,na.rm=T))){ ##SP carrying SNV and SP carrying CNV have same size, i.e. are identical:
      joinedFit=TRUE; ##LOH has to be associated with copy number variation --> SNV and CNV must be fit together
      snv=snvJ; 
      dm[k,"scenario"]=3;
    }else{ 
      if(maxP_S>=maxP_SbeforeC){
        snv=snvS;
        dm[k,"scenario"]=1;
      }else{
        snv=snvSbeforeC;
        dm[k,"scenario"]=4;
      }
    }
    debug = 0
    if(debug){
      wd = getwd()
      file=paste(wd,"/","debug","/","simulation_snv_",k,"PN_B",dm[k,"PN_B"],"AF",dm[k,"AF_Tumor"],"CN",dm[k,"CN_Estimate"],sep="")
      j_file = paste(file,precision,"_J.pdf",sep="")
      s_file = paste(file,precision,"_S.pdf",sep="")
      sb_file = paste(file,precision,"_SbeforeC.pdf",sep="")
      if(class(snvJ)!="try-error" && any(!is.na(snvJ$p))){
        best_p_index = localSum(snvJ$p)[1]
        maxp = max(snvJ$p)
        best_freq = freq[best_p_index]
        if(best_freq < 0.5){
        pdf(j_file)
        barplot(snvJ$p,names=freq,las=2,main=paste("PS",snvJ$p[best_p_index],"K",kurt_J,"Freq:",best_freq,dm[k,"scenario"],precision))
        dev.off()
        }
      }
      if(class(snvS)!="try-error" && any(!is.na(snvS$p))){
        best_p_index = localSum(snvS$p)[1]
        maxp = max(snvS$p)
        best_freq = freq[best_p_index]
        if(best_freq < 0.5){
        pdf(s_file)
        barplot(snvS$p,names=freq,las=2,main=paste("PS",snvS$p[best_p_index],"Freq:",best_freq,dm[k,"scenario"],precision))
        dev.off()
      }
      }
      if(class(snvSbeforeC)!="try-error" && any(!is.na(snvSbeforeC$p))){
        best_p_index = localSum(snvSbeforeC$p)[1]
        best_freq = freq[best_p_index]
        kurt = kurtosis(snvSbeforeC$p)
        maxp = max(snvSbeforeC$p)
        if(best_freq < 0.5){
        pdf(sb_file)
        barplot(snvSbeforeC$p,names=freq,las=2,main=paste("K",kurt_SC,"Freq:",best_freq,dm[k,"scenario"],precision))
        dev.off()
      }
      }
    }
    ##Save end result:
    best_p_index = localSum(snv$p,threshold=p_cutoff)[1]
    idx_new = which.min(abs(spFreq-freq[best_p_index]))
    idx=which.min(abs(spFreq-freq[which.max(snv$p)]))

    #print(paste("best_by_localS:",idx_new,"best_old:",idx))
    #dm[k,"SP"]=spFreq[idx];  
    dm[k,"SP"]=spFreq[idx_new]; 
    idx=which.min(abs(snv$fit[,"f"]-dm[k,"SP"]));

    
    dm[k,c("PM_B","PM")]=snv$fit[idx,c("PM_B","PM")]; ##(dm[k,"CN_Estimate"]*dm[k,"AF_Tumor"]-(1-dm[k,"SP"])*dm[k,"PN_B"])/dm[k,"SP"];  dm[k,"PM_B"]=max(1,dm[k,"PM_B"]);
    if (!is.na(dm[k,"PM"]) && dm[k,"PM"]<0){
      dm[k,"PM"]=NA; ##PM can be -1 if obtained with snv_cnv_flag=1; TODO --> get NA directy for jar and remove this. 
    }
    #dm[k,"%maxP"]=max(snv$p,na.rm=T); #snv$p[which.min(abs(freq-dm[k,"SP"]))];
    dm[k,"%maxP"] = snv$p[best_p_index]

    densities[k,]=snv$p;
    
    if(joinedFit){
      dm[k,"SP_cnv"]=dm[k,"SP"];
      dm[k,"PM_cnv"]=dm[k,"PM"];
    }else if (!is.na(pm) && pm==2){
      dm[k,"SP_cnv"]=dm[k,"SP"];
      dm[k,"PM_cnv"]=pm; dm[k,"PM"]=pm;
      #dm[k,"PM"]=pm;
    }else{
      dm[k,"SP_cnv"]=f_CNV;
      dm[k,"PM_cnv"]=pm; 

      ##PM of SP does not have to be 2, because SP may also be a descendant of SP_cnv, i.e. the clone that acquired the CNV
      ##IF SP is larger than SP_cnv, then SP cannot be descandant of SP_cnv and therefor cannot harbor the CNV --> ploidy = 2
      ##IF SP is smaller than SP_cnv than SP may descend from SP_cnv: 
      ##-->Reject descendant hypothesis (ploidy of SP = 2) --> if 2 >= PM_B > PM 
      ##-->Accept descendant hypothesis (ploidy of SP = PM) --> if PM >= PM_B > 2 
      ##-->Irresolvable otherwise (ploidy of SP cannot be assigned) 
      if(!is.na(dm[k,"SP"]) && !is.na(dm[k,"SP_cnv"])){
        if(dm[k,"SP"]>dm[k,"SP_cnv"]){
          dm[k,"PM"]=2;
        }else{
          if(dm[k,"PM_B"]>dm[k,"PM_cnv"] && dm[k,"PM_B"]<=2){
            dm[k,"PM"]=2;
          }else if(dm[k,"PM_B"]>2 && dm[k,"PM_B"]<=dm[k,"PM_cnv"]){
            dm[k,"PM"]=dm[k,"PM_cnv"];
          }else{
            dm[k,"PM"]=NA;
          }
        }
      }
    }

    success=success+1;
    if (mod(k,20)==0){
      print(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ",
                  success,"/",k))
    }
  }
  
  dm[dm[,"%maxP"]==0,"SP"]=NA;


  ##Remove SPs to which no mutations were assigned
  toRm=c();
  for (j in 1:size(finalSPs,1)){
    if(is.null(dim(finalSPs))){
      idx=which(dm[,"SP"]==finalSPs["Mean Weighted"]);
      finalSPs["nMutations"]=length(idx);   
    }else{
      idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
      finalSPs[j,"nMutations"]=length(idx);
      if(length(idx)==0){
        toRm=c(toRm,j);
      }
    }
  }
  if(length(toRm)>0){
    finalSPs=finalSPs[-1*toRm,];
  }

#calculate mean %maxP for mutations assigned to each SP
#for (j in 1:size(finalSPs,1)){
#    if(is.null(dim(finalSPs))){
#      idx=which(dm[,"SP"]==finalSPs["Mean Weighted"]);
#      finalSPs["nMutations"]=length(idx);   
#    }else{
#      idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
#      finalSPs[j,"nMutations"]=length(idx);
#      if(length(idx)==0){
#        toRm=c(toRm,j);
#      }
#    }
#  }


  output=list("dm"=dm,"finalSPs"=finalSPs);
  return(output);
}

localSum<-function(means,simple=TRUE,threshold = NULL){
  #toss any peak reagion that contains a probability > this. Low quality predictions result from these. 
  if(!missing(threshold)){
    above_thresh = which(means > threshold)
    means[above_thresh] = 0
  }
  peak_regions = sign(diff(means))
  
      ind = 1
      last_sign = 0
      last_ind = 0
      peakmax = c(0)
      peaksum = c(0)
      peakmax_ind = c(0)
      peaknum = 1
      peakvals=c()
      peak_kurtosis = c()
      for(i in peak_regions){
      
        if(i < 0){
        
          if(last_sign == -1){
            #same peak
            
            peaksum[peaknum] = peaksum[peaknum] + means[ind]
            if(peakmax[peaknum]<means[ind]){
              peakmax[peaknum] = means[ind]
              
              peakmax_ind[peaknum] = ind
            }
            peakvals=c(peakvals,means[ind])
            last_ind = ind
            ind = ind+1
          } else{
            #new peak
            peakmax[peaknum] = means[ind]
            peaksum[peaknum] = peaksum[peaknum] + means[ind]
            peakmax_ind[peaknum] = ind
            peakvals = c(means[ind])
            last_ind = ind
            ind = ind+1
          }
          last_sign = i
        } else{
          if(last_sign < 0){
            
            #print(paste("K",kurtosis_last_peak))
            if(!simple){
              kurtosis_last_peak = kurtosis(peakvals)
              peak_kurtosis[peaknum] = kurtosis_last_peak
            }
            peaknum=peaknum+1
            peaksum = c(peaksum,0)
            
          }
          last_sign = i
          last_ind = ind
          ind = ind+1
          
        }

      }
      if(!simple){
        kurtosis_last_peak = kurtosis(peakvals)
        peak_kurtosis[peaknum] = kurtosis_last_peak
      }
      #print(paste("K",kurtosis_last_peak))
      best_peak_idx = which.max(peaksum)
      best_peak_max_idx = peakmax_ind[best_peak_idx]
      max_val_best_peak = peakmax[best_peak_idx]
      #print(peaksum)
      if(simple){
          return(c(best_peak_max_idx,peaksum[best_peak_idx]))
      }
      return(c(best_peak_max_idx,peaksum[best_peak_idx],max_val_best_peak,peak_kurtosis[best_peak_idx]))
  }

.collapseSimilar <-function(SPs,precision){
  isNaNIdx=which(is.na(SPs[,"Mean Weighted"]));
  if (!isempty(isNaNIdx)){
    SPs=SPs[-isNaNIdx,];
  }
  if (size(SPs,1)<2){
    return(SPs);
  }
  spSize=unique(round(SPs[,"Mean Weighted"]*100)/100);
  for (n in 1:length(spSize)){
    idx=which(abs(SPs[,"Mean Weighted"]-spSize[n])<precision);
    if (length(idx)>1){
      ia=which.min(SPs[idx,"wilcTest"]);
      rmIdx=setdiff(idx,idx[ia]);
      SPs=SPs[-rmIdx,];
    }
    if (is.null(dim(SPs))){
      break;
    }
  }
  return(SPs);
}
weightedMean<-function(peakCl,freq){
#modified
  ##weighted mean
  #wMean=0;sumWeight=sum(apply(peakCl,1,na.rm=T,max));
  #does this make more sense?
  wMean=0;sumWeight=sum(apply(peakCl,1,na.rm=T,max));
  for (pI in 1:nrow(peakCl)){
    #maxIdx=which.max(peakCl[pI,]);
    maxIdx=localSum(peakCl[pI,])[1]
    wMean=wMean+(peakCl[pI,maxIdx]/sumWeight)*freq[maxIdx];
  }
  return(wMean);
}
.chooseRobustSPs <- function(allSPs,precision,min_CellFreq){
  ## input parameter SPs is a cell array with DataMatrix (DM) objects. Each row in
  ## DM contains the size and the p-value associated with a SP. DM is sorted
  ## in descending order of SP size.
  #count frequencies among predictions;
  freq=t(seq(min_CellFreq,1,by=precision));
  SPsizes=matrix(nrow = length(allSPs), ncol = length(freq),
                 dimnames = list(paste(c(1:length(allSPs))), freq))
  for (i in 1:length(allSPs)){
    SPs=allSPs[[i]];
    if(is.null(dim(SPs))){
#       if(SPs["Mean Weighted"]>1){ ##Should no longer be necessary
#         SPs["Mean Weighted"]=1; 
#       }
      idx=which.min(abs(SPs["Mean Weighted"]-freq));
      SPsizes[i,idx]=SPs["score"];
    }else{
      #SPs[SPs[,"Mean Weighted"]>1,"Mean Weighted"]=1; ##Should no longer be necessary
      sps=SPs[,"Mean Weighted"];
      for (j in 1:length(sps)){
        idx=which.min(abs(sps[j]-freq));
        SPsizes[i,idx]=SPs[j,"score"];
      }
    }
  }
  keep=which(apply(!is.na(SPsizes),2,sum)>length(allSPs)*0.5);
  ##keep only best among recurrent SPs
  finalSPs=c();
  for (i in 1:length(keep)){
    sp_freq=freq[keep[i]];
    spI=which(!is.na(SPsizes[,keep[i]]));
    ia=which.min(SPsizes[spI,keep[i]]);
    SPs=allSPs[[as.numeric(spI[ia])]];
    if(is.null(dim(SPs))){
      SPs["Mean Weighted"]=sp_freq; ##standardize sp size
      finalSPs=rbind(finalSPs,SPs);
    }else{
      ia=which.min(abs(SPs[,"Mean Weighted"]-sp_freq));
      SPs[ia,"Mean Weighted"]=sp_freq; ##standardize sp size
      finalSPs=rbind(finalSPs,SPs[ia,]);
    }
  }
  output=list("SPs"=finalSPs,"spGrid"=SPsizes);
  return(output);
}


