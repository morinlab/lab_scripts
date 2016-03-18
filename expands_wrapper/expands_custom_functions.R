  



newClusterCellFrequencies <- function(densities, precision, nrep=30, min_CellFreq=0.1, p_cutoff=0.8, local_sum = TRUE){ ##, plotF=0
    #local_sum forces the clustering to use the local maximum in the probability distribution rather than global maximum
    #p_cutoff is used to remove extremely high P estimates in the added localSum method
    #kurtosis feature is not used

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
      ia = NULL
      if(local_sum){
        ia=apply(clusterM,1,localSum,max_p_threshold=p_cutoff)
        
      } else {
        ia=apply(clusterM,1,which.max)
      }
      ia=order(freq[ia])
      clusterM=clusterM[ia,];
      
      ##Extend around maxima
      
      # replaced with weightedMean here too
      if(local_sum){
       peak=newWeightedMean(clusterM,freq);
       } else{
        peak=weightedMean(clusterM,freq);
       }
     
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
      if(local_sum){
        # Instead of global max, get weighted max (sum of P around peak)  
        ia = localSum(mean1,max_p_threshold=p_cutoff)
        ib = localSum(mean2,max_p_threshold=p_cutoff)
      } else{
        ia=which.max(apply(cl1,2,mean));
        ib=which.max(apply(cl2,2,mean));
      }
      
  

      meanCl=c(freq[ia],freq[ib]);

      peakIdx=which.min(abs(meanCl-peak));
      peakCl=densitiesOk[which(Tx==peakIdx),];
      meanClPeak=meanCl[peakIdx]; 
      if(local_sum){
       wMean=newWeightedMean(peakCl[,idx],freq[idx]);
      } else{
       wMean=weightedMean(peakCl[,idx],freq[idx]);
      }

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
  
  
  robSPs=.chooseRobustSPs(allSPs,precision,min_CellFreq);
  SPs=.collapseSimilar(robSPs$SPs,precision);
  
  outcols=c("Mean Weighted","score","precision","nMutations");##printlay only these columns
  if (is.null(dim(SPs))){
    SPs=SPs[outcols];
  }else{
    SPs=SPs[,outcols];
  }
  print("Done.");
  return(SPs);
}

.addColumn<-function(M,newCol,initVal){
  if (!any(colnames(M)==newCol)){
    M=matrix(cbind(M,matrix(initVal,nrow(M),1)),nrow=nrow(M),ncol=ncol(M)+1,
             dimnames = list(rownames(M), c(colnames(M),newCol)));
  }
  return(M);
}



newAssignMutations<-function( dm, finalSPs, max_PM=6, p_cutoff=0.8,prune_by_ploidy = 0,min_kurtosis=3,debug=0){
  #modified to use localSum method to select more appropriate fit based on local maximum in probability distributions
  #p_cutoff will cause high probabilities to be masked. these seem to be a symptom of a bad fit of the model.
    library(moments)
    if (is.null(dim(finalSPs))) {
        spFreq = finalSPs[ "Mean Weighted"]
        precision=finalSPs["precision"]
    }  else {
        spFreq = finalSPs[, "Mean Weighted"]
        precision=finalSPs[1,"precision"]
    }
    spFreq=sort(spFreq);
  
    added = matrix(nrow=length(dm[,1]),ncol=7)
    colnames(added) = c("%maxP","SP","PM_B","SP_cnv","PM","PM_cnv","scenario");
    dm=cbind(dm,added)
    if (!any(colnames(dm)=="f")){
        dm=.addF(dm,  max_PM);
    }
    if (!any(colnames(dm)=="AF_Tumor_Adjusted")){
        dm=.addColumn(dm,"AF_Tumor_Adjusted",NA);
    }
    dm[,"SP"]=NA;  dm[,"SP_cnv"]=NA; dm[,"PM_B"]=NA; dm[,"PM"]=NA; dm[,"PM_cnv"]=NA; dm[,"scenario"]=NA;##delete any potentially existing SP info
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
                #changed here to use localSum method
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
        
        #some variables to internally track the default values of maxP_J etc for comparison. These could all be removed because they are not actually used except for debugging purposes
        maxP_J=0; ##Maximum probability from joined fit
        maxP_J_original = 0; 
        maxP_S=0; ##Maximum probability from separate fit
        maxP_S_original=0;
        maxP_SbeforeC=0; ##Maximum probability from separate fit, under the assumption that CNV happened in descendant of SP with SNV
        maxP_SB_original = 0;
        if(class(snvJ)!="try-error" && any(!is.na(snvJ$p))){
            maxP_J_original=max(snvJ$p,na.rm=T)
            best_ind_original = which.max(snvJ$p)
            lsres = localSum(snvJ$p,max_p_threshold=p_cutoff,simple=FALSE)
            maxP_J = lsres[2]
            # contents: maxp_idx,maxP_J,maxP_best_ind,kurt_J
            changed = 0
            #the kurtosis cutoff being used here may be unnecessary. It could be removed in favor of relying just on localSum to see if the results are any different.
            if(is.na(lsres[4]) | lsres[4] < min_kurtosis){ #kurtosis is too low, ignoring this fit completely
                maxP_J = 0
                changed = 1
            }
            #can be removed, just here for debugging
            if(debug){
                if(changed || !(best_ind_original ==lsres[1])){
                    wd = getwd()
                    file=paste(wd,"/","debug","/","snv",k,"PN_B",dm[k,"PN_B"],"AF",dm[k,"AF_Tumor"],"CN",dm[k,"CN_Estimate"],"_J.pdf",sep="")
                    pdf(file)
                    barplot(snvJ$p,names=freq,las=2,main=paste("LSp",snvJ$p[lsres[1]],"K",signif(lsres[4],3),"Freq:",signif(freq[lsres[1]],3)))
                    dev.off()
                }
            }
        }
        if(class(snvS)!="try-error" && any(!is.na(snvS$p))){
            maxP_S_original=max(snvS$p,na.rm=T)
            best_ind_original = which.max(snvS$p)
            lsres = localSum(snvS$p,max_p_threshold=p_cutoff,simple=FALSE)
            maxP_S = lsres[2]
            # contents: maxp_idx,maxP_J,maxP_best_ind,kurt_J

            #can be removed, just here for debugging
            if(debug){
                if(!(best_ind_original ==lsres[1])){
                    wd = getwd()
                    file=paste(wd,"/","debug","/","snv",k,"PN_B",dm[k,"PN_B"],"AF",dm[k,"AF_Tumor"],"CN",dm[k,"CN_Estimate"],"_S.pdf",sep="")
                    pdf(file)
                    barplot(snvS$p,names=freq,las=2,main=paste("LSp",snvS$p[lsres[1]],"K",signif(lsres[4],3),"Freq:",signif(freq[lsres[1]],3)))
                    dev.off()
                }
          
            }
        }
        if(!is.null(snvSbeforeC) && class(snvSbeforeC)!="try-error" && any(!is.na(snvSbeforeC$p))){
            maxP_SB_original=max(snvS$p,na.rm=T)
            best_ind_original = which.max(snvS$p)
            lsres = localSum(snvS$p,max_p_threshold=p_cutoff,simple=FALSE)
            maxP_SBeforeC = lsres[2]
            # contents: maxp_idx,maxP_J,maxP_best_ind,kurt_J
            
            changed = 0
            if(is.na(lsres[4]) | lsres[4] < min_kurtosis){
                #kurtosis is too low, ignoring
                maxP_SBeforeC = 0
                
                changed = 1
            }
            #can be removed, just here for debugging
            if(debug){
                if(changed || !(best_ind_original ==lsres[1])){
                    
                    wd = getwd()
                    file=paste(wd,"/","debug","/","snv",k,"PN_B",dm[k,"PN_B"],"AF",dm[k,"AF_Tumor"],"CN",dm[k,"CN_Estimate"],"_SbeforeC.pdf",sep="")
                    pdf(file)
                    barplot(snvSbeforeC$p,names=freq,las=2,main=paste("LSp",snvSbeforeC$p[lsres[1]],"K",signif(lsres[4],3),"Freq:",signif(freq[lsres[1]],3)))
                    dev.off()
                }
            }
        }
    
        ##Skip if no solution found
        if (maxP_J==0 && maxP_S==0 && maxP_SbeforeC==0){
            dm[k,"SP"]=NA;
            
            next;
        }
        joinedFit=FALSE;
        joinedFitDefault =FALSE;
        if (class(snvJ)!="try-error" && (dm[k,"PN_B"]==1 || maxP_J >=max(maxP_S,maxP_SbeforeC,na.rm=T))){ ##SP carrying SNV and SP carrying CNV have same size, i.e. are identical:
            joinedFit=TRUE; ##LOH has to be associated with copy number variation --> SNV and CNV must be fit together
            snv=snvJ; 
            dm[k,"scenario"]=3;
        } else{ 
            if(maxP_S>=maxP_SbeforeC){
                snv=snvS;
                dm[k,"scenario"]=1;
            }else{
                snv=snvSbeforeC;
        
                dm[k,"scenario"]=4;
            }
        }
        message = ""
        #compare to default behaviour for debugging and evaluation purposes
        scenario_switch = 0
        
            scenario = 0
            if (class(snvJ)!="try-error" && (dm[k,"PN_B"]==1 || maxP_J_original >=max(maxP_S_original,maxP_SB_original,na.rm=T))){ ##SP carrying SNV and SP carrying CNV have same size, i.e. are identical:
                scenario =3
                snv_default = snvJ
                joinedFitDefault =TRUE;
            }else{ 
                if(maxP_S_original>=maxP_SB_original){
                    scenario =1
                    snv_default = snvS
          
                    }else{
                        scenario =4
                        snv_default = snvSbeforeC          
                    }
            }
            if(!(scenario == dm[k,"scenario"])){
                scenario_switch = 1
            }
        
    
    
        #new behaviour
        best_p_index = localSum(snv$p,max_p_threshold=p_cutoff)[1]
        idx = which.min(abs(spFreq-freq[best_p_index]))
        dm[k,"SP"] = spFreq[idx] 
        idx1=which(abs(snv$fit[,"f"]-dm[k,"SP"])<=precision/2); ##index of fits matching SP size
    
        idx_new=idx1[which.min(snv$fit[idx1,"dev"])] ##index of fit of matching SP size with minimal residual (dev)

        #default behaviour

        idx=which.min(abs(spFreq-freq[which.max(snv_default$p)]))
        SP = spFreq[idx]; 
        if(! (SP == dm[k,"SP"])){
            if(debug){
             print(paste("oldSP:",SP,"newSP:",dm[k,"SP"]),freq[which.max(snv_default$p)],freq[best_p_index])
            }
        }
    
        idx=which(abs(snv_default$fit[,"f"]-SP)<=precision/2); ##index of fits matching SP size
    
        idx=idx[which.min(snv_default$fit[idx,"dev"])] ##index of fit of matching SP size with minimal residual (dev)
    

        if(isempty(idx_new)){
            #no good fit found
            
            dm[k,"SP"]=NA
            next
        }

       
        dm[k,c("PM_B","PM")]=snv$fit[idx_new,c("PM_B","PM")]; ##(dm[k,"CN_Estimate"]*dm[k,"AF_Tumor"]-(1-dm[k,"SP"])*dm[k,"PN_B"])/dm[k,"SP"];  dm[k,"PM_B"]=max(1,dm[k,"PM_B"]);
        
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
        PM_B = 0
        PM = 0
        PM_cnv =0
        SP_cnv = 0
        #--> is all just for debugging purposes
        if(!isempty(idx)){
            PM_B = snv_default$fit[idx,"PM_B"]
            PM = snv_default$fit[idx,"PM"]
            if (!is.na(PM) && PM<0){
                PM = NA
            }
            if(joinedFitDefault){
                SP_cnv = SP
                PM_cnv=PM
            }else if (!is.na(pm) && pm==2){
                SP_cnv=SP
                PM_cnv=pm;
                PM=pm;
                #dm[k,"PM"]=pm;
            }else{
                SP_cnv=f_CNV;
                PM_cnv=pm 
                if(!is.na(SP) && !is.na(SP_cnv)){
                    if(SP>SP_cnv){
                        PM=2;
                    }else{
                        if(PM_B>PM_cnv && PM_B<=2){
                            PM=2;
                        }else if(PM_B>2 && PM_B<=PM_cnv){
                            PM=PM_cnv
                        }else{
                            PM=NA;
                        }
                    }
                }
            }
            if(scenario_switch && debug){
                message = paste("DefScen:",scenario,"NewScen:",dm[k,"scenario"],dm[k,"PM"],dm[k,"PM_B"],PM,PM_B)
                print(message)
            }
        }
        if(debug){
           print(paste(PM_B,dm[k,"PM_B"],PM,dm[k,"PM"], dm[k,"SP"],SP,sep="="))
         }
        #<-----
       
        
    success=success+1;
    if (mod(k,20)==0){
        print(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ", success,"/",k))
    }
  }
  
  dm[dm[,"%maxP"]==0,"SP"]=NA;

  #add CCF/adjucted VAF to each mutation based on results
  dm[,"AF_Tumor_Adjusted"]=(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]-dm[,"PN_B"])/(dm[,"PM_B"]-dm[,"PN_B"])

  #EXPERIMENTAL: assess the quality of SPs based on the fraction of variants with PM_B = 1 (should be the most common scenario in most cases because homozygous SNVs should be rare as they can only arise due to deletion or neutral LOH)
  #a more sensible extension here may be to only consider variants in regions of CN 2 or higher. Note, this change may be useless now given the improved clustering in v 1.7
  if(prune_by_ploidy){
    toRm=c();
    for (j in 1:size(finalSPs,1)){
      idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
      PM_B1=sum(dm[idx,"PM_B"]==1 | dm[idx,"PN_B"]==1)
      tot = length(idx)
      print(paste("For SP:",finalSPs[j,"Mean Weighted"],PM_B1,"out of",tot,"are ploidy 1 or LOH"))
      if(PM_B1 < tot/2 | tot < 10){
        #if a minority of SNVs in this cluster have PM_B 1, remove it and then rerun this to assign these to another SP
        toRm=c(toRm,j);
        print("removing")
      }
      #also adjust SP based on mean AF_Adjusted value, since it should improve the accuracay of the SP and subsequent assignment of variants to it
      adj_vaf_mean = mean(dm[idx,"AF_Tumor_Adjusted"])
      adj_vaf_ploidy1 = mean((dm[idx,"AF_Tumor"]*dm[idx,"CN_Estimate"]-dm[idx,"PN_B"])/(1-dm[idx,"PN_B"]))

      print(paste("SP:",finalSPs[j,"Mean Weighted"],'mean CCF:',adj_vaf_mean,"and assuming ploidy 1",adj_vaf_ploidy1))

      finalSPs[j,"Mean Weighted"] = adj_vaf_mean
      dm[idx,"SP"] = adj_vaf_mean
    }
    if(length(toRm)>0){
     finalSPs=finalSPs[-1*toRm,];
    }
  } else{
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
  }


  output=list("dm"=dm,"finalSPs"=finalSPs);
  return(output);
}

localSum<-function(probs,simple=TRUE,max_p_threshold = NULL){
  # This function was devised to handle the observation that the probability distributions seem to favour LOH events even when the copy number is normal and generally seem to
  # favour a higher ploidy value. The peaks for the higher ploidy in the probability distributions are sharp but there is a clear local maximum often for the alternate ploidy
  # This function attempts to find the local maximum by calculating the sum of all probabilities in individual peaks and returning the index of the peak with the highest sum rather
  # than the maximum point probability, which is used in the default behaviour of Expands.

  #another feature of this function is to assess the kurtosis of each peak. This is to handle another scenario that was observed in which some probability distributions contain
  #regions with very broad local maxima (blocky peaks). Such low kurtosis peaks seem to usually represent a poor quality fit and should ideally be removed. This function returns the kurtosis of the chosen peak
  #to help in determining if the fit is worth keeping.
  
  #"simple" mode skips kurtosis calculations

  # if threshold is supplied, toss any peak region that contains a probability > this. Lower quality predictions seem to result from fits of the model with such high values. 

  #one known (or suspected) limitation of this function is that probabilities associated with different SPs are not considered separately and almost certainly should be. It's not clear how serious this problem may be

  if(!missing(max_p_threshold)){
    above_thresh = which(probs > max_p_threshold)
    probs[above_thresh] = 0
  }
  peak_regions = sign(diff(probs))
  
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
            
            peaksum[peaknum] = peaksum[peaknum] + probs[ind]
            if(peakmax[peaknum]<probs[ind]){
              peakmax[peaknum] = probs[ind]
              
              peakmax_ind[peaknum] = ind
            }
            peakvals=c(peakvals,probs[ind])
            last_ind = ind
            ind = ind+1
          } else{
            #new peak
            peakmax[peaknum] = probs[ind]
            peaksum[peaknum] = peaksum[peaknum] + probs[ind]
            peakmax_ind[peaknum] = ind
            peakvals = c(probs[ind])
            last_ind = ind
            ind = ind+1
          }
          last_sign = i
        } else{
          if(last_sign < 0){
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
      best_peak_idx = which.max(peaksum)
      best_peak_max_idx = peakmax_ind[best_peak_idx]
      max_val_best_peak = peakmax[best_peak_idx]
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
  ##weighted mean
  wMean=0;sumWeight=sum(apply(peakCl,1,na.rm=T,max));
  for (pI in 1:nrow(peakCl)){
    maxIdx=which.max(peakCl[pI,]);
    wMean=wMean+(peakCl[pI,maxIdx]/sumWeight)*freq[maxIdx];
  }
  return(wMean);
}
newWeightedMean<-function(peakCl,freq){
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


