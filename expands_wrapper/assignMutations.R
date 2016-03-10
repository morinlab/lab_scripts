assignMutations<-function( dm, finalSPs, max_PM=6){
  
  if (is.null(dim(finalSPs))) {
    spFreq = finalSPs[ "Mean Weighted"]
    precision=finalSPs["precision"]
  }  else {
    spFreq = finalSPs[, "Mean Weighted"]
    precision=finalSPs[1,"precision"]
  }
  spFreq=sort(spFreq);
  
  ##PM_B is the ploidy of the B-allele in SP; PM is the total ploidy in SP_cnv
  addCols=c("%maxP","SP","PM_B","SP_cnv","PM","PM_cnv","scenario");
  for (k in 1:length(addCols)){
    dm=.addColumn(dm,addCols[k],NA);
  }
  if (!any(colnames(dm)=="f")){
    dm=.addF(dm,  max_PM);
  }
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
        idx=which.min(abs(spFreq-freq[which.max(cnv$p)]))
        f_CNV=spFreq[idx];
        idx=which.min(abs(cnv$fit[,"f"]-f_CNV));
        pm=cnv$fit[idx,"PM"]; ##pm=(dm[k,"CN_Estimate"]-(1-f_CNV)*2)/f_CNV;  pm=max(0,pm); pm=min(max_PM,pm,na.rm=T);
        ##Fit under the assumption that SNV happened before CNV
        snvSbeforeC=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=NA, snv_cnv_flag=4, SP_cnv=f_CNV, PM_cnv=pm),silent=TRUE);
      }
    }
    ##Max_PM is either 2 if SP with SNV is not a descendant of SP with CNV.... 
    snvS_noDesc=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=2, snv_cnv_flag=1),silent=TRUE)
    ##.. or Max_PM is PM of SP with CNV otherwise
    snvS_Desc=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq[freq<=f_CNV+precision/2], max_PM=pm, snv_cnv_flag=1),silent=TRUE)
    ##Choose better solution between the two
    snvS=snvS_noDesc;
    if(class(snvS_Desc)!="try-error" ){
      if(class(snvS_noDesc)=="try-error" || max(snvS_noDesc$p,na.rm=T)<max(snvS_Desc$p,na.rm=T)){
        tmp=matrix(0,length(freq),1);        
        tmp[freq<=f_CNV+precision/2]=snvS_Desc$p;  snvS_Desc$p=tmp;  ##Complement to cover entire spFreq space
        snvS=snvS_Desc;
      }
    }
  
    maxP_J=0; ##Maximum probability from joined fit
    maxP_S=0; ##Maximum probability from separate fit
    maxP_SbeforeC=0; ##Maximum probability from separate fit, under the assumption that CNV happened in descendant of SP with SNV
    
    if(class(snvJ)!="try-error" && any(!is.na(snvJ$p))){
      maxP_J=max(snvJ$p,na.rm=T)
    }
    if(class(snvS)!="try-error" && any(!is.na(snvS$p))){
      maxP_S=max(snvS$p,na.rm=T)
    }
    if(!is.null(snvSbeforeC) && class(snvSbeforeC)!="try-error" && any(!is.na(snvSbeforeC$p))){
      maxP_SbeforeC=max(snvSbeforeC$p,na.rm=T)
    }
    
    ##Skip if no solution found
    if (maxP_J==0 && maxP_S==0 && maxP_SbeforeC==0){
      dm[k,"SP"]=NA;
      next;
    }
    
    joinedFit=FALSE;
    if (class(snvJ)!="try-error" && (dm[k,"PN_B"]==1 || maxP_J>=max(maxP_S,maxP_SbeforeC,na.rm=T))){ ##SP carrying SNV and SP carrying CNV have same size, i.e. are identical:
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
    
    ##Save end result:
    idx=which.min(abs(spFreq-freq[which.max(snv$p)]))
    dm[k,"SP"]=spFreq[idx];  
    idx=which(abs(snv$fit[,"f"]-dm[k,"SP"])<=precision/2); ##index of fits matching SP size
    #print(paste(k,idx))
    idx=idx[which.min(snv$fit[idx,"dev"])] ##index of fit of matching SP size with minimal residual (dev)
    #print(paste(k,idx))
    if(!isempty(idx)){
      dm[k,c("PM_B","PM")]=snv$fit[idx,c("PM_B","PM")]; ##(dm[k,"CN_Estimate"]*dm[k,"AF_Tumor"]-(1-dm[k,"SP"])*dm[k,"PN_B"])/dm[k,"SP"];  dm[k,"PM_B"]=max(1,dm[k,"PM_B"]);
    } else{
      print("empty")
      dm[k,"PM_B"] = NA
      dm[k,"SP_cnv"] = NA
      dm[k,"SP"] = NA
      next
    }
  
    if (!is.na(dm[k,"PM"]) && dm[k,"PM"]<0){
      dm[k,"PM"]=NA; ##PM can be -1 if obtained with snv_cnv_flag=1; TODO --> get NA directy for jar and remove this. 
    }
    dm[k,"%maxP"]=max(snv$p,na.rm=T); #snv$p[which.min(abs(freq-dm[k,"SP"]))];
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
      print(paste(dm[k,"PM_B"],dm[k,"PM_cnv"],dm[k,"PM_B"]))
      print(k)
      print("-----------------------------------")
      if(!is.na(dm[k,"SP"]) && !is.na(dm[k,"SP_cnv"])){
        print(paste("HERE:",dm[k,"SP"],dm[k,"SP_cnv"]))
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
  output=list("dm"=dm,"finalSPs"=finalSPs);
  return(output);
}


.addF<-function (dm,  max_PM){
  #.jaddClassPath("ExPANydS.jar")
  .jinit(classpath="ExPANdS.jar")
  #javaImport(packages = "core.analysis.ngs.algorithms.*")
  dm=.addColumn(dm,"f",NA);
  snv_cnv_flag=3; ##co-occurrence assumption of SNV and CNV
  
  for (k in 1:nrow(dm)){
    expands <-try(.jnew("ExPANdS", as.double(dm[k,"AF_Tumor"]),as.double(dm[k,"CN_Estimate"]),
                        as.integer(dm[k,"PN_B"]),as.integer(max_PM)));
    if (class(expands)=="try-error"){
      print(expands);
      print(paste("At SNV ",k,": -->"));
      print(dm[k,]);
    }else{
      .jcall(expands,,"run",as.integer(snv_cnv_flag))
      dm[k,"f"]<-.jcall(expands,"D","getF");
    }
  }
  return(dm);
}
