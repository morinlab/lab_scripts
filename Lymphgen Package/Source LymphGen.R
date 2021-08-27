

Predict9=function(Flatfile, Sample.annot,Index.List,Genclass.lab,Full.Uber.2020,Fullmat.2020,Featord.2020, Arm.file=NULL,CGH.class=0, Has.Trunc=T,mut.genelist=unique(Full.Uber.2020[,3]),CGH.genelist=unique(Full.Uber.2020[,3]),Test.set=c("BN2","EZB","MCD","N1","ST2","A53"),originalpred,ForceCGH=c(6689,5966,1029,4853,7874),Has.CGH=rep(T,length(lab)),CalcWilcox=F,report.N=T,mincount=c(2,2,2,1,2,2),forcedrop=rep(F,dim(Full.Uber.2020)[1]),Version=2)   #CHANGE version2.0
{	# Flatfile is a dataframe with the following columns in this order:  SampleID, geneID, Feature.type.  Feature Type will only include the following values  {"GAIN","AMP","HETLOSS","HOMDEL","MUTATION","TRUNC","Synon","L265P"}, Mutation location in HG38 coordinates
  # Armfile is a flat file that has the following columns in this order:  Sample ID, Armname, Feature.type  Where Armname will be of the form (1q, 1p or 1Chrom)
  # Sample.annot is a data frame with the following columns in this order:  Sample ID, CGH= {0,1} BCL2 = {0,1,NA}, BCL6 = {0,1,NA}
  # CGH.class indicates which CGH types are allowed.  0 = full CGH  1=No CGH  2 = HOMDEL and AMP only 3=HETLOSS and GAIN only
  # mut.genelist and CGH.genelist is the list of genes to be included in the mutation and CGH features.  It defaults to the total list 
  # Test.set is the list of categories to be considered
  # Index.List is the feature type by abnormality index file
  # Genclass.lab is the result of the Genclass algorithm on the training data. In my database its Study.samp[,18]
  #  Full.Uber2020 and Full.Mat2020 are the data files for the training data
  # originalpred is the original prediction matrix that is used to evaluate submodel prediction ability = Final.paper.pred[,5:10]
  # Force CGH is a list of genes to be defaulted to within the CGH.  Leave as default unless told otherwise
  # Has.CGH is in an indicator of whether a training sample has CGH. By default it should be Study.samp[,7] 
  # CalcWilcox indicates whehter or not to calculate the Cross validated Wilcox test.
  
  Warn.set=NULL
  
  if((Has.Trunc)&(sum(Flatfile[,3]=="TRUNC")==0))
  {	Warn.set=c(Warn.set,"Truncations indicated as available but not included in Flat file")
  }
  
  if(CGH.class==1)					#if CGH is not used set all samples to no CGH
  {	Sample.annot[,2]=0
  Test.set=Test.set[Test.set!="A53"]
  }
  
  if((sum(Sample.annot[,2])==0)&(is.null(Arm.file)))    # if no samples have CGH set CGH to not used
  {	CGH.class=1
  Test.set=Test.set[Test.set!="A53"]
  }	
  
  if(!is.element(4851,mut.genelist))
  {	Test.set=Test.set[Test.set!="N1"]
  Warn.set=c(Warn.set,"NOTCH1 mutations not available.  Prediction will exclude N1 subtype")
  }		
  
  
  lst=c("BN2","EZB","MCD","N1","ST2","A53")
  mod.lst=c("Full","NOCGH","NoBCL2Fus", "NoBCL2CGH" ,"NoBCL6Fus", "NoBCL6CGH", "NoFus", "NoFusCGH")
  
  Include.type=getindex(is.element(lst,Test.set))
  
  #Genclass.lab[!is.element(Genclass.lab,lst[Include.type])]="Other"
  
  
  sampnam=as.character(Sample.annot[,1])  
  
  
  badset=rep(F,dim(Flatfile)[1]) 
  if(Has.Trunc)
  {	badset[(is.element(Flatfile[,2],c(4853,4851)))&(is.element(Flatfile[,3],c("Synon","MUTATION")))]=T   #NOTCH2 and NOTCH1 must be Truncation  
  
  }  
  
  if(dim(Flatfile)[2]==4)
  {	Flatfile[is.na(Flatfile[,4]),4]=-1
  Flatfile[is.element(Flatfile[,3],c("AMP","HETLOSS","HOMDEL","GAIN")),4]=-1
  
  set=(Flatfile[,2]==4615)&(Flatfile[,4]==38182641)
  Flatfile[set,3]="L265P"
  
  badset[(Flatfile[,2]==4853)&(Flatfile[,4]>120459150)]=T  	#Notch2
  badset[(Flatfile[,2]==4851)&(Flatfile[,4]>139391455)]=T	#NOTCH1
  badset[(Flatfile[,2]==2416)&(Flatfile[,4]>148508764)]=T	#EZH2 lower
  badset[(Flatfile[,2]==2416)&(Flatfile[,4]<148506238)]=T	#EZH2 upper
  badset[(Flatfile[,2]==2416)&(Flatfile[,4]==-1)]=F	#EZH2 non-missing
  badset[(Flatfile[,2]==974)&(Flatfile[,4]>63007172)]=T
  badset[(Flatfile[,2]==974)&(Flatfile[,2]=="TRUNC")&(Flatfile[,4]>63006800)]=T
  
  }
  
  Flatfile=Flatfile[!badset,]
  
  
  
  
  set2=((Flatfile[,3]=="L265P")&(Flatfile[,2]!=4615))
  if(sum(set2)>0)
  {	Warn.set=c(Warn.set,"L265P indicated as type for non-MYD88 gene")
  }
  if((sum(Flatfile[,3]=="L265P")==0)&(sum((Flatfile[,2]==4615)>0))&(dim(Flatfile)[2]==3))
  {	Warn.set=c(Warn.set,"No MYD88mutation indicated as L265P: check to make sure MYD88 designation is indicated")
  }
  
  
  #set up MYD88.L265P	
  set=(Flatfile[,3]=="L265P")&(Flatfile[,2]==4615)
  Flatfile[set,2]=1000000004
  Flatfile[set,3]="MUTATION"
  
  if(sum(set)>0)
  {	mut.genelist=c(mut.genelist,1000000004)
  0
  }
  
  if(dim(Flatfile)[2]==3)
  {    Flatfile=data.frame(Flatfile,-1)
  }
  
  
  
  
  nsamp=length(sampnam)
  Model.index=2-Sample.annot[,2]+2*(is.na(Sample.annot[,3]))+4*(is.na(Sample.annot[,4]))    #indicator of which sample is run on which model
  modset=unique(Model.index)																			# indicator of which models to use
  Sample.annot=data.frame(Samp.index=1:nsamp,Sample.annot)
  
  modelcheck=unique(Model.index)
  if(is.element("A53",Test.set))
  {	modelcheck=c(modelcheck,9)
  0
  }	
  
  #Change 12-13-19
  
  names(Flatfile)[1:4]=c("sampnam","ENTREZ.ID","Type","Loc")
  sampnam=data.frame(Samp.index=1:nsamp,sampnam)
  

  
  if(CalcWilcox)                                                                   #prepare mutation frequency for permutation test
  {	tst=Flatfile[is.element(Flatfile$Type,c("MUTATION","TRUNC","SYNON")),]
  tst=tst[!duplicated(tst$ENTREZ.ID+tst$Samp.index/(2*max(tst$Samp.index))),]
  prob.wt=tabulate(tst$Samp.index)
  0
  }
  
  Flatfile=merge(Flatfile,Index.Flat[,c(7,2)])
  Flatfile=merge(Flatfile,Full.Uber.2020[,c(1,3,6,22,23,8,9)])
  
  badset=((Flatfile$Loc<Flatfile$Feat.start)|(Flatfile$Loc>Flatfile$Feat.end))
  gdset=(Flatfile$Full.Cat.index>200)&(is.element(Flatfile$Type,c("MUTATION","TRUNC")))
  badset[gdset]=((Flatfile$Loc[gdset]<Flatfile$Gene.start[gdset])|(Flatfile$Loc[gdset]>Flatfile$Gene.end[gdset]))
                                       
  badset[Flatfile$Loc==-1]=F
  Flatfile=Flatfile[!badset,c(1,2,3,4,6)]
  
  #add arms into flatfile
  
  if((CGH.class!=1)&(!is.null(Arm.file)))
  {	names(Arm.file)[1:3]=c("sampnam","Gene.Symbol","Type")
    Arm.file[,2]=sub("X","23",Arm.file[,2],ignore.case=T)    #use chromsome 23 and 24 instead of X and Y
    Arm.file[,2]=sub("y","24",Arm.file[,2],ignore.case=T)
    Arm.file=merge(Arm.file,Index.Flat[,c(7,2)])
    Arm.file=expand.arm(Arm.file)                                #  make sure that Chrom  implies p and q arms and vis versa   
    Arm.file=merge(Arm.file,Full.Uber.2020[,c(1,2,6)])
    set=nchar(Arm.file[,1])>5    #Chrom features
    Arm.file[set,1]=substring(Arm.file[set,1],first=1,last=nchar(Arm.file[set,1])-5)
    Arm.file[!set,1]=substring(Arm.file[!set,1],first=1,last=nchar(Arm.file[!set,1])-1)
    Arm.file[,1]=-as.integer(Arm.file[,1])
    names(Arm.file)=names(Flatfile)
    Flatfile=rbind(Flatfile,Arm.file)
    CGH.genelist=c(CGH.genelist,-(1:24))
  0
  }
 else
 {  forcedrop[Full.Uber.2020[,3]<0]=T
 }

  #End change 12-13-19
  
  
  
  
  
  badsamp1=!is.element(Flatfile$sampnam,sampnam$sampnam)
  if(sum(badsamp1)>0)
  { valx=unique(Flatfile$sampnam[badsamp1])
    cur.warn="Samples found in Flatfile but not in sample annotation file were excluded:"
    for(i in length(valx))
    { cur.warn=paste(cur.warn,valx[i])
    }
    Warn.set=c(Warn.set,cur.warn)
    Flatfile=Flatfile[!badsamp1,]
  }
  
  badsamp1=!is.element(sampnam$sampnam,Flatfile$sampnam)
  if(sum(badsamp1)>0)
  { valx=unique(Flatfile$sampnam[badsamp1])
    cur.warn="Samples found in Flatfile but not in sample annotation file were excluded:"
    for(i in length(valx))
    { cur.warn=paste(cur.warn,valx[i])
    }
    Warn.set=c(Warn.set,cur.warn)
    
  }
  
  Flatfile=merge(sampnam,Flatfile)
  
  Flatfile=Flatfile[,-1]
  flat.sub=data.frame(Samp.index=Flatfile$Samp.index,Flatfile$Feature.Index)

  
  #add BCL2 and BCL6 into flat file	
  
  if(sum(Sample.annot[,4],na.rm=T)>0)
  {	bc2=Sample.annot[,4]==1
  bc2[is.na(bc2)]=F
  bc2=data.frame(getindex(bc2),137)   #CHANGE version2.0   
  names(bc2)=names(flat.sub)
  flat.sub=rbind(flat.sub,bc2)
  0
  }
  
  if(sum(Sample.annot[,5],na.rm=T)>0)
  {	bc6=Sample.annot[,5]==1
  bc6[is.na(bc6)]=F
  bc6=data.frame(getindex(bc6),157)    #CHANGE version2.0   
  names(bc6)=names(flat.sub)
  flat.sub=rbind(flat.sub,bc6)
  0
  }
  
  #If sample has no features, add ineffectual ones
  mis.set=!is.element(sampnam[,1],flat.sub[,1])
  if(sum(mis.set)>0)
  { add.null=data.frame(sampnam[mis.set,1],-100)
    names(add.null)=names(flat.sub)
    flat.sub=rbind(flat.sub,add.null)
  }
    
  #generate models
  model.9=Model9(lab=Genclass.lab,Full.Uber.2020,Fullmat.2020,ForceCGH,CGH.class,mut.genelist,CGH.genelist,Has.CGH,modelcheck=modelcheck,Has.Trunc=Has.Trunc,forcedrop=forcedrop,Featord.2020=Featord.2020,Version=Version)   #CHANGE version2.0
  
  outmod=Nmod=matrix(NA,nsamp,6)
  
  #Predict samples	
  for(i in modset)
  {	cat("modcheckrun=",i,"\n")
    checkset=getindex(Model.index==i)
    cur=predict.flatfile(flat.sub,model=model.9[[i]],mincount=mincount)
    chk=cur$pmat[checkset,]
    outmod[checkset,1:5]=cur$pmat[checkset,]
    Nmod[checkset,1:5]=cur$nhit[checkset,]
    
    if(CalcWilcox)
    {	cur1=log(abs(cur$pmat)/(1-abs(cur$pmat)))
    newcheck=c(i,1)
    if(is.element(i,c(4,6,8)))
    {	newcheck=c(newcheck,2)
    }	
    if(is.element(i,c(4,7,8)))
    {	newcheck=c(newcheck,3)
    }	
    if(is.element(i,c(6,7,8)))
    {	newcheck=c(newcheck,5)
    }
    if(i==8)
    {	newcheck=1:8	
    }
    checkset2=	getindex(is.element(Model.index,newcheck))
    
    wilcoxval=Wilcoxrun(flat.sub,mod=model.9[[i]],matx=cur1,checkset=checkset2,prob.wt=prob.wt)
    model.9[[i]]$Mod.type=data.frame(model.9[[i]]$Mod.type,wilcoxval)
    }	
    
    
  }	
  if(is.element("A53",Test.set))
  {	checkset=getindex(is.element(Model.index,c(1,3,5,7)))
  if(length(checkset)>0)
  {	cur=predict.flatfile(flat.sub[is.element(flat.sub[,1],checkset),],model=model.9[[9]],mincount=mincount[6])
  outmod[checkset,6]=cur$pmat[checkset,]
  Nmod[checkset,6]=cur$nhit[checkset,]
  
  cur1=log(outmod[,6]/(1-outmod[,6]))
  
  atcheck=checkset
  if(CalcWilcox)
  {	wilcoxval=Wilcoxrun(flat.sub[is.element(flat.sub[,1],checkset),],mod=model.9[[9]],as.matrix(cur1),checkset)
  model.9[[9]][[1]]=data.frame(model.9[[9]][[1]],Wilcox=wilcoxval)		
  }
  }
  }
  outmod[,-Include.type]=NA
  
  
  
  originalpred[,-Include.type]=0
  call=combine.res(outmod,lst)
  outmod=abs(outmod[,Include.type])	
  #outmod=abs(outmod)	
  outmod=data.frame(trunc(outmod*100)/100)
  names(outmod)=c(paste("Confidence",lst[Include.type],sep="."))
  newAnnot=Sample.annot[,-1]
  newAnnot[,2]=ifelse(newAnnot[,2]==1,"Available","Not Available")
  newAnnot[,3]=ifelse(newAnnot[,3]==1,"yes","no")
  newAnnot[,4]=ifelse(newAnnot[,4]==1,"yes","no")
  newAnnot[is.na(newAnnot[,3]),3]="Not Available"
  newAnnot[is.na(newAnnot[,4]),4]="Not Available"
  newAnnot[,1]=as.character(newAnnot[,1])
  names(newAnnot)=c("Sample.Name","Copy.Number","BCL2.Translocation","BCL6.Translocation")
  
  NewCall=data.frame(newAnnot,Model.Used=mod.lst[Model.index],outmod)
  if(report.N)	
  {	Nmod=data.frame(Nmod[,Include.type])
  names(Nmod)=paste(lst[Include.type],"Feature.Count",sep=".")
  NewCall=data.frame(NewCall,Nmod)
  
  }	
  NewCall=data.frame(NewCall,Subtype.Prediction=call)
  
  #Check performance
  outCheck=NULL	
  Try.A53=rep(0,574)
  if(is.element("A53",Test.set))
  {	curcheck=getindex(Has.CGH)
  for(i in 1:length(curcheck))
  {	samp.val=model.9[[9]]$startdat[,i]
  cur2=predict.samp(samp.val,model.9[[9]],droptype=NULL)
  val=cur2[[1]][1]
  val[cur2[[2]][1]<mincount[6]]=0
  Try.A53[curcheck[i]]=val
  
  
  }
  }	
  
  Fullout=array(NA,c(574,6,8))
  
  for(j in modset)
  {	outmat=matrix(0,574,6)
  cur.original=originalpred
  if(j%%2==1)
  {	outmat[,6]=Try.A53
  
  }	
  else
  {	cur.original[,6]=0
  }	
  for(i in 1:574)
  {	if((j%%2==1)&!Has.CGH[i])
  {	Model=model.9[[j+1]]
  }
    else
    {	Model=model.9[[j]]
    }	
    samp.val=Model$startdat[,i]
    cur2=predict.samp(samp.val,Model)
    val=cur2[[1]]
    nv=length(val)
    val[cur2[[2]][1:nv]<mincount[1:nv]]=0
    
    outmat[i,1:5]=val[1:5]
  }
  outmat=outmat/(1+outmat)
  Fullout[,,j]=outmat
  curMatch=checkMatch(outmat,cur.original,lst)
  curMatch=data.frame(Model=mod.lst[j],curMatch)
  outCheck=rbind(outCheck,curMatch)
  
  }
  outCheck=data.frame(outCheck,Sensitivity=outCheck[,6]/(outCheck[,4]+outCheck[,6]),Specificity=outCheck[,3]/(outCheck[,3]+outCheck[,5]), Precision=outCheck[,6]/(outCheck[,5]+outCheck[,6]))
  

  
  c1=list(Prediction=NewCall,Compare=outCheck,model.9=model.9,Fullout=Fullout,Warn.set=Warn.set,modset=modelcheck)
  
  
  
  
  
}



Model9=function(lab,Full.Uber.2020,Fullmat.2020,ForceCGH=c(6689,5966,1029,4853,7874),CGH.class=0,mut.genelist=unique(Full.Uber[,3]),CGH.genelist=unique(Full.Uber[,3]),Has.CGH=rep(T,length(lab)),modelcheck=1:9,Has.Trunc=T,forcedrop=rep(F,dim(Full.Uber.2020)[1]),Featord.2020=Featord.2020,Version=2)    #CHANGE version2.0
{	#lab is the class calls for the training set in my database its Study.samp[,18]
  #Full.Uber.2020 and Full.Mat.2020  is the training data and annotation
  #ForceCGH  is a list of geneID's to force inclusion
  #CGH.class indicates which type of CGH is available  0=Full CGH  1= No CGH  2=Only HOMDEL and AMP   3=Only Gain and Loss 
  #mut.genelist and CGH.genelist are the list of Gene ID's available on the data set
  #Has.CGH is in an indicator of whether a training sample has CGH. By default it should be Study.samp[,7]
  #modelcheck  is a list of which models to generate
  #lab=Genclass.lab
  
  lab1=lab
  lab1[lab1=="A53"]="Other"
  
  CGH.default=CGH.class
  all.genelist=mut.genelist[is.element(mut.genelist,CGH.genelist)]
  custom.Drop=rep(T,dim(Full.Uber.2020)[1])
  custom.Drop[is.element(Full.Uber.2020[,3],mut.genelist)&is.element(Full.Uber.2020[,19],c(5,9))]=F
  custom.Drop[is.element(Full.Uber.2020[,3],CGH.genelist)&is.element(Full.Uber.2020[,19],c(1:4))]=F
  custom.Drop[is.element(Full.Uber.2020[,3],all.genelist)]=F
  custom.Drop[forcedrop]=T
  
  
  if(!Has.Trunc)
  {	custom.Drop[is.element(Full.Uber.2020[,19],c(5:8,13))]=T
  }	
  Featord=Featord.2020[is.element(Featord.2020[,2],1:5),]  #CHANGE version2.0
  
  if(sum(is.element(modelcheck,1))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.CGH=CGH.default,custom.Drop=custom.Drop,Version=Version)#CHANGE version2.0
  Model.Full=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)  #CHANGE version2.0
  }
  else
  {	Model.Full=NULL
  }
  
  if(sum(is.element(modelcheck,c(1,2)))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.CGH=1,custom.Drop=custom.Drop,Version=Version)#CHANGE version2.0
  Model.NoCGH=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)   #CHANGE version2.0
  }
  else
  {	Model.NoCGH=NULL
  }
  
  
  if(sum(is.element(modelcheck,c(5)))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.BCL6fus=T,drop.CGH=CGH.default,custom.Drop=custom.Drop,Version=Version)#CHANGE version2.0
  Model.NoBCL6Fus=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)  #CHANGE version2.0
  }
  else
  {	Model.NoBCL6Fus=NULL
  }
  
  if(sum(is.element(modelcheck,c(5,6)))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.BCL6fus=T,drop.CGH=1,custom.Drop=custom.Drop,Version=Version)    #CHANGE version2.0
  Model.NoBCL6CGH=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)    #CHANGE version2.0
  }
  else
  {	Model.NoBCL6CGH=NULL
  }
  
  
  if(sum(is.element(modelcheck,3))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.BCL2fus=T,drop.CGH=CGH.default,custom.Drop=custom.Drop,Version=Version)    #CHANGE version2.0
  Model.NoBCL2Fus=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)    #CHANGE version2.0
  }
  else
  {	Model.NoBCL2Fus=NULL
  }
  
  if(sum(is.element(modelcheck,c(3,4)))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.BCL2fus=T,drop.CGH=1,custom.Drop=custom.Drop,Version=Version)    #CHANGE version2.0
  Model.NoBCL2CGH=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)    #CHANGE version2.0
  }
  else
  {	Model.NoBCL2CGH=NULL
  }
  
  
  if(sum(is.element(modelcheck,c(7)))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.BCL6fus=T,drop.BCL2fus=T,drop.CGH=CGH.default,custom.Drop=custom.Drop,Version=Version)    #CHANGE version2.0
  Model.NoFus=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)    #CHANGE version2.0
  }
  else
  {	Model.NoFus=NULL
  }
  
  if(sum(is.element(modelcheck,c(7,8)))>0)
  {	Keep.Mat=Get.Featureset(Full.Uber.2020=Full.Uber.2020,ForceCGH=ForceCGH,drop.BCL6fus=T,drop.BCL2fus=T,drop.CGH=1,custom.Drop=custom.Drop,Version=Version)    #CHANGE version2.0
  Model.NoFusCGH=getModel(Fullmat.2020=Fullmat.2020,lab=lab1,Full.Uber.2020=Full.Uber.2020,Keep.Mat,Featord=Featord,Version=Version)    #CHANGE version2.0
  }
  else
  {	Model.NoFusCGH=NULL
  }
  
  
  if(is.element(9,modelcheck))	
  {	lab2=lab[Has.CGH==1]
  lab2[lab2!="A53"]="Other"
  Featord=Featord.2020[is.element(Featord.2020[,2],6),]  #CHANGE version2.0
  Featord[,2]=1
  Keep.Mat=Get.Featureset(Full.Uber.2020,exclude.Other=F,ForceCGH=c(6689,5966,1029,4853,7874),Stage=2,custom.Drop=custom.Drop,Version=Version)	  #CHANGE version2.0
  Model.A53=getModel(Fullmat.2020[,Has.CGH],lab2,Full.Uber.2020,Keep.Mat,exclude.Other=F,Featord=Featord,Version=Version)    #CHANGE version2.0
  }	
  else
  {	Model.A53=NULL
  }	
  
  
  model.9=list(Model.Full, Model.NoCGH, Model.NoBCL2Fus, Model.NoBCL2CGH ,Model.NoBCL6Fus, Model.NoBCL6CGH, Model.NoFus, Model.NoFusCGH, Model.A53)
  
  model.9
  
}



predict.flatfile=function(flt,model,mincount=c(2,2,2,1,2))
{	mod1=model$Mod.scale
mod2=model$Mod.type

nsamp=max(flt[,1])
n1=max(mod2[,1])
mat=mcnt=matrix(0,nsamp,n1)
n2=dim(mod1)[2]
outmat=matrix(0,max(flt[,1]),n1)
for(j in 1:n1)
{	set=mod2[,1]==j
n=sum(set)
curmod=mod1[set,,drop=F]	
curidx=mod2[set,-(1:2),drop=F]

for(i in 1:n)
{	cat(i," ")
  cur=rep(curmod[i,1],nsamp)
  for(k in n2:2)
  {	if(!is.na(curidx[i,k-1]))
  {	set=flt[,2]==curidx[i,k-1]
  if(sum(set)>0)
  {	cur[flt[set,1]]=curmod[i,k]
  }	
  }
  }
  mat[,j]=mat[,j]+cur
  mcnt[,j]=mcnt[,j]+(cur>curmod[i,1])
}		 
}
outmat=exp(mat)/(1+exp(mat))
for(i in 1:n1)
{	outmat[(mcnt[,i]<mincount[i])&(outmat[,i]>0.5),i]=-outmat[(mcnt[,i]<mincount[i])&(outmat[,i]>0.5),i]
}

list(pmat=outmat,vmat=mat,nhit=mcnt)

}



combine.res=function(in.mat,lst)
{	in.mat[is.na(in.mat)]=0
  outmat=Quick.check(in.mat)
nsamp=dim(outmat)[1]
nlab=dim(outmat)[2]
call=rep("",nsamp)
for(i in 1:nlab)
{ 	set=outmat[,i]==1
call[set&(call!="")]=paste(call[set&(call!="")],"/",sep="")
call[set]=paste(call[set],lst[i],sep="")
}	
call[call==""]="Other"
call

}

Quick.check=function(outmat,cut1=0.5,cut2=0.9)
{	hit5=rowSums(outmat>cut1)
for(i in 1:dim(outmat)[2])
{	outmat[,i]=(outmat[,i]>cut2)|((outmat[,i]>cut1)&(hit5==1))
}	
outmat
}


checkMatch=function(new.call,original.call,lst)
{	
includeset=colSums(original.call)!=0
original.call=Quick.check(original.call)
new.call=Quick.check(new.call)
comp=1+original.call+2*new.call

comp=comp[,includeset]
lst=lst[includeset]
out=matrix(0,sum(includeset),4)
for(i in 1:4)
{	out[,i]=colSums(comp==i,na.rm=T)
}	
out=data.frame(out)
names(out)=c("true.neg","false.neg","false.pos","true.pos")
out=data.frame(Subtype=lst,out)
out

}	


Wilcoxrun=function(flt,mod,matx,checkset,nperm=10000,prob.wt=rep(1,nsamp))
{	mod1=mod$Mod.scale
mod2=mod$Mod.type

nsamp=dim(matx)[1]
flt=flt[is.element(flt[,1],checkset),]

ncheck=length(checkset)
prob.wt=prob.wt[checkset]

n1=max(mod2[,1])
n2=dim(mod1)[2]
outcheck=rep(0,dim(mod1)[1])

for(j in 1:n1)
{	opset=mod2[,1]==j
n=sum(opset)
curmod=mod1[opset,]	
curidx=mod2[opset,-(1:2)]
for(i in 1:n)
{	cat(i," ")
  cur=rep(curmod[i,1],nsamp)
  for(k in n2:2)
  {	if(!is.na(curidx[i,k-1]))
  {	set=flt[,2]==curidx[i,k-1]
  if(sum(set)>0)
  {	cur[flt[set,1]]=curmod[i,k]
  }	
  }
  }
  set=unique(flt[is.element(flt[,2],as.matrix(curidx[i,1:(n2-1)])),1])
  set=is.element(checkset,set)
  
  if(min(mean(set),mean(!set))>0)
  {	curval=matx[,j]-cur
  curval=round(curval,5)
  
  val.check=rank(curval[checkset])
  wilcoxval=sum(val.check[set])
  totbet=0
  for(i in 1:nperm)
  {	tstval=sum(val.check[sample(1:ncheck,sum(set),prob=prob.wt,minimal=F)])
  totbet=totbet+(tstval>=wilcoxval)
  }	
  
  outcheck[opset][i]=totbet/nperm
  }								
  else
  {	outcheck[opset][i]=NA
  }
  
}	

}
outcheck
}


expand.arm=function(Arm.file)
{	outarm=Arm.file
checkval=paste(Arm.file[,2],Arm.file[,4],sep="|")
for(i in 1:24)
{	set=Arm.file[,3]==paste(i,"Chrom",sep="")

if(sum(set)>0)
{	c1=Arm.file[set,]
c1[,3]=substring(c1[,3],first=1,last=nchar(c1[,3])-5)
c1a=c1b=c1
c1a[,3]=paste(c1a[,3],"p",sep="")
c1b[,3]=paste(c1b[,3],"q",sep="")
outarm=rbind(outarm,c1a,c1b)
}

setp=Arm.file[,3]==paste(i,"p",sep="")
setq=Arm.file[,3]==paste(i,"q",sep="")
matchset=setq&(is.element(checkval,checkval[setp]))
if(sum(matchset)>0)
{	c1=Arm.file[matchset,]
c1[,3]=paste(substring(c1[,3],first=1,last=nchar(c1[,3])-1),"Chrom",sep="")
outarm=rbind(outarm,c1)

}	

}

outarm	
}	


Get.Featureset=function(Full.Uber.2020,drop.CGH=F,drop.Fus=F,drop.Wes=F,exclude.Other=F,custom.Drop=rep(F,dim(Full.Uber.2020)[1]),ForceCGH=NULL,drop.BCL6fus=drop.Fus,drop.BCL2fus=drop.Fus,Stage=1,  Version=2)   #CHANGE version2.0  Substantial changes throughout function
{	# generates feature set from Feature annotation file and subtype indicators
	# custom.Drop allows a custom set of Features to be dropped
	# script is specific to exact format and order of Full.Uber.2020
	# ForceCGH indicates genesIDs which should be given preferential choice when generating the model 
	# Assumes Fullmat2020
	
	Full.Uber.2020=Full.Uber.2020[order(Full.Uber.2020[,1,drop=T]),]
	useset=!custom.Drop	
	nforce=length(ForceCGH)
	if((nforce>0)&(drop.CGH!=1))
	{	has.cgh=Full.Uber.2020[,12]
		for(i in 1:nforce)
		{  set1=(Full.Uber.2020[,3]==ForceCGH[i])&useset
		   if(sum(has.cgh&set1&!custom.Drop)>0)
	  	 { idx=getindex(set1)[1]
		  	 set2=Full.Uber.2020[,7]==Full.Uber.2020[idx,7]
		  	 set2=set2&(Full.Uber.2020[,8]<Full.Uber.2020[idx,9]+1000000)
		  	 set2=set2&(Full.Uber.2020[,9]>Full.Uber.2020[idx,8]-1000000)
		  	 set2[is.na(set2)]=F
		  	 useset[has.cgh&set2]=F
		  	 useset[has.cgh&set1&!custom.Drop]=T
		   }
		}
		
	}


		

	Keep.Mat=matrix(F,length(useset),5)								#Keepmat is matrix of logical variables indicating which features are used in which models.  
																			#rows indicate features columns indicate models in the order ("A53","BN2","EZB","MCD","N1","SDT")
																			
	if(Stage==1)																		
	{		useset[Full.Uber.2020[,3]<0]=F  #CGH arms	
	  
	  for(i in 1:5)
		{	Keep.Mat[,i]=useset
		}	
		Keep.Mat[is.element(Full.Uber.2020[,6],c(3,13,14,113,114,213,214,313,314)),]=F   #exclude single copy gain features from all but A53
    if(Version==1)
    {Keep.Mat[Full.Uber.2020[,6]==4,]=F     #orginal version drop HETLOSS
    }  
		
		set=Full.Uber.2020[,3]==604
		if(!drop.BCL6fus)
		{	
			Keep.Mat[set,]=F   #drop BCL6 non fusion features
		
		}
		Keep.Mat[set&(Full.Uber.2020[,6]==0),]=!drop.BCL6fus #include BCL6 Fusions 

		
    set=Full.Uber.2020[,3]==596
    if(!drop.BCL2fus)
    {	Keep.Mat[set,]=F   #drop BCL2 non fusion features

    }
    Keep.Mat[set&(Full.Uber.2020[,6]==0),]=!drop.BCL2fus #include BCL2 Fusions 
	
		if(drop.Wes)
		{	Keep.Mat[Full.Uber.2020[,13],]=F	
		}
	}
	else
	{	Keep.Mat=as.matrix(useset)
		Keep.Mat[Full.Uber.2020[,3]<0,1]=T			#include arms A53	
		if(sum(Full.Uber.2020[useset,3]==7157)>0)
	  {	Keep.Mat[(Full.Uber.2020[,7]==17)&(Full.Uber.2020[,12]==1),1]=F		#exclude non TP53 chrom 17 from A53	
		  Keep.Mat[Full.Uber.2020[,3]==7157,1]=T    #include TP53 in A53
		}
		Keep.Mat[custom.Drop,1]=F

	}	
	if(drop.CGH>0)
	{	if(drop.CGH==1)
		{	Keep.Mat[Full.Uber.2020[,12]==1,]=F
		}
		if(drop.CGH==2)  #only include Homdels and Amps
		{	Keep.Mat[is.element(Full.Uber.2020[,6],c(3,4,8,12,13,14,108,112,113,114,213,214,308,313,314)),]=F
		}
		if(drop.CGH==3)  #only include Loss and Gain
		{	Keep.Mat[is.element(Full.Uber.2020[,6],c(1,2,6,7,10,11,106,107,110,111,210,211,310,311)),]=F
		}

	}
	Keep.Mat
}



getModel=function(Fullmat.2020,lab,Full.Uber.2020,Keep.Mat,Featord,exclude.Other=F,drop.Dup=F,SD.adjust=1,alpha=0.001,minprev=0.2, Version=2)  #CHANGE version2.0
{ #Creates full model
	# Fullmat.2020 is a matrix of logicals indicating which features (rows) are included in which samples (columns)
	# lab is a vector indicating the class labels of the
	# Full.Uber.2020  feature annotation
	# Keep.Mat   indicator of which feature is useable in which model
	# exclude.Other   indicates whether "Other" samples are to be used in significance calculations
	# drop.Dup indicates whether to drop extra copies of genes that occur in more than 1 subtype model  
	# SD.adjust  number of standard deviations to adjust by   
	
  lab=as.character(lab)
  
  lst=unique(lab)                 #ordered list of subtypes with other as final element
  lst=lst[lst!="Other"]
  lst=lst[order(lst)]
  ntype=length(lst)
  lst=c(lst,"Other")
  
  if(exclude.Other)
  {	Subset=lab!="Other"
  0
  }	
  else
  {	Subset=rep(T,length(lab))
  0
  }	
  
  
  
  strtlab=lab
  
  dat=Fullmat.2020+0           #convert feature matrix from logical to integer
  
  ttlab=ttmat=NULL    # ttlab is matrix of model annotation   ttmat=matrix of model data
  
  #dat=dat[,goodset]
  
	for(i in 1:ntype)   #individually generate hierarchical models and data sets for each sbutype
	{	cat("Run",i,"\n")
		keepset=Full.Uber.2020[Keep.Mat[,i],1]
    Featord.use=Featord[Featord[,2]==i,]  #CHANGE version2.0
		lab=strtlab==lst[i]
		minnum=max(4,sum(lab)/10-1)
	 	Hier=getHier(dat,lab,Full.Uber.2020,keepset,minnum=minnum,droplossgain=i!=1 ,Subset=Subset,p1=alpha,minprev=minprev,Featord.use=Featord.use,Version=Version)         #Main script for individual subtype model generation   #CHANGE version2.0
		outlab=Hier$outlab
		outlab=data.frame(i,outlab)
		ttlab=rbind(ttlab,outlab)
		ttmat=rbind(ttmat,Hier$outmat)
	}
	


	#	Reformat model annotation 
	
	ttlab=data.frame(ttlab[,1:2],parseby(ttlab[,3],"|"))  
	for(i in 3:dim(ttlab)[2])
	{	ttlab[,i]=as.integer(ttlab[,i])
	}
	nindex=dim(ttlab)[2]-2
	tmp=ttlab[,3]
	tmp[tmp==0]=NA
	ttlab=data.frame(ttlab,Full.Uber.2020[tmp,c(2:4)])
	for(i in 4:(nindex+2))
	{	tmp=ttlab[,i]
	  tmp[tmp==0]=NA
	  ttlab=data.frame(ttlab,Full.Uber.2020[tmp,4])
	}	



	mxdpth=max(ttmat,na.rm=T)    # maximum number of submodels 

	Full.arry=array(0,c(dim(ttlab)[1],2,mxdpth+1))   #array of number of samples in various categories  
															 # First dimension is gene, Second dimesion is (in/out of subgroup) third dimension is subfeature within gene 

	

	for(i in 1:ntype)
	{	set=ttlab[,1]==i
		#lab=Study.samp[,15]==lst[i]
		lab=strtlab==lst[i]
	
  	for(j in 0:(mxdpth))
  	{	Full.arry[set,1,j+1]=rowSums(ttmat[set,lab,drop=F]==j,na.rm=T)
  		Full.arry[set,2,j+1]=rowSums(ttmat[set,Subset&!lab,drop=F]==j,na.rm=T)
  	}	
		
	}	

	set=(Full.arry[,1,]==0)&(Full.arry[,2,]>0)
	Full.arry[,1,][set]=0.25
	
	set=(Full.arry[,2,]==0)&(Full.arry[,1,]>0)
	Full.arry[,2,][set]=0.25
	

	#Modscale is the output scaling for each gene.
	Modscale=log(Full.arry[,1,]/Full.arry[,2,])-log(rowSums(data.frame(Full.arry[,1,]))/rowSums(data.frame(Full.arry[,2,,drop=F])))
	
	Modscale[Full.arry[,1,]+Full.arry[,2,]==0]=0
	
	qval=Full.arry
	qval[,1,]=	1-Full.arry[,1,]/rowSums(data.frame(Full.arry[,1,]))
	qval[,2,]=	1-Full.arry[,2,]/rowSums(data.frame(Full.arry[,2,]))
	adjst=sqrt(qval[,1,]/Full.arry[,1,]+qval[,2,]/Full.arry[,2,])    #1 standard deviation adjustment
	
	set=Modscale>0
	Modscale[set]=pmax(Modscale[set]-SD.adjust*adjst[set],0.0001)
	
	set=Modscale<0
	Modscale[set]=pmin(Modscale[set]+SD.adjust*adjst[set],-0.0001)
	
	
	
	
	

	
	if(drop.Dup)                                                #if mulitple models have the same gene, include the one that is most significanly positive in the core.
	{	chk=data.frame(1:dim(ttlab)[1],Full.Uber.2020[ttlab[,3],3],-Modscale[,2])
		chk=chk[order(chk[,3]),]
		chk=chk[!duplicated(chk[,2]),]
		chk=chk[order(chk[,1]),1]
		Modscale=Modscale[chk,]
		ttlab=ttlab[chk,]
		ttmat=ttmat[chk,]			
	}
			
		

	list(Mod.type=ttlab,Mod.scale=Modscale,startdat=ttmat,lst=lst)
}






getHier=function(dat,lab,ub,keepset,p1=0.001,pvcut2=0.05,minnum=3,droplossgain=F,Subset=rep(T,dim(dat)[2]),minprev=0.2,Featord.use, Version=2)   #Change version2.0: Substantial changes throughout function
{	#function for hierarchical model building on individual subtype 
	#  dat=data matrix  rows = features, columns= samples
	# lab binary indicator of whether the sample is inside/outside the current subset
	# ub feature annotation
	# keepset  list of indicies of useable features for this model
	# p1= minimal p-value for gene inclusion
	# pvcut2 = minimal pvalue for subFeature inclusion
	# minnum = minimum number of allowed instances for inclusion
	# droplossgain  = whehter to exclude single copy losses and gains as individual features unless accompanied by mutation feature
	# Subset over which to calculate significance  (e.g. whether "Other" samples should be excluded)
	
	#cat("drp",droplossgain,"\n")
	
	
	outlab=outmat=NULL	  #eventual storage place for Hierarchcial annotation and data
	
	

	
	val=Featord.use[,c(1,4:10)]                #ordered set of features from full data
	uby=ub[val[,1],]
	stillin=rep(T,dim(uby)[1])                           #indicators of which features are still potential candidates
	stillin[rowSums(dat[val[,1],],na.rm=T)<=minnum]=F 
	stillin[val[,7]<val[,6]]=F

	keepset2=is.element(val[,1],keepset)   # set of features useable in model 
	
	
	
	stillin[!keepset2]=F
  if(Version==2)
  {stillin[Featord.use[,11]]=T   #allow exemplar inclusion of original model genes
  }  
	
	run=0              #flags useful for error checking and in case of bug, keep script from getting out of hand.
	if(Version==1)
	{runn=33          #previous version capped at 34 genes for A53
	}  
	else
	{runn=100
	}  
	
	
	while((sum(stillin)>0)&(val[getindex(stillin)[1],8]<p1)&run<runn)   #while we still have feature with p<p1
	{	
		ubbest=getindex(stillin)[1]                #identify best feature
		
		set0=stillin&(uby[,5]==uby[ubbest,5])  #Features available as exemplars
		set=set0&keepset2	  # find all features available associated with that gene


		set[is.na(set)]=F
		check2=(val[set0,8]<p1)&(val[set0,7]>minprev)
		  if((sum(set)>0)&(sum(check2)>0))
  	  {curdat=dat[uby[set,1],]						  #data for features
		    curub=uby[set,]								  # annotation for features
		    if(sum(set)>1)
			  {	grab=grabHier(curdat,curub,lab,pvcut2,minnum,droplossgain=droplossgain,Subset=Subset)    # generate vector of heirarchical feature IDs
			  	cat("outub",grab$outub,"\n")
			    if((length(grab$outub)>0)&!is.na(grab$outub))			
				  {	outmat=rbind(outmat,grab$mat)
				    outlab=rbind(outlab,data.frame(run=run,idx=as.character(grab$outub),stringsAsFactors=F))
			 	    run=run+1   #CHANGE version2.0 
				   if(grab$hasCGH)  #CHANGE version2.0                                                                         #if CGH part of feature, exclude nearby CGH features
				   {	curub=uby[set0,][1,]
						  set2=(uby[,7]==curub[,7])&(uby[,8]<curub[,9]+15000000)&(uby[,9]>curub[,8]-15000000)&(is.element(uby[,19],c(1:4,6:8,10:14)))
					  	stillin[set2]=F
					 }
				 }
		  	}
  	 		else            
  		 	{ if((!is.element(curub[,6],c(3,4)))|!droplossgain)	     # If only 1 feature, skip hierarchy and just use it unless its a gain/loss and droplossgain==T
	  	 		{  outmat=rbind(outmat,dat[uby[set,1],,drop=F])
	  				addub=uby[set,,drop=F]
 	  				outlab=rbind(outlab,data.frame(run=run,idx=as.character(addub[,1]),stringsAsFactors=F))
		  			if(is.element(addub[,19],c(1:4,6:8,10:14)))
			  		{	set2=(uby[,7]==addub[,7])&(uby[,8]<addub[,9]+15000000)&(uby[,9]>addub[,8]-15000000)&(is.element(uby[,19],c(1:4,6:8,10:14)))
		  				stillin[set2]=F
	  				}
	  			}
		  	}
		  }
	  stillin[set0]=F
	  cat(run,ubbest,outlab[dim(outlab)[1],1],outlab[dim(outlab)[1],2],sum(set),sum(stillin),val[ubbest,8],"\n")
		#cat(ub[outlab[length(outlab)],1:5])
	
	}
	
	#fewer than 2 features are included bulk up model with non existent features to avoid dimension dropping errors
  if(is.null(outlab))
  {  outmat=matrix(0,2,dim(dat)[2])
	   outlab=data.frame(run=c(0,0),idx="0")
  }
	if(dim(outlab)[1]==1)
	{ outmat=rbind(outmat,t(rep(0,dim(dat)[2])))
	  outlab=rbind(outlab,data.frame(run=0,idx="0"))
	} 
	list(outlab=outlab,outmat=outmat)	           
	
}








grabHier=function(curdat,curub,lab,pvcut2,minnum=3,droplossgain=T,Subset)
{	#creates hierachy for single gene
	#assumes features are ordered according to statistical significance
	#curdat is data matrix for gene features
	#curub  is annotation for gene features
	# pvcut2 = minimal pvalue for subFeature inclusion
	# minnum = minimum number of allowed instances for inclusion
	# droplossgain  = whehter to exclude single copy losses and gains from possible feature components
	# Subset over which to calculate significance  (e.g. whether "Other" samples should be excluded)
	
	
	set1=!duplicated(curub[,19])          #include either sub on non-sub version of each feature based on statistical significance
	curub=curub[set1,]
	curdat=curdat[set1,,drop=F]
	stx=LabFish(curdat[,Subset,drop=F],lab[Subset])                  #stx= stats table including statistical (column 8) and effective (column 9) significance
	stx=data.frame(1:dim(stx)[1],stx,comboval(stx[,1:4]))
	
	pv=(stx[,8]<pvcut2)&(rowSums(curdat!=0,na.rm=T)>minnum)    # passes subfeature critera
	upset=is.element(curub[,19],c(1,3))&pv                     # upCGH features
	dnset=is.element(curub[,19],c(2,4))&pv						 # down CGH features
	if(min(sum(upset),sum(dnset))>0)								 # choose either up or down features based on single most statistically significant feature
	{	if(min(stx[upset,8])<min(stx[dnset,8]))
		{	cgset=upset
		}
		else
		{	cgset=dnset
		}	
	}
	else
	{	cgset=upset|dnset											#cgset  set of possible CGH subfeatures
	}	
	mutset=is.element(curub[,19],c(5,9))&pv				#mutset  set of possible Mutation subfeatures
	cmpset=is.element(curub[,19],c(6:8,10:14))&pv		#cmpset set of possible Mutation subfeatures

	
	use.bestcomp=F										#flag indicating whether to default to best comp
	
	#Several different cases are considered.  The results is the list object "out" wich contains a vector or feature indexs ordered from core to larger,
	# and a vector indicating the call for each sample 
	
	if((sum(cgset)>0)&(sum(mutset)>0))					#case where both CGH and Mutations are possibly included
	{	
		#order mutations features so that Truncs occur before Muations
		curub.mut=curub[mutset,]
		curdat.mut=curdat[mutset,,drop=F]
		stx.mut=stx[mutset,]
		ord=order(-curub.mut[,19])                
		curub.mut=curub.mut[ord,]
		curdat.mut=curdat.mut[ord,,drop=F]
		stx.mut=stx.mut[ord,]
		
		#order mutations features so that Amp/Hom.Dels  occur before  Gain/Hetloss
		curub.cg=curub[cgset,]
		curdat.cg=curdat[cgset,,drop=F]
		stx.cg=stx[cgset,]
		ord=order(-curub.cg[,19])
		curub.cg=curub.cg[ord,]
		curdat.cg=curdat.cg[ord,,drop=F]
		stx.cg=stx.cg[ord,]
		
		useboth=T                                  #indicate that both CGH and mutation heirarchies are to be used unless told otherwise
		
		#decide whether to use mutation or CGH as core based on largest effective signficance
		if(max(abs(stx.mut[,9]))>=max(abs(stx.cg[,9])))                      
		{	out1=submake(curub.mut,curdat.mut,lab,pvcut2,minnum=minnum,Subset=Subset)     #calculate mutation portion of heirarchy 
			setx=(is.na(out1$mat)|(out1$mat==0))								  # exclude mutation samples from CGH model
			stz=LabFish(curdat.cg[,setx&Subset	,drop=F],lab[setx&Subset])
			pvy=(stz[,7]<(pvcut2))&(rowSums(curdat.cg[,setx&Subset,drop=F],na.rm=T)>minnum)   #check to make sure CGH model still passes criteria
			if(sum(pvy)>0)
			{	out2=submake(curub.cg,curdat.cg[,setx,drop=F],lab=lab[setx],pvcut2,minnum=minnum,Subset=Subset[setx])
			}
			else    #proceed as if only Mutant side existed
			{	useboth=F
				cgset=rep(F,length(cgset))
			}		
		}
		else
		{	out1=submake(curub.cg,curdat.cg,lab,pvcut2,minnum=minnum,Subset=Subset)           #calculate mutation portion of heirarchy 
			setx=(is.na(out1$mat)|(out1$mat==0))  									 # exclude CGH samples from mutation model
			stz=LabFish(curdat.mut[,setx&Subset,drop=F],lab[setx&Subset])
			pvy=(stz[,7]<(pvcut2))&(rowSums(curdat.mut[,setx&Subset,drop=F],na.rm=T)>minnum)         #check to make sure CGH model still passes criteria
			if(sum(pvy)>0)
			{	out2=submake(curub.mut,curdat.mut[,setx,drop=F],lab[setx],pvcut2,minnum=minnum,Subset=Subset[setx])
			}
			else           #proceed as if only CGH side existed
			{	useboth=F    
				mutset=rep(F,length(mutset))
			}						
		}
		if(useboth)
		{	out=combsub(out1,out2,setx)   #combine two separate hierarchies	
			hascgh=T
		}
	}	
	
	
	if((sum(mutset)==0)&(sum(cgset)==0))  #if no individual subfeature passes criteria, use comp
	{	use.bestcomp=T
	}	
	
	if((sum(mutset)==0)&(sum(cgset)>0))    #if only CGH component passes
	{	hascgh=T
		if(sum(cmpset>0)&(min(stx[cmpset,8])<min(stx[cgset,8])))   #Check whether Comp feature is statistically better than all CGH features
		{	out1=submake(curub[cgset,],curdat[cgset,,drop=F],lab,pvcut2,minnum=minnum,Subset)      #make CGH model
			if(length(out1$outub)>1)                                                    #if CGH model includes two levels replace Het.loss/Gain with associated comp instead
		   {	setx=(is.na(out1$mat)|(out1$mat!=1))
				sety=cmpset&is.element(curub[,19],c(8,12:14))                           #only consider Het.loss or Gain comps
				if(sum(sety)>0)
				{	xmat=curdat[sety,setx,drop=F]
					xmat1=curdat[sety,setx&Subset,drop=F]
					sty=LabFish(xmat1,lab[setx&Subset])
					idx=getindex(sety)
					idx=idx[order(sty[,7])[1]]
					ub1=c(out1$outub[1],as.character(curub[idx,1]))
					mat=(out1$mat==1)+0
					mat[setx&curdat[idx,]&!is.na(curdat[idx,])]=2
					out=list(outub=ub1,mat=mat)
				}
				else
				{	use.bestcomp=T
				}	
			}
			else	
			{	use.bestcomp=T   #if only single CGHfeature use comp instead of single feature
			}
		}
		else
		{	out=submake(curub[cgset,],curdat[cgset,,drop=F],lab,pvcut2,minnum=minnum,Subset)  #if comp is statistically worse than CGH feature generate CGH feature model
			hascgh=T
		}
	}
	
	if((sum(cgset)==0)&(sum(mutset)>0))  #if only mutation component passes criteria
	{	if((sum(cmpset)>0)&(min(stx[cmpset,8])<min(stx[mutset,8])))        #Check whether Comp feature is statistically better than all CGH features
		{	out1=submake(curub[mutset,],curdat[mutset,,drop=F],lab,pvcut2,minnum=minnum,Subset)             #make mutation model
			if(length(out1$outub)>1)												 #if Mutation model includes Mutations and Truncs replace Mutaiton with associated comp instead
			{  setx=is.na(out1$mat)|(out1$mat!=1)
				sety=cmpset&is.element(curub[,19],c(10,11,12,14))				  #only consider mutation comps
				if(sum(sety)>0)
				{	xmat=curdat[sety,setx,drop=F]
					xmat1=curdat[sety,setx&Subset,drop=F]
					sty=LabFish(xmat1,lab[setx&Subset])
					idx=getindex(sety)
					idx=idx[order(sty[,7])[1]]
					ub1=c(out1$outub[1],as.character(curub[idx,1]))
					mat=(out1$mat==1)+0
					mat[setx&curdat[idx,]&!is.na(curdat[idx,])]=2
					out=list(outub=ub1,mat=mat)
					hascgh=T
				}
				else
				{	use.bestcomp=T
				}	
			}
			else	
			{	use.bestcomp=T  #if only single Mutationfeature use comp instead of single feature
			}
		}
		else
		{	out=submake(curub[mutset,],curdat[mutset,,drop=F],lab,pvcut2,minnum=minnum,Subset=Subset)
			hascgh=F	
		}
	}
	
	if(use.bestcomp)                              #If single comp is best use it
	{	curub=curub[cmpset,]
		curdat=curdat[cmpset,,drop=F]
		stx=stx[cmpset,]
		bst=order(stx[,8])[1]
		curub=curub[bst,]
		out=list(outub=as.character(curub[,1]),mat=curdat[bst,])
		hascgh=T
	}	
	
	ub1=out$outub
	mat=out$mat
	ubz=curub[is.element(curub[,1],ub1),]
	if(droplossgain&(sum(is.element(ubz[,19],c(5:14)))==0))              #if single copies are to be excluded and no mutation or comp is included
	{	ub1=ubz[!is.element(ubz[,19],c(3,4)),1]                           #  exclude single copy gains/losses from list 
		mat[!is.na(mat)&mat==2]=0												# in case of 2 level CGH exclude larger level
	}	
	if(length(ub1)>1)																# merge feature list into single string separated by "|"
	{	ub2=ub1[1]
		for(i in 2:length(ub1))
		{	ub2=paste(ub2,ub1[i],sep="|")
		}
		ub1=ub2	
	}	
	list(outub=ub1,mat=mat,hasCGH=hascgh)
}
	
	
	
#curub=curub.mut
#curdat=curdat.mut	
	
	
submake=function(curub,curdat,lab,pvcut2,minnum,Subset)
{	#make heirarchical Mutation or CGH list
	#curub is current data frame of features
	#curdat is current matrix of data calls
	#minnum is minimum allowed number of samples in subfeature
	#Subset is subest of samples used to calculate significance
	
	if(dim(curub)[1]==2)
	{	ord=order(curub[,19])   #order subfeatures so that core is first
		curub=curub[ord,]
		curdat=curdat[ord,]
		stx=LabFish(curdat[,Subset],lab[Subset])
		stx=data.frame(1:dim(stx)[1],stx,abs(comboval(stx[,1:4])))
		usebest=T
		if((stx[1,9]>stx[2,9]+.00001)&(stx[1,8]<pvcut2)&sum(curdat[1,Subset],na.rm=T)>minnum)  # require that core satisfies criteria and has more significant effect that higher order 
		{  set=!is.na(curdat[1,])&(curdat[1,]==0)	&Subset
			try=curdat[2,set,drop=F]
			sty=LabFish(try,lab[set])
			
			if((sty[,7]<pvcut2)&(sum(try,na.rm=T)>minnum))  #make sure that with core removed larger set still has sufficient size and signficance
			{	mat=curdat[1,]
				mat[curdat[2,]&set]=2
				out=list(outub=curub[,1],mat=mat)
				usebest=F
			}
		}
		if(usebest)   #if can't divide, use most statistically significant
		{	idx=order(stx[,8])[1]
			out=list(outub=curub[idx,1],mat=curdat[idx,])
		}
	
	}	
	else  #if only one feature use it.
	{	out=list(outub=as.character(curub[,1]),mat=curdat[1,])
	}	
	out
}

combsub=function(out1,out2,setx)
{	#combine Mutation and CGH submodels.
	#out1 is core, out2 is remaining samples
	
	mat=out1$mat
	n=max(mat,na.rm=T)
	mat2=out2$mat
	set=!is.na(mat2)&mat2>0
	mat2[set]=mat2[set]+n
	mat[setx][set]=mat2[set]
	list(outub=c(out1$outub,out2$outub),mat=mat)
}




comboval=function(x,freq,scale=1)
{	#calculate effective significance of feature list

	x[x==0]=0.25
	out=log(x[,1]*x[,4]/(x[,2]*x[,3]))	
	
	adj=sqrt(1/x[,1]+1/x[,2]+1/x[,3]+1/x[,4])  #adjust by one standard deviation
	set=out<0	
	out[set]=pmin(out[set]+scale*adj[set],0)
	set=out>0
	out[set]=pmax(out[set]-scale*adj[set],0)
	
	out	
}	




Fudge.arry=function(in.arry,frq,scale=0.25)
{	#adjusts frequency matrix according to beta-binomal model
	# with mean=freq and scale=alpha+beta
	in.arry[in.arry==0]=scale
	in.arry	
}
	
	


predict.samp=function(samp.val,Model,droptype=NULL)
{	# predict subtype for new sample
	
	Mod.type=Model$Mod.type[,1]
	Mod.scale=Model$Mod.scale
	
	badset=is.na(samp.val)           #exclude features not found on sample
	samp.val=samp.val[!badset]
	Mod.type=Mod.type[!badset]
	Mod.scale=Mod.scale[!badset,]
	maxmod=max(samp.val)
	ngen=length(samp.val)
	
	type.set=unique(Mod.type)              #list of subtypes
	type.set=type.set[order(type.set)]
	

	ntype=length(type.set)
	
	badtype=is.element(type.set,droptype)
	nhit=rep(0,ntype)
	outprob=rep(0,length(type.set))
	idx=rep(0,ngen)
	for(i in 0:maxmod)
	{	set=samp.val==i
		idx[set]=Mod.scale[set,i+1]
	}
	scr=rep(0,ntype)
	for(i in (1:ntype)[!badtype])
	{	set=Mod.type==type.set[i]
		#cat(sum(set),)
		scr[i]=exp(sum(idx[set]))
		nhit[i] = sum(samp.val[set] != 0)
	}
	scr=c(scr,1)  #add "Other" category
	#scr/(sum(scr))
	list(scr = scr, nhit = nhit)
	
}


Multifish.quk=function(dat)
{	#Performs a fishers exact test on a 2 by 2 table
	# dat is a matrix or data.frame with four columns.  Rows represent different tables to calculate
	# the columns represent the number of hits in the bins x11, x12, x21, x22
	# takes less memeory than Multifish, but only works for 2x2
	
	ngen=dim(dat)[1]
	if(!is.matrix(dat))
	{	dat=as.matrix(dat)
		0
	}
	N=rowSums(dat)
	revset=dat[,1]<dat[,4]
	if(sum(revset)>0)
	{	dat[revset,]=dat[revset,4:1]
		0
	}

	n00=dat[,1]
	n10=dat[,2]
	n01=dat[,3]
	n11=dat[,4]
	N1=n10+n11
	N2=n01+n11
	Nmin=pmin(N1,N2)
	Nmax=max(Nmin)
	tsub=ttot=rep(0,ngen)
	ttarg=lchoose(N,n10)+lchoose(N-n10,n01)+lchoose(N-n10-n01,n11)
	matchval=lchoose(N,N1)+lchoose(N,N2)

	for(m11 in 0:Nmax)
	{	#cat(m11,)
		useset=Nmin>=m11
		M=N[useset]
		m10=N1[useset]-m11
		m01=N2[useset]-m11
		ttst=lchoose(M,m10)+lchoose(M-m10,m01)+lchoose(M-m10-m01,m11)
		set=ttst-ttarg[useset]<10^-8
		ttst=exp(ttst-matchval[useset])
		tsub[useset][set]=tsub[useset][set]+ttst[set]
		0
	}	
	tsub
}


Multi3by2=function(dat)
{	#Performs a fishers exact test on a 2 by 2 table
	# dat is a matrix or data.frame with six columns.  Rows represent different tables to calculate
	# the columns represent the number of hits in the bins x11, x12, x13 x21, x22, x23
	# sided=2 indicated 2 sided test
	# sided=1 indicates a test for large values along the diagonal
	# sided=-1 indicates a test for large values along the off diagonal

	
	ngen=dim(dat)[1]
	dat=as.matrix(dat)
	N=rowSums(dat)
	
	N1=rowSums(dat[,c(1:3)])
	N2=rowSums(dat[,c(4:6)])
	rev=N1>N2
	tmp=dat[rev,1:3]
	dat[rev,1:3]=dat[rev,4:6]
	dat[rev,4:6]=tmp
	N0=rowSums(dat[,c(1:3)])
	N1=rowSums(dat[,c(4:6)])

	
	n00=dat[,1]
	n10=dat[,2]
	n20=dat[,3]
	n01=dat[,4]
	n11=dat[,5]
	n21=dat[,6]
	
	M0=n00+n01
	M1=n10+n11
	M2=n20+n21
	
	
		
	Mx1=min(max(M0),max(N1))
	rev1=M0>(pmax(M1,M2))

	{	t1=n00[rev1]
		t2=n01[rev1]
		n00[rev1]=n20[rev1]
		n01[rev1]=n21[rev1]
		n20[rev1]=t1
		n21[rev1]=t2
		
	}	
	
	rev1=M1>(pmax(M0,M2))
	{	t1=n10[rev1]
		t2=n11[rev1]
		n10[rev1]=n20[rev1]
		n11[rev1]=n21[rev1]
		n20[rev1]=t1
		n21[rev1]=t2
		
	}	
		

	M0=n00+n01
	M1=n10+n11
	M2=n20+n21
	

	
	tstval=lchoose(N,n00)+lchoose(N-n00,n10)+lchoose(N-n00-n10,n20)+lchoose(N1,n01)+lchoose(N1-n01,n11)
	Totval=lchoose(N,M0)+lchoose(M1+M2,M1)+lchoose(N,N0)

	downset=upset=rep(0,ngen)
	top=0
	for(m00 in 0:Mx1) 
	{	Mx2 = min(max(M1), max(N1 - m00))
		for(m10 in 0:Mx2) 
		{	cat(m00,m10,"\n")
			m01 = M0 - m00
			m11 = M1 - m10
			m20 = N0 - m10 - m00
			m21 = N1 - m11 - m01
			goodset = pmin(m00, m10, m20, m01, m11, m21) >= 0
			if(sum(goodset) > 0) 
			{	top=top+1
				goodset = getindex(goodset)
				val1=lchoose(N[goodset],m00)+lchoose((N-m00)[goodset],m10)+lchoose((N1+m20)[goodset],m20[goodset])+lchoose(N1[goodset],m01[goodset])+lchoose((m11+m21)[goodset],m11[goodset])				
				val2=val1-Totval[goodset]
				
				isup = (val1 - tstval[goodset])>10^-8

		
				if(sum(isup > 0)) 
				{	upset[goodset[isup]] = upset[goodset[isup]] + exp(val2[isup])
				}
				if(sum(!(isup > 0))) 
				{	downset[goodset[!isup]] = downset[goodset[!isup]] + exp(val2[!isup])
				}
			}
		}
	}		
	data.frame(downset,upset)
}



LabFish=function(data,lab,PV.only=F)
{	bad=is.na(lab)
	lab=lab[!bad]
	data=data[,!bad,drop=F]
	labset=unique(lab)
	labset=labset[order(labset)]
	nlab=length(labset)
	ngen=dim(data)[1]
	data=data+0
	if(is.element(nlab,c(2,3)))
	{	outmat=matrix(0,ngen,nlab*2)
		top=0
		for(i in 1:nlab)
		{	top=top+1
			outmat[,top]=rowSums(data[,lab==labset[i],drop=F]==1,na.rm=T)
		}	
		for(i in 1:nlab)
		{	top=top+1
			outmat[,top]=rowSums(data[,lab==labset[i],drop=F]==0,na.rm=T)
		}	
		if(nlab==2)
		{	Pv=Multifish.quk(outmat)
		}	
		else
		{	Pv=Multi3by2(outmat)
		}	
		if(!PV.only)
		{	rat1=outmat[,1:nlab,drop=F]
			rat2=outmat[,-(1:nlab),drop=F]
			rat=data.frame(rat1/(rat1+rat2))
			names(rat)=paste("prev",labset)
			outmat=data.frame(outmat)
			names(outmat)=paste(c(labset,labset),c(rep("Hit",nlab),rep("Miss",nlab)))
			Pv=data.frame(outmat,rat,Pv)
		}	
	}	
	else
	{	cat("Bad labnum",labset,"\n")
			Pv=rep(1,dim(data)[1])
			if(!PV.only)
			{	outmat=rowSums(data,na.rm=T)
				outmat=data.frame(outmat,0,rowSums(1-data,na.rm=T),0)
				prev=data.frame(outmat[,1]/(outmat[,1]+outmat[,3]),rep(NA,dim(data)[1]))
				names(outmat)=paste(c(labset,"NULL"),c(rep("Hit",2),rep("Miss",2)))
				Pv=data.frame(outmat,prev,Pv)
			}
	}	
	Pv
	
	
}


getindex=function(x)
{	#reports the positive indices of a logical vector
	(1:length(x))[x]
}



"parseby"=function(dat,reg,maxrun=100,fixed=T)
{	#creates multicolumn data frame out of character string by dividing columns according to regular or fixed expression delimater
  #number of columns is 1 plus the largest number of delimaters in any element of string vector
  # if vector elements don't contain identical numbers of delimaters, those with fewer delimaters will have missing values in later collumns.
  #dat is character vector
  #reg is regular expression
  
  
  dat=as.character(dat)
  options(stringsAsFactors=FALSE)
  
  n=length(dat)
  done=rep(F,n)
  ln=nchar(reg)
  run=0
  while((sum(!done)>0)&(run<maxrun))
  {	run=run+1
  c1=regexpr(reg,dat,fixed=fixed)
  set=c1>0
  set[is.na(set)]=F
  ot1=rep(NA,n)
  if(sum(set)>0)
  {	ot1=rep(NA,n)
  ot1[set]=substring(dat[set],first=1,last=c1[set]-1)
  dat[set]=substring(dat[set],first=c1[set]+ln)
  
  }
  set2=(!set)&(!done)
  ot1[set2]=dat[set2]
  if(run==1)
  {	out=ot1
  }	
  else
  {	out=data.frame(out,ot1)
  }
  done[!set]=T		
  cat(run,sum(set),sum(!done),"\n")
  }
  out	
}






"lchoose"=function(m,n)
{	# log of binomial coefficient
	lgamma(m+1)-lgamma(n+1)-lgamma(m-n+1)
}

makewarn=function(base,set,maxn)
{  if(length(set)>maxn)
   { set=c(set[1:maxn],"etc.")
   }
  
  for(i in 1:length(set))
  { base=paste(base,set[i])
  }  
  base
  
}
  
  
  
