options(stringsAsFactors = FALSE)
Warnset=NULL

#flags that should be set at website
CGH.class=0       #CGH.class indicates which type of CGH is available  0=Full CGH  1= No CGH  2=Only HOMDEL and AMP   3=Only Gain and Loss 
Test.set=c("BN2","EZB","MCD","N1","ST2","A53")   #Types to predict
HasL265P=T                                      # are MYD88 L265P locations specified in flat file
Has.Trunc=T                                     # are Truncations indentified in flat file
Version=1   # set to 1 to use the original version, set to 2 for the improved version which may work better on incomplete data

#the names and locations of the user input files should be here
mutflat.loc="Mutation flat.txt"
copyflat.loc="Copy.flat.txt"
mutlist.loc="Mut gene list.txt"
copylist.loc="copy genelist"
samp.annot.loc="samp annot.txt"
Armfile.loc="Arm.txt"


Mutation.flat=read.table(mutflat.loc,header=T,sep="\t")    
if(dim(Mutation.flat)[2]==3)
{   Mutation.flat=data.frame(Mutation.flat,Location=-1)
}  
Copy.flat=read.table(copyflat.loc,header=T,sep="\t")   

Copy.flat=data.frame(Copy.flat,Location=-1)
names(Mutation.flat)=names(Copy.flat)
Flatfile=rbind(Copy.flat,Mutation.flat)



if(dim(Mutation.flat)[2]==3)
{     Mutation.flat=data.frame(Mutation.flat,-1)
}  
mut.genelist=read.table(mutlist.loc,header=T,sep="\t")               #Change to User file location
mut.genelist=mut.genelist[,1]




Sample.annot=read.table(samp.annot.loc,header=T,sep="\t")               #Change to User file location




subtyp=unique(Flatfile[,3])
subtyp=subtyp[!is.element(subtyp,c("AMP","GAIN","HETLOSS","HOMDEL","TRUNC","Synon","L265P","MUTATION"))]

if(length(subtyp)>0)
{ Warnset=c(Warnset,makewarn("Unknown mutation types included:",subtyp,6))  
}
mutset=is.element(Flatfile[,3],c("TRUNC","Synon","L265P","MUTATION"))
mutset=unique(Flatfile[mutset,2])
mutset=mutset[!is.element(mutset,mut.genelist)]
if(length(mutset)>0)
{ Warnset=c(Warnset,makewarn("Genes with mutations included in flat file not found in mut.genelist:",mutset,6))
}



cgset=is.element(Flatfile[,3],c("AMP","GAIN","HETLOSS","HOMDEL"))
if(CGH.class!=1)
{ CGH.genelist=mut.genelist
  if(file.exists(Armfile.loc))                                                    #Change to User file location
  { Arm.file=read.table(Armfile.loc,header=T,sep="\t")                                             #Change to User file location
  }
  else
  {Arm.file=NULL
  } 
  cgset=unique(Flatfile[cgset,2])
  cgset=cgset[!is.element(cgset,CGH.genelist)]
  if(length(mutset)>0)
  { Warnset=c(Warnset,makewarn("Genes with copy number changes included in flat file not found in CGH.genelist:",cgset,6))  
  }
}
if(CGH.class==1)
{ Arm.file=NULL
  if(sum(cgset>0))
  {  Warnset=c(Warnset,"Copy number changes found in flat file but not included in analysis")}

}  



#  Inputs that are static for from the training data
#Fullmat.2020=Fullmat.save

Fullmat.2020=read.table("Fullmat.2020.txt",header=T,sep="\t")
Fullmat.2020=as.matrix(Fullmat.2020==1)

Index.Flat=read.table("Index.Flat.txt",header=T,sep="\t")
Study.samp=read.table("Study.samp.txt",header=T,sep="\t")

if(Version==1)
{  originalpred=read.table("originalpred.txt",header=T,sep="\t")
}

if(Version==2)
{  originalpred=read.table("originalpredV2.txt",header=T,sep="\t")
}


  
Full.Uber.2020=read.table("Full.Uber.2020.txt",header=T,sep="\t")
Genclass.lab=Study.samp[,18]
Genclass.lab[Genclass.lab=="SDT"]="ST2"
Featord.2020=read.table("Featord.2020.txt",header=T,sep="\t")      #CHANGE version2.0
Featord.2020[,11]=Featord.2020[,11]==1            #CHANGE version2.0

mut.genlist=CGH.genelist=unique(Full.Uber.2020[,3])
#mut.genelist=phx[,2]
# Run predictor
result=Predict9(Flatfile,Sample.annot,Index.Flat,Genclass.lab,CGH.class=CGH.class,Full.Uber.2020,Fullmat.2020,originalpred=originalpred,Has.CGH=Study.samp[,7]==1,Arm.file=Arm.file,CGH.genelist=CGH.genelist,mut.genelist=mut.genelist,CalcWilcox=F,Test.set=Test.set,Has.Trunc=Has.Trunc,Featord.2020=Featord.2020,Version=Version)  #CHANGE version2.0


#output results
write.table(result$Prediction,"Result.txt",quote=F,sep="\t",row.names=F,na="")
write.table(result$Compare,"Compare.txt",quote=F,sep="\t",row.names=F,na="")

Warnset=c(Warnset,result$Warn.set)
if(length(Warnset)==0)
{  Warnset="None"
}  
write.table(Warnset,"warn.txt",quote=F,sep="\t",row.names=F,col.names=F)



