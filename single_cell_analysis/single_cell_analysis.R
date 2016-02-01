
#setwd("~/Dropbox/Research and collaborations/PT_255 and plasma/R/single_cell")
setwd("/projects/rmorin/pt255_single_cell/analysis/")
#change the above line to wherever your matrix files are

#first batch, all T2
 batch1<-read.csv("AA_2015-12-16-AKJBB_VAFs_matrix",sep="\t",row.names=1)
 v=seq(1, 96, 2)
 t=seq(2, 96, 2)
 all.vaf=batch1[,v]
masked.vaf=all.vaf
all.cov=batch1[,t]
low.cov=all.cov<500
#play around with this  line to change how deep you need coverage to be
#low.cov=all.cov<50
masked.vaf[low.cov]=NA
odds=seq(1, 46, 2)
batch1.mat<-matrix(nrow=length(masked.vaf[1,]),ncol=23)
colnames(batch1.mat)=colnames(masked.vaf[,odds])
rownames(batch1.mat)=rownames(masked.vaf)

#add means to matrix
pos=1
for(i in odds){
pair=c(i,i+1)
means=apply(masked.vaf[,pair], 1, function(x){mean(x,na.rm=TRUE)})
batch1.mat[,pos] = means
pos=pos+1
}
batch1.nonan=batch1.mat
batch1.nonan[ is.nan(batch1.nonan) ] <- 0
batch1.round = apply(batch1.nonan,c(1,2),function(x){round(x*4)/4})

heatmap(batch1.round,scale='none',col=c("white",'yellow','orange',"red","brown"))

#second batch of data, also T2

t3.data=read.csv("AA_2015-12-04-AGV0P_VAFs_matrix",sep="\t",row.names=1)

t3.vaf=t3.data[,v]
t3.cov=t3.data[,t]

#low.cov=t3.cov<500
low.cov=t3.cov<100
T3.masked.vaf=t3.vaf
T3.masked.vaf[low.cov]=NA

odds=seq(1, 46, 2)

T3.mat<-matrix(nrow=length(T3.masked.vaf[,1]),ncol=23)
colnames(T3.mat)=colnames(T3.masked.vaf[,odds])
rownames(T3.mat)=rownames(T3.masked.vaf)
pos=1
for(i in odds){
	pair=c(i,i+1)
	means=apply(T3.masked.vaf[,pair], 1, function(x){mean(x,na.rm=TRUE)})
	T3.mat[,pos] = means
	pos=pos+1
}

T3.nonan=T3.mat
T3.nonan[ is.nan(T3.nonan) ] <- 0
T3.round = apply(T3.nonan,c(1,2),function(x){round(x*4)/4})


 heatmap(T3.round,scale='none',col=c("white",'yellow','orange',"red","brown"))


#three single cells, T3

three<-read.csv("AA_2015-11-09-AHFB6_VAFs_matrix",sep="\t",row.names=1)
v=seq(1, 12, 2)
t=seq(2, 12, 2)
all.vaf=three[,v]
masked.vaf=all.vaf
all.cov=three[,t]
low.cov=all.cov<500
Three.masked.vaf=all.vaf
Three.masked.vaf[low.cov]=NA

odds=seq(1, 6, 2)

Three.mat<-matrix(nrow=length(Three.masked.vaf[,1]),ncol=3)
colnames(Three.mat)=colnames(Three.masked.vaf[,odds])
rownames(Three.mat)=rownames(Three.masked.vaf)
pos=1
for(i in odds){
	pair=c(i,i+1)
	means=apply(Three.masked.vaf[,pair], 1, function(x){mean(x,na.rm=TRUE)})
	Three.mat[,pos] = means
pos=pos+1
}

Three.nonan=Three.mat
Three.nonan[ is.nan(Three.nonan) ] <- 0
Three.round = apply(Three.nonan,c(1,2),function(x){round(x*4)/4})

#untested, needs editing from here 
#sarah's fourth batch (all T3 cells)
t2.data<-read.csv("AA_2015-12-24-ALF73_VAFs_matrix",sep="\t",row.names=1) #sarah to edit
v=seq(1, 96, 2)
t=seq(2, 96, 2)
t2.vaf=t2.data[,v]
t2.cov=t2.data[,t]

low.cov=t2.cov<500
#low.cov=t3.cov<50
T2.masked.vaf=t2.vaf
T2.masked.vaf[low.cov]=NA

odds=seq(1, 46, 2)

T2.mat<-matrix(nrow=length(T2.masked.vaf[,1]),ncol=23)
colnames(T2.mat)=colnames(T2.masked.vaf[,odds])
rownames(T2.mat)=rownames(T2.masked.vaf)
pos=1
for(i in odds){
	pair=c(i,i+1)
	means=apply(T2.masked.vaf[,pair], 1, function(x){mean(x,na.rm=TRUE)})
	T2.mat[,pos] = means
	pos=pos+1
}

T2.nonan=T2.mat
T2.nonan[ is.nan(T2.nonan) ] <- 0
T2.round = apply(T2.nonan,c(1,2),function(x){round(x*4)/4})

#any genes you want to exclude from the figures
bad_genes=c("ARID1A","YEATS","OXTR","DND1","EME1","GXYLT1","HIST1H3H","IL1R1","KRT3","MN1","PCDH10","MYH4","OCA2")
bad_mask=!rownames(T3.round)%in%bad_genes

all.cells=cbind(batch1.round,T3.round,Three.round) #currently works, but need to add T2.round to include final batch of cells

#add fourth data set:
batch4<-read.csv("AA_2016-01-ALF4J_VAFs_matrix",sep="\t",row.names=1)

 v=seq(1, 96, 2)
 t=seq(2, 96, 2)
 all.vaf=batch4[,v]
masked.vaf=all.vaf
all.cov=batch4[,t]
low.cov=all.cov<500
#play around with this  line to change how deep you need coverage to be
#low.cov=all.cov<50
masked.vaf[low.cov]=NA
odds=seq(1, 46, 2)
batch4.mat<-matrix(nrow=length(masked.vaf[1,]),ncol=23)
colnames(batch4.mat)=colnames(masked.vaf[,odds])
rownames(batch4.mat)=rownames(masked.vaf)

#add means to matrix
pos=1
for(i in odds){
pair=c(i,i+1)
means=apply(masked.vaf[,pair], 1, function(x){mean(x,na.rm=TRUE)})
batch4.mat[,pos] = means
pos=pos+1
}
batch4.nonan=batch4.mat
batch4.nonan[ is.nan(batch4.nonan) ] <- 0
batch4.round = apply(batch4.nonan,c(1,2),function(x){round(x*4)/4})



#this needs to be updated to reflect the number of cells from each time point (azure for T2, blue for T3)
#cols=c(rep("azure2",12),rep("azure4",11),rep("azure2",12),rep("azure4",11),rep("blue",3),rep("azure2",3),rep("azure4",3),rep("blue",17))
#azure2 = CD20-, azure4 = CD20+
cols=c(rep("azure2",12),rep("azure4",11),rep("azure2",12),rep("azure4",11),rep("blue",3),rep("azure2",3),rep("azure4",3),rep("blue",17),rep("azure2",20),rep("azure4",3))
# heatmap(all.cells,scale='none',col=c("white",'yellow','orange',"red","brown"),ColSideColors=cols)


trunk.genes=c("TMEM218","IGF2BP2","TUSC5","DHDH","RTN1","TP53_1","TP53_2")
germ.genes=c("KRTA5-5","ZSCAN4","GPAM")

keep.genes=rownames(T3.round)%in%c(trunk.genes,germ.genes)
#trunk/germline only
#heatmap(all.cells[keep.genes,],scale='none',col=c("white",'yellow','orange',"red","brown"),ColSideColors=cols)


#all.vafs<-cbind(batch1.nonan,T3.nonan,Three.nonan,T2.nonan)
all.vafs<-cbind(batch1.nonan,T3.nonan,Three.nonan,T2.nonan,batch4.nonan)


boxplot(all.vafs[keep.genes,],las=2)

pruned.vaf=apply(all.vafs,2,cell_locut)
rownames(pruned.vaf)=rownames(all.vafs)

genewise.vaf.state = t(apply(pruned.vaf,1,cutoff))
colnames(genewise.vaf.state)=colnames(all.vafs)

keep =! colnames(genewise.vaf.state) == "AA.CD20.B8.01_S4"
#removed single cell with no data 
cols = cols[keep]
#heatmap(all.cells[keep,],scale='none',col=c("white",'yellow','orange',"red","brown"),ColSideColors=cols)


gene.order=c("ZSCAN4","GPAM","ID3","IGF2BP2","DERA","TUSC5","DHDH","TP53_1","NUP210","CP","SGK494","DPPA2","DDR2","URB2","AKAP9","ZNF45","GOLSYN/SYBU","LRP1","ZNF7","PPP1CA","UNC93B1","MS4A1_INS","TP53_2","MS4A1_SNV","SCMH1","SBK1","INPP4B","RELN","ARHGAP35","NR3C1")

#heatmap(genewise.vaf.state[gene.order,],scale='none',col=c("white","orange","red"),Rowv=NA)

#heatmap(genewise.vaf.state[gene.order,],scale='none',col=c("white","orange","red"),Rowv=NA,ColSideColors=cols,distfun = function(c) dist(c,method="canberra"),main="canberra")

o1=order(genewise.vaf.state["TP53_2",])
o2=order(genewise.vaf.state["TP53_1",])
 o3=order(genewise.vaf.state["MS4A1_SNV",])
 oa = o1 * -o2 * 1.5*o3

some.genes=c("IGFBP2","NR3C1","MS4A1_INS","MS4A1_SNV","TP53_1","TP53_2")
 #plot just the genes of interest
 heatmap(genewise.vaf.state[some.genes,],scale='none',col=c("white","orange","red"),Rowv=NA,ColSideColors=cols,distfun = function(c) dist(c,method="canberra"),main="canberra",Colv=oa)



#use actual mean VAFs
all.vafs=cbind(batch1.nonan,T3.nonan,Three.nonan,T2.nonan)


#to apply to all samples, to determine threshold to ignore VAFs below (completely), account for cell-specific noise/bleed etc

heatmap(genewise.vaf.state[gene.order,keep],scale='none',col=c("white","orange","red"),Rowv=NA,ColSideColors=cols,distfun = function(c) dist(c,method="binary"))

some.genes=c("ID3","IGF2BP2","TP53_2","TP53_1","MS4A1_INS","AKAP9","MS4A1_SNV","NR3C1")

t=rownames(genewise.vaf.state) %in% some.genes

#o1=order(genewise.vaf.state["TP53_2",keep])
#o2=order(genewise.vaf.state["TP53_1",keep])
#o3=order(genewise.vaf.state["MS4A1_SNV",keep])
#oa = o1 * -o2 * 1.5*o3
heatmap(genewise.vaf.state[t,keep],scale='none',col=c("white","orange","red"),Rowv=NA,ColSideColors=cols,distfun = function(c) dist(c,method="binary"))


#count cells positive for each mutation in each of the two samples
t2.cells = cols=="azure2" | cols =="azure4"
t3.cells = cols=="blue"

t2.mut.cells=genewise.vaf.state[,t2.cells]>0
t3.mut.cells=genewise.vaf.state[,t3.cells]>0

t3.genesums=apply(t3.mut.cells,1,sum)
t2.genesums=apply(t2.mut.cells,1,sum)
t3.percent = 100*t3.genesums/sum(t3.cells)
t2.percent = 100*t2.genesums/sum(t2.cells)

t2.order = order(t2.percent)
t3.order = order(t3.percent)
all.order = t2.order * -t3.order
heatmap(genewise.vaf.state[all.order,keep],scale='none',col=c("white","orange","red"),Rowv=NA,ColSideColors=cols,distfun = function(c) dist(c,method="binary"))

#this is based on manual curation of the phylogenetic tree, root, then Clone from T2 then clone from T3
gene.order=c("ZSCAN4","GPAM","ID3","IGF2BP2","DERA","TUSC5","TP53_2","TP53_1",
             "TMEM218","KRTAP5-5","DHDH","MS4A1_INS","AKAP9","ZNF45","CP","NUP210","DPPA2","GOLSYN/SYBU",
             "LRP1","PPP1CA","ZNF7","UNC93B1","URB2","DDR2","ARHGAP35",
             "MS4A1_SNV","SCMH1","RELN","INPP4B","SBK1","NR3C1")


for(i in c(1:length(gene.order)-1)){
  gene1 = gene.order[i]
  gene2 = gene.order[i+1]
  o = order_genes(gene1,gene2)
  print(paste(gene1,gene2,o))
}
#redone by comparing the mutations unique to adjacent gene pairs
gene.order=c("ZSCAN4","GPAM","ID3","IGF2BP2","DERA","TUSC5","TP53_2","TP53_1","DHDH",
             "TMEM218","KRTAP5-5","MS4A1_INS","AKAP9","ZNF45","CP","NUP210","PPP1CA","DPPA2","URB2","UNC93B1","GOLSYN/SYBU",
             "LRP1","ZNF7","DDR2","ARHGAP35",
             "MS4A1_SNV","SCMH1","RELN","INPP4B","SBK1","NR3C1")

gene.order=c("ZSCAN4","GPAM","ID3","IGF2BP2","DERA","TP53_2","TP53_1","UNC93B1","URB2","DDR2","ARHGAP35",
             "MS4A1_SNV","SCMH1","RELN","INPP4B","SBK1","NR3C1")


heatmap(genewise.vaf.state[rev(gene.order),keep],scale='none',col=c("white","orange","red"),Rowv=NA,ColSideColors=cols,distfun = function(c) dist(c,method="binary"))

cellcounts = rbind(t2.percent,t3.percent)
barplot(cellcounts[,gene.order],beside=TRUE,las=2,cex.names=0.7)

clone2.genes=c("MS4A1_INS","AKAP9","ZNF45","CP","NUP210","PPP1CA","DPPA2","URB2","UNC93B1","GOLSYN/SYBU","LRP1","ZNF7","DDR2")
clone3.genes=c("ARHGAP35","MS4A1_SNV","SCMH1","RELN","INPP4B","SBK1","NR3C1")

table(colSums(t2.mut.cells[clone2.genes,])>1)
clone2.cells = which(colSums(t2.mut.cells[clone2.genes,])>1)

table(colSums(t2.mut.cells[clone3.genes,])>1)
clone3.cells = which(colSums(t2.mut.cells[clone3.genes,])>1)

keep.genewise.vaf.state= genewise.vaf.state[,keep]
cell_locut = function(x){
	locut = mean(x[x>0.05 ],na.rm=TRUE) - sd(x[x>0.05 ],na.rm=TRUE)
	if(is.na(locut)){
		locut = 0.33
	}
	if (locut > 0.33){
		locut =0.33
	}
	col=rep(0,length(x))
	for(i in c(1:length(x))){
		#print(i)
		#print(locut)
		if(is.na(x[i])){
			col[i]=0
			next
		}
		if(x[i]< locut){
			col[i] = 0
		}
		else{
			col[i] = x[i]
		}
	}
	return(col)
}

cutoff = function(x){
 hicut = mean(x[x>0.05 ]) + sd(x[x>0.05] )
 locut = mean(x[x>0.05 ]) - sd(x[x>0.05 ])
 col = rep(0,length(x))
 if (is.na(hicut)){
 	return(col)
 }
 if (hicut > 0.8){
 	hicut = mean(x[x>0.05 & x < 0.9]) + sd(x[x>0.05 & x < 0.9] )
 	locut = mean(x[x>0.05 & x < 0.9]) - sd(x[x>0.05 & x < 0.9])
 	if(is.na(hicut)){
 		hicut = 0.66
 		locut = 0.33
 	}
 }
 #hard threshold for min locut
 if (locut < 0.2){
 	locut = 0.2
 }
 
col = rep(0,length(x))
for(i in c(1:length(x))){
	if(x[i]>= locut & x[i]< hicut){
		col[i] = 1
	}
	else{
		if(x[i] >= hicut){
			col[i] = 2
		}
	}
}
return(col)
}
function(gene1,gene2){
  one.only =sum( t2.mut.cells[gene1,] & !t2.mut.cells[gene2,])
  two.only = sum(!t2.mut.cells[gene1,] & t2.mut.cells[gene2,])
  if(two.only> one.only){return(2)}
  else{return(1)}
}

