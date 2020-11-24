suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(snow))
#devtools::install_github("wilkelab/ungeviz")
suppressPackageStartupMessages(library(ungeviz))
#script containing various helper functions

#return row max of dataframe
maxRow=function(x){
  maxes=c()
  for(i in 1:nrow(x)){
    maxes=c(maxes,max(x[i,]))
  }
  return(maxes)
}

#return row min of dataframe
minRow=function(x){
  mins=c()
  for(i in 1:nrow(x)){
    mins=c(mins,min(x[i,]))
  }
  return(mins)
}

#make rownames from col 1
rowNamer=function(x){
  rownames(x)=x[,1]
  x=x[,-1]
  return(x)
}

#split a string a get a given index back
splitDex=function(vals,splitter,index){
return(unlist(lapply(strsplit(vals,splitter),"[[",index)))
}

#function to map values to colour gradient
colorGradient=function(x, colors=colpal, colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

#collapse duplicate columns as means or as sums
mergDupCol=function(x,mean){
  res=matrix(nrow=nrow(x),ncol=length(unique(colnames(x))))
  for(i in 1:length(unique(colnames(x)))){
    subset=data.frame(x[,which(colnames(x)==unique(colnames(x))[i])])
    if(mean==TRUE){
      res[,i]=rowMeans(subset)
    }
    else{
     res[,i]=rowSums(subset) 
    }
  }
  res=data.frame(res)
  colnames(res)=unique(colnames(x))
  rownames(res)=rownames(x)
  return(res)
}


#function to remove observations with zero values across a given fraction of all samples
zeroFilt=function(datatab,fraction){
  zeros=rowSums(datatab==0)
  return(datatab[zeros<(fraction*ncol(datatab)),])
}

#convert counts to relative abundances (cols sum to 1)
#option to additionally log2(x+10e-6) transform
relAbund=function(table,trans){
  if(trans==FALSE){
  return(t(t(table)/colSums(table)))}
  else{
    table=t(t(table)/colSums(table))
    return(log2(table+10E-6))
  }
}

#function to merge all but top x rows (based on row means)
topX=function(table,n){
  table=table[order(rowMeans(table),decreasing = T),]
  Other=colSums(table[(n+1):nrow(table),])
  table=table[1:n,]
  table=rbind(table,Other)
  rownames(table)[n+1]="Other"
  return(table)
}

#function to plot a dendrogram of samples with colorbars of metadata
metaDend=function(obsdata,metadata,textsize,dendsize,margl,margr){
  if(all(rownames(obsdata)!=rownames(metadata))){
    warning("Make sure that rownames in the observations match the metadata.")
  }
  clustering=hclust(dist(obsdata))
  plots=list()
  dend=ggdendrogram(clustering)+theme_void()+theme(plot.margin=margin(t=0,r=margr,b=0,l=margl,unit="pt"))
  metadata=metadata[clustering$order,]
  continuous_variable=unlist(lapply(metadata, is.numeric))
  metadata=metadata[,c(which(continuous_variable),which(!continuous_variable))]
  for(i in 1:ncol(metadata)){
    df=data.frame(x=c(1:nrow(metadata)),fill=metadata[,i])
    bar=ggplot(df,aes(x,y=as.factor(1),fill=fill))+geom_tile()+
      scale_y_discrete(expand=c(0,0),labels=c("1"=colnames(metadata)[i]))+
      theme_void()+
      theme(axis.title=element_blank(),
            axis.ticks=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size = textsize),
            axis.line = element_blank(),
            legend.position="none",
            plot.margin = margin(0,0,0,0))
      plots[[paste0("bar",i)]]=bar
  }
  bargrid=plot_grid(plotlist=plots,nrow=length(plots),align="v")
  return(plot_grid(dend,bargrid,nrow=2,rel_heights =c(dendsize,1-dendsize),axis="l"))
}

#function to plot simple heatmaps from correlation matrices
simpHeat=function(cordat,pdat,pal){
  roword=hclust(dist(cordat))
  colord=hclust(dist(t(cordat)))
  cordat=cordat[roword$order,colord$order]
  pdat=pdat[roword$order,colord$order]
  melt=melt(cordat)
  pmelt=melt(pdat)
  pmelt$value[pmelt$value==0]=NA
  colnames(melt)=c("Var1","Var2","Correlation")
  ggplot(melt,aes(x=Var2,y=Var1,fill=Correlation))+geom_tile()+
    geom_rect(data=pmelt, size=0.5, fill=NA, aes(color=as.factor(pmelt$value), xmin=as.numeric(pmelt$Var2)-0.5,xmax=as.numeric(pmelt$Var2)+0.5,
              ymin=as.numeric(pmelt$Var1)-0.5,ymax=as.numeric(pmelt$Var1)+0.5))+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
      scale_fill_gradientn(colors=pal)+scale_color_manual(values=c("black",NA),labels=NULL,guide=F)
}

#function to plot simple heatmaps from correlation matrices, with 2 levels of p sig
simpHeatNomFdr=function(cordat,pdat,fdrdat,pal){
  roword=hclust(dist(cordat))
  colord=hclust(dist(t(cordat)))
  cordat=cordat[roword$order,colord$order]
  pdat=pdat[roword$order,colord$order]
  fdrdat=fdrdat[roword$order,colord$order]
  melt=melt(cordat)
  pmelt=melt(pdat)
  fmelt=melt(fdrdat)
  pmelt$value[pmelt$value==0]=NA
  fmelt$value[fmelt$value==0]=NA
  colnames(melt)=c("Var1","Var2","Correlation")
  ggplot(melt,aes(x=Var2,y=Var1,fill=Correlation))+geom_tile()+
    geom_rect(data=pmelt, size=0.5, fill=NA, aes(col=as.factor(pmelt$value), xmin=as.numeric(pmelt$Var2)-0.5,xmax=as.numeric(pmelt$Var2)+0.5,
                                                 ymin=as.numeric(pmelt$Var1)-0.5,ymax=as.numeric(pmelt$Var1)+0.5))+
    geom_point(data=fmelt, size=2, fill=NA, aes(col=as.factor(fmelt$value), xmin=as.numeric(fmelt$Var2)-0.5,xmax=as.numeric(fmelt$Var2)+0.5,
                                                 ymin=as.numeric(fmelt$Var1)-0.5,ymax=as.numeric(fmelt$Var1)+0.5))+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
    scale_fill_gradientn(colors=pal)+scale_color_manual(values=c("black",NA),labels=NULL,guide=F)
}

#function that performs DESeq2 LRT Multiple comparisons from raw counts, assuming a single grouping with no covariates
#returns number of signficant results as specified, ordered by p-val
DESeqLrt=function(counts,groups,nres){
  cdat=data.frame(Group=groups)
  dds=DESeqDataSetFromMatrix(countData = counts, colData = cdat,  design= ~ Group)
  dds=DESeq(dds, test="LRT", reduced=~1, parallel = T)
  res=results(dds)
  sig=data.frame(Variable=rownames(counts),pval=res$pvalue,padj=res$padj)
  sig=sig[sig$padj<0.05&!is.na(sig$padj),]
  sig=sig[order(sig$padj),]
  if(nres>nrow(sig)){
    return(sig)
  }
  else{
    return(sig[1:nres,])
  }
}

multiAnova=function(counts,groups,nres){
  cl=makeSOCKcluster("localhost")
  res=parLapply(cl,data.frame(t(counts),check.names = F),groups=groups,function(x,groups){r=anova(lm(x~groups));return(r$`Pr(>F)`[1])})
  res=unlist(res)
  res=res[order(res)]
  res=data.frame(Variable=names(res),p=res,check.names = F)
  return(res[1:nres,])
}






