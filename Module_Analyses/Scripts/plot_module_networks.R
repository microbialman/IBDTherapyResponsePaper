#Script to generate plots of the networks within each module

#use the igraph and qgraph libraries to plot networks
library(igraph)
library(qgraph)
library(RColorBrewer)
library(here)

#load in our data from the WGCNA analysis
#module assignments per gene
mods=read.csv(here("Data","Generated","discovery_WGCNAgeneInfo.csv"))
colnames(mods)[4]="GS.weight"
tpm=readRDS(here("Data","Generated","discovery_filtered_tpm.Rds"))

#TOM matrix
load(here("Data","Generated","TOM-block.1.RData"))
TOM=as.matrix(TOM)
colnames(TOM)=colnames(tpm)
rownames(TOM)=colnames(tpm)

#gene data to get names in place of Ensemble IDs
genenames=read.csv(here("Data","Input","discovery_tpm_counts.csv"))
genenames=genenames[,c(1,2)]

#function to plot an igraph network of a given module
netPlot=function(tom,genes,name,genenames,thresh,mods){
  #subset tom to just the genes in the module
  indices=which(colnames(tom)%in%genes)
  set=tom[indices,indices]
  
  #rename the genes to names from ids
  names=genenames$gene_name[match(colnames(set),genenames$gene_id)]
  colnames(set)=names
  rownames(set)=names
  
  #set threshold to only keep top edges
  set[abs(set)<thresh]=0
  
  #initiate igraph
  ig=graph_from_adjacency_matrix(set,mode="undirected",weighted=TRUE, diag=FALSE, add.colnames = "GeneName")
  #calculate node centrality scores
  ig=set_vertex_attr(ig,"Degree",value=degree(ig))
  #remove genes with no edges at threshold
  ig=delete.vertices(ig,V(ig)$Degree==0)
  #keep only top 100 most connected
  if(length(V(ig)$Degree)>100){
  ig=delete.vertices(ig,which(!rank(1/V(ig)$Degree,ties.method="max")<=100))}
  if(length(V(ig)>1)){
  #plot the network
  e=get.edgelist(ig)
  colRam=colorRampPalette(brewer.pal(9,"YlOrRd"))
  #specify plot location
  name=here("Plots","Networks",name)
  #color by correlation with inflammation
  corv=mods$GS.weight[match(V(ig)$GeneName,mods$geneName)]
  cols=colRam(length(corv))[rank(corv)]
  l=qgraph.layout.fruchtermanreingold(e,weights=E(ig)$weight,vcount = vcount(ig),area=20*vcount(ig)^2, repulse.rad=40+vcount(ig)^3)
  pdf(paste0(name,".pdf"),10,10)
  plot(ig,layout=l,vertex.size=5*log(degree(ig)),vertex.label=V(ig)$GeneName,vertex.label.cex=0.8,edge.width=5*log(1+E(ig)$weight),
       vertex.label.family="Helvetica",vertex.label.color="black",vertex.label.font=2,vertex.frame.color=NA,
       vertex.color = cols)
  dev.off()}
}

dir.create(here("Plots","Networks"))
#run the function for each module at thresholds 0.3, 0.2 and 0.1
for(i in unique(mods$moduleColor)){
  if(i!="grey"){
  dir.create(here("Plots","Networks",i))
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.3"),genenames,0.3,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.25"),genenames,0.25,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.2"),genenames,0.2,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.15"),genenames,0.15,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.1"),genenames,0.1,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.05"),genenames,0.05,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.025"),genenames,0.025,mods)
  netPlot(TOM,mods$geneid[which(mods$moduleColor==i)],paste0(i,"/",i,"_threshold_0.01"),genenames,0.01,mods)
    }
}

