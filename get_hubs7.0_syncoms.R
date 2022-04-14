### Creates and plots correlation networks by syncom treatment
### Creates random correlation networks by sub samplig to validate hub OTUs 
### Calculates and plots network metrics of random networks 

#Packages
library(Hmisc)
library(igraph)
library(rgexf)
library(ggplot2)
library(GGally)
library(extrafont)
library(tidyr)
library(Rmisc)
library(FSA)


#Gaphics
tema=theme(axis.text.x = element_text(color="black",size=7, angle=90,hjust=0.5,vjust=0,family = "Arial"),
           axis.text.y = element_text(color="black",size=7,family = "Arial"),
           axis.title = element_text(color="black",size=7,family = "Arial"),
           legend.text = element_text(color = "black",size=7,family = "Arial"),
           legend.key.size = unit(0.3,"cm"),
           panel.border =element_rect(color = "white", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=7, color="black",face="bold",family = "Arial"),
           strip.text.y = element_text(size=7, color="black",face="bold",family = "Arial"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")


colores=c("black","#29854E","#74A557","#B8C466","#FFE081", "#A02355","#C55DA1","#DB99EE", "#234EA0","#008CCF","#00C6D9")
colnet=c("#009E8E","sandybrown","gray25","gray75")

# Set PARAMETERS
setwd("C:/Users/victorfn/Desktop/Results_syncoms")
amplicon=c("16s","ITS")[1]
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
base=c("rdp")

s16=("C:/Users/victorfn/Desktop/Results_syncoms/16s")
its=("C:/Users/victorfn/Desktop/Results_syncoms/ITS")  
result_folder="redes_extra"
result.dir=paste(getwd(),result_folder,sep = "/")
dir.create(result.dir,showWarnings=F)

#Load Metadata
met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"

#Load data
#16S
potu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
colnames(potu)=gsub("X","",colnames(potu))
rownames(met)=met$AIM.sample.id
ptax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))

#ITS
fotu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
colnames(fotu)=gsub("X","",colnames(fotu))
rownames(met)=met$AIM.sample.id
ftax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))


#Rearrange taxonomy
rownames(ftax)=gsub("OTU","FOTU",rownames(ftax))
rownames(ptax)=gsub("OTU","POTU",rownames(ptax))
ftax$otu.id=rownames(ftax)
ptax$otu.id=rownames(ptax)
tax=rbind(ptax,ftax)

#Make sure to have the same samples in each data set
rownames(fotu)=gsub("OTU","FOTU",rownames(fotu))
rownames(potu)=gsub("OTU","POTU",rownames(potu))

potu=potu[,colnames(potu) %in% colnames(fotu)]
fotu=fotu[,colnames(fotu) %in% colnames(potu)]

otu=rbind(potu,fotu)

#colnames(potu)==colnames(fotu)

#generate metadata
met.2=met[colnames(potu),]
limit=3 ; thr=0.75

#Select samples
se=c("nov19","feb20","jul20","nov20")
tr=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")

#Construct networks and random networks 
result.metrics=data.frame()

for ( t in tr){
 
#empty objects for metrics   
dens=c();diam=c();modu=c();tran=c();mdeg=c();aspa=c();node=c();edge=c();bf=c();bfp=c();n.hub=c();bb=c();bbp=c();ff=c();ffp=c()

#Select samples
samples=rownames(met.2[met.2$season %in% se & met.2$treatment %in% t,])
otus=otu[,samples]
otus=otus[rowSums(otus >= limit) >= length(samples)*thr,]

#Generate correlations
cor=rcorr(t(otus), type = "spearman")
rom=cor$r
pvm=cor$P

#Filter data based on ro 
rom[rom > -0.6 & rom < 0.6] <- 0

#Remove NaN
rom[is.nan(rom)] <- 0

#Filter data based on P-value
for (i in 1:dim(rom)[1]){
  logic=!(pvm[i,]<=0.01)
  rom[i, logic]<-0
}

#generate igraph object
g=graph.adjacency(rom, mode = "undirected", weighted = TRUE, diag = FALSE)

#delate isolated vertices
g=delete.vertices(g, V(g)[degree(g) == 0])

#save
save(g, file = paste(result.dir,paste(t,"network","igraph",sep = "."),sep = "/"))


#Select the most connected OTUS
grado=degree(g, v = V(g), mode = c("all"),loops = F, normalized = FALSE)
grado=grado[grado>summary(grado)[5]]

#Select the most central OTUS
central=betweenness(g, v = V(g), directed = F, weights = abs(E(g)$weight),nobigint = F, normalized = FALSE)
central=central[central>summary(central)[5]]

#Select the hub OTUS
hub=names(grado[names(grado) %in% names(central)])

#Save previous data
grafica=data.frame(degree=degree(g, v = V(g), mode = c("all"),loops = F, normalized = FALSE),
                   betweenesscentrality=betweenness(g, v = V(g), directed = F, weights = abs(E(g)$weight),nobigint = F, normalized = FALSE))
grafica=cbind(grafica,tax[rownames(grafica),])

grafica$hub=0 ; grafica$central=0 ; grafica$grado=0
grafica[hub,"hub"]=1 ; grafica[names(central),"central"]=1 ; grafica[names(grado),"grado"]=1

#Prepare data frame for random networks
result=data.frame(row.names = rownames(grafica), PRS=rep(0, times=dim(grafica)[1]), HRS=rep(0, times=dim(grafica)[1]),
                  GRS=rep(0, times=dim(grafica)[1]),CRS=rep(0, times=dim(grafica)[1]),
                  NRS=rep(0, times=dim(grafica)[1]))

#Plot the network
vertex_attr(g, name="kingdom", index = V(g)) <- tax[V(g)$name,1]
edge_attr(g, name="peso",index=E(g)) <- ifelse(E(g)$weight>0, "grey50", "red")
a=ggnet2(g, size = 0, node.color="kingdom", edge.color = "peso",
       mode="kamadakawai")+
  geom_point(size=2,aes(fill=color),shape=21)+
  scale_fill_manual(values =  colnet)+ggtitle(t)

png(paste(result.dir,paste(t,"net","png",sep = "."),sep = "/"), 
    width = 900, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(a)
dev.off()


#RANDOM NETWORKS (100)

initial.seed=1587571314
the.seed=71314

for (seed in seq(the.seed, initial.seed, 15880000 )){ print((seed))

#Result data frame   
random=data.frame(row.names = rownames(grafica), PRS=rep(0, times=dim(grafica)[1]), HRS=rep(0, times=dim(grafica)[1]),
                  GRS=rep(0, times=dim(grafica)[1]),CRS=rep(0, times=dim(grafica)[1]),
                  NRS=rep(0, times=dim(grafica)[1]))
#Select samples
otus2=otu[,samples]
set.seed(seed)
otus2=otus2[sample(x=rownames(otus2),size= 1290),]
rotus=otus2[rowSums(otus2 >= limit) >= length(samples)*thr,]

#Generate correlations
cor=rcorr(t(rotus), type = "spearman")
rom=cor$r
pvm=cor$P

#Filter data based on ro 
rom[rom > -0.6 & rom < 0.6] <- 0

#Remove NaN
rom[is.nan(rom)] <- 0

#Filter data based on P-value
for (i in 1:dim(rom)[1]){
  logic=!(pvm[i,]<=0.01)
  rom[i, logic]<-0
}

#generate igraph object
r.g=graph.adjacency(rom, mode = "undirected", weighted = TRUE, diag = FALSE)
#delate isolated vertices
r.g=delete.vertices(r.g, V(r.g)[degree(r.g) == 0])

#Select the most connected OTUS
r.grado=degree(r.g, v = V(r.g), mode = c("all"),loops = F, normalized = FALSE)
r.grado=r.grado[r.grado>summary(r.grado)[5]]
#Select the most central OTUS
r.central=betweenness(r.g, v = V(r.g), directed = F, weights = abs(E(r.g)$weight),nobigint = F, normalized = FALSE)
r.central=r.central[r.central>summary(r.central)[5]]
#Select the hub OTUS
r.hub=names(r.grado[names(r.grado) %in% names(r.central)])

#PRS times present in random networks
random[rownames(random) %in% rownames(otus)[rownames(otus) %in% rownames(otus2)],"PRS"]=1
#HRS times is a hub in random networks
random[rownames(random) %in% rownames(otus2)[rownames(otus2) %in% r.hub],"HRS"]=1
#GRS times ins highly connected in random networks
random[rownames(random) %in% rownames(otus2)[rownames(otus2) %in% names(r.grado)],"GRS"]=1
#CRS times is highly central in random networks
random[rownames(random) %in% rownames(otus2)[rownames(otus2) %in% names(r.central)],"CRS"]=1
#NRS times is present but not a hub in random samples
random[random$PRS==1 & random$HRS==0, "NRS"]=1

#merge
result=result+random

#Calculate metrics
dens=c(dens,edge_density(r.g, loops = FALSE)) 
diam=c(diam,diameter(r.g, directed=F, weights=NA))
modu=c(modu,modularity(r.g, membership(cluster_fast_greedy(r.g,weights = abs(E(r.g)$weight))), weights = abs(E(r.g)$weight)))
tran=c(tran,transitivity(r.g, type = "undirected", vids = NULL, weights = abs(E(r.g)$weight), isolates = c("NaN")))
mdeg=c(mdeg,mean(degree(r.g, v = V(r.g), mode = c("all"),loops = F, normalized = FALSE)))
aspa=c(aspa,mean_distance(r.g, directed = TRUE, unconnected = TRUE))
node=c(node,length(V(r.g)))
edge=c(edge,length(E(r.g)))
n.hub=c(n.hub,length(r.hub))

#Determine bacterial-fungal edges
ver=as_data_frame(r.g, what = c("edges"))
one=ver[grep("POTU_",ver$from),]
two=dim(one[grep("FOTU_",one$to),])[1]
three=ver[grep("FOTU_",ver$from),]
four=dim(three[grep("POTU_",three$to),])[1]

bf=c(bf,c(two+four))
bfp=c(bfp,c((two+four)/dim(ver)[1]*100))

#Determine bacterial-bacteria edges
one=ver[grep("POTU_",ver$from),]
two=dim(one[grep("POTU_",one$to),])[1]
bb=c(bb,c(two))
bbp=c(bbp,c((two)/dim(ver)[1]*100))

#Determine fungal-fungal edges
one=ver[grep("FOTU_",ver$from),]
two=dim(one[grep("FOTU_",one$to),])[1]
ff=c(ff,c(two))
ffp=c(ffp,c((two)/dim(ver)[1]*100))


}

#Calculates % of hub, central and connected OTUs
result$HF=round(result$HRS/result$PRS*100,digits=1)
result$GF=round(result$GRS/result$PRS*100,digits=1)
result$CF=round(result$CRS/result$PRS*100,digits=1)

grafica$random.hub=grafica$hub
grafica$random.grado=grafica$grado
grafica$random.central=grafica$central
#sum(rownames(grafica)==rownames(result))==dim(grafica)[1]


#Filter OTUs that appears as hub, central or connected in more than 50% of the random networks
grafica[result$HF<50,"random.hub"]=0
grafica[result$GF<50,"random.grado"]=0
grafica[result$CF<50,"random.central"]=0

write.table(grafica, file = paste(result.dir,paste(t,"R_nodes","txt",sep = "."),sep = "/"),
            row.names = F, quote = F, sep = "\t" ,col.names = T)

write.table(result, file = paste(result.dir,paste(t,"random_networks","txt",sep = "."),sep = "/"),
            row.names = T, quote = F, sep = "\t", col.names = T)

#merge metrics
metrics=data.frame(filter=c(limit), network=t,
                   density=dens,diameter=diam,modularity=modu,clustering=tran,
                   av.degree=mdeg,av.short.phat=aspa,vertex=node,edges=edge,
                   p.e.edges=bfp, no.hubs=n.hub,p.p.edges=bbp,e.e.edges=ffp)

result.metrics=rbind(result.metrics,metrics)
}

write.table(result.metrics, file = paste(result.dir,paste("metrics","txt",sep = "."),sep = "/"),
            row.names = T, quote = F, sep = "\t", col.names = T)



#Alternative plot network for Figure 4
hub=data.frame()
cen=data.frame()
con=data.frame()
for (t in tr){
  
  file = paste(result.dir,paste(t,"network","igraph",sep = "."),sep = "/")
  
  print(t)
  load(file = paste(result.dir,paste(t,"network","igraph",sep = "."),sep = "/"))
  nodos = read.delim(paste(result.dir,paste(t,"R_nodes","txt",sep = "."),sep = "/"))
  nodos$treatment=t
  
  vertex_attr(g, name="kingdom", index = V(g)) <- tax[V(g)$name,1]
  vertex_attr(g, name="hub", index = V(g)) <- ifelse(names(V(g)) %in% nodos[nodos$random.hub==1,"otu.id"],"hub","no_hub")
  colnet=c("#009E8E","sandybrown","gray25","gray75")
  edge_attr(g, name="peso",index=E(g)) <- ifelse(E(g)$weight>0, "grey50", "red")
  
  
  a=ggnet2(g, size = 0, node.color="kingdom", edge.color = "peso", node.size = "hub",
           mode="fruchtermanreingold")+
    geom_point(aes(fill=color,size=size),shape=21)+
    scale_fill_manual(values =  colnet)+
    scale_size_manual(values = c(3.5,2))+
    ggtitle(t)
    
    png(paste(result.dir,paste(t,"net_2","png",sep = "."),sep = "/"), 
        width = 900, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(a)
    dev.off()
    
    
    
  hub=rbind(hub,nodos[nodos$random.hub==1,])
  cen=rbind(hub,cen[cen$random.hub==1,])
  con=rbind(hub,con[con$random.hub==1,])
  
}

write.table(hub, file = paste(result.dir,paste("hubs","txt",sep = "."),sep = "/"),
            row.names = F, quote = F, sep = "\t", col.names = T)



#Plot network metrics

niveles=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
result.metrics$network=factor(result.metrics$network, levels = niveles)

mean=summarySE(result.metrics, "vertex", "network")
ggplot(mean, aes(x=network, y=vertex, fill=network))+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=vertex-sd, ymax=vertex+sd), width=0.2,size=0.7)+
  scale_y_continuous(limits = c(0,120))+ylab("no. of OTUs")
write.table(mean,file = paste(result.dir,paste("vertex.metric.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

mean=summarySE(result.metrics, "edges", "network")
ggplot(mean, aes(x=network, y=edges, fill=network))+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=edges-sd, ymax=edges+sd), width=0.2,size=0.7)+
  scale_y_continuous(limits = c(0,600))+ylab("no. of edges")
write.table(mean,file = paste(result.dir,paste("edges.metric.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

mean=summarySE(result.metrics, "no.hubs", "network")
ggplot(mean, aes(x=network, y=no.hubs, fill=network))+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=no.hubs-sd, ymax=no.hubs+sd), width=0.2,size=0.7)+
  scale_y_continuous(limits = c(0,20))+ylab("no. of hubs")
write.table(mean,file = paste(result.dir,paste("no.hubs.metric.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

mean=summarySE(result.metrics, "p.e.edges", "network")
ggplot(mean, aes(x=network, y=p.e.edges, fill=network))+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=p.e.edges-sd, ymax=p.e.edges+sd),width=0.2,size=0.7)+
  scale_y_continuous(limits = c(0,50))+ylab("% interdomain edges")
write.table(mean,file = paste(result.dir,paste("p.e.edges.metric.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

mean=summarySE(result.metrics, "p.p.edges", "network")
ggplot(mean, aes(x=network, y=p.p.edges, fill=network))+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=p.p.edges-sd, ymax=p.p.edges+sd), width=0.2,size=0.7)+
  scale_y_continuous(limits = c(0,100))+ylab("% prok-prok edges")
write.table(mean,file = paste(result.dir,paste("p.p.edges.metric.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

mean=summarySE(result.metrics, "e.e.edges", "network")
ggplot(mean, aes(x=network, y=e.e.edges, fill=network))+tema+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=e.e.edges-sd, ymax=e.e.edges+sd), width=0.2,size=0.7)+
  scale_y_continuous(limits = c(0,60))+ylab("% euk-euk edges")
write.table(mean,file = paste(result.dir,paste("e.e.edges.metric.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

mean=summarySE(result.metrics, "density", "network")
ggplot(mean, aes(x=network, y=density, fill=network))+tema+
  geom_bar(stat = "identity",size=0.7,colour="black")+scale_fill_manual(values = colores)+
  geom_errorbar(aes(ymin=density-sd, ymax=density+sd), width=0.2,size=0.7)+
  scale_y_continuous()+ylab("density")
write.table(mean,file = paste(result.dir,paste("density.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)



##Perform multiple comparison of metrics

test=dunnTest(vertex ~ network, data = result.metrics, method = "bh")$res
write.table(test[grep("MOCK",test$Comparison),],
            file = paste(result.dir,paste("vertex.kruskal.test.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

test=dunnTest(edges ~ network, data = result.metrics, method = "bh")$res
write.table(test[grep("MOCK",test$Comparison),],
            file = paste(result.dir,paste("edges.kruskal.test.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

test=dunnTest(no.hubs ~ network, data = result.metrics, method = "bh")$res
write.table(test[grep("MOCK",test$Comparison),],
            file = paste(result.dir,paste("no.hubs.kruskal.test.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

test=dunnTest(p.e.edges ~ network, data = result.metrics, method = "bh")$res
write.table(test[grep("MOCK",test$Comparison),],
            file = paste(result.dir,paste("p.e.edges.kruskal.test.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

test=dunnTest(p.p.edges ~ network, data = result.metrics, method = "bh")$res
write.table(test[grep("MOCK",test$Comparison),],
            file = paste(result.dir,paste("p.p.edges.kruskal.test.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)

test=dunnTest(e.e.edges ~ network, data = result.metrics, method = "bh")$res
write.table(test[grep("MOCK",test$Comparison),],
            file = paste(result.dir,paste("e.e.edges.kruskal.test.txt",sep = "."),sep = "/"),
            sep = "\t",quote = F, row.names = F, col.names = T)



for (t in tr){
  print(t)
  print(length(rownames(met.2[met.2$season %in% se & met.2$treatment %in% t,])))
}
    
    
#Generates supplementary table of hubs 

hub=data.frame()

for (t in tr){
  print(t)
  nodos = read.delim(paste(result.dir,paste(t,"R_nodes","txt",sep = "."),sep = "/"))
  nodos$treatment=t
  hub=rbind(hub,nodos[nodos$random.hub==1,])
  
}

hub.c=data.frame(row.names=unique(hub$otu.id))
hub$otu.id=factor(hub$otu.id, levels = unique(hub$otu.id))
for (t in tr){
count=summary(hub[hub$treatment==t,"otu.id"])
hub.c=cbind(hub.c,as.data.frame(count))
} 

colnames(hub.c)=tr
hub.c=cbind(hub.c, tax[rownames(hub.c),])
write.table(hub.c, "hubs.syncoms.txt", sep="\t", row.names = F)


#Generates a heatmap of metrics (not in manuscript)


result=data.frame(row.names = niveles)
for (j in colnames(result.metrics)[3:14]){
  
  mean=summarySE(data=result.metrics, measurevar = j, groupvars = c("network"),na.rm = T)
  a=log2(mean[,j]/mean[mean$network=="MOCK",j])
  result[,j]=a
  
}

#annotation data
anot=data.frame(row.names = niveles,
                treatment = niveles)
col=list(treatment=c(MOCK=colores[1],
                     PFC=colores[2],
                     PFCF=colores[3],
                     PFCS=colores[4],
                     PFCT=colores[5],
                     AAP=colores[6],
                     CC=colores[7],
                     BC=colores[8],
                     RC1=colores[9],
                     RC2=colores[10],
                     RC3=colores[11]))


breaksList = seq(-2, 2, by = 0.02)
der=colorRampPalette(c("#004B40","#009E8E" ,"#89D9CF","white"))
izq=colorRampPalette(c("white","peachpuff2","sandybrown","chocolate4"))
my_pal=c(rev(izq(100)),rev(der(101)))


pheatmap(result[niveles,c(7,8,9,11,12,10,1,5,4,3,2,6)],cluster_rows = F, cluster_cols = F, color = my_pal,
         annotation_row = anot, annotation_colors = col,
         border_color = "white", breaks = breaksList, gaps_col = c(2,5,6,9), gaps_row = c(1,5,8))

