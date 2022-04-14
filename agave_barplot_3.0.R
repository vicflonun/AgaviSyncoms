### Generates relative abundance barplots by treatment and/or season for Figure 1 ###

#Packages
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(phyloseq)
library(tidyr)
library(dplyr)
library(FSA)
library(VennDiagram)
library(extrafont)

#Directories
setwd("C:/Users/victorfn/Desktop")  
result="C:/Users/victorfn/Desktop/Results_syncoms/Bars" #result directories
dir.create(result)

s16=("C:/Users/victorfn/Desktop/Results_syncoms/16s") #data directories
its=("C:/Users/victorfn/Desktop/Results_syncoms/ITS") 

#Graphics
colores.t=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
nivel=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
etiquetas=c("Mock","Functional","Functional-F","Functional-As","Functional-At","Photo-anox","Photo-oxi","Biofilm","RC1","RC2","RC3")

colores=c("#A74476","#F56AB0",
          "#804795","#DF87FF",
          "#524CA1","#9A94EF",
          "#2C79AE","#6EC4FF",
          "#369DA4","#6EF6FF",
          "#29A876","#43FFB5",
          "#55A054","#A5F77A",
          "#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A",
          "#A12F2F","#FF4848")

tema=theme(axis.text.x = element_text(color="black",size=7, angle=0,hjust=0.5,vjust=0,family = "Arial"),
           axis.text.y = element_text(color="black",size=7,family = "Arial"),
           axis.title = element_text(color="black",size=7,family = "Arial"),
           legend.text = element_text(color = "black",size=7,family = "Arial"),
           legend.key.size = unit(0.2,"cm"),
           panel.border =element_rect(color = "white", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=7, color="black",face="bold",family = "Arial"),
           strip.text.y = element_text(size=7, color="black",face="bold",family = "Arial"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")

# Set PARAMETERS
amplicon=c("16s","ITS")[2]
soil=F
bootstrap=0.80 #not implemented yet
filter=1000 #minimum read count for a valid sample
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
base=c("rdp")

met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"


#Bar plot by treatment 

#Load data
if (amplicon=="16s"){
  #16S
otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
colnames(otu)=gsub("X","",colnames(otu))
tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
#match=read.delim(file=paste("C:/Users/victorfn/Desktop/Results_syncoms/16s",paste("isolate.match","txt",sep = "."),sep = "/"), sep = " ")
} else if (amplicon=="ITS"){
  #ITS2
otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
colnames(otu)=gsub("X","",colnames(otu))
tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
}

#Remove not necessary samples
if (soil==T){
met.2=met[colnames(otu),]
met.2=met.2[met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]
} else if (soil==F){
  met.2=met[colnames(otu),]
  met.2=met.2[met.2$treatment!="VOC",]
  met.2=met.2[met.2$plant.compartment!="soil",]
  otu.2=otu[,colnames(otu) %in% rownames(met.2)]
}

#Select taxon 
t=c("domain","phylum","class","order","family","genus","otu.id")[2]

#Sub sample reads
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2))))

#Organize data
otu.s=as.data.frame(cbind(otu.s,tax))
otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu.2))
otu.t[is.na(otu.t[,t]),t]="Unclassified"
for (i in met.2$AIM.sample.id){
  otu.t[otu.t$sample == i, c( "block","square","lane","season","treatment","plant.compartment" )] = met.2[i,c( "block","square","lane","season","treatment","plant.compartment" )]
}

#Calculate relative abundance of each taxon in each treatment 
otu.t0 = otu.t %>% group_by(treatment,plant.compartment) %>% summarise(absab = sum(abs))
otu.t1 = otu.t %>% group_by(otu.t[,t],treatment,plant.compartment) %>% summarise(absab = sum(abs)) ; colnames(otu.t1)[1]=c("rank")
sum(otu.t1$treatment==otu.t0$treatment)==dim(otu.t1)[1]
sum(paste(otu.t1$treatment,otu.t1$plant.compartment,otu.t1$rank)==paste(otu.t0$treatment,otu.t0$plant.compartment,otu.t1$rank))==dim(otu.t1)[1]
otu.t1$relab=otu.t1$absab/otu.t0$absab*100
#for (i in unique(met$treatment)){print(sum(otu.t1[otu.t1$treatment==i,"relab"]))}

#Remove low abundant taxa for better plotting
#Prokaryota : order=0.8 phylum=0.1 
#Eukaryota : order=0.5 phylum=0 

orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
low = orden[orden$relab < 0.5, "rank"]
otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "Low abundant"
n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 

#Reorder taxa if needed
orden=orden[order(orden$relab,decreasing = T),]
niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "Low abundant"))
otu.t1$treatment=factor(otu.t1$treatment, level = nivel)

#Plot
ggplot(data=otu.t1, aes(y=relab, x=treatment, fill=rank))+
  geom_bar(stat="identity",width = 0.95,size=0.5)+
  scale_fill_manual(values = c(colores[c(1:(n-2))],"grey50","black"), name=t)+tema+
  facet_wrap(~plant.compartment,scales = "free")+ylab("Relative abundance (%)")


#Bar plot by season 

#Load data
amplicon=c("16s","ITS")[2]


if (amplicon=="16s"){
  #16S
  otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
  #match=read.delim(file=paste("C:/Users/victorfn/Desktop/Results_syncoms/16s",paste("isolate.match","txt",sep = "."),sep = "/"), sep = " ")
} else if (amplicon=="ITS"){
  #ITS2
  otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
}

#Remove not necessary samples
if (soil==T){
  met.2=met[colnames(otu),]
  met.2=met.2[met.2$treatment!="VOC",]
  otu.2=otu[,colnames(otu) %in% rownames(met.2)]
} else if (soil==F){
  met.2=met[colnames(otu),]
  met.2=met.2[met.2$treatment!="VOC",]
  met.2=met.2[met.2$plant.compartment!="soil",]
  otu.2=otu[,colnames(otu) %in% rownames(met.2)]
}


#Select taxon 
t=c("domain","phylum","class","order","family","genus","otu.id")[2]

#Sub sample reads
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2))))
#sum(rownames(otu.s)==row.names(tax))==dim(otu.s)[1]

#Organize data 
otu.s=as.data.frame(cbind(otu.s,tax))
otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu.2))
otu.t[is.na(otu.t[,t]),t]="Unclassified"
for (i in met.2$AIM.sample.id){
  otu.t[otu.t$sample == i, c( "block","square","lane","season","treatment","plant.compartment" )] = met.2[i,c( "block","square","lane","season","treatment","plant.compartment" )]
}

#block C will not be plotted
met.2=met.2[met.2$block!="C",]
otu.t=otu.t[otu.t$block!="C",]

#Calculate relative abundance of each taxon in each season 
otu.t0 = otu.t %>% group_by(season,plant.compartment) %>% summarise(absab = sum(abs))
otu.t1 = otu.t %>% group_by(otu.t[,t],season,plant.compartment) %>% summarise(absab = sum(abs)) ; colnames(otu.t1)[1]=c("rank")
sum(otu.t1$season==otu.t0$season)==dim(otu.t1)[1]
sum(paste(otu.t1$season,otu.t1$plant.compartment,otu.t1$rank)==paste(otu.t0$season,otu.t0$plant.compartment,otu.t1$rank))==dim(otu.t1)[1]
otu.t1$relab=otu.t1$absab/otu.t0$absab*100
#for (i in unique(met$treatment)){print(sum(otu.t1[otu.t1$treatment==i,"relab"]))}


#Remove low abundant taxa for better plotting
#Prokaryota : order=0.8 phylum=0.1 
#Eukaryota : order=0.5 phylum=0.0 
orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
low = orden[orden$relab < 0.0, "rank"]
otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "Low abundant"
n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 

#Reorder taxa if needed
orden=orden[order(orden$relab,decreasing = T),]
niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "Low abundant"))
otu.t1$season=factor(otu.t1$season, level = c("nov19","feb20","jul20","nov20"))

#Plot
tema$legend.position="none"
plot=ggplot(data=otu.t1, aes(y=relab, x=season, fill=rank))+
  geom_bar(stat="identity",width = 0.95,size=0.5)+
  #scale_fill_manual(values = c(pal(n-2),"grey50","black"), name=t)+tema+
  scale_fill_manual(values = c(colores[c(1:(n-2))],"grey50","black"), name=t)+tema+
  #facet_wrap(~plant.compartment,scales = "free")+
  ylab("Relative abundance (%)")+xlab("Season")+
  scale_x_discrete(label=c("t0","t3","t8","t12"))


png(paste(result,paste(amplicon,t,"barplot.png",sep = "."),sep = "/"), 
    width = 500, height = 700, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()

#Plot legend
tema$legend.position="right"
plot=ggplot(data=otu.t1, aes(y=relab, x=season, fill=rank))+
  geom_bar(stat="identity",width = 0.95,size=0.5)+
  #scale_fill_manual(values = c(pal(n-2),"grey50","black"), name=t)+tema+
  scale_fill_manual(values = c(colores[c(1:(n-2))],"grey50","black"), name=t, 
                    label=gsub("[a-z]:","",levels(otu.t1$rank)))+
  tema+
  #facet_wrap(~plant.compartment,scales = "free")+
  ylab("Relative abundance (%)")+xlab("Season")+
  scale_x_discrete(label=c("t0","t3","t8","t12"))+
  guides(fill=guide_legend(ncol=1))


png(paste(result,paste(amplicon,t,"barplot.legend.png",sep = "."),sep = "/"), 
    width = 2000, height = 2000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()

#Save the plotted taxa for enrichment figure
if(amplicon=="16s"){
  save(niveles, file = paste(s16, paste(t,"names",sep = "."),sep = "/"))
}else if(amplicon=="ITS"){
  save(niveles, file = paste(its, paste(t,"names",sep = "."),sep = "/"))
}

