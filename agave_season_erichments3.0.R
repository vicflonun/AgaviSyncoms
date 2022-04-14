### Performs OTU enrichment analysis by season and plots Triplots from Figure 1 ###

#Packages
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(Rmisc)
library(tidyr)
library(dplyr)
library(FSA)
library(ggtern)
library(extrafont)

#Graphics
colores=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
nivel=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
etiquetas=c("Mock","Functional","Functional-F","Functional-As","Functional-At","Photo-anox","Photo-oxi","Biofilm","RC1","RC2","RC3")

tema=theme(axis.text= element_text(color="gray30",size=7, angle=0,hjust=0.5,vjust=0,family = "Arial"),
           axis.title = element_text(color="black",size=7,family = "Arial"),
           legend.text = element_text(color = "black",size=7,family = "Arial"),
           legend.key.size = unit(0.8,"cm"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=7, color="black",face="bold",family = "Arial"),
           strip.text.y = element_text(size=7, color="black",face="bold",family = "Arial"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")


# Set parameters
amplicon=c("16s","ITS")[2]
bootstrap=0.80 #not implemented yet
filter=1000 #minimum read count for a valid sample
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
base=c("rdp")
t=c("domain","phylum","class","order","family","genus","otu.id")[4] #colors by
s=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")[] #Analyze only one treatment or all

# Directories
setwd("C:/Users/victorfn/Desktop")
s16=("C:/Users/victorfn/Desktop/Results_syncoms/16s")
its=("C:/Users/victorfn/Desktop/Results_syncoms/ITS")   

# Metadata
met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"

#Loading data 
if (amplicon=="16s"){
  #16S
  otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
  #match=read.delim(file=paste("C:/Users/victorfn/Desktop/Results_syncoms/16s",paste("isolate.match","txt",sep = "."),sep = "/"), sep = " ")
  load(file = paste(s16, paste(t,"names",sep = "."),sep = "/"))
  
  } else if (amplicon=="ITS"){
  #ITS2
  otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
  load(file = paste(its, paste(t,"names",sep = "."),sep = "/"))
}

# Removing unwanted samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

# Sub sample reads and calculate relative abundance  
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) #rarefication 
otu.r=t(t(otu.s)/colSums(otu.s)*100)
otu.r=cbind(otu.r,tax)
otu.t=gather(data=otu.r, key = "sample", value = "relab", colnames(otu.s))

for (i in met$AIM.sample.id){
  otu.t[otu.t$sample == i, c( "block","square","lane","season","plant.compartment","treatment" )] = met.2[i,c( "block","square","lane","season","plant.compartment","treatment" )]
} 


#Subset the data 
otu.t1=otu.t[otu.t$treatment %in% s,]
otu.t1=otu.t1[otu.t1$season!="nov19",]

#Test
kruskal=data.frame(row.names = rownames(otu.r), otu.id=rownames(otu.r), p.value=rep(1, dim(otu.r)[1]))

for (id in rownames(otu.r)){
  sub1=otu.t1[otu.t1$otu.id==id,]
  kruskal[id,"p.value"]=kruskal.test(relab ~ season, data = sub1)$p.value
} #Kruskal test

if (sum(kruskal$p.value<=0.05, na.rm = T)>0){
  
  dunn=data.frame(row.names = rownames(otu.r))
  
  for (h in rownames(otu.r)){
    sub1=otu.t1[otu.t1$otu.id==h,]
    
    if (sum(sub1$relab)>0){
      test1=dunnTest(relab ~ season, data = sub1, method = "bh")
      p.value=test1$res$P.adj
      dunn[h,1:choose(length(unique(sub1$season)),2)]=p.value
    } else if (sum(sub1$relab)==0){
      p.value=rep(1,choose(length(unique(sub1$season)),2))
      dunn[h,1:choose(length(unique(sub1$season)),2)]=p.value
    }
    
  } #Dunn multiple comparison 
}

colnames(dunn)=test1$res$Comparison
dunn[is.na(dunn)]=1
pv=dunn
  
##Plot the enrichment

#Calculate log fold
otu.d=summarySE(otu.t1, measurevar = "relab", groupvars = c("otu.id","season")) #mean abundance of each OTU
en=data.frame(row.names=unique(otu.d$otu.id))

for (i in colnames(pv)){
  j=str_split(i, " - ")[[1]]
  logfold=log2(otu.d[otu.d$season==j[1],"relab"]/otu.d[otu.d$season==j[2],"relab"])
  en=cbind(en, logfold)
}

colnames(en)=colnames(pv)
en=en[rownames(pv), colnames(pv)]
en[is.infinite(as.matrix(en))]=0
en[is.nan(as.matrix(en))]=0

#Generate plot data
plot.data=spread(otu.d[,c(1,2,4)], key = c("season"), value = "relab") #transform data from tidy to untidy in order to have the 3 axis for the triplot
plot.data=cbind(plot.data, tax[plot.data$otu.id,])
plot.data=plot.data[,-1]

#Select the differential OTUs in each season p<= 0.05 LF>=1.5 

df=unique(c(rownames(pv[pv$`feb20 - jul20`<=0.05 & en$`feb20 - jul20`>= 1.5,]),
            rownames(pv[pv$`feb20 - nov20`<=0.05 & en$`feb20 - nov20`>= 1.5,])))

dj=unique(c(rownames(pv[pv$`feb20 - jul20`<=0.05 & en$`feb20 - jul20`<=-1.5,]),
            rownames(pv[pv$`jul20 - nov20`<=0.05 & en$`jul20 - nov20`>= 1.5,])))

dn=unique(c(rownames(pv[pv$`feb20 - nov20`<=0.05 & en$`feb20 - nov20`<=-1.5,]),
            rownames(pv[pv$`jul20 - nov20`<=0.05 & en$`jul20 - nov20`<=-1.5,])))

sig=unique(c(df,dj,dn))
print(paste("differential OTUs", length(sig)))

#Generate labels for OTUS
plot.data$enrichment="not_enriched"
plot.data[plot.data$otu.id %in% sig,"enrichment"] ="enriched"

plot.data$label=plot.data[,t]
plot.data[is.na(plot.data$label), "label"]="Unclassified" 
plot.data[!plot.data$label %in% niveles, "label"]="Low.abundant"
plot.data[plot.data$enrichment=="not_enriched", "label"]="no.differential" 

#Generate colors per each taxa in bar plots 
niveles=niveles[niveles!="Unclassified"]
colores=c("#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848")
col=data.frame(color=c(colores[1:length(niveles)],"grey50","black"),taxa=c(niveles,"Unclassified","Low.abundant" ))

#Subset only the differential
over=plot.data[plot.data$label!="no.differential",]
over$label=factor(over$label, levels = c(niveles,"Unclassified","Low.abundant" ))
over2=over[!over$label %in% c("Unclassified","Low.abundant"),]

#Plot 
tema$legend.position="none"
plot=ggtern(data=over, aes(x=feb20, y=jul20, z=nov20, color=label))+
  geom_mask()+
  geom_point(alpha=1, shape=19,size=1)+
  geom_point(data=over2, aes(x=feb20, y=jul20, z=nov20, color=label),alpha=1, shape=19,size=1)+
  tema+
  scale_color_manual(values=c(col[col$taxa %in% unique(over$label),"color"]))+
  #scale_shape_manual(values=c(15:18,8))+
  #scale_size_continuous(range = c(1,4))+
  xlab("t3")+ylab("t8")+zlab("t12")

result="C:/Users/victorfn/Desktop/Results_syncoms/Enrichment"
dir.create(result)
result=paste("C:/Users/victorfn/Desktop/Results_syncoms/Enrichment","season",sep = "/")
dir.create(result)

png(paste(result,paste(amplicon,t,"triangle.png",sep = "."),sep = "/"), 
    width = 700, height = 700, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()


#Prepare data for saving
pv$otu.id=rownames(pv)
write.table(pv,paste(result,paste(amplicon,"p.value.differential.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)

en$otu.id=rownames(en)
write.table(en,paste(result,paste(amplicon,"log.fold.differential.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)


#Generates a data frame with the Differential OTUs in each season (Supplementary files)
for (i in colnames(pv)[1:3]){
  j=str_split(i, " - ")[[1]]
  
  print(j)
  
  p=rownames(pv[pv[,i]<=0.05,])
  e1=rownames(en[en[,i]>=1.5,])
  e2=rownames(en[en[,i]<=-1.5,])
  
  print(j[1]); print(sum(p%in%e1))
  print(j[2]); print(sum(p%in%e2))
  
  pv[p[p%in%e1],i]=j[1]
  pv[p[p%in%e2],i]=j[2]
  
}

res=data.frame(pv,tax[rownames(pv),])

write.table(res,paste(result,paste(amplicon,"enriched.differential.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)

#Check the number of differential OTUs

pv=dunn
#Differential OTUs in each season
for (i in colnames(pv)[1:3]){
  j=str_split(i, " - ")[[1]]
  
  print(j)
  
  p=rownames(pv[pv[,i]<=0.05,])
  e1=rownames(en[en[,i]>=1.5,])
  e2=rownames(en[en[,i]<=-1.5,])
  
  print(paste(j[1],sum(p%in%e1)))
  print(sort(summary(factor(tax[rownames(pv[p[p%in%e1],]),"order"]), maxsum = length(unique(tax$order))),decreasing = T))
  
  print(paste(j[2],sum(p%in%e2)))
  print(sort(summary(factor(tax[rownames(pv[p[p%in%e2],]),"order"]), maxsum = length(unique(tax$order))),decreasing = T))
  
}












