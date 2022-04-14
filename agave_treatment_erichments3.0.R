### Performs OTU enrichment analysis between treatments and plots ###

#Directories
setwd("C:/Users/victorfn/Desktop")
s16=("C:/Users/victorfn/Desktop/Results_syncoms/16s")
its=("C:/Users/victorfn/Desktop/Results_syncoms/ITS")   

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
library(ggrepel)

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
t=c("domain","phylum","class","order","family","genus","otu.id")[7] #taxa level analysis
s=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")[] #Analyze only one treatment or all
m=c("nov19","feb20","jul20","nov20")[]
test=c("kruskal")[1]

# Load metadata
met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"

#Load data 

if (amplicon=="16s"){
  #16S
  otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
  #match=read.delim(file=paste("C:/Users/victorfn/Desktop/Results_syncoms/16s",paste("isolate.match","txt",sep = "."),sep = "/"), sep = " ")
  load(file = paste(s16, paste("order","names",sep = "."),sep = "/"))
  
} else if (amplicon=="ITS"){
  #ITS2
  otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
  load(file = paste(its, paste("order","names",sep = "."),sep = "/"))
}

# Removing unwanted samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

# Sub sample reads 
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) 
otu.r=cbind(otu.s,tax)
otu.t=gather(data=otu.r, key = "sample", value = "abs", colnames(otu.s))

for (i in met$AIM.sample.id){
  otu.t[otu.t$sample == i, c( "block","square","lane","season","plant.compartment","treatment" )] = met.2[i,c( "block","square","lane","season","plant.compartment","treatment" )]
} 

otu.t[is.na(otu.t[,t]),t]="Unclassified"

#Calculates relative abundance for each taxa (OTU) in each treatment
otu.t0 = otu.t %>% group_by(sample,treatment) %>% summarise(absab = sum(abs))
otu.t1 = otu.t %>% group_by(otu.t[,t],sample,treatment) %>% summarise(absab = sum(abs)) 
colnames(otu.t1)[1]=c("rank")

#Sanity verification 
sum(paste(otu.t1$sample,otu.t1$treatment,otu.t1$rank)==paste(otu.t0$sample,otu.t0$treatment,otu.t1$rank))==dim(otu.t1)[1]
otu.t1$relab=otu.t1$absab/otu.t0$absab*100

#Subset the data 
otu.t1=otu.t1[otu.t1$sample %in% met.2[met.2$season %in% m, "AIM.sample.id"],]
otu.t1=otu.t1[otu.t1$sample %in% met.2[met.2$season != "nov19", "AIM.sample.id"],]
taxon=sort(unique(otu.t1$rank))

#Test 
if(test=="kruskal"){

kruskal=data.frame(row.names = taxon, rank=taxon, p.value=rep(1, length(taxon)[1]))

for (id in taxon){
  
  sub1=otu.t1[otu.t1$rank==id,]
  kruskal[id,"p.value"]=kruskal.test(relab ~ treatment, data = sub1)$p.value
}

if (sum(kruskal$p.value<=0.05, na.rm = T)>0){
  
  dunn=data.frame(row.names = taxon)
  
  for (h in taxon){
    sub1=otu.t1[otu.t1$rank==h,]
    
    if (sum(sub1$relab)>0){
      test1=dunnTest(relab ~ treatment, data = sub1, method = "bh")
      p.value=test1$res$P.adj
      #p.value=test1$res$P.unadj
      dunn[h,1:choose(length(unique(sub1$treatment)),2)]=p.value
    } else if (sum(sub1$relab)==0){
      p.value=rep(1,choose(length(unique(sub1$treatment)),2))
      dunn[h,1:choose(length(unique(sub1$treatment)),2)]=p.value
    }
    
  }
}

colnames(dunn)=test1$res$Comparison
dunn[is.na(dunn)]=1
pv=dunn[,grep("MOCK",colnames(dunn))]
colnames(pv)=gsub("MOCK - ","",colnames(pv))
colnames(pv)=gsub(" - MOCK","",colnames(pv))

for (h in colnames(pv)){
  print(h)
  print(sum(pv[,h]<=0.05))
}


} 

#Calculate log fold for plotting 
otu.d=summarySE(otu.t1, measurevar = "relab", groupvars = c("rank","treatment")) #mean abundance of each OTU
en=data.frame(row.names=taxon)
for (i in colnames(pv)){
  logfold=log2(otu.d[otu.d$treatment==i,"relab"]/otu.d[otu.d$treatment=="MOCK","relab"])
  en=cbind(en, logfold)
}
colnames(en)=colnames(pv)
en=en[rownames(pv), colnames(pv)]
en[is.infinite(as.matrix(en))]=0
en[is.nan(as.matrix(en))]=0

#Select any colors 
colores=c("#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848")


niveles=niveles[niveles!="Unclassified"]

# Generates a plot per each comparison 
for (w in colnames(pv) ){

  # threshold for enrichment 
thf=1.5
thp=0.05
lab="genus"

ab=otu.d[otu.d$treatment == w,]
plot.data=data.frame(row.names = rownames(pv), pvalue=pv[,w], logf=en[,w], relab=ab$relab)


for (i in rownames(plot.data)){
  plot.data[i, c("domain","phylum","class","order","family","genus","otu.id")]=unique(tax[tax[,t] %in% i, c("domain","phylum","class","order","family","genus","otu.id")])
}

plot.data$order[is.na(plot.data$order)]="Unclassified"
plot.data[!plot.data$order %in% c(niveles,"Unclassified"),"order"]="Low.abundant"
plot.data$order=factor(plot.data$order, levels = c(niveles,"Unclassified","Low.abundant"))
col=data.frame(color=c(colores[1:length(niveles)],"grey50","black"),taxa=c(niveles,"Unclassified","Low.abundant" ))


plot.data$label=plot.data[,lab]
plot.data[plot.data$pvalue>thp | abs(plot.data$logf)<thf,"label"]=""

plot=ggplot(data=plot.data, aes(x=logf, y=-log10(pvalue), color=order))+
  geom_point(size=1.5)+
  scale_color_manual(values=c(col[col$taxa %in% unique(plot.data$order),"color"]))+
  geom_text_repel(aes(label=label),size=1.5, color="black", max.overlaps = 40, force_pull = 2, force = 2)+
  geom_hline(yintercept = -log10(thp), linetype=2)+
  geom_vline(xintercept = c(-thf,thf), linetype=2)+tema+#+ggtitle(paste("MOCK","vs",w))
  scale_x_continuous(limits = c(-7.5,7.5), breaks = c(seq(-7.5,7.5,2.5)))


result=paste("C:/Users/victorfn/Desktop/Results_syncoms/Enrichment","treatment",sep = "/")
dir.create(result)

res=plot.data[plot.data$pvalue<=thp & abs(plot.data$logf)>=thf,]
res[res$logf>=thf,"label"]=w
res[res$logf<=-thf,"label"]="MOCK"


png(paste(result,paste(amplicon,w,"volcano.plot.png",sep = "."),sep = "/"), 
    width = 400, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot)
dev.off()

write.table(res,paste(result,paste(amplicon,w,"differential.plot.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)

}

