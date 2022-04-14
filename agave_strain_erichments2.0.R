# Directories
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
library(extrafont)
library(ggrepel)

#Graphics
colores=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
nivel=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
etiquetas=c("Mock","Functional","Functional-F","Functional-As","Functional-At","Photo-anox","Photo-oxi","Biofilm","RC1","RC2","RC3")

tema=theme(axis.text= element_text(color="gray30",size=7, angle=0,hjust=0.5,vjust=0,family = "Arial"),
           axis.title = element_text(color="black",size=7,family = "Arial"),
           legend.text = element_text(color = "black",size=7,family = "Arial"),
           legend.key.size = unit(0.7,"cm"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=7, color="black",face="bold",family = "Arial"),
           strip.text.y = element_text(size=7, color="black",face="bold",family = "Arial"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")

#Set parameters
amplicon=c("16s","ITS")[1]
bootstrap=0.80 #not implemented yet
filter=1000 #minimum read count for a valid sample
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
base=c("rdp")
t=c("domain","phylum","class","order","family","genus","otu.id")[4] #colors by
s=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")[] #Analyze only one treatment or all

#Load metadata
met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"

#Load data 
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

#Removing unwanted samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

#Sub sample reads
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2))))  
otu.r=t(t(otu.s)/colSums(otu.s)*100) #Relative abundance 
otu.r=cbind(otu.r,tax)
otu.t=gather(data=otu.r, key = "sample", value = "relab", colnames(otu.s))

for (i in met$AIM.sample.id){
  otu.t[otu.t$sample == i, c( "block","square","lane","season","plant.compartment","treatment" )] = met.2[i,c( "block","square","lane","season","plant.compartment","treatment" )]
} 


#Check matching OTUs
match=read.delim(file=paste("C:/Users/victorfn/Desktop/Results_syncoms",paste("isolate.match.16s.consensus","txt",sep = "."),sep = "/"), sep = "\t")
match$isolate.match=gsub("_[a-z]+_","_",match$isolate.match)
print("did not pass the thresholds") ; print(match[is.na(match$reads.my.data),"isolate.match"])
print("Double hit") ; print(match[duplicated(match$consensus.match),"otu.id"]) ; print(match[duplicated(match$consensus.match),"isolate.match"])
match=match[!is.na(match$reads.my.data),]
match=match[!duplicated(match$consensus.match),]
rownames(match)=match$consensus.match

#Relative abundance of matching OTUs
otu.r=otu.r[match$consensus.match,]
otu.t=otu.t[otu.t$otu.id %in% match$consensus.match,]
print("pct de reads") ; sum(rowSums(otu[match$consensus.match,])) / sum(rowSums(otu.2)) * 100

#Enrichment between TREATMENTS by season

se=c("nov19","feb20","jul20","nov20")[1] #Choose 1 or more seasons

#Group treatments 
otu.n=otu.t[otu.t$season %in% se,]
otu.n[otu.n$treatment %in% nivel[c(2,3,4,5)],"treatment"]="Network"
otu.n[otu.n$treatment %in% nivel[c(9:11)],"treatment"]="Random"
otu.n[otu.n$treatment %in% nivel[c(6:8)],"treatment"]="Metagenome"
#otu.n=otu.n[!otu.n$treatment %in% nivel[c(7,8)],]

#Test 
kruskal=data.frame(row.names = rownames(otu.r), otu.id=rownames(otu.r), p.value=rep(1, dim(otu.r)[1]))
for (id in rownames(otu.r)){
  sub1=otu.n[otu.n$otu.id==id,]
  kruskal[id,"p.value"]=kruskal.test(relab ~ treatment, data = sub1)$p.value
}

if (sum(kruskal$p.value<=0.05)>0){
  dunn=data.frame(row.names = rownames(otu.r))
  for (h in rownames(otu.r)){
    sub1=otu.n[otu.n$otu.id==h,]
    test1=dunnTest(relab ~ treatment, data = sub1, method = "bh")
    p.value=test1$res$P.adj
    dunn[h,1:choose(length(unique(otu.n$treatment)),2)]=p.value
    
  }
} else {
  
  print("No differences per treatment")
}

colnames(dunn)=test1$res$Comparison
pv=cbind(dunn,strain=match[rownames(dunn),"isolate.match"])

#Add strain name
for (i in match$consensus.match){
  otu.n[otu.n$otu.id==i,"isolate.match"]=match[match$consensus.match==i,"isolate.match"]}

#Generate data frame of mean abundances 
mean.n=summarySE(otu.n, measurevar = "relab", groupvars = c("isolate.match","treatment","season"))
mean.n$treatment=factor(mean.n$treatment, levels = c("MOCK","Network","Random","Metagenome"))
mean.n$season=factor(mean.n$season, levels = se)

#Calculate log fold
en=data.frame(row.names=match$isolate.match)
for (i in colnames(pv)[1:6]){
  j=str_split(i, " - ")[[1]]
  logfold=log2(mean.n[mean.n$treatment==j[1],"relab"]/mean.n[mean.n$treatment==j[2],"relab"])
  en=cbind(en, logfold)
}

colnames(en)=colnames(pv)[1:6]
#pv$strain=rownames(pv)
en$strain=rownames(en)

#Save results
result.dir=paste("C:/Users/victorfn/Desktop/Results_syncoms/Enrichment","strains",sep = "/")
dir.create(result.dir)
write.table(pv,paste(result.dir,paste(se,amplicon,"p_value_by_treatment.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)

write.table(en,paste(result.dir,paste(se,amplicon,"logfold_by_treatment.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)


#Make an abundance plot

#Group treatments 
se=c("nov19","feb20","jul20","nov20")[]
otu.n=otu.t[otu.t$season %in% se ,]
otu.n[otu.n$treatment %in% nivel[c(2,3,4,5)],"treatment"]="Network"
otu.n[otu.n$treatment %in% nivel[c(9:11)],"treatment"]="Random"
otu.n[otu.n$treatment %in% nivel[c(6:8)],"treatment"]="Metagenome"
#otu.n=otu.n[!otu.n$treatment %in% nivel[c(7,8)],]

#Generate dataframe 
for (i in match$consensus.match){
  otu.n[otu.n$otu.id==i,"isolate.match"]=match[match$consensus.match==i,"isolate.match"]}
mean.n=summarySE(otu.n, measurevar = "relab", groupvars = c("isolate.match","treatment","season"))
mean.n$treatment=factor(mean.n$treatment, levels = c("MOCK","Network","Metagenome","Random"))
mean.n$season=factor(mean.n$season, levels = se)

#Plot
tema$legend.position="none"
plot=ggplot(mean.n,(aes(x=season,y=relab, fill=treatment)))+
  geom_bar(stat="identity", position=position_dodge(),colour="black",size=0.4)+tema+
  geom_errorbar(aes(ymin=relab, ymax=relab+sd), width=0.5,size=0.4,
                position=position_dodge(0.9),stat="identity")+
  scale_fill_manual(values = colores[c(1,3,6,10)])+ylab("relative abundance (%)")+
  facet_wrap(~isolate.match, scales = "free", nrow = 6, ncol = 4)+#+scale_y_continuous(limits = c(0,37))+
  geom_hline(yintercept = 1, linetype=2,size=0.4)
  #geom_hline(yintercept = 0.1, linetype=2)

png(paste(result.dir,paste(amplicon,"strainsplot.png",sep = "."),sep = "/"), 
    width = 1700, height = 2000, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot)
dev.off()

write.table(mean.n,paste(result.dir,paste(amplicon,"strain.abundances.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)


# Enrichment per SEASON by treatment 

nt=c("MOCK","Network","Metagenome","Random")

#Group treatments 
otu.n=otu.t
otu.n[otu.n$treatment %in% nivel[c(2,3,4,5)],"treatment"]="Network"
otu.n[otu.n$treatment %in% nivel[c(9:11)],"treatment"]="Random"
otu.n[otu.n$treatment %in% nivel[c(6:8)],"treatment"]="Metagenome"


#Test enrichment per season in each group of treatments
pv.treat=list()
en.treat=list()
for (s in nt){ 
  otu.t1=otu.n[otu.n$treatment==s,]
  #test per compartment
  kruskal=data.frame(row.names = rownames(otu.r), otu.id=rownames(otu.r), p.value=rep(1, dim(otu.r)[1]))
  for (id in rownames(otu.r)){
    sub1=otu.t1[otu.t1$otu.id==id,]
    kruskal[id,"p.value"]=kruskal.test(relab ~ season, data = sub1)$p.value
  }
  
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
      
    }
  }
  
  colnames(dunn)=test1$res$Comparison
  dunn[is.na(dunn)]=1
  pv.treat[[s]]=dunn
  pv=dunn
  
  #Generata dat frame of mean relatiev abundances 
  otu.d=summarySE(otu.n, measurevar = "relab", groupvars = c("otu.id","season"))
  en=data.frame(row.names=unique(otu.d$otu.id))
  comp=colnames(pv)
  
  #Calculate logFold change 
  for (i in comp){
    
    j=str_split(i, " - ")[[1]]
    
    logfold=log2(otu.d[otu.d$season==j[1],"relab"]/otu.d[otu.d$season==j[2],"relab"])
    en=cbind(en, logfold)
  }
  colnames(en)=comp
  en=en[rownames(pv), colnames(pv)]
  
  en.treat[[s]]=en
}


#How many differential are there?
comp=data.frame()
for (se in nt){
  print(se)
  pv=pv.treat[[se]]
  r=colSums(pv<=0.05)
  comp=rbind(comp,r)
}
colnames(comp)=colnames(pv)
comp$treatment=nt


####


#Generate a dot plot of seasonal enrichments (Not used in manuscript) 
result=data.frame()
for (t in nt){
  
  pv=pv.treat[[t]] ; pv$otu.id=rownames(pv)
  en=en.treat[[t]] ; en$otu.id=rownames(en)
  
  pv=pv[,c("feb20 - nov19","feb20 - jul20","jul20 - nov20","otu.id")]
  en=en[,c("feb20 - nov19","feb20 - jul20","jul20 - nov20","otu.id")]
  
  pv.t=gather(pv, key = "comparison", value = "p.value", c("feb20 - nov19","feb20 - jul20","jul20 - nov20"))
  en.t=gather(en, key = "comparison", value = "logfold", c("feb20 - nov19","feb20 - jul20","jul20 - nov20"))
  
  pv.t$logfold=en.t$logfold
  pv.t$treatment=rep(t,dim(pv.t)[1])
  result=rbind(result,pv.t)
  
}

for (i in match$consensus.match){
  result[result$otu.id==i,"isolate.match"]=match[match$consensus.match==i,"isolate.match"]}


result$label=result$isolate.match
result[result$p.value>0.05,"label"]="not.significant"
result$label=factor(result$label,levels = c(match$isolate.match,"not.significant"))
result$treatment=factor(result$treatment, levels = nt)
result$comparison=factor(result$comparison, levels = c("feb20 - nov19","feb20 - jul20","jul20 - nov20"))


write.table(result,paste(result.dir,paste(amplicon,"strains_p_value_logfold_byseason.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)



colores=c("#A12F2F","#FF4848","#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","grey75","grey75")

result[result$comparison=="jul20 - nov20","logfold"]=result[result$comparison=="jul20 - nov20","logfold"]*-1
result$label2=result$isolate.match
result[!result$label2 %in% c(match$isolate.match)[],"label2"]="" #c(3,10:12,20,18,19,7)
result[result$p.value>0.05,"label2"]=""
result$label2=substr(result$label2,1,3)

tema$legend.position="right"
plot=ggplot(result,aes(x=comparison, y=logfold, size=-log10(p.value), colour=label))+
  geom_jitter(width = 0.35, alpha=1)+facet_wrap(~treatment,ncol = 4)+
  scale_size_continuous(range=c(2,6))+tema+
  geom_hline(yintercept = c(1,-1,0), linetype=3)+
  scale_color_manual(values = c(colores[1:16],"grey50"))+
  geom_text_repel(aes(label=label2), color="black", size=3, force=0.5, max.overlaps = 10)


png(paste(result.dir,paste(amplicon,"strainsplot_byseason.png",sep = "."),sep = "/"), 
    width = 2000, height = 1500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot)
dev.off()


write.table(result,paste(result.dir,paste(amplicon,"strains_p_value_logfold_byseason.txt",sep = "."),sep = "/"),
            sep = "\t", quote = F, row.names = F, col.names = T)



