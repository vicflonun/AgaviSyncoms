### Generates plots and statistical analyses for alpha and beta diversity ###

#Packages
library(vegan) 
library(tidyr) 
library(dplyr) 
library(ggplot2)
library(MASS)
library(Rmisc)
library(extrafont)
library(Rmisc)
library(FSA)
library(colorspace)
library(pheatmap)


#Directories
setwd("C:/Users/victorfn/Desktop")  
result="C:/Users/victorfn/Desktop/Results_syncoms/Diversity" #result directories
dir.create(result)

s16=("C:/Users/victorfn/Desktop/Results_syncoms/16s") #data directories
its=("C:/Users/victorfn/Desktop/Results_syncoms/ITS") 


#Graphics
colores.t=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
niveles=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
etiquetas=c("Mock","Functional","Functional-F","Functional-As","Functional-At","Photo-anox","Photo-oxi","Biofilm","RC1","RC2","RC3")
#colores.s=c("#BE9E80","#886B4F","#447690","#594635")
colores.s=c("#FFC5C7","#E6959F","#6091AC","#CB627B")
scaleFUN <- function(x) sprintf("%.2f", x) #Put 2 extra decimal zeros

tema=theme(axis.text.x = element_text(color="black",size=7, angle=0,hjust=0.5,vjust=0,family = "Arial"),
           axis.text.y = element_text(color="black",size=7,family = "Arial"),
           axis.title = element_text(color="black",size=7,family = "Arial"),
           legend.text = element_text(color = "black",size=7,family = "Arial"),
           legend.key.size = unit(0.3,"cm"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=7, color="black",face="bold",family = "Arial"),
           strip.text.y = element_text(size=7, color="black",face="bold",family = "Arial"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")


#### FIGURE 2 & supplementary

# Set PARAMETERS
amplicon=c("16s","ITS")[2]
bootstrap=0.80 #not implemented yet
filter=1000 #minimum read count for a valid sample
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
base=c("rdp")

# Load metadata
met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"

# Load data
if (amplicon=="16s"){
  #16s
  otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
} else if (amplicon=="ITS"){
  #ITS
  otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
}


# Subset samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",] #Remove samples that are not part of the analysis
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

#Subsample reads
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) #Subsampling 
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#Calculate alpha-diversity
ric=rarefy(t(otu.2), sample = min(colSums(otu.2)))
shan=diversity(t(otu.s),index = "shannon")

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan )
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c( "block","square","lane","season","treatment","AIM.sample.id","plate")])
alpha$season=factor(alpha$season, level = c("nov19","feb20","jul20","nov20"))
alpha$treatment=factor(alpha$treatment, level = niveles)

#block C will not be plotted
met.2=met.2[met.2$block!="C",]
otu.r=otu.r[,rownames(met.2)]
alpha=alpha[alpha$AIM.sample.id %in% rownames(met.2),]

# Calculate beta-diversity
#otu.n=log10(otu.n+1)
#otu.n=sqrt(otu.n)
set.seed(23171341)
scaling=vegdist(t(otu.r), method = "bray", binary = T) #calculate distance
scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
#scaling2=metaMDS(otu.n, distance = "bray")
#scaling2=monoMDS(scaling)
scaling3=data.frame(scaling2$points) #select cordenates
scaling3=cbind(scaling3,alpha)


### plots by SEASON

# alpha
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "season")
plot=ggplot(data=mean.r, aes(x=season, y=richness,fill=season))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  #scale_fill_manual(values =colores[c(2,4,5,8)])+tema+
  scale_fill_manual(values =colores.s)+tema+
  scale_y_continuous(limits = c(0,500),breaks = c(seq(0,500,100)))+ylab("OTU richness")+
  scale_x_discrete(labels=c("t0","t3","t8","t12"))+xlab("months post inoculation")

png(paste(result,paste(amplicon,"season.richness.png",sep = "."),sep = "/"), 
    width = 320, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()

mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "season")
plot=ggplot(data=mean.s, aes(x=season, y=round(shannon,1), fill=season))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #scale_fill_manual(values = colores[c(2,4,5,8)])+tema+
  scale_fill_manual(values =colores.s)+tema+
  scale_y_continuous(limits = c(0,8),breaks = c(seq(0,8,2)), labels = scaleFUN)+ylab("Shannon index")+
  scale_x_discrete(labels=c("t0","t3","t8","t12"))+xlab("months post inoculation")

png(paste(result,paste(amplicon,"season.shannon.png",sep = "."),sep = "/"), 
    width = 320, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()

#NMDS
plot=ggplot(data=scaling3, aes(x=X1, y=X2, fill=season, shape=block))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(alpha=1,size=2)+tema+
  scale_fill_manual(values=colores.s)+ylab("NMDS2")+xlab("NMDS1")+
  scale_shape_manual(values=c(21:26))

png(paste(result,paste(amplicon,"season.NMDS.png",sep = "."),sep = "/"), 
    width = 430, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()


#beta
mod = with(met.2, betadisper(scaling, season))
alpha$centroid=mod$distances
mean.c=summarySE(data = alpha, measurevar = "centroid", groupvars = "season")
plot=ggplot(data=mean.c, aes(x=season, y=centroid, fill=season))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=centroid-sd, ymax=centroid+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=1,dotsize=0.5,aes(fill=season))+
  scale_fill_manual(values =colores.s)+tema+
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,1,0.2),labels = scaleFUN)+ylab("Distance to centroid")+
  scale_x_discrete(labels=c("t0","t3","t8","t12"))+xlab("months post inoculation")

png(paste(result,paste(amplicon,"season.distance.png",sep = "."),sep = "/"), 
    width = 310, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()


#test
kruskal.test(richness~season, data = alpha) ; r=dunnTest(richness~season, data = alpha, method = "bh")$res
write.table(r,file = paste(result,paste(amplicon,"season.richness.Kruskal-Dunn.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)
kruskal.test(shannon~season, data = alpha) ; r=dunnTest(shannon~season, data = alpha, method = "bh")$res
write.table(r,file = paste(result,paste(amplicon,"season.shannon.Kruskal-Dunn.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)
kruskal.test(centroid~season, data = alpha) ; r=dunnTest(centroid~season, data = alpha, method = "bh")$res
write.table(r,file = paste(result,paste(amplicon,"season.distance.Kruskal-Dunn.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)


### plot by BLOCK

#alpha
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "block")
plot=ggplot(data=mean.r, aes(x=block, y=richness,fill=block))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values =c("grey75","grey25"))+tema+
  scale_y_continuous(limits = c(0,500),breaks = c(seq(0,500,100)))+
  ylab("OTU richness")+
  xlab("block")

png(paste(result,paste(amplicon,"block.richness.png",sep = "."),sep = "/"), 
    width = 320, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot); dev.off()

mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "block")
plot=ggplot(data=mean.s, aes(x=block, y=shannon, fill=block))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values =c("grey75","grey25"))+tema+
  scale_y_continuous(limits = c(0,8),breaks = c(seq(0,8,2)),labels = scaleFUN)+
  ylab("Shannon index")+xlab("block")

png(paste(result,paste(amplicon,"block.shannon.png",sep = "."),sep = "/"), 
    width = 320, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot); dev.off()

#NMDS
plot=ggplot(data=scaling3, aes(x=X1, y=X2, fill=block, shape=block))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(alpha=1,size=2)+tema+ylab("NMDS2")+xlab("NMDS1")+
  scale_fill_manual(values = c("grey75","grey25"))+
  scale_shape_manual(values=c(21:26))

png(paste(result,paste(amplicon,"block.NMDS.png",sep = "."),sep = "/"), 
    width = 430, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot)
dev.off()
 
#beta
mod = with(met.2, betadisper(scaling, block))
alpha$centroid=mod$distances
mean.c=summarySE(data = alpha, measurevar = "centroid", groupvars = "block")
plot=ggplot(data=mean.c, aes(x=block, y=centroid, fill=block))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=centroid-sd, ymax=centroid+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=1.5, aes(fill=season))+
  scale_fill_manual(values = c("grey75","grey25"))+tema+
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2),labels = scaleFUN)+ylab("Distance to centroid")+
  xlab("block")

png(paste(result,paste(amplicon,"block.distance.png",sep = "."),sep = "/"), 
    width = 320, height = 400, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot); dev.off()


#test
r=print(wilcox.test(alpha[alpha$block=="A","richness"],alpha[alpha$block=="B","richness"]))
s=print(wilcox.test(alpha[alpha$block=="A","shannon"],alpha[alpha$block=="B","shannon"]))
c=print(wilcox.test(alpha[alpha$block=="A","centroid"],alpha[alpha$block=="B","centroid"]))

r=data.frame(index=c("richness","shannon","distance.to.centroid"), W=c(r$statistic,s$statistic,c$statistic),p=c(r$p.value,s$p.value,c$p.value))
write.table(r,file = paste(result,paste(amplicon,"block.wilcox.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)


# plot by TREATMENT

#alpha
tema$axis.text.x$angle=90
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "treatment")
plot=ggplot(data=mean.r, aes(x=treatment, y=richness,fill=treatment))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values =colores.t)+tema+
  scale_y_continuous(limits = c(0,500),breaks = c(seq(0,500,100)))+
  ylab("OTU richness")

png(paste(result,paste(amplicon,"treatment.richness.png",sep = "."),sep = "/"), 
    width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()

mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "treatment")
plot=ggplot(data=mean.s, aes(x=treatment, y=shannon, fill=treatment))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values = colores.t)+tema+
  scale_y_continuous(limits = c(0,8),breaks = c(seq(0,8,2)),labels = scaleFUN)+
  ylab("Shannon index")

png(paste(result,paste(amplicon,"treatment.shannon.png",sep = "."),sep = "/"), 
    width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()

#NMDS
tema$axis.text.x$angle=0
plot=ggplot(data=scaling3, aes(x=X1, y=X2, fill=treatment, shape=block))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(alpha=1,size=2)+tema+ylab("NMDS2")+xlab("NMDS1")+
  scale_fill_manual(values =colores.t)+
  scale_shape_manual(values=c(21:25))+
  scale_size_continuous(range = c(1,3))

png(paste(result,paste(amplicon,"treatment.NMDS.png",sep = "."),sep = "/"), 
    width = 510, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()

#beta
tema$axis.text.x$angle=90
mod = with(met.2, betadisper(scaling, treatment))
alpha$centroid=mod$distances
mean.c=summarySE(data = alpha, measurevar = "centroid", groupvars = "treatment")
plot=ggplot(data=mean.c, aes(x=treatment, y=centroid, fill=treatment))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=centroid-sd, ymax=centroid+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=1.5, aes(fill=season))+
  scale_fill_manual(values = colores.t)+tema+
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2),labels = scaleFUN)+
  ylab("Distance to centroid")


png(paste(result,paste(amplicon,"treatment.distance.png",sep = "."),sep = "/"), 
    width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot);dev.off()

#test
r=c()
for (i in niveles[-1]){
  r[i]=wilcox.test(alpha[alpha$treatment=="MOCK","richness"],alpha[alpha$treatment==i,"richness"])$p.value}

s=c()
for (i in niveles[-1]){
  s[i]=wilcox.test(alpha[alpha$treatment=="MOCK","shannon"],alpha[alpha$treatment==i,"shannon"])$p.value}

c=c()
for (i in niveles[-1]){
  c[i]=wilcox.test(alpha[alpha$treatment=="MOCK","centroid"],alpha[alpha$treatment==i,"centroid"])$p.value}


r=data.frame(treatment=niveles[-1],richness.p=r,shannon.p=s,distance.p=c)
write.table(r,file = paste(result,paste(amplicon,"treatment.wilcox-MOCK.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)


# plot TREATMENT by SEASON 

#Load data

  if (amplicon=="16s"){
    #16s
    otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
  } else if (amplicon=="ITS"){
    #ITS
    otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
  }
  
  
#Subset samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

#Subsample otus
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) #rarefication 
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#calculate alpha-diversity
ric=rarefy(t(otu.2), sample = min(colSums(otu.2)))
shan=diversity(t(otu.s),index = "shannon")

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan )
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c( "block","square","lane","season","treatment","AIM.sample.id","plate")])
alpha$season=factor(alpha$season, level = c("nov19","feb20","jul20","nov20"))
alpha$treatment=factor(alpha$treatment, level = niveles)

#block C will not be plotted
met.2=met.2[met.2$block!="C",]
otu.r=otu.r[,rownames(met.2)]
alpha=alpha[alpha$AIM.sample.id %in% rownames(met.2),]

#plot
mes=c("feb20","jul20","nov20")

for (se in mes){

  #subset data
  met.s=met.2[met.2$season==se,]
  alpha.s=alpha[alpha$season==se,]
  
  #Calculate beta-diversity 
  set.seed(23171341)
  scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T)
  scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
  scaling3=data.frame(scaling2$points) #select coordenates
  scaling3=cbind(scaling3,alpha.s)
  
  #alpha
  tema$axis.text.x$angle=90
  mean.r=summarySE(data = alpha.s, measurevar = "richness", groupvars = "treatment")
  plot=ggplot(data=mean.r, aes(x=treatment, y=richness,fill=treatment))+
    geom_bar(stat="identity",colour="black",size=0.5)+
    geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
    scale_fill_manual(values =colores.t)+tema+
    scale_y_continuous(limits = c(0,500),breaks = c(seq(0,500,100)))+ylab("OTU richness")
  
  png(paste(result,paste(amplicon,se,"treatment.richness.png",sep = "."),sep = "/"), 
      width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot); dev.off()
  
  mean.s=summarySE(data = alpha.s, measurevar = "shannon", groupvars = "treatment")
  plot=ggplot(data=mean.s, aes(x=treatment, y=shannon, fill=treatment))+
    geom_bar(stat="identity",colour="black",size=0.5)+
    geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
    scale_fill_manual(values = colores.t)+tema+
    scale_y_continuous(limits = c(0,8),breaks = c(seq(0,8,2)),labels = scaleFUN)+ylab("Shannon index")
  
  png(paste(result,paste(amplicon,se,"treatment.shannon.png",sep = "."),sep = "/"), 
      width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot); dev.off()
  
  #NMDS
  tema$axis.text.x$angle=0
  plot=ggplot(data=scaling3, aes(x=X1, y=X2, fill=treatment, shape=block))+
    geom_hline(yintercept = 0, linetype=2)+
    geom_vline(xintercept = 0,linetype=2)+
    geom_point(alpha=1,size=2)+tema+ylab("NMDS2")+xlab("NMDS1")+
    scale_fill_manual(values =colores.t)+
    scale_shape_manual(values=c(21:25))

  png(paste(result,paste(amplicon,se,"treatment.NMDS.png",sep = "."),sep = "/"), 
      width = 510, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot) ; dev.off()

  #beta
  mod = with(met.s, betadisper(scaling, treatment))
  alpha.s$centroid=mod$distances
  mean.c=summarySE(data = alpha.s, measurevar = "centroid", groupvars = "treatment")
  plot=ggplot(data=mean.c, aes(x=treatment, y=centroid, fill=treatment))+
    geom_bar(stat="identity",colour="black",size=0.5)+
    geom_errorbar(size=0.5,aes(ymin=centroid-sd, ymax=centroid+sd),width=0.2,position=position_dodge(1.5))+
    #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=1.5, aes(fill=season))+
    scale_fill_manual(values = colores.t)+tema+
    scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2))+
    ylab("Distance to centroid")+coord_flip()
  
  
  png(paste(result,paste(amplicon,se,"treatment.distance.png",sep = "."),sep = "/"), 
      width = 300, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot); dev.off()

  #test
  r=c()
  for (i in niveles[-1]){
    r[i]=wilcox.test(alpha.s[alpha.s$treatment=="MOCK","richness"],alpha.s[alpha.s$treatment==i,"richness"])$p.value}
  s=c()
  for (i in niveles[-1]){
    s[i]=wilcox.test(alpha.s[alpha.s$treatment=="MOCK","shannon"],alpha.s[alpha.s$treatment==i,"shannon"])$p.value}
  c=c()
  for (i in niveles[-1]){
    c[i]=wilcox.test(alpha.s[alpha.s$treatment=="MOCK","centroid"],alpha.s[alpha.s$treatment==i,"centroid"])$p.value}
  
  
  r=data.frame(treatment=niveles[-1],richness.p=r,shannon.p=s,distance.p=c)
  write.table(r,file = paste(result,paste(amplicon,se,"treatment.wilcox-MOCK.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)
}


  
# plot TREATMENT by  BLOCK 
  
  if (amplicon=="16s"){
    #16s
    otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
  } else if (amplicon=="ITS"){
    #ITS
    otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
  }
  
  
#Subset samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

#Subsample otus
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) #rarefication 
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#calculate alpha-diversity
ric=rarefy(t(otu.2), sample = min(colSums(otu.2)))
shan=diversity(t(otu.s),index = "shannon")

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan )
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c( "block","square","lane","season","treatment","AIM.sample.id","plate")])
alpha$season=factor(alpha$season, level = c("nov19","feb20","jul20","nov20"))
alpha$treatment=factor(alpha$treatment, level = niveles)

#block C will not be plotted
met.2=met.2[met.2$block!="C",]
otu.r=otu.r[,rownames(met.2)]
alpha=alpha[alpha$AIM.sample.id %in% rownames(met.2),]

block=c("A","B")
  
for (se in block){
    
    #Subset data
    met.s=met.2[met.2$block==se,]
    alpha.s=alpha[alpha$block==se,]
    
    #Calculate beta-diversity
    set.seed(23171341)
    scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T)
    scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
    scaling3=data.frame(scaling2$points) #select cordenates
    scaling3=cbind(scaling3,alpha.s)
    
    #alpha
    tema$axis.text.x$angle=90
    mean.r=summarySE(data = alpha.s, measurevar = "richness", groupvars = "treatment")
    plot=ggplot(data=mean.r, aes(x=treatment, y=richness,fill=treatment))+
      geom_bar(stat="identity",colour="black",size=0.5)+
      geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
      scale_fill_manual(values =colores.t)+tema+
      scale_y_continuous(limits = c(0,500),breaks = c(seq(0,500,100)))+ylab("OTU richness")
    
    png(paste(result,paste(amplicon,se,"treatment.richness.png",sep = "."),sep = "/"), 
        width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(plot); dev.off()
    
    mean.s=summarySE(data = alpha.s, measurevar = "shannon", groupvars = "treatment")
    plot=ggplot(data=mean.s, aes(x=treatment, y=shannon, fill=treatment))+
      geom_bar(stat="identity",colour="black",size=0.5)+
      geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
      scale_fill_manual(values = colores.t)+tema+
      scale_y_continuous(limits = c(0,8),breaks = c(seq(0,8,2)),labels = scaleFUN)+ylab("Shannon index")
    
    png(paste(result,paste(amplicon,se,"treatment.shannon.png",sep = "."),sep = "/"), 
        width = 500, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(plot); dev.off()
    
    #NMDS
    tema$axis.text.x$angle=0
    plot=ggplot(data=scaling3, aes(x=X1, y=X2, fill=treatment, shape=block))+
      geom_hline(yintercept = 0, linetype=2)+
      geom_vline(xintercept = 0,linetype=2)+
      geom_point(alpha=1,size=2)+tema+ylab("NMDS2")+xlab("NMDS1")+
      scale_fill_manual(values =colores.t)+
      scale_shape_manual(values=c(21:25))
    
    png(paste(result,paste(amplicon,se,"treatment.NMDS.png",sep = "."),sep = "/"), 
        width = 510, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(plot); dev.off()
    
    #beta
    mod = with(met.s, betadisper(scaling, treatment))
    alpha.s$centroid=mod$distances
    mean.c=summarySE(data = alpha.s, measurevar = "centroid", groupvars = "treatment")
    plot=ggplot(data=mean.c, aes(x=treatment, y=centroid, fill=treatment))+
      geom_bar(stat="identity",colour="black",size=0.5)+
      geom_errorbar(size=0.5,aes(ymin=centroid-sd, ymax=centroid+sd),width=0.2,position=position_dodge(1.5))+
      #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=1.5, aes(fill=season))+
      scale_fill_manual(values = colores.t)+tema+
      scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2))+
      ylab("Distance to centroid")+coord_flip()
    
    
    png(paste(result,paste(amplicon,se,"treatment.distance.png",sep = "."),sep = "/"), 
        width = 300, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
    print(plot); dev.off()
    
    #test
    r=c()
    for (i in niveles[-1]){
      r[i]=t.test(alpha.s[alpha.s$treatment=="MOCK","richness"],alpha.s[alpha.s$treatment==i,"richness"])$p.value}
    s=c()
    for (i in niveles[-1]){
      s[i]=t.test(alpha.s[alpha.s$treatment=="MOCK","shannon"],alpha.s[alpha.s$treatment==i,"shannon"])$p.value}
    c=c()
    for (i in niveles[-1]){
      c[i]=t.test(alpha.s[alpha.s$treatment=="MOCK","centroid"],alpha.s[alpha.s$treatment==i,"centroid"])$p.value}
    
    r=data.frame(treatment=niveles[-1],richness.p=r,shannon.p=s,distance.p=c)
    write.table(r,file = paste(result,paste(amplicon,se,"treatment.wilcox-MOCK.txt",sep = "."), sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)
  }
  
  
  
# plot BLOCK by SEASON


#plot
mes=c("feb20","jul20","nov20")

r=data.frame(row.names = mes)
for (se in mes){
  
  #subset data
  met.s=met.2[met.2$season==se,]
  alpha.s=alpha[alpha$season==se,]
  
  #Calculate beta-diversity 
  set.seed(23171341)
  scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T)
  scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
  scaling3=data.frame(scaling2$points) #select coordenates
  scaling3=cbind(scaling3,alpha.s)
  
  #alpha
  tema$axis.text.x$angle=0
  mean.r=summarySE(data = alpha.s, measurevar = "richness", groupvars = "block")
  plot=ggplot(data=mean.r, aes(x=block, y=richness,fill=block))+
    geom_bar(stat="identity",colour="black",size=0.5)+
    geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
    scale_fill_manual(values =c("grey75","grey25"))+tema+
    scale_y_continuous(limits = c(0,500),breaks = c(seq(0,500,100)))+ylab("OTU richness")
  
  png(paste(result,paste(amplicon,se,"block.richness.png",sep = "."),sep = "/"), 
      width = 220, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot); dev.off()
  
  mean.s=summarySE(data = alpha.s, measurevar = "shannon", groupvars = "block")
  plot=ggplot(data=mean.s, aes(x=block, y=shannon, fill=block))+
    geom_bar(stat="identity",colour="black",size=0.5)+
    geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
    scale_fill_manual(values = c("grey75","grey25"))+tema+
    scale_y_continuous(limits = c(0,8),breaks = c(seq(0,8,2)),labels = scaleFUN)+ylab("Shannon index")
  
  png(paste(result,paste(amplicon,se,"block.shannon.png",sep = "."),sep = "/"), 
      width = 220, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot); dev.off()
  
  #NMDS
  plot=ggplot(data=scaling3, aes(x=X1, y=X2, fill=block, shape=block))+
    geom_hline(yintercept = 0, linetype=2)+
    geom_vline(xintercept = 0,linetype=2)+
    geom_point(alpha=1,size=2)+tema+ylab("NMDS2")+xlab("NMDS1")+
    scale_fill_manual(values =c("grey75","grey25"))+
    scale_shape_manual(values=c(21:25))
  
  png(paste(result,paste(amplicon,se,"block.NMDS.png",sep = "."),sep = "/"), 
      width = 330, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot) ; dev.off()
  
  #beta
  mod = with(met.s, betadisper(scaling, block))
  alpha.s$centroid=mod$distances
  mean.c=summarySE(data = alpha.s, measurevar = "centroid", groupvars = "block")
  plot=ggplot(data=mean.c, aes(x=block, y=centroid, fill=block))+
    geom_bar(stat="identity",colour="black",size=0.5)+
    geom_errorbar(size=0.5,aes(ymin=centroid-sd, ymax=centroid+sd),width=0.2,position=position_dodge(1.5))+
    #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=1.5, aes(fill=season))+
    scale_fill_manual(values = c("grey75","grey25"))+tema+
    scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2))+
    ylab("Distance to centroid")
  
  
  png(paste(result,paste(amplicon,se,"block.distance.png",sep = "."),sep = "/"), 
      width = 220, height = 300, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  print(plot); dev.off()
  
  #test

  
    r[se,"richness"]=wilcox.test(alpha.s[alpha.s$block=="A","richness"],alpha.s[alpha.s$block=="B","richness"])$p.value
    r[se,"shannon"]=wilcox.test(alpha.s[alpha.s$block=="A","shannon"],alpha.s[alpha.s$block=="B","shannon"])$p.value
    r[se,"distance"]=wilcox.test(alpha.s[alpha.s$block=="A","centroid"],alpha.s[alpha.s$block=="B","centroid"])$p.value
  
}

r$season=mes
write.table(r,file = paste(result,paste(amplicon,"season","block.wilcox-MOCK.txt",sep = "."),sep = "/"),sep = "\t",quote = F, row.names = F, col.names = T)



#### FIGURE 3 & supplementary

# Set PARAMETERS
amplicon=c("16s","ITS")[1]

# Load data
if (amplicon=="16s"){
  #16s
  otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
} else if (amplicon=="ITS"){
  #ITS
  otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
}

# Subset samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",] #Remove samples that are not part of the analysis
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

#Subsample reads
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) #Subsampling 
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#Calculate alpha-diversity
ric=rarefy(t(otu.2), sample = min(colSums(otu.2)))
shan=diversity(t(otu.s),index = "shannon")

#Generate plot data
alpha=data.frame(richness=ric, shannon=shan )
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c( "block","square","lane","season","treatment","AIM.sample.id","plate")])
alpha$season=factor(alpha$season, level = c("nov19","feb20","jul20","nov20"))
alpha$treatment=factor(alpha$treatment, level = niveles)

#block C will not be plotted
met.2=met.2[met.2$block!="C",]
otu.r=otu.r[,rownames(met.2)]
alpha=alpha[alpha$AIM.sample.id %in% rownames(met.2),]


#alpha
tema$axis.text.x$angle=90
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = "treatment")
plot=ggplot(data=mean.r, aes(x=treatment, y=richness,fill=treatment))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values =colores.t)+tema+
  scale_y_continuous(limits = c(0,400),breaks = c(seq(0,400,100)))+
  ylab("OTU richness")

png(paste(result,paste(amplicon,"treatment.richness.lowscale.png",sep = "."),sep = "/"), 
    width = 500, height = 250, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()

mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = "treatment")
plot=ggplot(data=mean.s, aes(x=treatment, y=shannon, fill=treatment))+
  geom_bar(stat="identity",colour="black",size=0.5)+
  geom_errorbar(size=0.5,aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values = colores.t)+tema+
  scale_y_continuous(limits = c(0,6),breaks = c(seq(0,6,2)),labels = scaleFUN)+
  ylab("Shannon index")

png(paste(result,paste(amplicon,"treatment.shannon.lowscale.png",sep = "."),sep = "/"), 
    width = 500, height = 250, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
print(plot) ; dev.off()


#heatmap

dat=amplicon
res=data.frame(row.names = niveles)
index=c("richness","shannon")[2]

for (s in c("feb20","jul20","nov20")) {
  met.s=met.2[met.2$season==s,]
  alpha.s=alpha[alpha$season==s,]
  mean.r=summarySE(data = alpha.s, measurevar = index, groupvars = "treatment")
  a=log2(mean.r[,index]/mean.r[mean.r$treatment=="MOCK",index])
  res[,s]=a
  
}


anot=data.frame(row.names = niveles,
                treatment = niveles)

col=list(treatment=c(MOCK=colores.t[1],
                     PFC=colores.t[2],
                     PFCF=colores.t[3],
                     PFCS=colores.t[4],
                     PFCT=colores.t[5],
                     AAP=colores.t[6],
                     CC=colores.t[7],
                     BC=colores.t[8],
                     RC1=colores.t[9],
                     RC2=colores.t[10],
                     RC3=colores.t[11]))

#Color palette
breaksList = seq(-1, 1, by = 0.01)
der=colorRampPalette(c("#004B40","#009E8E" ,"#89D9CF","white"))
izq=colorRampPalette(c("white","peachpuff2","sandybrown","chocolate4"))
my_pal=c(rev(izq(100)),rev(der(101)))



png(paste(result,paste(amplicon,index,"heatmap.png",sep = "."),sep = "/"), 
    width = 715, height = 275, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")

pheatmap(t(res[niveles,]),cluster_rows = F, cluster_cols = F, color = my_pal,
         annotation_col = anot, annotation_colors = col, fontsize = 7,
         breaks = breaksList, border_color = "black",labels_row = c("t3","t8","t12"), 
         main = paste(index,dat))

dev.off()

png(paste(result,paste(amplicon, index,"heatmap.legend.png",sep = "."),sep = "/"), 
    width = 1000, height = 500, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")

pheatmap(t(res[niveles,]),cluster_rows = F, cluster_cols = F, color = my_pal,
         annotation_col = anot, annotation_colors = col, fontsize = 7, 
         breaks = breaksList, border_color = "black",labels_row = c("t3","t8","t12"), 
         main = paste(index,dat))

dev.off()







  
  
  
  
  
  
  