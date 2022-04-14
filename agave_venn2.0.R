
### Generates Venn diagrams for Figure 1 ###

#Directories
setwd("C:/Users/victorfn/Desktop")  
result="C:/Users/victorfn/Desktop/Results_syncoms/Bars" #result directories
dir.create(result)

s16=("C:/Users/victorfn/Desktop/Results_syncoms/16s")
its=("C:/Users/victorfn/Desktop/Results_syncoms/ITS")     


#Packages
library(vegan)
library(tidyr)
library(dplyr) 
library(ggplot2)
library(MASS) 
library(Rmisc)
library(FSA)
library(VennDiagram)
library(ggVennDiagram)

#Graphics
colores=c("black", "#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
niveles=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
etiquetas=c("Mock","Functional","Functional-F","Functional-As","Functional-At","Photo-anox","Photo-oxi","Biofilm","RC1","RC2","RC3")

tema=theme(axis.text.x = element_text(color="black",size=10, angle=0,hjust=0.5,vjust=0),
           axis.text.y = element_text(color="black",size=10),
           axis.title = element_text(color="black",size=10),
           legend.text = element_text(color = "black",size=10),
           legend.key.size = unit(0.3,"cm"),
           panel.border =element_rect(color = "black", fill = NA) ,#element_blank(),
           strip.text.x = element_text(size=10, color="black",face="bold"),
           strip.text.y = element_text(size=10, color="black",face="bold"),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")

# Set PARAMETERS
amplicon=c("16s","ITS")[1]
bootstrap=0.80 #not implemented yet
filter=1000 #minimum read count for a valid sample
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
base=c("rdp")

# Load metadata
met=read.delim(file=paste(s16,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
met[is.na(met$treatment),"treatment"]="none"

# Load data
if(amplicon=="16s"){
  otu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(s16,paste("taxa.table",data,"rdp","txt",sep = "."),sep = "/"))
} else if (amplicon=="ITS"){
  otu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  tax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
}

#Remove unecessary samples
met.2=met[colnames(otu),]
met.2=met.2[met.2$plant.compartment!="soil" & met.2$treatment!="VOC",]
otu.2=otu[,colnames(otu) %in% rownames(met.2)]

#Sub sample reads
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2)))) #rarefication 
otu.r=t(t(otu.s)/colSums(otu.s)*100)
#otu.r=otu.s

#block C will not be plotted
met.2=met.2[met.2$block!="C",]
otu.r=otu.r[,rownames(met.2)]


p=0.5 #prevalence  or frequency 

n19=rownames(otu.s[rowSums(otu.s[,met.2[met.2$season=="nov19","AIM.sample.id"]]>0)>(length(met.2[met.2$season=="nov19","AIM.sample.id"])*p),])
f20=rownames(otu.s[rowSums(otu.s[,met.2[met.2$season=="feb20","AIM.sample.id"]]>0)>(length(met.2[met.2$season=="feb20","AIM.sample.id"])*p),])
j20=rownames(otu.s[rowSums(otu.s[,met.2[met.2$season=="jul20","AIM.sample.id"]]>0)>(length(met.2[met.2$season=="jul20","AIM.sample.id"])*p),])
n20=rownames(otu.s[rowSums(otu.s[,met.2[met.2$season=="nov20","AIM.sample.id"]]>0)>(length(met.2[met.2$season=="nov20","AIM.sample.id"])*p),])

s=list(Nov19=n19,Feb20=f20,Jul20=j20,Nov29=n20)
#so=calculate.overlap(s)

venn=Venn(s)
plotdata=process_data(venn)

#ggVennDiagram(s, label_alpha = 0, category.names = c("No19","Feb20","Jul20","Nov20"),edge_color="white")

plot=ggplot() +
  geom_sf(aes(fill=count), data = venn_region(plotdata)) +
  scale_fill_gradient(low = "gray90", high = "grey60")+
  geom_sf(size = 1, lty = "solid", color = "black", data = venn_setedge(plotdata), show.legend = F) +
  geom_sf_text(aes(label = c("","","","")), data = venn_setlabel(plotdata),size=5) +
  geom_sf_label(aes(label=count), alpha=0,label.size = 0, size=4.5,family="Arial", data = venn_region(plotdata))+
  theme_void()

png(paste(result,paste(amplicon,"venn.png",sep = "."),sep = "/"), 
    width = 1000, height = 1000, units = "px", pointsize = 15, bg = "white", res=300, type="cairo")
print(plot)
dev.off()

#Relative abundance of core
all=discern_overlap(venn)
mean(colSums(otu.r[all,]))
mean(colSums(otu.r[all,met.2[met.2$season=="nov19","AIM.sample.id"]]))
mean(colSums(otu.r[all,met.2[met.2$season=="feb20","AIM.sample.id"]]))
mean(colSums(otu.r[all,met.2[met.2$season=="jul20","AIM.sample.id"]]))
mean(colSums(otu.r[all,met.2[met.2$season=="nov20","AIM.sample.id"]]))

s.core=tax[tax$otu.id %in% all,]
write.table(s.core,"seasonal.core2.txt",sep = "\t")

