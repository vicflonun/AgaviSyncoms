### Performs a PERMANOVA analysis ###

#Packages
library(vegan) 
library(tidyr) 
library(dplyr) 
library(ggplot2) 
library(MASS)
library(Rmisc)
library(extrafont)

#Set directories
setwd("C:/Users/victorfn/Desktop")
result="C:/Users/victorfn/Desktop/Results_syncoms/PERMANOVA"
dir.create(result)

#Graphics
colores.t=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
niveles=c("MOCK","PFC","PFCF","PFCS","PFCT","AAP","CC","BC","RC1","RC2","RC3")
etiquetas=c("Mock","Functional","Functional-F","Functional-As","Functional-At","Photo-anox","Photo-oxi","Biofilm","RC1","RC2","RC3")
colores.s=c("#BE9E80","#886B4F","#447690","#594635")

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

#Sub sample reads
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

# PERMANOVA all
set.seed(23171341)
scaling=vegdist(t(otu.r), method = "bray", binary = T) 
#adonis2(scaling~season, data = met.2, permutations = 1000)
#adonis2(scaling~block, data = met.2, permutations = 1000)
#adonis2(scaling~treatment, data = met.2, permutations = 1000)
per=adonis2(scaling~season*block*treatment, data = met.2, permutations = 1000)

write.table(per, paste(result,paste(amplicon,"permanova.global.txt",sep="."),sep="/"),sep = "\t")


# PERMANOVA nov19
met.s=met.2[met.2$season=="nov19" & met.2$block!="C" ,]
set.seed(23171341)
scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T) 
per=adonis2(scaling~block, data = met.s, permutations = 1000)
write.table(per, paste(result,paste(amplicon,"permanova.nov19.block.txt",sep="."),sep="/"),sep = "\t")

per=adonis2(scaling~treatment, data = met.s, permutations = 1000)
write.table(per, paste(result,paste(amplicon,"permanova.nov19.treatment.txt",sep="."),sep="/"),sep = "\t")

# PERMANOVA feb20
met.s=met.2[met.2$season=="feb20",]
set.seed(23171341)
scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T) 
per=adonis2(scaling~block*treatment, data = met.s, permutations = 1000)

write.table(per, paste(result,paste(amplicon,"permanova.feb20.txt",sep="."),sep="/"),sep = "\t")

# PERMANOVA jul20
met.s=met.2[met.2$season=="jul20",]
set.seed(23171341)
scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T) 
per=adonis2(scaling~block*treatment, data = met.s, permutations = 1000)
write.table(per, paste(result,paste(amplicon,"permanova.jul20.txt",sep="."),sep="/"),sep = "\t")

# PERMANOVA nov20
met.s=met.2[met.2$season=="nov20",]
set.seed(23171341)
scaling=vegdist(t(otu.r[,rownames(met.s)]), method = "bray", binary = T) 
per=adonis2(scaling~block*treatment, data = met.s, permutations = 1000)
write.table(per, paste(result,paste(amplicon,"permanova.nov20.txt",sep="."),sep="/"),sep = "\t")











