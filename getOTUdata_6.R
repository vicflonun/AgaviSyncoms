##### Generates filtered OTU, taxa and metadata tables for upstream analysis #####

    # 1. Removes not processed libraries
    # 2. Removes control libraries 
    # 3. Removes contaminant OTUs
    # 4. Removes not microbial OTUs
    # 5. Removes low frequent and abundant OTUs
    # 6. Sub sample the libraries
    # Optional plots

library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Biostrings)

  # Select the files
table=c("otutab.txt")[1]
metadata=c("vicdani_metadata.txt")[1] #this can be a different name but must be the complete 309 sample metadata
fasta=c("otus.fasta")[1]

  # Set PARAMETERS
amplicon=c("16s","ITS")[1]
bootstrap=0.80 #not implemented yet
filter=1000 #minimum read count for a valid sample
data=c("syncoms","daniel", "all")[1]  #what sub sample would you need
rm.na=FALSE #remove unclassified sequences, not implemented yet but could be an option

if ( amplicon == "16s"){
  setwd("~/Documents/agave16s/result_def")
  taxa=c(rdp="otus.sintax.rdp",silva="otus.sintax.silva",green="otus.sintax.green")[1]
  mr=7 #minimum read count for abundant OTUs
  ms=5 #minim sample frequency for abundant OTUs 
  nr=7 #maximum read count for contaminant OTUs
  ns=2 #maximum sample frequency for contaminant OTUs 
} else if ( amplicon == "ITS" ){ 
  setwd("~/Documents/agaveITS/result_dani") 
  taxa=c(rdp="otus.sintax.rdp", unite="otus.sintax.unite", uniteall="otus.sintax.uniteall")[3]
  mr=2 #minimum read count for abundant OTUs
  ms=5 #minim sample frequency for abundant OTUs 
  nr=2 #maximum read count for contaminant OTUs
  ns=2 #maximum sample frequency for contaminant OTUs 
}

  # Load OTU table
otu=read.delim(paste(getwd(),table,sep = "/"), header = T) #load
otu$X.OTU.ID=as.numeric(gsub("OTU","",otu$X.OTU.ID)) #change OTU names and sort
otu=otu[order(otu$X.OTU.ID),]
otu$X.OTU.ID=paste("OTU",otu$X.OTU.ID,sep="_")
rownames(otu)=otu$X.OTU.ID ; otu=otu[,-1]

  # Load TAXA classification
ranks=c("otu.id","domain","phylum","class","order","family","genus")
sintax=read.delim(paste(getwd(),taxa,sep = "/"), header = F)
sintax=sintax[,-c(2,3,5)] #Format taxa
sintax=splitstackshape::cSplit(indt = sintax, splitCols = "V4", sep = ",", type.convert = F)
sintax=data.frame(sintax[,c(1:7)])
colnames(sintax)=ranks 

sintax$otu.id=as.numeric(gsub("OTU","",sintax$otu.id)) #change OTU names and sort
sintax=sintax[order(sintax$otu.id),]
sintax$otu.id=paste("OTU",sintax$otu.id,sep="_")
rownames(sintax)=sintax$otu.id ; sintax=sintax[,c(2:7,1)]

  # Load METADATA 
meta=read.delim(paste("~/Documents/agaveITS",metadata,sep = "/"))
meta$AIM.sample.id=gsub("-","_",meta$AIM.sample.id)
rownames(meta)=meta$AIM.sample.id
meta=meta[order(meta$AIM.sample.id),] #order samples
meta$plate=c(rep("E116",80),rep("E117",80),rep("E118",80),rep("E119",69)) #add plate information 

  # Load SEQUENCES
seq=readDNAStringSet(paste(getwd(),fasta,sep = "/"))
names(seq)=gsub("OTU","OTU_",names(seq)) #change names

  # FILTERING SAMPLES
total=colSums(otu) #calculate total read per sample 

 # Samples that were no processed
missing=c() 
m=meta$AIM.sample.id[!meta$AIM.sample.id %in% str_extract(colnames(otu),"20210219\\D+\\d+")] #identify missing samples
print("Not processed") ; print(m)
missing=c(missing,m)

 # Identify  OTUS in negative controls 
otu.n=otu[,grep("neg",colnames(otu))] #subset the negative controls
otu.n=otu.n[rowSums(otu.n>nr)>ns,] #subset the contaminant OTUs
n=rownames(otu.n)
print(paste("Contaminants", length(n), "OTUs"))
sum(rowSums(otu[n,]))/ sum(total) *100 # Percentage reads contaminants 
                                       #**might be very high because a plant OTU resulted a contaminant

 ## Identify  non microbial sequences
    #For 16s is easier to identify them in any database for ITS only using the UNITE all Eukaryotas database

if ( amplicon == "16s"){
  e=unique(c(sintax[grep("Eukaryota",sintax$domain),"otu.id"],sintax[grep("Mitochondria",sintax$family), "otu.id"]))
  c=sintax[grep("Chloroplast",sintax$class),"otu.id"]

  print(paste("Chloroplast", length(c), "OTUs"))
  print(sum(rowSums(otu[c,])) / sum(total) *100)
  
  print(paste("Eukaryota", length(e), "OTUs"))
  print(sum(rowSums(otu[e,])) / sum(total) * 100)
  
  u=sintax[is.na(sintax$domain),"otu.id"]
  print(paste("Unclassified", length(u), "OTUs"))
  print(sum(rowSums(otu[u,])) / sum(total) * 100)
  
} else if (amplicon == "ITS"){
  
  e=unique(sintax[grep("Viridiplantae",sintax$domain),"otu.id"])
  print(paste("Viridiplantae", length(e), "OTUs"))
  print(sum(rowSums(otu[e,])) / sum(total) * 100)
  
  u=unique(sintax[grep("unidentified",sintax$domain),"otu.id"])
  print(paste("unidentified", length(u), "OTUs"))
  print(sum(rowSums(otu[u,])) / sum(total) * 100)
  
  u=sintax[is.na(sintax$domain),"otu.id"]
  print(paste("Unclassified", length(u), "OTUs"))
  print(sum(rowSums(otu[u,])) / sum(total) * 100)
}


## Remove not desired sequences and negative samples
otu.a=otu[!rownames(otu) %in% unique(c(n,e,c)),grep("neg",colnames(otu), invert = T)] 

## Change the name of the samples 
colnames(otu.a)=str_extract(colnames(otu.a),"20210219\\D+\\d+")

## Recalculate total 
total.a=colSums(otu.a)



## Check low read samples
read=data.frame(total=total.a,
                plate=meta[colnames(otu.a),"plate"],
                name=str_extract(names(total.a),"C_\\d+"), 
                row.names =names(total.a))

ggplot(data=read[,], aes(x=plate, y=total, fill=plate))+
  geom_boxplot()+scale_y_continuous(breaks = seq(0,max(total.a),5000), limits =c(0,max(total.a)) )
ggplot(data=read, aes(x=total, fill=plate))+
  geom_histogram()+ scale_x_continuous(breaks = seq(0,80000,5000))
ggplot(data=read, aes(x=log10(total), fill=plate))+
  geom_histogram()+ scale_x_continuous()

print("Losing samples due to low read count")
summary(as.factor(read[read$total<filter,"plate"]))

m=str_extract(names(total.a[total.a<filter]),"20210219\\D+\\d+")
missing=c(missing,m)

 ## Remove low read samples

total.b=total.a[total.a>=filter]
otu.b=otu.a[,names(total.b)]
read.b=read[colnames(otu.b),]

## Subset data if needed, (note that we are not removing samples from metadata at this step)

if (data == "syncoms"){ 
    metas=meta[meta$sample.by=="Victor" | meta$sample.by=="David",]
    otus=otu.b[,colnames(otu.b) %in% sort(metas[,"AIM.sample.id"])]
    reads=read.b[colnames(otus),]
  } else if (data == "daniel"){
    metas=meta[meta$sample.by=="Daniel",]
    otus=otu.b[,colnames(otu.b) %in% sort(metas[,"AIM.sample.id"])]
    reads=read.b[colnames(otus),]
  } else {
    otus=otu.b
    reads=read.b
    metas=meta
  }

  ## Removing OTUs with cero reads (important when subseting)
otus=otus[rowSums(otus==0)<dim(otus)[2],]

    ## Optional: Compare the read number against the richness

rarecurve(t(otus), step = 2000, xlab = "Sample Size", ylab = "Species", label = F)
reads$richness=colSums(otus>0)
ggplot(data=reads,aes(x=total,y=richness))+
  geom_point(aes(color=plate), size=2)+
  geom_smooth(method = "lm",color="black")+
  scale_x_continuous(breaks = seq(0,max(read$total),5000))

## Removing OTUs with low reads 

otus=otus[rowSums(otus>=mr)>=ms,]
#otus=otus[rowSums(otus>=mr)>=(dim(otus)[2]*.20),]

    ## Optional: Compare the read number against the richness

rarecurve(t(otus), step = 2000, xlab = "Sample Size", ylab = "Species", label = F)
reads$richness=colSums(otus>0)
reads$total=colSums(otus)
ggplot(data=reads,aes(x=total,y=richness))+
  geom_point(aes(color=plate))+
  geom_smooth(method = "lm",color="black")+
  scale_x_continuous(breaks = seq(0,max(read$total),5000))


  ## Generate metadata

metas$done="yes"
metas[metas$AIM.sample.id %in% missing, "done"] = "no"

  ## Extract sequences

seqs=seq[rownames(otus)]

  ## Extract taxa


taxas=sintax[rownames(otus),]


 ## Save the data

write.table(otus,
            file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")


write.table(metas,
            file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")


write.table(taxas,
            file=paste(getwd(),paste("taxa.table",data,names(taxa),"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")


writeXStringSet(seqs, 
                file=paste(getwd(),paste("otus",data,"fasta",sep = "."),sep = "/"),
                format = "fasta")

