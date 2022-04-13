#!/bin/bash

#Version 2.5
#Script must be run in the main_dir were data is storaged
#main_dir should contain a data folder
#TRIMMOMATIC, FLASH, VSEARCH AND R NEEDED

echo "###### ============================ ######"
echo "###### Setting programs and folders ######"
echo "###### ============================ ######"


#Remove out files so they dont concatenate the same result
rm -r ./stats

#Seting directories
mkdir -p ./trimmed
mkdir -p ./merged
mkdir -p ./process
mkdir -p ./result
mkdir -p ./stats

#set variables #avoid $o, $f, $n
#Trimmomatic
q=12 #minimum quality 
l=30 #minimun length

#Flash 
m=10 #minimum overlap
M=303 #maximum overlap

#Filtering #more can be added if needed
sr=18 #strip bases from right 
sl=39 #strip bases from left
mi=296 #minumum length
ma=394 #maximum length
ee=1 #expected error

#clustering
i=0.95 #threshold OTU for clustering

#otu table
b=0.80 #minimum boostrap value

#Create variables for the tools 
export trimmomatic=/home/astrid/Documents/Trimmomatic-0.39/trimmomatic-0.39.jar
export flash=/home/astrid/Documents/FLASH-1.2.11/flash
export vsearch=/home/astrid/Documents/vsearch-2.18.0-linux-x86_64/bin/vsearch

echo "###### ============================ ######"
echo "###### Q trimming and merging pairs ######"
echo "###### ============================ ######"

#TRIMM low quality bases from raw read
for f in $(find ./data/* -type f | grep _R1_001.fastq.gz); do
o=./trimmed/$(awk -F '/|_L001' '{print $3}' <<< "$f").fastq.gz
java -jar $trimmomatic PE -basein $f -baseout $o -summary ./trimmed/out.txt LEADING:$q TRAILING:$q SLIDINGWINDOW:4:$q MINLEN:$l AVGQUAL:$q
echo $o >> ./stats/trimmed.out.txt
cat ./trimmed/out.txt >> ./stats/trimmed.out.txt
done

#MERGE paired reads
for f in $(find ./trimmed/* -type f | grep _1P.fastq.gz) ; do 
o=./merged/$(awk -F '/|_1P' '{print $3}' <<< "$f")
$flash $f ${f/1P/2P} -o $o -m $m -M $M -z >> ./stats/merged.out.txt 
done

echo "###### ============================ ######"
echo "###### Q filtering and cropping     ######"
echo "###### ============================ ######"

#RENAME, CROP to size and FILTER merged pairs
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_quality.fastq.gz 
$vsearch --fastx_filter $f --relabel $n. --fastq_stripright $sr --fastq_stripleft $sl --fastq_maxlen $ma --fastq_minlen $mi --fastq_maxee $ee --fastqout $o --log ./process/out.txt 
cat ./process/out.txt >> ./stats/filter.out.txt
done

echo "###### ============================ ######"
echo "###### Dereplicating sequences      ######"
echo "###### ============================ ######"

#DEREPLICATE by sample 
for f in $(find ./process/* -type f | grep _quality.fastq.gz) ; do 
o=./process/$(awk -F '/|_quality' '{print $3}' <<< "$f")_derep.fasta
$vsearch --derep_fulllength $f --strand plus --output  $o --sizeout --fasta_width 0 --log ./process/out.txt 
cat ./process/out.txt >> ./stats/derep.sample.out.txt
done 

#DEREPLICATE across samples
cat ./process/*derep.fasta > ./result/all_derep.fasta
$vsearch --derep_fulllength ./result/all_derep.fasta --strand plus --output ./result/all_unique.fasta --sizein --sizeout --fasta_width 0 --uc ./result/all-unique.uc --log ./stats/derep.all.out.txt 

echo "###### ============================ ######"
echo "###### Clustering OTUS              ######"
echo "###### ============================ ######"

# OTU CLUSTERING 
$vsearch --cluster_size ./result/all_unique.fasta --id $i --strand plus --sizein  --sizeout  --fasta_width 0 --centroids ./result/centroids.fasta --log ./result/out.txt
echo Clusters generated: $(grep -c "^>" ./result/centroids.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#Remove SINGLETONS
$vsearch --sortbysize ./result/centroids.fasta --sizein --sizeout  --fasta_width 0 --minsize 2 --output ./result/sorted.fasta --log ./result/out.txt
echo Clusters after removal of singletons: $(grep -c "^>" ./result/sorted.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#Remove CHIMERAS de novo 
$vsearch --uchime_denovo ./result/sorted.fasta --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras ./result/denovo.nonchimeras.fasta --log ./result/out.txt
echo Clusters after removal of chimeras de novo: $(grep -c "^>" ./result/denovo.nonchimeras.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#remove CHIMERAS base on reference data base 16S only
$vsearch --uchime_ref ./result/denovo.nonchimeras.fasta --db ~/Documents/databases/gold.db --sizein --sizeout --fasta_width 0 --qmask none -dbmask none --nonchimeras ./result/nonchimeras.fasta --log ./result/out.txt
echo Clusters after removal of chimeras db: $(grep -c "^>" ./result/nonchimeras.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#RENAME OTUs
$vsearch --fastx_filter ./result/nonchimeras.fasta --fasta_width 0 --relabel OTU --fastaout ./result/otus.fasta

echo "###### ============================ ######"
echo "###### Mapping reads and OTU table  ######"
echo "###### ============================ ######"

#Construct concatenated SEMIRAW reads (merged reads before quality filtering)
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_renamed.fasta 
$vsearch --fastx_filter $f --fastaout $o --relabel $n.
done

cat ./process/*renamed.fasta > ./result/all_semiraw.fasta
rm  ./process/*renamed.fasta

#Create OTU TABLE based on semiraw merged pairs
$vsearch --usearch_global ./result/all_semiraw.fasta --db ./result/otus.fasta --id $i --strand plus --sizein --sizeout --fasta_width 0  --qmask none --dbmask none --otutabout ./result/otutab.txt --log ./stats/table.out.txt

#Create OTU TABLE based on dereplicated and filtered merged pairs
#$vsearch --usearch_global ./result/all_derep.fasta --db ./result/otus.fasta --id $i --strand plus --sizein --sizeout --fasta_width 0  --qmask none --dbmask none --otutabout ./result/otutab.txt --#log ./stats/table.out.txt

echo "###### ============================ ######"
echo "###### Assign taxonomy to OTUs      ######"
echo "###### ============================ ######"

#OTU CLASSIFICATION use the appropriate database for 16S and ITS2 
$vsearch --sintax ./result/otus.fasta --db ~/Documents/databases/rdp_its.udb --tabbedout ./result/otus.sintax.rdp --strand both --sintax_cutoff $b --log ./stats/taxa.out.txt
#$vsearch --sintax ./result/otus.fasta --db ~/Documents/databases/unite_its.udb --tabbedout ./result/otus.sintax.unite --strand both --sintax_cutoff $b --log ./stats/taxa.out.txt

#Past OTU classification. 
#$vsearch --sintax ./result/otus.fasta --db ~/Documents/databases/potu_db.udb --tabbedout ./result/potus.sintax --strand both --sintax_cutoff $b --log ./stats/ptaxa.out.txt

echo "###### ============================ ######"
echo "###### Stats                        ######"
echo "###### ============================ ######"

#Process STAT data
Rscript ../parse.output2_0_ITS.R 
#Rscript ../parse.output2.R 



