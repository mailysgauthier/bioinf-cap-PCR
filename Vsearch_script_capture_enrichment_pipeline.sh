#!/bin/bash
#Analysis pipeline for capture raw data 

###########################################################################
#I.PROCESSING OF RAW SAMPLES
###########################################################################

mkdir /0.Raw
#Place capture sequencing raw fastq.qz files in this directory
raw=/0.Raw

#I.1.VSEARCH PIPELINE : MERGE + TRIMMING + QUALITY FILTERING 
#Defining variables
mkdir /1.Process.raw.samples.vsearch/
mkdir /1.Process.raw.samples.vsearch/1.1.Merging
mkdir /1.Process.raw.samples.vsearch/1.2.Trimming
mkdir /1.Process.raw.samples.vsearch/1.3.Filtering

source=/1.Process.raw.samples.vsearch/
merge=/1.1.Merging
trim=/1.2.Trimming
clean=/1.3.Filtering

VSEARCH=$(which vsearch)
THREADS=5

#Copy raw fastq to raw
cd $raw
gunzip *.fastq.gz 

#Process samples 
for f in *_R1_*.fastq; do
	cd $raw
    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(cut -d_ -f1 <<< "$f")

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================

    $VSEARCH --threads $THREADS \
             --fastq_mergepairs $f \
             --reverse $r \
             --fastq_minovlen 10 \
             --fastq_maxdiffs 5 \
             --fastqout $merge/$s.merged.fastq \
             --fastq_eeout \
             --eetabbedout $merge/$s.error.stats.mergepairs.txt \
             --fastqout_noMC2erged_fwd $merge/$s.R1.unmerged.fastq \
             --fastqout_noMC2erged_rev $merge/$s.R2.unmerged.fastq  &>> $merge/$s.vsearch.merge.log                               

	cd $merge
    for m in $s*.fastq; do

    n=$(cut -d. -f-2 <<< "$m")

	echo
    echo Trimming of adaptors	  
    fastq-mcf -o $trim/$n.trimmed.fastq -q 30 -x 20 -S $source/adaptors.fasta $m &>> $trim/$n.fastqmcf.trimming.log

    echo
    echo Calculate quality statistics
    $VSEARCH --threads $THREADS \
             --fastq_eestats $trim/$n.trimmed.fastq \
             --output $trim/$n.stats

    echo
    echo Quality filtering

    $VSEARCH --threads $THREADS \
             --fastq_filter $trim/$n.trimmed.fastq \
             --fastq_maxee 1 \
             --fastq_minlen 150 \
             --fastq_maxns 50 \
             --fastq_stripleft 1 \
             --fastqout $clean/$n.filtered.fastq &>> $clean/$n.vsearch.quality.log 

     done
done

#Concatenate output of merging
cd $merge
tail -vn +1 *vsearch.merge.log > All.stat.merge.log
grep 'Merged' -A 1  *vsearch.merge.log > All.stat.merge.count.log

#Concatenate output of quality filter
cd $clean
tail -vn +1 *vsearch.quality.log > All.stat.quality.log
fastqc *.fastq
multiqc .

#Count the number of unique filtered sequences 
cd $clean
for f in *R1.filtered.fastq ; do

    r=$(sed -e "s/.R1./.R2./" <<< "$f")
    s=$(cut -d. -f1 <<< "$f")

	grep '^@M' < $f | sort -k1,1 > $s.r1.name
	grep '^@M' < $r | sort -k1,1 > $s.r2.name
	cat $s.r1.name $s.r2.name | sort -u -k1,1 | wc -l > $s.filtered.count.uniq.r1.r2.out
done

tail -v  *.filtered.count.uniq.r1.r2.out > All.filtered.count.uniq.r1-r2.out
rm *.filtered.count*.r1.r2.out
rm *.name 



###########################################################################
#II.ASSIGNATION - VSEARCH PIPELINE : ALIGNMENT + ABUNDANCE QUANTIFICATION
###########################################################################
#First on 16S + mt genes -> remove of the sequences in the filtered fasta files -> db COI db -> abundance quantification
#Personal haplotype db for 10 species mock communities (hereafter MC1) samples = COI.10sps 
#Personal haplotype db for 52 taxa mock communities (hereafter MC2) samples = MC2.haplotypes

############################################################################################
###### Reference databases

####16S + mt genes (non COI) database####
####Database containing NAD1, NAD4, NAD5, CYTB and ATP6 complete genes of 
#Asellus aquaticus ADA69754.1
#Daphnia pulex AAD33231.1
#Dinocras cephalotes AGZ03516.1
#Gammarus fossarum YP_009379680.1
#Physella acuta YP_008994230.1
#Radix balthica HQ330989.1
#Sericostoma personatum AJR19241.1
#Thremma gallicum AJR19254.1
####plus 16S of the MC1 (GenBank Accession Numbers MK584525:MK584534)
####plus 16S corresponding to the MC2 (downloaded from NCBI)
#(NC_023927, FJ263197.1, EF472580.1, KX087288.1, NC_029246, NC_016439, EF623244.1, KX453695.1, FJ443044.1, NC_016167.1, KT876896.1, KY568891.1, KF855860.1, NC_026219.1, EF623195.1,
#LN868660.1, KJ493406.1, EF623176.1, KT970065.1, GU073059.1, AF428198.1, AY749775.1, EF623221.1, FJ002820.1, DQ062646.1, FN806889.1, AY620141.1, AF117817, EF623194.1, HQ330989,
#HM637016.1, KY261229.1, AY577466.1, FJ443047.1, GQ118289.1, EF623192.1, KP455291, EU215181.1, EU215174.1, KP455290, FN806898.1, DQ198627.1, AY648139.1, EU477688.1, GU130252,
#FJ471666.1, EF623193.1, HG810353.1, HG810354.1, HG810363.1, HG810360.1, FJ002822.1, EU005435.1, NC_016203, KT251040.1, NC_016129, NC_023927, NC_014683)
mkdir /Database/DB_16S_mt_genes
#Put sequences of the database in this folder
cd /Database/DB_16S_mt_genes
cat *.fas > all_nonCOI_genes.fas
makeblastdb -in all_nonCOI_genes.fas -dbtype nucl
MC_db_16Sco=/Database/DB_16S_mt_genes/all_nonCOI_genes.fas

####COI sequences barcoded for the MC1 (MK584516:MK584524)
mkdir /Database/DB_COI_10sps/
#Put sequences of the database in this folder
cat /Database/DB_COI_10sps/*.fasta > COI_database_10sps.fas
MC1_db_COI=/Database/DB_COI_10sps/COI_database_10sps.fas
makeblastdb -in $MC1_db_COI -dbtype nucl

####Contaminants database: sequences used in Francois et al. (in review)
mkdir /Database/DB_contam/
#Put sequences of the database in this folder
cat /Database/DB_contam/.fasta > /Database/DB_contam/customDB_arthropod_db.fa
diamond makedb --in /Database/DB_contam/customDB_arthropod_db.fa -d nr
MC_db_conta_prot=/Database/DB_contam/customDB_arthropod_db.fa.dmnd
list_protostomia = list_protostomia #Available in github folder Resources 

############################################################################################
####MC1 samples####
clean=/1.Process.raw.samples.vsearch/1.3.Filtering #fastq files
mkdir /2.Assignation.vsearch/
mkdir /2.Assignation.vsearch/2.1.MC1
MC_root=/2.Assignation.vsearch.2019/2.1.MC1

cd $MC1_root 

mkdir $MC1_root/2.1.1.Cat.filtered
mkdir $MC1_root/2.1.2.16S.co
mkdir $MC1_root/2.1.3.nonCOI.filtered.out
mkdir $MC1_root/2.1.4.COI.10sps
mkdir $MC1_root/2.1.5.COI.blastx
mkdir $MC1_root/2.1.6.Conta
mkdir $MC1_root/2.1.7.Clean.table

MC_cat=$MC1_root/2.1.1.Cat.filtered
MC_16Sco=$MC1_root/2.1.2.16S.co
MC_fas=$MC1_root/2.1.3.nonCOI.filtered.out
MC_COI=$MC1_root/2.1.4.COI.10sps
MC_blastx=$MC1_root/2.1.5.COI.blastx
MC_conta=$MC1_root/2.1.6.Conta
MC_clean=$MC1_root/2.1.7.Clean.table

#Fastq to fasta
cd $clean
for f in MC*filtered.fastq; do
s=$(cut -d. -f-2 <<< "$f")
awk 'NR % 4 == 1 || NR%4 == 2' $f | sed -e 's/@/>/' > $MC1_root/$s.filtered.fasta
done

#Cat of the filtered files
cd $MC1_root
for f in $MC1_root/*.merged.*.fasta; do
R1=$(sed 's/.merged./.R1./' <<< "$f")
R2=$(sed 's/.merged./.R2./' <<< "$f")
s=$(cut -d. -f1 <<< "$f")
cat $R1 $R2 $f > $MC1_cat/$s.all.filtered.fasta
done

######################################################################################################
cd $MC1_cat
for f in *.fasta; do

s=$(cut -d. -f1 <<< "$f")
grep '>' < $f | sed 's/>//g' | cut -f1 -d' '| sort -k1,1 > $MC1_cat/$s.seq.names

echo
echo ==========================
echo  Processing $s 
echo ==========================

#Blastn short on the 16S and non COI mt genes (sequences to remove)
echo 16s and co alignment
 
blastn  -db  $MC1_db_16Sco \
        -task 'blastn-short' -query $f \
        -outfmt 6 -penalty -1 -dust no \
        -evalue 1E-10 -num_threads 30 > $MC1_16Sco/$s.blastn.16Sco.out
sort -u -k1,1 $MC1_16Sco/$s.blastn.16Sco.out | cut -f1 | sort > $MC1_16Sco/$s.blastn.16Sco.hits


#Remove 16S&co + chimera reads from the fasta files
comm -13 $MC1_16Sco/$s.blastn.16Sco.hits $MC1_cat/$s.seq.names > $MC1_fas/$s.seq2conserve
select_sequences.pl $f $MC1_fas/$s.seq2conserve $MC1_fas/$s.filt.out.fas #perl script available here : 

#10 species COI db alignment 
echo COI 10sps alignment

blastn -db $MC1_db_COI -task 'blastn' -query $MC1_fas/$s.filt.out.fas -outfmt 6 -evalue 1E-10 > $MC1_COI/$s.aln.COI.blastn.out 
awk '$3>=97' < $MC1_COI/$s.aln.COI.blastn.out | awk '$4>=250' | sort -k1,1 -k3nr,3nr | sort -u -k1,1 > $MC1_COI/$s.aln.COI.blastn.uniq.id97.qc250.out
cut -f1 $MC1_COI/$s.aln.COI.blastn.uniq.id97.qc250.out > $MC1_COI/$s.aln.COI.blastn.uniq.id97.qc250.hits

#contamination db alignment
diamond blastx -d $MC1_db_conta_prot -q $MC1_fas/$s.filt.out.fas --evalue 1E-10 --query-gencode 5 --threads 20 --more-sensitive --out $MC1_conta/$s.conta.blastx.msens.out >> $MC1_conta/blastx.conta.out.log
sort -u -k1,1 $MC1_conta/$s.conta.blastx.msens.out > $MC1_conta/$s.conta.blastx.uniq.msens.out	
grep -f $list_protostomia $MC1_conta/$s.conta.blastx.uniq.msens.out | awk '{print $1, "protostomia"}' > $MC1_conta/$s.conta.blastx.uniq.msens.proto.out
grep -v -f $list_protostomia $MC1_conta/$s.conta.blastx.uniq.msens.out | sed 's/|/\t/g' | sed 's/\t/ /g' | cut -f1,3 -d' ' >> $MC1_conta/$s.conta.blastx.uniq.msens.proto.out
for cat in `cat $categ_conta`; do
echo $s $cat `cut -f2 -d' ' $MC1_conta/$s.conta.blastx.uniq.msens.proto.out | grep -c "$cat"` >> $MC1_clean/quanti.conta.blastx.msens.proto
done

done

#Count of the abundance per taxa for COI 10sps - Script available in R_scripts folder
Rscript 0.Rscript.affiliation.count.blastn.paper.R $MC1_COI $MC1_clean


############################################################################################
####MC2 samples####
clean=/1.Process.raw.samples.vsearch/1.3.Filtering #fastq files
mkdir /2.Assignation.vsearch/2.2.MC2
MC2_root=/2.Assignation.vsearch/2.2.MC2

cd $MC2_root 

mkdir $MC2_root/2.2.1.Cat.filtered
mkdir $MC2_root/2.2.2.16S.co
mkdir $MC2_root/2.2.3.nonCOI.filtered.out
mkdir $MC2_root/2.2.4.COI.10sps
mkdir $MC2_root/2.2.5.COI.blastx
mkdir $MC2_root/2.2.6.Conta
mkdir $MC2_root/2.2.7.Clean.table

MC2_cat=$MC2_root/2.2.1.Cat.filtered
MC2_16Sco=$MC2_root/2.2.2.16S.co
MC2_fas=$MC2_root/2.2.3.nonCOI.filtered.out
MC2_COI=$MC2_root/2.2.4.COI.10sps
MC2_blastx=$MC2_root/2.2.5.COI.blastx
MC2_conta=$MC2_root/2.2.6.Conta
MC2_clean=$MC2_root/2.2.7.Clean.table

############################################################################################
###### Reference databases

####16S + mt genes (non COI) database####
####Database containing NAD1, NAD4, NAD5, CYTB and ATP6 complete genes of 
#Asellus aquaticus ADA69754.1
#Daphnia pulex AAD33231.1
#Dinocras cephalotes AGZ03516.1
#Gammarus fossarum YP_009379680.1
#Physella acuta YP_008994230.1
#Radix balthica HQ330989.1
#Sericostoma personatum AJR19241.1
#Thremma gallicum AJR19254.1
####plus 16S of the 10 species mock communities (hereafter MC1) (GenBank Accession Numbers MK584525:MK584534)
####plus 16S corresponding to the 52 taxa mock communities (hereafter MC2) (downloaded from NCBI)
#(NC_023927, FJ263197.1, EF472580.1, KX087288.1, NC_029246, NC_016439, EF623244.1, KX453695.1, FJ443044.1, NC_016167.1, KT876896.1, KY568891.1, KF855860.1, NC_026219.1, EF623195.1,
#LN868660.1, KJ493406.1, EF623176.1, KT970065.1, GU073059.1, AF428198.1, AY749775.1, EF623221.1, FJ002820.1, DQ062646.1, FN806889.1, AY620141.1, AF117817, EF623194.1, HQ330989,
#HM637016.1, KY261229.1, AY577466.1, FJ443047.1, GQ118289.1, EF623192.1, KP455291, EU215181.1, EU215174.1, KP455290, FN806898.1, DQ198627.1, AY648139.1, EU477688.1, GU130252,
#FJ471666.1, EF623193.1, HG810353.1, HG810354.1, HG810363.1, HG810360.1, FJ002822.1, EU005435.1, NC_016203, KT251040.1, NC_016129, NC_023927, NC_014683)
mkdir /Database/DB_16S_mt_genes
#Put sequences of the database in this folder
cd /Database/DB_16S_mt_genes
cat *.fas > all_nonCOI_genes.fas
makeblastdb -in all_nonCOI_genes.fas -dbtype nucl
MC2_db_16Sco=/Database/DB_16S_mt_genes/all_nonCOI_genes.fas

####COI sequences barcoded for the MC2 
mkdir /Database/DB_MC2_haplotypes
#Put sequences of the database in this folder
MC2_db_COI=/Database/DB_MC2_haplotypes
MC2_db_COI_A=/Database/DB_MC2_haplotypes/TierMix-A.COI.ref.fasta
MC2_db_COI_B=/Database/DB_MC2_haplotypes/TierMix-B.COI.ref.fasta
MC2_db_COI_C=s/Database/DB_MC2_haplotypes/TierMix-C.COI.ref.fasta
MC2_db_COI_D=/Database/DB_MC2_haplotypes/TierMix-D.COI.ref.fasta
MC2_db_COI_E=/Database/DB_MC2_haplotypes/TierMix-E.COI.ref.fasta
MC2_db_COI_F=/Database/DB_MC2_haplotypes/TierMix-F.COI.ref.fasta
MC2_db_COI_G=/Database/DB_MC2_haplotypes/TierMix-G.COI.ref.fasta
MC2_db_COI_H=/Database/DB_MC2_haplotypes/TierMix-H.COI.ref.fasta
MC2_db_COI_i=/Database/DB_MC2_haplotypes/TierMix-i.COI.ref.fasta
MC2_db_COI_j=/Database/DB_MC2_haplotypes/TierMix-j.COI.ref.fasta
#Files made from sequences of Elbrecht and Leese, 2015 
#Available in github folder Resources 
for f in /Database/DB_MC2_haplotypes/*COI.ref.fasta; do
	makeblastdb -in $f -dbtype nucl
done

####Contaminants database: sequences used in Francois et al. (in review)
mkdir /Database/DB_contam/
#Put sequences of the database in this folder
cat .fasta > customDB_arthropod_db.fa
diamond makedb --in customDB_arthropod_db.fa -d nr
MC2_db_conta_prot=/Database/customDB_arthropod_db.fa.dmnd
list_protostomia = list_protostomia #Available in github folder Resources 

#Fastq to fasta
cd $clean
for f in TierMix*filtered.fastq; do
s=$(cut -d. -f-2 <<< "$f")
	awk 'NR % 4 == 1 || NR%4 == 2' $f | sed -e 's/@/>/' > $MC2_root/$s.filtered.fasta
done

#Cat of the filtered files
cd $MC2_root
for f in *.merged.*.fasta; do
	R1=$(sed 's/.merged./.R1./' <<< "$f")
	R2=$(sed 's/.merged./.R2./' <<< "$f")
	s=$(cut -d. -f1 <<< "$f")
	cat $R1 $R2 $f > $MC2_cat/$s.all.filtered.fasta
done


######################################################################################################
cd $MC2_cat
for f in *.fasta; do

s=$(cut -d. -f1 <<< "$f")
grep '>' < $f | sed 's/>//g' | cut -f1 -d' '| sort -k1,1 > $MC2_cat/$s.seq.names

echo
echo ==========================
echo  Processing $s 
echo ==========================

#Blastn short on the 16S and non COI mt genes (sequences to remove)
echo 16s and co alignment
 
blastn  -db  $MC2_db_16Sco \
        -task 'blastn-short' -query $f \
        -outfmt 6 -penalty -1 -dust no \
        -evalue 1E-10 -num_threads 30 > $MC2_16Sco/$s.blastn.16Sco.out
sort -u -k1,1 $MC2_16Sco/$s.blastn.16Sco.out | cut -f1 | sort > $MC2_16Sco/$s.blastn.16Sco.hits


#Remove 16S&co + chimera reads from the fasta files
comm -13 $MC2_16Sco/$s.blastn.16Sco.hits $MC2_cat/$s.seq.names > $MC2_fas/$s.seq2conserve
select_sequences.pl $f $MC2_fas/$s.seq2conserve $MC2_fas/$s.filt.out.fas #perl script available here: 

#10 species COI db alignment 
echo COI 10sps alignment
blastn -db $MC2_db_COI/$s.COI.ref.fasta  -task 'blastn' -query $MC2_fas/$s.filt.out.fas -outfmt 6 -evalue 1E-10 > $MC2_COI/$s.aln.COI.blastn.out 
awk '$3>=97' < $MC2_COI/$s.aln.COI.blastn.out | awk '$4>=250' | sort -k1,1 -k3nr,3nr | sort -u -k1,1 > $MC2_COI/$s.aln.COI.blastn.uniq.id97.qc250.out
cut -f1 $MC2_COI/$s.aln.COI.blastn.uniq.id97.qc250.out > $MC2_COI/$s.aln.COI.blastn.uniq.id97.qc250.hits
sed 's/_[a-zA-Z]_Assembly//g'  $MC2_COI/$s.aln.COI.blastn.uniq.id97.qc250.out | sed 's/_j_BF1_BR2_66_Bivalvia//g' | sed 's/_g_BF1_BR2_102_Bivalvia//g' | sed 's/_Contig_[0-9]//g' > $MC2_COI/$s.aln.COI.blastn.uniq.id97.qc250.clean.out

#contamination db alignment
diamond blastx -d $MC2_db_conta_prot -q $MC2_fas/$s.filt.out.fas --evalue 1E-10 --query-gencode 5 --threads 20 --more-sensitive --out $MC2_conta/$s.conta.blastx.msens.out >> $MC2_conta/blastx.conta.out.log
sort -u -k1,1 $MC2_conta/$s.conta.blastx.msens.out > $MC2_conta/$s.conta.blastx.uniq.msens.out	
grep -f $list_protostomia $MC2_conta/$s.conta.blastx.uniq.msens.out | awk '{print $1, "protostomia"}' > $MC2_conta/$s.conta.blastx.uniq.msens.proto.out
grep -v -f $list_protostomia $MC2_conta/$s.conta.blastx.uniq.msens.out | sed 's/|/\t/g' | sed 's/\t/ /g' | cut -f1,3 -d' ' >> $MC2_conta/$s.conta.blastx.uniq.msens.proto.out
for cat in `cat $categ_conta`; do
echo $s $cat `cut -f2 -d' ' $MC2_conta/$s.conta.blastx.uniq.msens.proto.out | grep -c "$cat"` >> $MC2_clean/quanti.conta.blastx.msens.proto
done

done

#Count of the abundance per taxa for COI 52 taxa - Script available in R_scripts folder
Rscript /0.Rscript.affiliation.count.blastn.MC2.paper.R $MC2_COI $MC2_clean
