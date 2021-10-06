#!/bin/bash
#Analysis pipeline for PCR enrichment raw data 

###########################################################################
#I.PROCESSING OF RAW SAMPLES
###########################################################################

mkdir /0.Raw
raw=/0.Raw
#Place metabarcoding sequencing raw fastq.qz files in this directory

mkdir /1.0.Merging
mkdir /1.1.Trimming
mkdir /1.2.Filtering
mkdir /1.3.Dereplication
mkdir /1.4.Preclustering
mkdir /1.5.Chimera
mkdir /1.6.Clustering
mkdir /1.7.Assignement

merge=/1.0.Merging
trim=/1.1.Trimming
clean=/1.2.Filtering
derep=/1.3.Dereplication
preclust=/1.4.Preclustering
chim=/1.5.Chimera
clust=/1.6.Clustering
ass=/1.7.Assignement

VSEARCH=$(which vsearch)
PERL=$(which perl)
THREADS=5

#Rename sequence to avoid problems in dowstream analysis
rename -e 's/-//g' -f -v $raw/*

cd $raw
#Process samples amplified with BF2/BR2 primer pair
for f in *bfr*_R1_001.fastq; do

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
             --fastqout_notmerged_fwd $merge/$s.R1.unmerged.fastq \
             --fastqout_notmerged_rev $merge/$s.R2.unmerged.fastq  &>> $merge/$s.vsearch.merge.log                               

    echo Trimming of primers BFR
    cutadapt3 -g ^GCHCCHGAYATRGCHTTYCC --discard-untrimmed  -o $trim/temp.fastq $merge/$s.merged.fastq &>> $trim/$s.vsearch.trim.fw.log
    cutadapt3 -a TGRTTYTTYGGNCAYCCHGA$  -o $trim/$s.trimmed.fastq $trim/temp.fastq &>> $trim/$s.vsearch.trim.fw.rv.log
    rm $trim/temp.fastq

	echo
    echo Calculate quality statistics
    $VSEARCH --threads $THREADS \
             --fastq_eestats $trim/$s.trimmed.fastq \
             --output $trim/$s.stats

    echo
    echo Quality filtering

    $VSEARCH --threads $THREADS \
             --fastq_filter $trim/$s.trimmed.fastq \
             --fastq_maxee 1 \
             --fastq_minlen 200 \
             --fastq_maxns 0 \
             --fastqout $clean/$s.filtered.fastq &>> $clean/$s.vsearch.quality.log 

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --threads $THREADS \
        --derep_fulllength $clean/$s.filtered.fastq \
        --strand plus \
        --output $derep/$s.derep.fasta \
        --sizeout \
        --uc $derep/$s.derep.uc \
        --relabel $s. \
        --relabel_keep \
        --fasta_width 0 &>> $derep/$s.vsearch.derep.log
done

#Process samples amplified with Fwh1 primer pair
for f in *Fwh*_R1_001.fastq; do

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
             --fastqout_notmerged_fwd $merge/$s.R1.unmerged.fastq \
             --fastqout_notmerged_rev $merge/$s.R2.unmerged.fastq  &>> $merge/$s.vsearch.merge.log                                

    echo Trimming of primers Fwh
    cutadapt3 -g ^YTCHACWAAYCAYAARGAYATYGG --discard-untrimmed  -o $trim/temp.fastq $merge/$s.merged.fastq &>> $trim/$s.vsearch.trim.fw.log
    cutadapt3 -a GGDGGDTTYGGWAAYTGAYT$  -o $trim/$s.trimmed.fastq $trim/temp.fastq &>> $trim/$s.vsearch.trim.fw.rv.log
    rm $trim/temp.fastq
	
	echo
    echo Calculate quality statistics
    $VSEARCH --threads $THREADS \
             --fastq_eestats $trim/$s.trimmed.fastq \
             --output $trim/$s.stats

    echo
    echo Quality filtering

    $VSEARCH --threads $THREADS \
             --fastq_filter $trim/$s.trimmed.fastq  \
             --fastq_maxee 1 \
             --fastq_minlen 150 \
             --fastq_maxlen 250 \
             --fastq_maxns 0 \
             --fastqout $clean/$s.filtered.fastq &>> $clean/$s.vsearch.quality.log 

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --threads $THREADS \
        --derep_fulllength $clean/$s.filtered.fastq \
        --strand plus \
        --output $derep/$s.derep.fasta \
        --sizeout \
        --uc $derep/$s.derep.uc \
        --relabel $s. \
        --relabel_keep \
        --fasta_width 0 &>> $derep/$s.vsearch.derep.log
done

#Concatenate software outputs of merging
cd $merge
tail -vn +1 *vsearch.merge.log > All.stat.merge.log
grep 'Merged' -A 1  *vsearch.merge.log > All.stat.merge.count.log

#Concatenate software outputs of trimming
cd $trim
tail -vn +1 *.vsearch.trim.fw.log > All.stat.trim.fw.log
tail -vn +1 *.vsearch.trim.fw.rv.log > All.stat.trim.fw.rv.log

#Concatenate software outputs of quality filter
cd $clean
tail -vn +1 *vsearch.quality.log > All.stat.quality.log

#Concatenate software outputs of dereplication
cd $derep
tail -vn +1 *vsearch.derep.log > All.stat.derep.log
grep -c '^>' *derep.fasta > uniq.seq.each.spl.txt


echo
echo ==============================================
echo Fwh and BFR samples independently processed
echo For each - concatenate all samples in one file
echo ==============================================

echo
echo Merge all samples

rm -f $derep/*all.derep.fasta $chim/*all.nonchimeras.derep.fasta

cat $derep/*Fwh*.derep.fasta > $derep/fwh.all.fasta
grep '>' < $derep/fwh.all.fasta | cut -d. -f1 | sort -u

cat $derep/*bfr*.derep.fasta > $derep/bfr.all.fasta
grep '^>' < $derep/bfr.all.fasta | cut -d. -f1 | sort -u

for i in $(ls | grep -v "Fwh.derep.fasta"); do 
echo "$i" $(grep -c "Fwh" "$i"); done


target='fwh bfr'

for t in $target; do
echo
echo $t Dereplicate across samples and remove singletons 

$VSEARCH --threads $THREADS \
    --derep_fulllength $derep/$t.all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc $derep/$t.all.derep.uc \
    --output $derep/$t.all.derep.fasta

echo $t Unique non-singleton sequences: $(grep -c "^>" $derep/$t.all.derep.fasta)

echo
echo $t Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
    --cluster_size $derep/$t.all.derep.fasta \
    --id 0.98 \
    --iddef 4 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc $preclust/$t.all.preclustered.uc \
    --centroids $preclust/$t.all.preclustered.fasta

echo $t Unique sequences after preclustering: $(grep -c "^>" $preclust/$t.all.preclustered.fasta)

echo
echo De novo chimera detection
$VSEARCH --threads $THREADS \
    --uchime_denovo $preclust/$t.all.preclustered.fasta \
    --abskew 10 \
    --sizein \
    --sizeout \
    --xsize \
    --fasta_score \
    --fasta_width 0 \
    --uchimealns $chim/$t.all.denovo.algns.fasta \
    --chimeras $chim/$t.all.denovo.chimeras.fasta \
    --nonchimeras $chim/$t.all.denovo.nonchimeras.fasta \
    --uchimeout $chim/$t.all.denovo.chimeras.table 

echo Unique sequences after denovo chimera detection: $(grep -c "^>" $chim/$t.all.denovo.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated
$PERL $source/map.pl $derep/$t.all.derep.fasta $preclust/$t.all.preclustered.uc $chim/$t.all.denovo.nonchimeras.fasta > $clust/$t.all.nonchimeras.derep.fasta
echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" $clust/$t.all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample
$PERL $source/map.pl $derep/$t.all.fasta $derep/$t.all.derep.uc $clust/$t.all.nonchimeras.derep.fasta > $clust/$t.all.nonchimeras.each.smple.fasta
echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" $clust/$t.all.nonchimeras.each.smple.fasta)

echo
echo $t Cluster at 97% and relabel with OTU_n, generate OTU table

$VSEARCH --threads $THREADS \
    --cluster_size $clust/$t.all.nonchimeras.each.smple.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --sizeorder \
    --xsize \
    --fasta_width 0 \
    --uc $clust/$t.all.clustered.uc \
    --relabel OTU_ \
    --centroids $clust/$t.all.otus.fasta \
    --otutabout $clust/$t.all.otutab.txt

echo
echo Number of OTUs: $(grep -c "^>" $clust/$t.all.otus.fasta)

done

grep '^>' $clust/bfr.all.otus.fasta > $clust/bfr.seqname.abd
grep '^>' $clust/fwh.all.otus.fasta > $clust/fwh.seqname.abd

sed -i 's/;.*;//g' $clust/bfr.all.otus.fasta
sed -i 's/;.*;//g' $clust/fwh.all.otus.fasta

Rscript 0.Rscript.amplicon.vsearch.bfr.MC2.MC1.subtab.R bfr.all.otutab.txt $clust

#OTUs in common for the 10 species mock communities (hereafter MC1) and the 52 taxa mock communities (hereafter MC2)
sort $clust/bfr.MC1.seqnames > $clust/bfr.MC1.seqnames.sorted
sort $clust/bfr.MC2.seqnames > $clust/bfr.MC2.seqnames.sorted
comm -12 $clust/bfr.MC2.seqnames.sorted $clust/bfr.MC1.seqnames.sorted > $clust/comm.otu.MC1.MC2
wc -l $clust/comm.otu.MC1.MC2 

sed 's/"//g' $clust/bfr.MC1.seqnames > $clust/bfr.MC1.seqnames2
sed 's/"//g' $clust/bfr.MC2.seqnames > $clust/bfr.MC2.seqnames2

select_sequences.pl $clust/bfr.all.otus.fasta $clust/bfr.MC1.seqnames2 $clust/bfr.MC1.otus.fasta 
select_sequences.pl $clust/bfr.all.otus.fasta $clust/bfr.MC2.seqnames2 $clust/bfr.MC2.otus.fasta
#perl script available here: 

###########################################################################
#II.ASSIGNATION - VSEARCH PIPELINE : ALIGNMENT + ABUNDANCE QUANTIFICATION
###########################################################################

############################################################################################
###### Reference databases

####COI sequences barcoded for the MC1 (available on NCBI Accession Numbers = MK584516:MK584524)
mkdir ./Database/
mkdir ./Database/DB_COI_10sps/
#Put sequences of the database in this folder
cat ./Database/DB_COI_10sps/*.fasta > COI_database_10sps.fas
MC1_db_COI=./Database/DB_COI_10sps/COI_database_10sps.fas
makeblastdb -in $MC1_db_COI -dbtype nucl

####COI sequences barcoded for the MC2 
mkdir ./Database/DB_MC2_haplotypes
#Put sequences of the database in this folder
MC2_db_COI=./Database/DB_MC2_haplotypes
MC2_db_COI_A=./Database/DB_MC2_haplotypes/TierMix-A.COI.ref.fasta
MC2_db_COI_B=./Database/DB_MC2_haplotypes/TierMix-B.COI.ref.fasta
MC2_db_COI_C=./Database/DB_MC2_haplotypes/TierMix-C.COI.ref.fasta
MC2_db_COI_D=./Database/DB_MC2_haplotypes/TierMix-D.COI.ref.fasta
MC2_db_COI_E=./Database/DB_MC2_haplotypes/TierMix-E.COI.ref.fasta
MC2_db_COI_F=./Database/DB_MC2_haplotypes/TierMix-F.COI.ref.fasta
MC2_db_COI_G=./Database/DB_MC2_haplotypes/TierMix-G.COI.ref.fasta
MC2_db_COI_H=./Database/DB_MC2_haplotypes/TierMix-H.COI.ref.fasta
MC2_db_COI_i=./Database/DB_MC2_haplotypes/TierMix-i.COI.ref.fasta
MC2_db_COI_j=./Database/DB_MC2_haplotypes/TierMix-j.COI.ref.fasta
#Files made from sequences of Elbrecht and Leese, 2015 
#Available in github folder Resources 
for f in ./Database/DB_MC2_haplotypes/*COI.ref.fasta; do
	makeblastdb -in $f -dbtype nucl
done

############################################################################################
####Fwh1 samples - MC1 samples only ####
blastn -db ./$MC1_db_COI -task 'blastn' -query $clust/fwh.all.otus.fasta -outfmt 6 -evalue 1E-10 > $ass/fwh.all.otus.nofilter.blastn.out
awk '$3>=97' < $ass/fwh.all.otus.nofilter.blastn.out | awk '$4>=90' > $ass/fwh.all.otus.nofilter.blastn.id97.qc90.out

Rscript 0.Rscript.amplicon.vsearch.quantification.blastn.nofilter.R $clust/fwh.all.otutab.txt $ass/fwh.all.otus.nofilter.blastn.id97.qc90.out $ass &>> $ass/fwh.quanti.blastn.nofilter.log

####BF2/BR2 samples - MC1 samples ####
blastn -db .s/$MC1_db_COI -task 'blastn' -query $clust/bfr.MC1.otus.nofilter.fasta -outfmt 6 -evalue 1E-10 > $ass/bfr.MC1.otus.nofilter.blastn.out
awk '$3>=97' < $ass/bfr.MC1.otus.nofilter.blastn.out | awk '$4>=200' > $ass/bfr.MC1.otus.nofilter.blastn.id97.qc200.out

Rscript 0.Rscript.amplicon.vsearch.quantification.blastn.nofilter.R $clust/bfr.MC1.otutab.nofilter.txt $ass/bfr.MC1.otus.nofilter.blastn.id97.qc200.out $ass &>> $ass/bfr.MC1.quanti.blastn.nofilter.log

####BF2/BR2 samples - MC2 samples ####
#NB: one alignment table per reference data
MC2="A B C D E F G H i j"

for a in $MC2; do
echo
echo Processing TierMix-$a
echo
db0="MC2_db_COI_$a"
db=${!db0}
blastn -db $db -task 'blastn' -query $clust/bfr.MC2.otus.nofilter.fasta -outfmt 6 -evalue 1E-10 > $ass/MC2.$a.bfr.otus.nofilter.blastn.out
awk '$3>=97' < $ass/MC2.$a.bfr.otus.nofilter.blastn.out | awk '$4>=200' > $ass/MC2.$a.bfr.otus.nofilter.blastn.id97.qc200.out
done 

Rscript $source/0.Rscript.amplicon.vsearch.quantification.MC2.spl.blast.nofilter.R &>> $ass/MC2.bfr.quanti.nofilter.log


