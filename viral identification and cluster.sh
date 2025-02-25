###1 viral contigs idendification
##1.1 Virsorter2
virsorter run \
-i S1_assembled_contigs.fa \
-w S1.vs2 \
-d vs2_db \
--min-length 5000 \
--keep-original-seq \
-j 60 all

#Screening viral contigs that meet the criteria
cd S1.vs2
awk '($6>=5000)&&($4>=0.9){print $1}' final-viral-score.tsv > S1.vs2_good_id.txt
awk '($6>=5000)&&(0.9>$4)&&(0.7<=$4)&&($7>=1){print $1}' final-viral-score.tsv >> S1.vs2_good_id.txt
sed -i '1d' S1.vs2_good_id.txt
seqkit grep --pattern-file S1.vs2_good_id.txt final-viral-combined.fa -o S1.vs2_good.fa



##1.2 DeepVirFinder
python dvf.py \
-i S1_assembled_contigs.fa \
-o S1.dvf \
-c 30 \
-l 5000

#Screening viral contigs that meet the criteria
cd S1.dvf
awk -F'\t' '($3>=0.9)&&($4<=0.01){print $1}' S1_assembled_contigs.fa_gt5000bp_dvfpred.txt > S1.dvf_good_id.txt
seqkit grep -n --pattern-file S1.dvf_good_id.txt S1_assembled_contigs.fa -o S1.dvf_good.fa

##1.3 VIBRANT
#Extract the sequence required for running VIBRANT
#sequences with scores ≥ 0.5 in Virsorter2 and scores ≥ 0.7 and p ≤ 0.05 in DeepVirFinder
awk -F'\t' '($3<0.9)&&($3>=0.7)&&($4<=0.05){print $1}' S1_assembled_contigs.fa_gt5000bp_dvfpred.txt >> vibrant_id_1.txt
awk -F'\t' '($3>=0.9)&&($4>0.01)&&($4<=0.05){print $1}' S1_assembled_contigs.fa_gt5000bp_dvfpred.txt >> vibrant_id_1.txt
awk '($6>=5000)&&($4<0.9)&&($4>=0.7)&&($7<1){print $1}' final-viral-score.tsv >> vibrant_id_2.txt
awk '($6>=5000)&&($4<0.7)&&($4>=0.5){print $1}' final-viral-score.tsv >> vibrant_id_2.txt
cut -d'|' -f 1 vibrant_id_2.txt >> vibrant_id_3.txt
cat vibrant_id_1.txt vibrant_id_3.txt | sort |uniq > S1.vibrant_id.txt
seqkit grep -n --pattern-file S1.vibrant_id.txt S1_assembled_contigs.fa -o S1.vibrant.fa;

#run vibrant
VIBRANT_run.py -i S1.vibrant.fa -f nucl -t 25


##1.4 Merge all the viral sequences and De-duplication
cat S1.vs2_good.fa S1.dvf_good.fa S1.vibrant.phages_combined.fna >> S1.vir_merged.fa
seqkit rmdup -s S1.vir_merged.fa -o S1.vir_rep.fa



###2 Clustering to generate vOTUs
##2.1 build an index
makeblastdb \
 -dbtype nucl \
 -in all_viral_contigs.fasta \
 -input_type fasta \
 -parse_seqids \
 -out blastdb

##2.2 blast
blastn \
  -query all_viral_contigs.fasta \
  -db blastdb \
  -out blast.tsv \
  -outfmt '6 std qlen slen' \
  -perc_identity 90 \
  -num_threads 50
 
##2.3 cluster using a python script (https://github.com/snayfach/MGV).
python blastani.py -i blast.tsv -o ani.tsv
python cluster.py \
--fna all_viral_contigs.fasta \
--ani ani.tsv \
--out clusters.tsv \
--min_ani 95 \
--min_qcov 0 \
--min_tcov 85

###3 Quality assessment using CheckV
checkv end_to_end vOTU.fa checkv -d db/checkv-db-v1.5 -t 50 


###4 Abundance calculation
##Calculate the relative abundance of metagenome-derived vOTUs and virome-derived vOTUs separately.
coverm contig \
--reference vOTU.fa \
--interleaved hq_reads.fastq \
--mapper bwa-mem \
--min-read-percent-identity 95 \
--min-read-aligned-percent 90 \
--methods trimmed_mean \
--bam-file-cache-directory out_bam \
--output-file coverage.tsv \
-t 50
