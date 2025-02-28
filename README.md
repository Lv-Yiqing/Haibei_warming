# About

**This is the main script used for viral contigs extraction and analysis in the following study (We take sample S1 as an example):**

Enhanced viral activity and opposite diversity responses in surface and deep soils under experimental warming. *Unpublished.* (2025) Yiqing Lv, Xiaomin Fan, Xiangyu Fan, Wenkuan Qin, Xudong Wang, Xin Sun, Biao Zhu, Linwei Wu


# 1 metagenomic/viromic sequencing process
## 1.1 quality control
```
bbduk.sh \
in1=S1.R1.fq.gz \
in2=S1.R2.fq.gz \
out=S1.interleave.clean.fastq \
ref='PATH/to/adapters.fa,'PATH/to/phix_adapters.fa.gz \
ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=20 maq=20 -t=4 tbo
```
## 1.2 fastqc
```
fastqc  S1.interleave.clean.fastq
```
## 1.3 Assembly Using Megahit
```
megahit --12 S1.interleave.clean.fastq \
-t 20 \
--k-list 31,51,71,91,111,131 \
--kmin-1pass \
--min-contig-len 500 \
--continue \
--out-prefix S1 \
-o S1
```

# 2 viral contigs idendification
## 2.1 Virsorter2
```
virsorter run \
-i S1.contigs.fa \
-w S1.vs2 \
-d vs2_db \
--min-length 5000 \
--keep-original-seq \
-j 60 all
```

### Screening viral contigs that meet the criteria
```
cd S1.vs2
awk '($6>=5000)&&($4>=0.9){print $1}' final-viral-score.tsv > S1.vs2_good_id.txt
awk '($6>=5000)&&(0.9>$4)&&(0.7<=$4)&&($7>=1){print $1}' final-viral-score.tsv >> S1.vs2_good_id.txt
sed -i '1d' S1.vs2_good_id.txt
seqkit grep --pattern-file S1.vs2_good_id.txt final-viral-combined.fa -o S1.vs2_good.fa
```

## 2.2 DeepVirFinder
```
python dvf.py \
-i S1.contigs.fa \
-o S1.dvf \
-c 30 \
-l 5000
```
### Screening viral contigs that meet the criteria
```
cd S1.dvf
awk -F'\t' '($3>=0.9)&&($4<=0.01){print $1}' S1_assembled_contigs.fa_gt5000bp_dvfpred.txt > S1.dvf_good_id.txt
seqkit grep -n --pattern-file S1.dvf_good_id.txt S1.contigs.fa -o S1.dvf_good.fa
```

## 2.3 VIBRANT
### Extract the sequence required for running VIBRANT(sequences with scores ≥ 0.5 in Virsorter2 and scores ≥ 0.7 and p ≤ 0.05 in DeepVirFinder)
```
awk -F'\t' '($3<0.9)&&($3>=0.7)&&($4<=0.05){print $1}' S1.contigs.fa_gt5000bp_dvfpred.txt >> vibrant_id_1.txt
awk -F'\t' '($3>=0.9)&&($4>0.01)&&($4<=0.05){print $1}' S1.contigs.fa_gt5000bp_dvfpred.txt >> vibrant_id_1.txt
awk '($6>=5000)&&($4<0.9)&&($4>=0.7)&&($7<1){print $1}' final-viral-score.tsv >> vibrant_id_2.txt
awk '($6>=5000)&&($4<0.7)&&($4>=0.5){print $1}' final-viral-score.tsv >> vibrant_id_2.txt
cut -d'|' -f 1 vibrant_id_2.txt >> vibrant_id_3.txt
cat vibrant_id_1.txt vibrant_id_3.txt | sort |uniq > S1.vibrant_id.txt
seqkit grep -n --pattern-file S1.vibrant_id.txt S1.contigs.fa -o S1.vibrant.fa;
```

### run vibrant
```
VIBRANT_run.py -i S1.vibrant.fa -f nucl -t 25
```

## 2.4 Merge all the viral sequences and De-duplication
```
cat S1.vs2_good.fa S1.dvf_good.fa S1.vibrant.phages_combined.fna >> S1.vir_merged.fa
seqkit rmdup -s S1.vir_merged.fa -o S1.vir_rep.fa
```

# 3 Clustering to generate vOTUs
## 3.1 build an index
```
makeblastdb \
 -dbtype nucl \
 -in all_viral_contigs.fasta \
 -input_type fasta \
 -parse_seqids \
 -out blastdb
```

## 3.2 blastn
```
blastn \
  -query all_viral_contigs.fasta \
  -db blastdb \
  -out blast.tsv \
  -outfmt '6 std qlen slen' \
  -perc_identity 90 \
  -num_threads 50
```
 
## 3.3 cluster using a python script (https://github.com/snayfach/MGV).
```
python blastani.py -i blast.tsv -o ani.tsv
python cluster.py \
--fna all_viral_contigs.fasta \
--ani ani.tsv \
--out clusters.tsv \
--min_ani 95 \
--min_qcov 0 \
--min_tcov 85
```
#### Extract vOTU sequences based on clustering results.

# 4 Quality assessment using CheckV
```
checkv end_to_end vOTU.fa checkv -d db/checkv-db-v1.5 -t 50 
```

# 5 Abundance calculation
### Calculate the relative abundance of metagenome-derived vOTUs and virome-derived vOTUs separately.(We only show metagenome-derived vOTU as an example)
```coverm contig \
--reference metagenome_derived_vOTU.fa \
--interleaved path/to/genomic_hq_reads.fastq \
--mapper bwa-mem \
--min-read-percent-identity 95 \
--min-read-aligned-percent 90 \
--methods trimmed_mean \
--bam-file-cache-directory out_bam \
--output-file coverage.tsv \
-t 50
```
#### The following viral community diversity and composition analyses were performed in R4.3.1. Please refer to the corresponding R code for details.
