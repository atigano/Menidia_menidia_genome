# Menidia_menidia_genome
Scripts accompanying manuscript on sequence and structural variation in the Atlantic silverside (Menidia menidia)

## Genome assembly with 10X data

```
supernova run --id supernova_run_output_270mreads --fastqs /workdir/chromium/H7NF5BCX2_10388115_supernova_mkfastq_output_24Jan2018/outs/fastq_path/ --maxreads 270000000
supernova mkoutput --asmdir supernova_run_output_270mreads/outs/assembly --outprefix meme_270mreads_1_pseudohap --style pseudohap
```
## Minor fixes to Dovetail Genome
a. Filter scaffolds smaller than 1000 bp
```
python /workdir/Menidia_Genome/newassembly_at/scripts_assembly/filter_scaffolds.py 500 menidia_menidia_30Jul2018_eLAsD.fasta 
python /workdir/Menidia_Genome/newassembly_at/scripts_assembly/filter_scaffolds.py 1000 menidia_menidia_30Jul2018_eLAsD.fasta 
```
b. Sort by length (from longest to shortest)
```
seqkit sort -l -r menidia_menidia_30Jul2018_eLAsD.filter1000.fasta > menidia_menidia_30Jul2018_eLAsD_sorted.filter1000.fasta
```
c. Rename scaffolds
```
awk '/^>/{print ">" ++i; next}{print}' < menidia_menidia_30Jul2018_eLAsD_sorted.filter1000.fasta > menidia_menidia_tidy.filter1000.fasta
```
d. Reduce genome to longest 27 scaffolds
```
cd /workdir/anna/dovetail
samtools faidx menidia_menidia_tidy.filter1000.fasta {1..27} > menidia_menidia_1-27.fasta
```


## BUSCO
Assessment of completeness of each of the intermediate and final assemblies
```
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_1-27.fasta -o busco_dovetail_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome 
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_tidy.filter1000.fasta -o busco_dovetail_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_26Jul2018_hZUSG.fasta -o busco_dovetail_chicago_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_30Jul2018_eLAsD.fasta -o busco_dovetail_hic_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/chromium/meme270m500.fasta -o busco_10x_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
```
## K-mer analysis on 10x data from Georgia
a. Pre-process 10x data 
```
longranger basic \
    --id=meme_10xreads \
    --fastqs=/workdir/anna/chromium/H7NF5BCX2_10388115_supernova_mkfastq_output_24Jan2018/outs/fastq_path/Project_9189_10561/Sample_65612_H7NF5BCX2_MEME1/decompressed \
    --localcores 24 
```
b. Deinterleave fastq data
```
/workdir/anna/scripts/deinterleave_fastq.sh < barcoded.fastq barcoded1.fastq barcoded2.fastq
```
c. Cut all reads to the same length
```
cd /workdir/anna/chromium/10x_reads_basic
cutadapt -l 129 -o barcoded_trimmed1.fastq -p barcoded_trimmed2.fastq /workdir/anna/chromium/10x_reads_basic/meme_10xreads/outs/barcoded1.fastq /workdir/anna/chromium/10x_reads_basic/meme_10xreads/outs/barcoded2.fastq
```
d. Run Jellyfish
```
cd /workdir/anna/chromium/10x_reads_basic
cutadapt -l 129 -o barcoded_trimmed1.fastq -p barcoded_trimmed2.fastq /workdir/anna/chromium/10x_reads_basic/meme_10xreads/outs/barcoded1.fastq /workdir/anna/chromium/10x_reads_basic/meme_10xreads/outs/barcoded2.fastq
```
## K-mer analysis on shotgun data from Connecticut
Run Jellyfish
```
jellyfish count -t 12 -C -m 25 -s 5G -o meme_25mer menidia_R*.fastq
jellyfish histo -o meme_25mer.histo meme_25mer
```

## Synteny with medaka
a. Create database for medaka genome
```
lastdb -P12 -uMAM8 -cR11 medakadb medaka_assembly.fasta
```
b. Align Atlantic silverside scaffolds to medaka reference genome
```
lastal -P12 -m100 -E0.05 medakadb /workdir/anna/dovetail/menidia_menidia_tidy.filter1000.fasta | last-split -m1 > last_meme2medaka.maf
```

## Variant calling for direct estimates of heterozygosity
Align 10X data from Georgia to Georgia reference genome
```
bwa mem -t 12 menidia_menidia_1-27.fasta \
/workdir/anna/chromium/10x_reads_basic/barcoded_trimmed1.fastq \
/workdir/anna/chromium/10x_reads_basic/barcoded_trimmed2.fastq \
| samblaster --removeDups | samtools view -h -b -F 4 -@6 - | samtools sort -o 10x2dovetail_1-27.bam
```
Align shotgun data from Connecticut to Georgia reference genome
```
bwa mem -t 12 menidia_menidia_1-27.fasta \
/workdir/Menidia_Genome/newassembly_at/menidia_trimmed_q30_R1.fq \
/workdir/Menidia_Genome/newassembly_at/menidia_trimmed_q30_R2.fq \
| samblaster --removeDups | samtools view -h -b -F 4 -@6 - | samtools sort -o shotgun2dovetail_1-27.bam
```
Estimate genome coverage for each dataset for downstream variant filtering
```
genomeCoverageBed -ibam 10x2dovetail_1-27.bam -max 300 > coverage_10x2dovetail.txt ###13865250 had no coverage, should be removed from heterozygosity calculations
genomeCoverageBed -ibam shotgun2dovetail_1-27.bam -max 300 > coverage_shotgun2dovetail.txt
```
Call variants from 10X data from Georgia
```
bcftools mpileup -Ou --threads 12 -f menidia_menidia_1-27.fasta 10x2dovetail_1-27.bam | \
    bcftools call -Ou -mv --threads 12| \
    bcftools filter --threads 12 -s LowQual -e '%QUAL<20 || DP>190' > 10x2dovetail_1-27.vcf
```
Call variants from shotgun data from Connecticut
```
bcftools mpileup -Ou --threads 12 -f menidia_menidia_1-27.fasta shotgun2dovetail_1-27.bam | \
    bcftools call -Ou -mv --threads 12| \
    bcftools filter --threads 12 -s LowQual -e '%QUAL<20 || DP>148' > shotgun2dovetail_1-27.vcf 
```
Second round of filtering for low qual, DP < 20, and homozygous sites
```
bcftools filter --threads 12 -e'DP<20 || FILTER="LowQual"' shotgun2dovetail_1-27.vcf > shotgun2dovetail_1-27_filt.vcf 
###to select only SNPs
bcftools view -v snps shotgun2dovetail_1-27_filt.vcf > shotgun2dovetail_1-27_filt_snp.vcf
###to retain only heterozygous sites
bcftools view -g het shotgun2dovetail_1-27_filt_snp.vcf > shotgun2dovetail_1-27_filt_snp_het.vcf
```
To calculate heterozygosity, we divided the number of heterozygous sites surviving the filters by the total number of bases that had depth of coverage between 20 and twice the mode the sequencing depth for each library
```
shotgun Connecticut --> 422616334 --> 6160030/422616334 = 0.01458 --> 1.46%
10x Georgia --> 439762099 --> 5823111/439762099 = 0.01324 --> 1.32%
```
Estimate heterozygity in coding regions only
```
Heterozygosity coding vs. non-coding

grep "CDS" mme_annotation_draftgenome_clean.noseq.gff > mme_annotation_CDS.gff 
sort -k1,1n -k4,4n mme_annotation_CDS.gff
awk '{if($1 < 28){print}}' mme_annotation_CDS.gff | sort -k1,1n -k4,4n > mme_annotation_CDS_1-27.gff
bedtools map -a mme_annotation_CDS_1-27.gff -b ../dovetail/10x2dovetail_1-27_filt_snp_het.vcf -o count > 10x2dovetail_1-27_filt_snp_het_CDS_1-27.txt
Number of variants in CDS is sum of last column in 10x2dovetail_1-27_filt_snp_het_CDS_1-27.txt
awk '{SUM+=$10}END{print SUM}' 10x2dovetail_1-27_filt_snp_het_CDS_1-27.txt ###236098 variant sites in CDS

#Heterozygosity 0.679% Georgia


bedtools map -a mme_annotation_CDS_1-27.gff -b ../dovetail/shotgun2dovetail_1-27_filt_snp_het.vcf -o count > shotgun2dovetail_1-27_filt_snp_het_CDS_1-27.txt
Number of variants in CDS is sum of last column in 10x2dovetail_1-27_filt_snp_het_CDS_1-27.txt
awk '{SUM+=$10}END{print SUM}' shotgun2dovetail_1-27_filt_snp_het_CDS_1-27.txt ###243376 variant sites in CDS

#Heterozygosity 0.700% Connecticut
```
## Structural variant calling
Call structural variants with Delly2 from shotgun data from Connecticut
```
samtools index shotgun2dovetail_1-27.bam
delly call -g menidia_menidia_1-27.fasta -o SV_shotgun.bcf shotgun2dovetail_1-27.bam
bcftools view SV_shotgun.bcf > SV_shotgun.vcf
```
Filter low quality calls
```
grep -v "LowQual" SV_shotgun.vcf > SV_shotgun_hiq.vcf
```
Further filtering was done in Excel and combined with SVs merging, for each SV type separately, to collaps redundant calls
```
bedtools merge -i ins_20-100.bed > ins_20-100_reduced.bed
bedtools merge -i del_20-100.bed > del_20-100_reduced.bed
bedtools merge -i dup_20-100.bed > dup_20-100_reduced.bed
```

## Manhattan plot for levels of heterozygosity in each of the two genome
Calculate number of variants in each window
```
cd /workdir/anna/dovetail
### Make 50-kb sliding windows
bedtools makewindows -g menidia_menidia_tidy.filter1000.fasta.fai -w 50000 | awk '$3 ~ /0000$/' | sed 's/ /\t/g' > genome_windows_50k_4GC.bed
### Count heterozygous sites for each window
bedtools map -a genome_windows_50k_4GC.bed -b 10x2dovetail_1-27_filt_snp_het.vcf -o count > 10x_variants_windows.txt
bedtools map -a genome_windows_50k_4GC.bed -b shotgun2dovetail_1-27_filt_snp_het.vcf -o count > shotgun_variants_windows.txt
```


