# cb-RNA-seq
Scripts used in the analysis on the chromatin-bound RNA-seq (cb-RNA-seq) data. The related article is https://doi.org/10.1186/s13007-022-00967-y .

The calculation method of the splicing ratio of intron at 5'splice site (SS5) and 3' splice site (SS3) was same as the method mentioned in the Zhu et al., (2020, MP).

1. Preparatory work. Download the scripts and necessary files from GitHub repository: https://github.com/flzh628/cb-RNA-seq. TAIR10 genome sequence and annotation files were downloaded from the Gramene website (https://ensembl.gramene.org/Arabidopsis_thaliana/Info/Index).

2. Quality control of raw data. And the reads shorter than 36 bp were also removed.

<code>$ java -jar /path/to/Trimmomatic/trimmomatic-0.39.jar PE -phred33 -threads 64 Read1.fastq.gz Read2.fastq.gz Read1.paired.fq.gz Read1.unpaired.fq.gz Read2.paired.fq.gz Read2.unpaired.fq.gz ILLUMINACLIP:/path/to/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36</code>

3. Read mapping. The reads with more than one reported alignment were excluded. 

<code>$ hisat2 -p 64 --dta -x /path/to/TAIR10/index/TAIR10 -1 Read1.paired.fq.gz -2 Read2.paired.fq.gz -S Sample.sam</code><br/>
<code>$ less Sample.sam | grep '\bNH:i:1\b' | samtools view -T /path/to/TAIR10/TAIR10.fa -h - | samtools view -bS -@ 32 - > Sample.uniq.bam</code><br/>
<code>$ samtools sort -@ 32 -O bam -o Sample.sort.bam -T temp Sample.uniq.bam</code><br/>
<code>$ rm Sample.uniq.bam && mv Sample.sort.bam Sample.uniq.bam</code><br/>
<code>$ samtools index Sample.uniq.bam</code>

4. Quantification of gene expression.

<code>$ stringtie -G TAIR10.gtf -o Sample.uniq.gtf -p 64 -e -b ballgown_out_ Sample -A Sample.uniq.gene.abundance Sample.uniq.bam</code><br/>
<code>$ echo "LOC_name		Sample.uniq.bam" > Sample.uniq.gene.FPKM</code><br/>
<code>$ tail -n +2 Sample.uniq.gene.abundance | awk '{print $1"\t"$(NF-1)}' >> Sample.uniq.gene.FPKM</code>

5. The calculation of SS ratio. Genes with FPKM ≥ 1 were used for the calculation of 5’ SS or 3’ SS ratio. The structure of the longest transcript was used as representative gene structure. 

For one sample:

<code>$ python IntronSpliceRatio_final.py --inputs Sample.uniq.bam --intron TAIR10.LongRNA.intron --output Sample.intron.SS-ratio --ex_chr ChrC,ChrM --geneflt Sample.uniq.gene.FPKM</code>

For Multiple samples:

<code>$ python IntronSpliceRatio_final.py --inputs Sample1.uniq.bam,Sample2.uniq.bam --intron TAIR10.LongRNA.intron --output Sample1_Sample2.intron.SS-ratio --ex_chr ChrC,ChrM --geneflt Sample1_Sample2.uniq.gene.FPKM</code>

6. Downstream analyses. Downstream analyses included identifications of differentially cleaved introns, genes, and correlation analysis between SS ratio and intron number of gene. And the average SS ratio of all the introns represents the gene’s SS ratio.

<code>$ python intron_SS_ratio_diff.py --ss_res Sample1_Sample2.intron.SS-ratio --tot_Col 10000000(# total read number of Sample1) --tot_mut 11000000(# total read number of Sample2) --names Sample1,Sample2</code>

<code>$ python intronSS_2_GeneLevel.py --inSS Sample1_Sample2.intron.SS-ratio --prefix Sample1_Sample2</code><br/>
<code>$ python gene_SS_ratio_diff.py Sample1_Sample2.gene.SS-ratio --tot_Col 10000000(# total read number of Sample1) --tot_mut 11000000(# total read number of Sample2) --names Sample1,Sample2</code><br/>

<code>$ awk '{print $1"\t"$7}' Sample1_Sample2.gene.SS-ratio >Sample1.gene.SS5-ratio</code><br/>
<code>$ python ssRatio_intronNum.py --pcg TAIR10.protein_coding_genes --intronNum TAIR10.LongRNA.intron_num --ratio Sample1.gene.SS5-ratio --prefix Sample1 --mark SS5</code>
