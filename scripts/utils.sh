# Get Mean Coverage
mosdepth --by last_exon_final.bed mosdepth/sample1 Sample1.bam

#Merge samples for mean coverage
paste <(zcat mosdepth/sample1.regions.bed.gz |awk '{print $4"\t"$5}') <(zcat mosdepth/sample2.regions.bed.gz |awk '{print $5}') <(zcat mosdepth/sample3.regions.bed.gz |awk '{print $5}') <(zcat mosdepth/sample4.regions.bed.gz |awk '{print $5}') <(zcat mosdepth/sample5.regions.bed.gz |awk '{print $5}') <(zcat mosdepth/sample6.regions.bed.gz |awk '{print $5}') > lastExon_MeanCov_All

###First filter // MeanCov > 10 in one condition
cat lastExon_MeanCov_All |awk '{if (($2+$3+$4)/3 >= 10 || ($5+$6+$7)/3 >= 10) print $1}' > lastExon_MeanSup10
grep -f lastExon_MeanSup10 last_exon_final.bed > lastExon_MeanSup10Final.bed


### Select non-overlaping genes

bedtools intersect -a lastExon_MeanSup10Final.bed -b gencode.bed -wb > intersect.signif
cat intersect.signif |awk '{print $4"\t"$10}'|sort -k1,2|uniq > uni_intersect
cat uni_intersect |awk '{print $1}'|uniq -c |awk '{if ($1 == 1) print $2}' > non_overlap_exons

###Filter last Exon with at least 2 PolyA sites
bedtools intersect -a data/lastExon_MeanSup10Final.bed -b data/PolyA_Site.bed -s -wa|sort -V |uniq -c |awk '{if ($1 > 1) print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > data/lastExon_MeanSup10Final_PolyAplus.bed


###Extract last exon fasta sequence from hg38
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed 's/>/>chr/g' > ref.fa
rm Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
grep -f non_overlap_exons last_exon_final.bed > signif_last_exon.bed

bedtools getfasta -s -fi ref.fa -bed signif_last_exon.bed -fo signif_last_exon.fa
