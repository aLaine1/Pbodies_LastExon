### Select non-overlaping genes
#grep_tmp ==> transcript names of significant last exon
grep -f grep_tmp last_exon_final.bed > signif_last_exon.bed
bedtools intersect -a signif_last_exon.bed -b gencode.bed -wb > intersect.signif
cat intersect.signif |awk '{print $4"\t"$10}'|sort -k1,2|uniq > uni_intersect
cat uni_intersect |awk '{print $1}'|uniq -c |awk '{if ($1 == 1) print $2}' > non_overlap_exons

###Extract last exon fasta sequence from hg38
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed 's/>/>chr/g' > ref.fa
rm Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
grep -f non_overlap_exons last_exon_final.bed > signif_last_exon.bed

bedtools getfasta -s -fi ref.fa -bed signif_last_exon.bed -fo signif_last_exon.fa
