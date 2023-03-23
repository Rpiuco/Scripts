

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq; print seq}' 16751776.fasta > concat.16751776.temp
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq; print seq}' 16751777.fasta > concat.16751777.temp
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq; print seq}' 21989093.fasta > concat.21989093.temp
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq; print seq}' 22313525.fasta > concat.22313525.temp

cat concat.16751776.temp concat.16751776.temp concat.16751776.temp concat.16751776.temp > concat.fasta

sort -k 1 concat.fasta | uniq > concat.fasta.unique

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq; print seq}' hsa_pirna.fa > hsa_pirna.fa.unique

cat concat.fasta.unique hsa_pirna.fa.unique | sort -k 1 | uniq -c | awk '{if($1==1){print$2}}' | awk 'BEGIN{num=1} {print ">seq"num"\n"$1; num++;}' > pirnaNoAlign.fa


scp -P 1717 pirnaNoAlign.fa rpiuco@www.bioinfo.mochsl.org.br:/home/projects2/local_projects/pirna/pirnaNoAlign

ln -s /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa 
ln -s /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa 


bwa aln -t 8 -n 1 /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa pirnaNoAlign.fa 2> bwa/pirnaNoAlign.hg38.n1.bwalog | bwa samse -n 9999999 /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa - pirnaNoAlign.fa | samtools view -bhS - -o bwa/pirnaNoAlign.hg38.n1.bam &
bwa aln -t 8 -n 2 /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa pirnaNoAlign.fa 2> bwa/pirnaNoAlign.hg38.n2.bwalog | bwa samse -n 9999999 /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa - pirnaNoAlign.fa | samtools view -bhS - -o bwa/pirnaNoAlign.hg38.n2.bam &
bwa aln -t 8 -n 1 /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa pirnaNoAlign.fa 2> bwa/pirnaNoAlign.hg19.n1.bwalog | bwa samse -n 9999999 /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa - pirnaNoAlign.fa | samtools view -bhS - -o bwa/pirnaNoAlign.hg19.n1.bam &
bwa aln -t 8 -n 2 /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa pirnaNoAlign.fa 2> bwa/pirnaNoAlign.hg19.n2.bwalog | bwa samse -n 9999999 /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa - pirnaNoAlign.fa | samtools view -bhS - -o bwa/pirnaNoAlign.hg19.n2.bam &

bwa aln -t 8 -n 0 /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa pirnaNoAlign.fa 2> bwa/pirnaNoAlign.hg19.n0.bwalog | bwa samse -n 9999999 /home/genomes/Homo_sapiens/hg19/toBWA.25chrs/hg19.fa - pirnaNoAlign.fa | samtools view -bhS - -o bwa/pirnaNoAlign.hg19.n0.bam &
bwa aln -t 8 -n 0 /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa pirnaNoAlign.fa 2> bwa/pirnaNoAlign.hg38.n0.bwalog | bwa samse -n 9999999 /home/genomes/Homo_sapiens/hg38/toBWA.25chrs/hg38.fa - pirnaNoAlign.fa | samtools view -bhS - -o bwa/pirnaNoAlign.hg38.n0.bam &


samtools flagstat bwa/pirnaNoAlign.hg19.n0.bam > bwa/flagstat.hg19.n0.txt
samtools flagstat bwa/pirnaNoAlign.hg19.n1.bam > bwa/flagstat.hg19.n1.txt
samtools flagstat bwa/pirnaNoAlign.hg19.n2.bam > bwa/flagstat.hg19.n2.txt
samtools flagstat bwa/pirnaNoAlign.hg38.n0.bam > bwa/flagstat.hg38.n0.txt
samtools flagstat bwa/pirnaNoAlign.hg38.n1.bam > bwa/flagstat.hg38.n1.txt
samtools flagstat bwa/pirnaNoAlign.hg38.n2.bam > bwa/flagstat.hg38.n2.txt


samtools view bwa/pirnaNoAlign.hg19.n0.bam | awk '{if($2!=4){print$1}}' > alignment/alignNoMissMatches.hg19.txt

samtools view bwa/pirnaNoAlign.hg19.n1.bam | awk '{if($13!="NM:i:1" && $13!="NM:i:0"){print$0}}' > alignment/notAligned.hg19.n1.txt
samtools view bwa/pirnaNoAlign.hg19.n1.bam | awk '{if($13=="NM:i:1"){print$0}}' > alignment/aligned.hg19.n1.txt
samtools view bwa/pirnaNoAlign.hg19.n1.bam | awk '{if($13=="NM:i:1"){print$19}}' > alignment/cigar.hg19.n1.txt
samtools view bwa/pirnaNoAlign.hg19.n2.bam | awk '{if($13!="NM:i:2" && $13!="NM:i:0" && $13!="NM:i:1"){print$0}}' > alignment/notAligned.hg19.n2.txt
samtools view bwa/pirnaNoAlign.hg19.n2.bam | awk '{if($13=="NM:i:2"){print$0}}' > alignment/aligned.hg19.n2.txt
samtools view bwa/pirnaNoAlign.hg19.n2.bam | awk '{if($13=="NM:i:2"){print$19}}' > alignment/cigar.hg19.n2.txt

samtools view bwa/pirnaNoAlign.hg38.n1.bam | awk '{if($13!="NM:i:1" && $13!="NM:i:0"){print$0}}' > alignment/notAligned.hg38.n1.txt
samtools view bwa/pirnaNoAlign.hg38.n1.bam | awk '{if($13=="NM:i:1"){print$0}}' > alignment/aligned.hg38.n1.txt
samtools view bwa/pirnaNoAlign.hg38.n1.bam | awk '{if($13=="NM:i:1"){print$19}}' > alignment/cigar.hg38.n1.txt
samtools view bwa/pirnaNoAlign.hg38.n2.bam | awk '{if($13!="NM:i:2" && $13!="NM:i:0" && $13!="NM:i:1"){print$0}}' > alignment/notAligned.hg38.n2.txt
samtools view bwa/pirnaNoAlign.hg38.n2.bam | awk '{if($13=="NM:i:2"){print$0}}' > alignment/aligned.hg38.n2.txt
samtools view bwa/pirnaNoAlign.hg38.n2.bam | awk '{if($13=="NM:i:2"){print$19}}' > alignment/cigar.hg38.n2.txt



cat > testeAbundanceProtein/concat.kallisto.norxhigh.header.proteinCoding
	Sample_140	Sample_237	Sample_119	Sample_138	Sample_142	Sample_130	Sample_169	Sample_231	Sample_123	Sample_167	Sample_229	Sample_110	Sample_104	Sample_125	Sample_248	Sample_102	Sample_113	Sample_242	Sample_128	Sample_225	Sample_114	Sample_163	Sample_116	Sample_117

cat > testeAbundanceProtein/concat.kallisto.norxlow.header.proteinCoding
	Sample_140	Sample_237	Sample_119	Sample_138	Sample_105	Sample_112	Sample_249	Sample_247	Sample_124	Sample_145	Sample_151	Sample_bn27	Sample_133	Sample_147	Sample_120

cat > testeAbundanceProtein/concat.kallisto.lowxhigh.header.proteinCoding
	Sample_105	Sample_112	Sample_249	Sample_247	Sample_124	Sample_145	Sample_151	Sample_bn27	Sample_133	Sample_147	Sample_120	Sample_142	Sample_130	Sample_169	Sample_231	Sample_123	Sample_167	Sample_229	Sample_110	Sample_104	Sample_125	Sample_248	Sample_102	Sample_113	Sample_242	Sample_128	Sample_225	Sample_114	Sample_163	Sample_116	Sample_117




