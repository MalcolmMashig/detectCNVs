
# SECNVs Simulator ----------

# python SECNVs/SECNVs.py -G <CHROMOSOME-1-FASTA> -T inputs/chr1-exons.bed -e_chr 100 -o_chr 10 -o test_human_nn -rn human -ssr -sb -f 1000 -ol 100 -tf 50 -clr 100 -sc -eN all -pr -n 20 -q 33 -s_r 0.001 -i_r 0.00001 -nr 1000000 -picard <ABSOLUTE-PATH-TO-detectCNVs>/picard -GATK <ABSOLUTE-PATH-TO-detectCNVs>/gatk

# OR...

# python SECNVs/SECNVs.py -o test-output -G <CHROMOSOME-1-FASTA> -T inputs/target-regions-chr1-long.bed -e_chr 10 -sc -sb -ssr

# CNV-Sim Simulator --------

# python CNV-Sim/cnv-sim.py -o test-output -n 1000000 exome <CHROMOSOME-1-FASTA> inputs/target-regions-chr1-long.bed

# NOTE: adjust paths if SECNVs was utilized instead of CNV-Sim

# make SAMs
# ./CNV-Sim/minimap2/minimap2 -t 8 -a -x sr <CHROMOSOME-1-FASTA> test-output/cnv_1.fastq test-output/cnv_2.fastq -o test-output/chr1-long.sam
# ./CNV-Sim/minimap2/minimap2 -t 8 -a -x sr <CHROMOSOME-1-FASTA> test-output/control_1.fastq test-output/control_2.fastq -o test-output/chr1-control-long.sam

# make BAMs
# ./samtools-1.14/samtools fixmate -O bam,level=1 -m test-output/chr1-long.sam test-output/chr1-long.bam
# ./samtools-1.14/samtools fixmate -O bam,level=1 -m test-output/chr1-control-long.sam test-output/chr1-control-long.bam

# sort BAMs
# ./samtools-1.14/samtools sort -l 1 -@8 -o test-output/chr1-long-srt.bam -T /tmp/example_prefix test-output/chr1-long.bam
# ./samtools-1.14/samtools sort -l 1 -@8 -o test-output/chr1-control-long-srt.bam -T /tmp/example_prefix test-output/chr1-control-long.bam

# mark dups in BAMs
# ./samtools-1.14/samtools markdup -O bam,level=1 test-output/chr1-long-srt.bam test-output/chr1-long-md.bam
# ./samtools-1.14/samtools markdup -O bam,level=1 test-output/chr1-control-long-srt.bam test-output/chr1-control-long-md.bam

# produce final BAMs
# ./samtools-1.14/samtools view -@8 test-output/chr1-long-md.bam -o test-output/chr1-long-final.bam
# ./samtools-1.14/samtools view -@8 test-output/chr1-control-long-md.bam -o test-output/chr1-control-long-final.bam

# index BAMs
# ./samtools-1.14/samtools index test-output/chr1-long-final.bam
# ./samtools-1.14/samtools index test-output/chr1-control-long-final.bam

# generate BAM counts per region (needed to detect CNVs)
# library(ExomeDepth)
# bamCounts <- getBamCounts(
#   bed.file = "test_human_nn/human1.cnv.overlap_target.bed",
#   bam.files = c(
#     "test_human_nn/control.final.bam",
#     "test_human_nn/human1.final.bam"
#   ),
#   referenceFasta = "test_human_nn/control.fa",
#   read.width = 100
# )
