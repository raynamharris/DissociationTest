# Mapping with STAR

STAR '‘Spliced Transcripts Alignment to a Reference" is a faster alternative to tophat for splice-aware read alignment. Unlike the multi step mapping process used by tophat, STAR can align the non-contiguous sequences directly to the genome. The STAR algorithm consists of two major steps: seed searching step and clustering/stitching/scoring step. STAR is more memory intensive (30 gb of RAM required for human genome as compared to ~4 gb required by tophat), but it is significantly faster.


## Set Stampede environment  variables

For each new project or new batch of samples, we can reset these variables and then all the code will work below, with out having to recode all the project specific file names.

~~~ {.bash}
## set the enviornment variables 
RNAseqProject=BehavEphyRNAseq
RNAseqJob=JA17009
~~~ 

## Download the STAR reference GENOME index: Only do this once!

Navigate to the `refs` directory and download mouse genome and the genome annotation file from a place on the internet with some resources for STAR: http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/GENCODE/

This was gonna take a hour on the home node, which isn't great etiquette, so I created a job. 
~~~ {.bash}
cd $SCRATCH/$RNAseqProject/refs
mkdir STAR
cd STAR
# mouse genome
echo 'curl -O http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/GENCODE/GRCh38_Gencode25/GRCh38.primary_assembly.genome.fa' > 05_STAR_download.cmds
# mouse genome annotation file
echo 'curl -O http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/GENCODE/GRCh38_Gencode25/gencode.v25.chr_patch_hapl_scaff.annotation.gtf' >> 05_STAR_download.cmds
launcher_creator.py -t 24:00:00 -j 05_STAR_download.cmds -n 05_STAR_download -l 05_STAR_download.slurm -e rayna.harris@utexas.edu -q normal
sbatch 05_STAR_download.slurm
~~~

Now, that we have the genome and the .gtf file, we need to build a STAR-specific index. Let's run this as a job just incase.

~~~ {.bash}
echo 'STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v25.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 74 --genomeChrBinNbits 14' > 05_STAR_index.cmds
launcher_creator.py -t 12:00:00 -j 05_STAR_index.cmds -n 05_STAR_index -l 05_STAR_index.slurm -e rayna.harris@utexas.edu -q normal
sbatch 05_STAR_index.slurm
~~~

## Mapping reads with STAR

STAR is currently not a stampede module, but it is installed in Dhivya's work directory, which should be in your path.

Create the commands file. 

~~~ {.bash}
rm 05_STAR_mapping.cmds
for R1 in *R1_001.filtrim.fastq.gz
do
    R2=$(basename $R1 R1_001.filtrim.fastq.gz)R2_001.filtrim.fastq.gz
    samp=$(basename $R1 _L002_R1_001.filtrim.fastq.gz)
    echo $R1 $R2 $samp
    echo "STAR --runThreadN 1 --outFileNamePrefix $samp --outSAMstrandField intronMotif --genomeDir $SCRATCH/BehavEphyRNAseq/refs/STAR --outSAMtype BAM Unsorted --readFilesIn $R1 $R2" >> 05_STAR_mapping.cmds
done
~~~

Now, run the job on Stampede.

~~~ {.bash}
launcher_creator.py -t 08:00:00 -j 05_STAR_mapping.cmds -n 05_STAR_mapping -l 05_STAR_mapping.slurm -e rayna.harris@utexas.edu -q normal
sbatch 05_STAR_mapping.slurm
~~~

Now, move the outut to a new directory
~~~ {.bash}
mkdir ../05_STAR
mv STAR* ../05_STAR
cd ../05_STAR
~~~

## Assessing alignment



Load the program Samtools. Then sort and index newly created BAM file.

~~~ {.bash}
module load samtools
samtools sort STARAligned.out.bam STARAligned.out.sorted.bam
samtools index STARAligned.out.sorted.bam
~~~