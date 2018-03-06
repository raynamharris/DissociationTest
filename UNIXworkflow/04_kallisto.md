# Kallisto: Quantifying RNA Transcript Abundances

To quantify transcripts, we first need to have a reference genome or transciptome to which the reads can be mapped. 

## Building a refernece index: Only do this once!

### Download a reference transcriptome

Download mouse transcriptome from https://www.gencodegenes.org/mouse_releases/current.html

~~~ {.bash}
# make a directory the reference transcriptome and index
mkdir $SCRATCH/refs
cd $SCRATCH/refs
# download mouse transcriptome, version M11, from Gencode
curl -O ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.pc_transcripts.fa.gz
~~~

### 04_kallistoindex

The kallisto index only needs to be built once for each. The manual for running kallisto can be found [here](https://pachterlab.github.io/kallisto/manual). I like to keep the same really long prefix "gencode.vM11.pc_transcripts" for the name so that I always know where it came from, rather than shortening it to something like "mouse" or "M11" becuase the full name is more informative. Then, I add "_kallisto.idx" to the end because this tells me that the index is specifically for kallisto, rather than any other alignment/mapper program.

~~~ {.bash}
# do this on an idev node
idev -m 60
kallisto index -i gencode.vM11.pc_transcripts_kallisto.idx gencode.vM11.pc_transcripts.fa.gz
~~~

## Pseudoalignment with Kallisto

I'm a big fan of the kallisto program because its super fast and easy to use! Its also becoming more widely used and trusted.

### 04_kallistoquant

Navigate to the directory with the processed reads and make a directory where the output can be stored. 

~~~ {.bash}
cd $SCRATCH/DissociationTest/02_filtrimmedreads
mkdir ../04_kallistoquant
~~~

Now, we will use the `kallistoquant` function to quantify reads! Again, we use a for loop to create the commands file. The output for each pair of samples will be stored in a subdirectory.  

The launcher creator command that I was using on Stampede 1 doesn't work on Stampede2, so I've modified my workflow a bit. Rather than use the launcher, I created a slurm file from scratch (see 04_kallistoquant.slurm) and put all the commands in a .exe file. 

~~~ {.bash}
for R1 in *R1_001.filtrim.fastq.gz
do
    R2=$(basename $R1 R1_001.filtrim.fastq.gz)R2_001.filtrim.fastq.gz
    samp=$(basename $R1 _L002_R1_001.filtrim.fastq.gz)
    echo $R1 $R2 $samp
    echo "kallisto quant -b 100 -t 16 -i ../../refs/gencode.vM11.pc_transcripts_kallisto.idx  -o ../04_kallistoquant/${samp} $R1 $R2 &> ${samp}.errout.txt" >> 04_kallistoquant.exe
done

chmod u+a 04_kallistoquant.exe

sbatch 04_kallistoquant.slurm
~~~

### checking status

I'm impatient, so I like to check the status of my files while they are being processed. I use this for loop to print the name of each file and the number of reads process and pseudoaligned.

~~~ {.bash}
for file in *errout.txt
do
echo $file
cat $file | grep 'processed'
echo " " 
done
~~~

### job summary
Note: this took 45 minutes running 14 samples on 4 cores. Next time, I'll ask for 14 cores, but I wasn't quite sure who to do that. 

### Calculating percent mapped and stdev

I can do this with multi-QC (see next step), but I can also modify the for loop a bit  to calculate the number of reads pseudoaligned in R.  

~~~ {.bash}
for file in *errout.txt
do
echo $file
cat $file | grep 'processed' >> pseudoaligned.txt
echo " " 
done
~~~

Using `scp`, I saved this file to my results directory, then I cleaned it up using find and replace to save only the number of reads processed and reads pseudoaligned. Then, I used `05_readcounts.Rmd` to calculate the percent of reads mapped. 



## Renaming files

I wanted to remove the uninformative bits of the "sample name" so they match up with the actual sample name that I use for downstream analysis. To do this, I first remove the parts added by the sequencing facility (_S*) and then I replace underscores with dashes. 

~~~ {.bash}
for file in *
do
    sample=${file//_S*/}
    echo $file $sample
    mv $file $sample
done
~~~

Then, replace the `_` with `-`

~~~ {.bash}
for file in *
do
    sample=${file//_/-}
    echo $file $sample
    mv $file $sample
done
~~~


## References
- Kallisto: https://pachterlab.github.io/kallisto/
- Kallisto on Stampede: https://wikis.utexas.edu/display/bioiteam/Pseudomapping+with+Kallisto

