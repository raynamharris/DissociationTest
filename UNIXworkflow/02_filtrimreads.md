# Cutadapt: Filter and Trim Low Quality Reads

Lets remove adapters (trim) and low quality reads (filter) with the program cutadapt.

## Navigate to your working directory on stampede, 

If necessary, rest the environment variables.

~~~ {.bash}
## set the enviornment variables 
RNAseqProject=<nameofproject>
RNAseqJob=<jobnumber>
~~~

Navigate to your working directory on stampede.

~~~ {.bash}
ssh <username>@stampede.tacc.utexas.edu
cd $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
~~~

## Using cutadapt on TACC

The cutadapt program is not maintained by TACC. We need to use the program stored in a colleagues working directory. To access the program, set the path. 

~~~ {.bash}
PATH=/work/01184/daras/bin/cutadapt-1.3/bin:$PATH
~~~

## Part 1

~~~ {.bash}
rm 02_filtrimmedreads.cmds
for R1 in *R1_001.fastq.gz
do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    R1filtrim=$(basename $R1 fastq.gz)filtrim.fastq.gz
    R2filtrim=$(basename $R2 fastq.gz)filtrim.fastq.gz
    tempR1=$(basename $R1 fastq.gz)temp.fastq.gz
    tempR2=$(basename $R2 fastq.gz)temp.fastq.gz
    echo $R1 $R2 $R1filtrim $R2filtrim
    echo "cutadapt -q 10 -a GATCGGAAGAGCACACGTCTGAACTCCA  -m 22 --paired-output $tempR2 -o $tempR1 $R1 $R2" >> 02_filtrimmedreads.cmds
done 
~~   
    

~~~ {.bash}
launcher_creator.py -t 4:00:00 -n 02_filtrimmedreads -j 02_filtrimmedreads.cmds -l 02_filtrimmedreads.slurm -A NeuroEthoEvoDevo -q 'normal'  
sbatch 02_filtrimmedreads.slurm
~~~

## Part 2

~~~ {.bash}
rm 02_filtrimmedreads.cmds
for R1 in *R1_001.fastq.gz
do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    R1filtrim=$(basename $R1 fastq.gz)filtrim.fastq.gz
    R2filtrim=$(basename $R2 fastq.gz)filtrim.fastq.gz
    tempR1=$(basename $R1 fastq.gz)temp.fastq.gz
    tempR2=$(basename $R2 fastq.gz)temp.fastq.gz
    echo "cutadapt -q 15 -a ATCGTCGGACTGTAGAACTCTGAACGTG -o $R2filtrim -p $R1filtrim $tempR2 $tempR1" >> 02_filtrimmedreads.cmds
done 
~~  

~~~ {.bash}
launcher_creator.py -t 4:00:00 -n 02_filtrimmedreads -j 02_filtrimmedreads.cmds -l 02_filtrimmedreads.slurm -A NeuroEthoEvoDevo -q 'normal'  
sbatch 02_filtrimmedreads.slurm
~~~ 

## Part 3 : remove temp files

~~~ {.bash}
rm 02_filtrimmedreads.sh
for R1 in *R1_001.fastq.gz
do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    R1filtrim=$(basename $R1 fastq.gz)filtrim.fastq.gz
    R2filtrim=$(basename $R2 fastq.gz)filtrim.fastq.gz
    tempR1=$(basename $R1 fastq.gz)temp.fastq.gz
    tempR2=$(basename $R2 fastq.gz)temp.fastq.gz
    echo $tempR2 $tempR1
    echo "rm $tempR2 $tempR1" >> 02_filtrimmedreads.sh
done
chmod a+x 02_filtrimmedreads.sh
bash 02_filtrimmedreads.sh
~~~

## Clean up

Now, let's make our processed reads read only so we don't accidentally modify them. 

~~~ {.bash}
chmod a-w *filtrim.fastq.gz 
~~~

Now, move the processed reads to a new file.

~~~ {.bash}
mkdir ../02_filtrimmedreads
mv *filtrim.fastq.gz ../02_filtrimmedreads
~~~




## References
Cutadapt: http://cutadapt.readthedocs.io/en/stable/guide.html
BioITeam Launcher Creator: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py
Dhivya's Cutadapt: 
