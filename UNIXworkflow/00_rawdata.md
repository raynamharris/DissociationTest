# Get, Store, and/or Retrieve Raw Data

Includes the process and commands needed to:
* download the data to scratch on Stampede with `00_gsaf_download`
* save a copy on Corral with `00_storeoncorral`
* retrieve the data from coral with `00_getfromcorral` 

The first two steps are down by the person who submitted the RNA samples for sequencing. In theory, this only has to be done once, but it could be done again if data are lost or compromised. 

The third step will be performed by any collaborator who wants to copy the raw data to their own working directory. 

## 00_gsaf_download 

The first step in the bioinformatics pipeline is to download the data from the sequencing facility, in this case GSAF. The GSAF provides the script for downloading the data from there amazon web server. 

1. Setup project and data directories on Stampede.
2. Copy script to download data: Save a script to download the data (called "00_gsaf_download.sh") in each subdirector  
3. Download the data using TACC: For this three-part step I created a commands file (called "00_gsaf_download.cmd") with the command to execute the download script,  created a launcher script to execute the commands file to execute the launcher script (I know, sounds like a lot of steps, but this is so I use TACC's compute power not my own), and launched the job on TACC  
4. Repeat steps 2 and 3 for all RNAseq jobs 

### Setup project and RNAseq job directories 

Login to TACC using ssh with your password and authentication credentials. Replace "<username>" with your TACC user name. 

~~~ {.bash}
ssh <username>@stampede.tacc.utexas.edu
~~~

Raw data from separate jobs need to be processed as separate jobs. Later, the read counts can be combined, but original processing must be done by job. So, let's create a an environment three environment variables for the project, job, and path to data so that the scripts can easily be co-opted for each new RNAseq job. 

~~~ {.bash}
## set the enviornment variables 
RNAseqProject=<name of project directory>
RNAseqJob=<name of sequencing job>
AmazonAddress=<"website to data key">
~~~

On scratch, create the project directory (SingleNeuronSeq), with a subdirectory for each job (in this case JA15597 and JA16033) and subsubdirectory called 00_rawdata. The argument `-p` will create the parent and subdirectories if they do not already exist.

~~~ {.bash}
mkdir -p $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
~~~

### Copy script to download data 

~~~ {.bash}
## Copy the text for the gsaf download script found here:  https://wikis.utexas.edu/display/GSAF/How+to+download+your+data 
~~~ 

Navigate to one of the subjectories. Use the program nano to open a text file.  I use the program nano to open a new text file. Paste the script and save it as `00_gsaf_download.sh`.

~~~ {.bash}
cd $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
nano
~~~ 

Now, you should have one file called `00_gsaf_download.sh`. Check with `ls`

~~~ {.bash}
ls
~~~ 

### Download the data using TACC
Technically, you can download the data with one command using your own comptuer's compute power, but I prefer to have TACC do it. So, instead of type the command which will do just that in the command line, I will save it to a script, like so. (Note: the webaddress provided was sent by secure email to the person who submitted the samples to GSAF). Then, you must make the bash script executable with `chmod`.

~~~ {.bash}
echo '00_gsaf_download.sh $AmazonAddress' > 00_gsaf_download.cmds
chmod a+x 00_gsaf_download.sh
~~~

Now, I use `launcher_creator.py` to create a launcher script that will tell how to launch this job on TACC. The arguments are defined clearly on this website: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py. Then I will use `sbatch 00_gsaf_download.slurm` to launch the job.

~~~ {.bash}
launcher_creator.py -t 12:00:00 -n 00_gsaf_download -j 00_gsaf_download.cmds -l 00_gsaf_download.slurm -A NeuroEthoEvoDevo -q normal
sbatch 00_gsaf_download.slurm
~~~

### Make all fastq.gz files READ ONLY!

To protect our data from accidental deletion or overriding, let's make all the files read only!

~~~ {.bash}
chmod a-w *fastq.gz
ls -l
~~~

### Clean up a little

Create a new directory called `wget_log` and move all the output files associated with downloading the data.

~~~ {.bash}
mkdir wget_log
mv *.wget.log wget_log
mv files.html wget_log
mv md5.txt wget_log
ls
~~~

### Repeat for all jobs. 
The great thing about using TACC for this is that you can go about doing other things while the files are download. Reset the environment variables above and repeat the process for all other RNAseq jobs.


## 00_storeoncorral

### Things to do while working on Corral

Open a second terminal window and login to corral, where we store data long-term

~~~ {.bash}
ssh <username>@corral.tacc.utexas.edu
~~~

Set environment variables on corral.

~~~ {.bash}
RNAseqProject=SingleNeuronSeq
RNAseqJob=JA15597
~~~

Create the directories on coral (which is meant for longterm storage) in the Hofmann lab repository called `NeuroEthoEvoDevo`. 


~~~ {.bash}
cd /corral-tacc/utexas/NeuroEthoEvoDevo
mkdir -p $RNAseqProject/$RNAseqJob/00_rawdata
~~~ 

### Things to do while working on Stampede

Return to your scratch directory on Stampede!

~~~ {.bash}
cd $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
~~~

Confirm the enviornment variables with `env` or reset them if necessary.

~~~ {.bash}
## set the enviornment variables 
RNAseqProject=SingleNeuronSeq
RNAseqJob=JA16033
~~~

Copy the files to corral

~~~ {.bash}
cd $SCRATCH/$RNAseqProject/$RNAseqJob/00_rawdata
scp *.fastq.gz <username>@login1.corral.tacc.utexas.edu:/corral-tacc/utexas/NeuroEthoEvoDevo/$RNAseqProject/$RNAseqJob/00_rawdata
~~~

Repeat for other datasets.

## 00_getfromcorral

This is for people who want to copy data from corral to their own working directory.

Login to TACC using ssh with your password and authentication credentials. Replace "<username>" with your TACC user name. 

~~~ {.bash}
ssh <username>@stampede.tacc.utexas.edu
~~~

Create a two environment variables for the project and the job. By setting variables, we can easily repeat this for new projects. In this case "<nameofproject>" should be replaced with "SingleNeuronSeq" and "<jobnumber>" should be repalced with either "JA15597" or "JA16033"

~~~ {.bash}
## set the enviornment variables 
RNAseqProject=SingleNeuronSeq
RNAseqJob=JA16033
~~~

Now, make directories on scratch where the data can be stored. This utilizes the variable names we just created. The `-p` argument says create parent and sub directories, if necessary.  Then, navigate to the new directory

~~~ {.bash}
mkdir -p $SCRATCH/$RNAseqProject/$RNAseqJob
cd $SCRATCH/$RNAseqProject/$RNAseqJob
~~~

Now, copy the data from corral with the `scp` command. Th.e last period meas "copy here"

~~~ {.bash}
scp -r <username>@login1.corral.tacc.utexas.edu:/corral-tacc/utexas/NeuroEthoEvoDevo$RNAseqProject/$RNAseqJob/00_rawdata .
~~~

To repeat for more directories, modify the RNAseqJob variable and execute again.
