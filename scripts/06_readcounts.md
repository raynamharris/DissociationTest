count read in fastq.gz files
----------------------------

The files in my directory are zipped. Each R1 and R2 have the same
number of reads, so I just count the R1 reads.

I use zcat to unzip the files and then I pipe it to wc and then I send
the output to a file.

    for file in *R1_001.fastq.gz
    do
    echo $file
    zcat $file | echo $((`wc -l`/4)) >> readcounts.txt
    done 

I did this on TACC then saved the file locally with scp.

Now, I can calculate the average numbe of reads and standard deviation.

    reads <- read.table("../results/readcounts.txt")
    summary(reads)

    ##        V1         
    ##  Min.   : 831222  
    ##  1st Qu.:3447633  
    ##  Median :4547688  
    ##  Mean   :4918879  
    ##  3rd Qu.:6600756  
    ##  Max.   :9496400

    mean(reads$V1)/1000000 # to show in millions

    ## [1] 4.918879

    sd(reads$V1)/1000000 # to show in millions

    ## [1] 2.606619

On average, my samples yielded 4.9 +/- 2.6 million reads.

    aligned <- read.table("../results/pseudoaligned_clean.txt", header = T)

    aligned$percent <- (aligned$pseudoaligned / aligned$processed) * 100
    summary(aligned)

    ##    processed       pseudoaligned        percent      
    ##  Min.   : 313015   Min.   :  65987   Min.   : 5.419  
    ##  1st Qu.:1479506   1st Qu.: 748285   1st Qu.:49.000  
    ##  Median :2958364   Median :2128611   Median :70.112  
    ##  Mean   :3302250   Mean   :2324577   Mean   :61.233  
    ##  3rd Qu.:3434580   3rd Qu.:2464398   3rd Qu.:74.323  
    ##  Max.   :8513137   Max.   :6653688   Max.   :79.540

    mean(aligned$percent)

    ## [1] 61.23341

    sd(aligned$percent)

    ## [1] 20.84956

On average, 61.2% +/1 20.8% of the trimmed reads were psuedoaligned to
the transcriptome.

MultiQC
-------

Here are the results QC before and after filtering and trimming reads

Before ![before](../figures/06_readcounts/before.png)

After ![after](../figures/06_readcounts/after.png)
