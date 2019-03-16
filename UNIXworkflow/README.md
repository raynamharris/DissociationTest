# RNAseq workflow on TACC
Last modified October 25, 2016

## The workflow
* **00_rawdata:** Includes the commands need to download the data to scratch on Stampede with `00_gsaf_download`.
* **01_fastqc:** Includes the commands needed to evaluate the quality of the reads using the program FastQC.
* **02_filtrimreads:** Includes the commands needed to filter low quality reads and trim adapters using the program cutadapt.
* **03_fastqc:** Assess the quality of the processed reads
* **04_kallisto:** Quantify transcript expression.
* **04_kallistoquant.slurm:** An example slurm file