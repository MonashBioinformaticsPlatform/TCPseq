TCPseq data analysis pipeline.

UPDATE 6/9/16:
For steps involving STAR, issue 143 may affect running. The work-around is to delete the --readFilesCommand zcat part and use --readFilesIn <(zcat filename.gz) as per https://github.com/alexdobin/STAR/issues/143

The pipeline below requires software dependencies that are all included in the virtual machine http://bioinformatics.erc.monash.edu/home/powell/TCPseq.box. For further information on how this VM was generated see vagrant/README.txt and vagrant/setup.sh.

The TCPseq pipeline is divided into two parts:

1) Downloading and setting up indices and other data sets for mapping against. This is only required if the default data sets, viewable in the setup_config.sh script, are not suitable. These default references and the required directory structure are all included in the abovementioned virtual machine. To change these, see the code in scripts/setup_vbox.sh 

2) Converting fastq sequencing data from a TCP-seq experiment into FP start/end coordinates within each annotated transcript. This pipeline is included in the manuscript and in the file pipeline_summary.sh; updates in this git repository should be included in the latter.   

