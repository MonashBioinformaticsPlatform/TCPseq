TCPseq data analysis pipeline.

UPDATE 27/1/17 (2):
Installed PANDAS:
$ sudo pip install pandas
# pandas (0.19.2) 
$ sudo pip install docopt
# docopt (0.6.2)

----------------------------------------------------------------------

UPDATE 27/1/17 (1):
This development branch (not fully tested) uses R3.3.2 and some extra R packages not installed on the VM by default.
$ sudo vim /etc/apt/sources.list
# in vim append the following line to the .list file: 
deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu/ trusty/
$ sudo apt-get update
$ sudo apt-get install r-base 
# as of 27/1/17 this is installs R 3.3.2
$ sudo R -e 'install.packages("plyr", repos = "http://cran.us.r-project.org")' # ENSURE this is v1.8.4
$ sudo R -e 'install.packages("dplyr", repos = "http://cran.us.r-project.org")'
$ sudo R -e 'install.packages("ggplot2", repos = "http://cran.us.r-project.org")'
$ sudo R -e 'install.packages("docopt", repos = "http://cran.us.r-project.org")'
# The following is the sessionInfo() in R after loading these 4 libraries:
R version 3.3.2 (2016-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.5 LTS

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] docopt_0.4.5  ggplot2_2.2.1 dplyr_0.5.0   plyr_1.8.4

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.9      assertthat_0.1   grid_3.3.2       R6_2.2.0
 [5] gtable_0.2.0     DBI_0.5-1        magrittr_1.5     scales_0.4.1
 [9] stringi_1.1.2    lazyeval_0.2.0   tools_3.3.2      stringr_1.1.0
[13] munsell_0.4.3    colorspace_1.3-2 tibble_1.2

Example of how to run the genewise bootstrapping script and the plotting script on the 'tbl3_' bed files (100 resamples) :
mkdir mtgn
python scripts/metagene_agg.py -b tabulated_data/tbl3_169_100Kreads.bed.gz -g sacCer3/annotation/genetable_dist_to_nextgene.txt -o mtgn/sa169 -p 100
Rscript scripts/dens_plt.R -i mtgn/sa169 -l -30 -r 20

----------------------------------------------------------------------

UPDATE 6/9/16:
For steps involving STAR, issue 143 may affect running. The work-around is to delete the --readFilesCommand zcat part and use --readFilesIn <(zcat filename.gz) as per https://github.com/alexdobin/STAR/issues/143

----------------------------------------------------------------------

The pipeline below requires software dependencies that are all included in the virtual machine http://bioinformatics.erc.monash.edu/home/powell/TCPseq.box. For further information on how this VM was generated see vagrant/README.txt and vagrant/setup.sh.

The TCPseq pipeline is divided into two parts:

1) Downloading and setting up indices and other data sets for mapping against. This is only required if the default data sets, viewable in the setup_config.sh script, are not suitable. These default references and the required directory structure are all included in the abovementioned virtual machine. To change these, see the code in scripts/setup_vbox.sh 

2) Converting fastq sequencing data from a TCP-seq experiment into FP start/end coordinates within each annotated transcript. This pipeline is included in the manuscript and in the file pipeline_summary.sh; updates in this git repository should be included in the latter.   

