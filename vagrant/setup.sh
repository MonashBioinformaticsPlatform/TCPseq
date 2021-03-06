#!/bin/sh

# Script to provision a base ubuntu/trusty64 image to run TCPseq 

# Add repository for R 3.3.2
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
sudo add-apt-repository 'deb http://cran.cnr.berkeley.edu/bin/linux/ubuntu/ trusty/'

# apt-get packages
sudo apt-key update
sudo apt-get update
sudo apt-get install -y build-essential curl git \
                        python-setuptools python-pip python-dev \
                        language-pack-en \
                        openjdk-7-jre-headless
sudo apt-get install -y samtools r-base-core=3.3.2-1trusty0 trimmomatic r-cran-reshape2 r-cran-plyr bowtie2

# R packages
sudo R -e 'install.packages(c("plyr","dplyr","ggplot2","docopt"), repos = "http://cran.us.r-project.org");'

# STAR
mkdir STAR-tmp && cd STAR-tmp \
 && curl -sL https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz > star.tgz \
 && tar xzf star.tgz \
 && sudo mv */bin/Linux_x86_64_static/* /usr/local/bin \
 && sudo mv /usr/local/bin/STAR /usr/local/bin/STAR.orig \
 && cd ~ && rm -rf STAR-tmp

## Trimmomatic
curl -sL http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip > trimmomatic.zip \
  && sudo unzip trimmomatic.zip -d /usr/local
sudo sh -c 'cat > /usr/local/bin/trimmomatic' <<EOM
#!/bin/bash
exec java -Xmx3g -jar "/usr/local/Trimmomatic-0.33/trimmomatic-0.33.jar" \$@
EOM
sudo chmod +x /usr/local/bin/trimmomatic
rm trimmomatic.zip

# BioPerl
sudo apt-get install -y --no-install-recommends bioperl

# BioPython
sudo apt-get install -y python-biopython

# PIP libraries.
sudo pip install intermine 
sudo pip install cutadapt
sudo pip install pandas
sudo pip install docopt

# Cleanup
sudo apt-get clean -y
sudo apt-get autoclean -y

# On login, copy into place.  And move to shared directory on login
cat >>~/.bash_profile <<EOM
if [ ! -e /vagrant/TCPseq ];
then
  echo "########## Copying initial TCPseq files into shared directory"
  rsync -a ~/TCPseq /vagrant
fi
cd /vagrant/TCPseq
EOM

# Clean the MOTD
sudo chmod -x /etc/update-motd.d/*
sudo rm /run/motd.dynamic
sudo sh -c 'cat > /etc/motd' <<EOM
-----------------------------------------
Welcome to the TCPseq vagrant box

EOM

# Checkout the repo
git clone https://github.com/MonashBioinformaticsPlatform/TCPseq.git

# Put the STAR wrapper in place
sudo cp ~/TCPseq/vagrant/STAR-wrapper.sh /usr/local/bin/STAR
sudo chmod +x /usr/local/bin/STAR

# Now build all the TCPseq refs
(cd ~/TCPseq ; bash ./scripts/setup_vbox.sh)

