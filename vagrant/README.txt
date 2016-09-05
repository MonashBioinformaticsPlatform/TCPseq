# To build the TCPseq vagrant box from scratch

# Get the trusty64 box:
vagrant box add ubuntu/trusty64

# Create and provision
vagrant up

# Optional: ssh in and look around
vagrant ssh


# Create the box
vagrant package --output TCPseq.box






#####################################
# To use on a new system

vagrant init http://bioinformatics.erc.monash.edu/home/powell/TCPseq.box
vagrant up
vagrant ssh

# Then use. eg:
source TCPseq_config.sh
cd trimmed
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}.fq.gz
  trimmomatic SE -threads $NCORES ../$INPUT_DATADIR/${FQ_FNS[$i]} tr_$bn SLIDINGWINDOW:7:24 MINLEN:26
  cutadapt --discard-untrimmed -a AAAAAAAAAAAA -e 0.0 -O 12 tr_$bn | gzip > pA_tr_$bn
done
