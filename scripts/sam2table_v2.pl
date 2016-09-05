#!/usr/bin/perl -w

#
# WARNING: not tested for Out By One errors in reverse ori
#

# example: 
# cd /Users/sarcher/VMs/bioinf/shared_data/projects/sarcher/43S
# export pm=./tabulate_mRNA_al/sam2table_v2.pl
# for f in ro_data/mR*.bam;  do  if [[ $f =~ mR900f_(str.+).bam ]]; then samtools view -h $f | perl $pm > ./tabulate_mRNA_al/tbl_${BASH_REMATCH[1]}\.txt;  fi; done

#for f in *Aligned.out.sam; do 
#    echo $f; 
#    cat $f | perl ../../tabulate_mRNA_al/sam2table_v2.pl -min 17 \
#        -out ../../int_processing/remap_160915/${f%Aligned.out.sam}_min17.sam \
#        > ../../tabulate_mRNA_al/remap_160915/${f%Aligned.out.sam}_m17al.txt ; 
#done

# TEST:
# cd 43S
# samtools view -hS ./int_processing/FAKEDdata/FAKEAligned_YOL082W.out.sam | perl ./tabulate_mRNA_al/sam2table_v2.pl -out  ./int_processing/FAKEDdata/_min17.sam > ./int_processing/FAKEDdata/FAKE_al.txt


use strict;
use Getopt::Long;


my $minlen = 17;
my $outsamfile='';
my $count_leading_Ss = '0'; #<--- FALSE: will now require no correction in Rscript Process_temp_before_plotting.R

GetOptions ( 'min=s' => \$minlen,
            'out=s' => \$outsamfile,
            'countleadS=s' => \$count_leading_Ss);

my $lcount = 0;
my $excount = 0;
my $FOUT;
if ($outsamfile ne '') {
    open($FOUT, '>' , $outsamfile) or die ('Cant open file for writing sam format.');
}


while(<STDIN>){
    chomp;
    my $line = $_;
    s/\'//g;
    my @splut = split/\t/;
    if (scalar(@splut) > 10) {
        my $flag = $splut[1];
        my $ori = "1";
        my $cig = $splut[5];
        my $readseq = $splut[9];
        $lcount ++;
        
        # regex to filter out tRNAs, snRNAs and RDNs:
        #if ($splut[2] =~ /t.\([ACGTU][ACGTU][ACGTU]\).+/ or $splut[2] =~ /^snR\d.+/ or $splut[2] =~ /^RDN\d/) {
        #    $len = 0;
        #}
        
        if(  $flag & hex("0x200") or $flag & hex("0x4") or $flag & hex("0x10") or length($splut[9]) < 12 ) {
            $excount ++;
            #if ($flag & 0x4  ) {print "segment unmapped"}  # this is unexpected as reference consisted of forward orientation mRNA sequences
            #if ($flag & 0x10 ) {print "SEQ being reverse complemented"}
            #if ($flag & 0x200 ) {print "not passing quality controls"}
        } else {
            my $secondary = "pri";
            if ($flag & 0x100 ) {
                $secondary = "sec"
            }
            my @NH_tag = grep {/NH:i/} @splut; 
            my $NH = 'NA';
            if ($NH_tag[0] =~ /NH:i:(\d+)/){
                $NH=$1;
            }
            
            # skip reads with large (>3 nt) gaps in alignment (with downstream M's: 3' terminal SNIDs we will fix below)
            my $mm = 0;  # mm = max conseq. mismatches
            while ($cig =~ /(\d+)[SNID]/g){
                my $current = $1;
                if (  $mm < $current ) {
                    my $rest_of_cig = substr $cig, $-[0]+1;
                    if ($rest_of_cig =~ /M/ ){ # if it's internal mismatch (ie there's a downstream M in cig)
                        $mm = $current;
                    }
                } 
            }
            if ($mm > 3) {
                $excount ++;
                next;
            }
            
            my $trimmedlen="NA";
            # trim of AAAAAAAG/T read 3' ends not mapping to anything
            if ($cig =~ /(\d+)[SNID]$/) {
                my $terminal_s = $1;
                if ($terminal_s > 3) {
                    if ($readseq =~ /([A]+[GT])$/) {
                        if (length($1)> 2) {  # if the AAAAAG/T length is nearly as long as (or longer than) number of terminal Ss in cig
                            # i.e if alignment has 4 or more 3' nts soft-clipped AND 3' sequence is 3 or more A's followed by a G/T nucleotide, THEN:
                            #print $cig."\t".$readseq."\t".$splut[10]."\n";
                            $readseq =~ s/([A]+[GTC])$//;   # trim the AAAAA[GTC]
                            $trimmedlen = length($1);
                            # repair cigar to reflect trimmed read. This will leave 3' S/N/Is intact if that's where trimming ended!
                            my $i = 0;
                            while ($i < $trimmedlen){  #trimmedlen is the number of nucleotides we want to trim off
                                if ($cig =~ s/(\d+)([SNIDM])$//) {  #chops off last number-letter pair from cigar and loads into $1,$2
                                    my $num = $1;
                                    my $letter = $2;
                                    if ($letter ne 'D') { # Intervening D's (deletions in the read) are trimmed for free from cigar. (eg. AA-AAG = 5 nt trim length but 6 elements must be trimmed from cigar)    
                                        $i ++;
                                    }
                                    $num -= 1;
                                    if ($num > 0) {
                                        $cig=$cig.$num.$letter
                                    }     
                                } else {
                                    warn "Warning: cigar stripped from $splut[5] to nothing... read: $readseq \n";  # this should never happen
                                    continue;
                                }
                            }
                            
                            #print $cig."\t".$readseq."\t".$splut[10]."\n";
                        }
                    }
                }
            }
            if ($cig =~ /(.+M)/) {
                $cig = $1  #  <<<< with greedy matching, ensures last element in cig is a match ('M'). Gets rid of 3' soft-clipped bases prior to ciglen calculation below >>>>
            } else {
                warn "Warning: rest of cigar $splut[5] had no matches. Read: $readseq \n";
                next;
            }
            my $ciglen = crunch_cigar($cig);
            my $readlen = length($readseq);
            
            my $initial_nonmatching = '0';
            if ($cig =~ /^(\d+?)[SI]/) {    # look for leading soft-clipped nucleotides and move the 5' end back accordingly
                #warn "$cig   ...  $1 \n";
                $initial_nonmatching = $1;
            }
            if ($count_leading_Ss) {
                $splut[3] = $splut[3] - $initial_nonmatching; # This was done by default before
            } else {
                $ciglen = $ciglen - $initial_nonmatching;
            }
              
            @splut = @splut[0,2,3,5];
            push @splut, ($readlen, $cig, $ciglen, $secondary, $NH, $initial_nonmatching);
            print join "\t", @splut;
            print "\n";
            
            if ($outsamfile ne '' && $ciglen >= $minlen) {
                print $FOUT "$line\n";
            }
            
            
        }
    } elsif(scalar(@splut) < 5 && $line =~ /^@/ && $outsamfile ne '') {  # print sam header to filtered sam file if specified
                print $FOUT "$line\n"; 
    }
}


warn "all lines: $lcount . excluded:  $excount \n\n";

#### crunch cigar
sub crunch_cigar {
    my $cigar = shift;
    my %cig_sum;
    my %symb = ('M' => '1', 'S' => '1', 'I' => '0', 'D' => '1');
    # reminder: I's are insertions in the READ relative to the REF (AKA deletions in REF rel READ). D's are DELETIONS in the READ relative the REF (AKA insertions in REF rel READ)
    #NOTE: S->1 is different from how I did it with Ion Torrent data
    my $adj_len = 0;

    foreach(keys(%symb)){
        while($cigar =~ m/(\d+)$_/g){
            $cig_sum{$_} += $1;
        }  
    }

    $adj_len += $symb{$_}*$cig_sum{$_} foreach(keys(%cig_sum));
    return $adj_len;
}
