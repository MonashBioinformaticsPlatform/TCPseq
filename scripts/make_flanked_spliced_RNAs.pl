#!/usr/bin/perl -w

use strict;
use Getopt::Long;

# usage example:
# export PERL5LIB=/shared_data/dev/bioperl-live/
# perl make_flanked_orfs.pl -dir /shared_data/common_data/genomes/sacCer3/S_cer_genome_gbk_April2013 -out /shared_data/common_data/genomes/sacCer3/orfs_800flank.fa -5p 800 -3p 800

# for structural RNA filter: (@terms = ('misc_RNA', 'ncRNA', 'tRNA', 'snoRNA'))
# perl make_flanked_orfs.pl -5p 50 -3p 50 -t misc_RNA ncRNA tRNA snoRNA -dir /shared_data/common_data/genomes/sacCer3/S_cer_genome_gbk_April2013 -out /shared_data/common_data/genomes/sacCer3/miscR_50flank.fa

# for no-dubious + overlapping orfs:
# perl make_flanked_orfs.pl -dir ~/VMs/bioinf/shared_data/common_data/genomes/sacCer3/S_cer_genome_gbk_April2013 \
# -out ~/VMs/bioinf/shared_data/common_data/genomes/sacCer3/custom_yeast/mRNAs_1000flank.fa -5p 1000 -3p 1000 -exdub yes

# 16/9/15:
# perl make_flanked_orfs.pl -dir ~/VMs/bioinf/shared_data/common_data/genomes/sacCer3/S_cer_genome_gbk_April2013 \
# -blacklist ~/VMs/bioinf/shared_data/projects/sarcher/43S/int_processing/Dubious_ol_orfs2.txt \
# -out ~/VMs/bioinf/shared_data/common_data/genomes/sacCer3/custom_yeast/mRNAs_exdub_1000flank.fa -5p 1000 -3p 1000 -exdub yes


#use lib '/shared_data/dev/bioperl-live';
use lib '/Users/sarcher/VMs/bioinf/shared_data/dev/bioperl-live/';

use Bio::SeqIO;
use Bio::Location::Simple;
use Bio::SeqFeature::Primer;

$/ = "\n"; # set newline character to linefeed

my $dir = "";
my $UTRs="";
my $fpdf="0"; # 5' default
my $tpdf="0"; # 3' default
my $format="genbank";
my $out="Flanked_ORFs_output.txt";
my $genome_file_term = ('.dat');
my $exclude_dubious = "no";
my $external_blacklist = "0";  # <---- FALSE
my $external_whitelist = "0";  # <---- FALSE
my @terms;
#my @terms = ('rRNA','misc_RNA');
#@terms = ('mRNA', 'misc_RNA', 'ncRNA', 'transcript', 'snoRNA');  # may change this to user input in later version....     -t <transcript types (multiple possible)> -ex <exclusion list file>   -in <inclusion list file> 
#my @terms = ('mRNA');  # <------- to make mRNAs_1000flank.fa
#@terms = ('misc_RNA', 'ncRNA', 'tRNA', 'snoRNA');  # may change this to user input in later version....     -t <transcript types (multiple possible)> -ex <exclusion list file>   -in <inclusion list file> 
my %keywords = ('RDN18-1' => 'y', 'RDN25-1' => 'y');


########## Parse input ############

GetOptions ( 'dir=s' => \$dir,
             'u=s' => \$UTRs,
             '5p=s' => \$fpdf,
             '3p=s' => \$tpdf,
             'f=s' => \$format,
             'out=s' => \$out,             
	     't=s{0,100}' => \@terms,
	     'k=s' => \$genome_file_term,
	     'exdub=s' => \$exclude_dubious,
	     'blacklist=s' => \$external_blacklist,
		 'whitelist=s' => \$external_whitelist);

my @errors;
push @errors, "folder '$dir' not found in current directory " unless (-e $dir);
push @errors, "no folder given for search sequences" if ($dir eq "");
push @errors, "Default 5-prime length nominated, \"$fpdf\", is not integer value" if ($fpdf =~/\D/);
push @errors, "Default 3-prime length nominated, \"$tpdf\", is not integer value" if ($tpdf =~/\D/);
push @errors, 'format (-f) must be either "genbank" or "ensembl" or left blank.' unless ($format eq "genbank" or $format eq "ensembl");
push @errors, "file $UTRs not found" unless ($UTRs eq "" or -e $UTRs);

my $seq_out = Bio::SeqIO->new(-file => ">$out", -format => 'fasta');

unless ($UTRs eq ""){ open FINU, "<$UTRs" or push @errors, "Couldn't open UTRs file $UTRs" } #Gene name \t 5p UTR len \t 3p UTR len

if (exists $errors[0]){         #print all input errors and exit
    print "\n$_" foreach @errors;
	print "\n";
    die "";
}
                         
########### load UTR values into hash of arrays %gene_UTRs: key== gene name, [0]==5' UTR length, [1]==3' UTR length.###########

my %gene_UTRs;
unless ($UTRs eq ""){
    $/ = "\n";
    while (<FINU>){
        chomp;
        my @splut = split /\t/, $_;
        $gene_UTRs{$splut[0]}[0]=$splut[1];
	$gene_UTRs{$splut[0]}[1]=$splut[2];
    }
}
$/ = "\n";


########################### open genome sequence files ##############################

opendir(FOLDER, "$dir") or die "Can't open folder $dir \n";
my @genome_files = grep {/$genome_file_term/} readdir(FOLDER);
closedir(FOLDER);

print "$_\n" foreach(@genome_files);
# example for 1 chromosome: @genome_files = ('Homo_sapiens.GRCh37.71.chromosome.22.dat', 'mm_alt_Mm_Celera_chr10.gbk');


########################### get blacklisted genes if specified ##############################

my %blacklisted_genes;
if ($external_blacklist) {
    my $blin;
    open($blin, '<', $external_blacklist) or die ("couldn\'t find file $external_blacklist. Quitting.");
    while (<$blin>) {
		chomp;
		$blacklisted_genes{$_}='y';
    }
}

my @whitelisted_genes;
if ($external_whitelist) {
    my $wlin;
    open($wlin, '<', $external_whitelist) or die ("couldn\'t find file $external_whitelist. Quitting.");
    while (<$wlin>) {
		chomp;
		push @whitelisted_genes, $_
    }
}

#################################### main loop #####################################  
foreach my $seqfile (@genome_files) {
    print "\nGetting transcript information from file: $seqfile\n";
    my $all = Bio::SeqIO->new(    -format => $format,
                              -file   => "$dir/$seqfile");
    while (my $contig = $all -> next_seq){
	my $contig_name = $contig -> display_id;
	my %contig_matches; # keys: {probe}{position}{strand} -> [0]-coord relative to start codon [1]-length [2]-sequence [3] tm1(highest at pos) [4] tm2  [5]-transcript ID [6] mRNA ID
			    # Note "position" is gDNA position where the 5' end of the oligo would sit on the transcript if the match had aligned all the way to the 5' end of the oligo. This coord is then mapped back to a genomic coordinate. The match itself may start downstream.
	my %gene_info;
	print "New sequence within in $seqfile: \"$contig_name\"\n\n";
	my $contig_seq = $contig ->seq; #hold sequence as string
	my $counter=0;
	my $UTR_count=0;
	my @sf = $contig->get_SeqFeatures;
	my @inc_sf;
	my %excluded_feats;
	my %included_feats;
	
	# Pre-process: get 'gene' info for all genes
	foreach my $feat_ob (@sf){
	    if ($feat_ob -> primary_tag eq 'gene') {
		if ($feat_ob->has_tag('gene')){
		    my @name = $feat_ob->get_tag_values('gene');
		    #print "\nGENE Nametag: $_ \t" foreach(@name);
		    if ( $feat_ob -> has_tag('note') ){
		        my @notes = $feat_ob -> get_tag_values('note');
		        @{$gene_info{$name[0]}} = @notes
		    }
		} else {
		    print "Warning: found feature object with primary tag \'gene\' but no secondary \'gene\' tag! Ignoring."
		}
	    }
	}
	
	# get transcripts with a primary tag as listed in @terms (this does not include 'gene') 
	for my $t (@terms){
	    push @inc_sf, grep { $_->primary_tag eq $t } @sf;
	}
	
	# add transcripts in whitelist
	my %already_added;
	foreach my $feat_ob (@sf){
		if ($feat_ob->has_tag('gene')){
			my @name = $feat_ob->get_tag_values('gene');
			my @transcript;
			if ($feat_ob -> has_tag('transcript_id')){
				@transcript = $feat_ob -> get_tag_values('transcript_id')  
			} elsif ($feat_ob->has_tag('note')) {
				foreach ($feat_ob->get_tag_values('note')){
				    if ($_ =~ /transcript_id\=(.+)/){push @transcript, $1};  # sometimes the transcript name is hidden away as a "note" 
				}
			}
			if($transcript[0]){ 
				for my $wl (@whitelisted_genes){
					if ($name[0] eq $wl & (! $already_added{$name[0]}) ){
						push @inc_sf, $feat_ob;
						$already_added{$name[0]} = 'y';
					}
				}
			}
		}
	}
	
	# iterate thru transcripts 
    foreach my $feat_ob (@inc_sf){  # iterate through mRNAs in contig, put in $feat_ob
	    my @name;
	    my $exclude_this_feat = ' ';
	    if ($feat_ob->has_tag('gene')){
		@name = $feat_ob->get_tag_values('gene');
	    } else{
		#print "\n Warning: encountered feature with no associated gene ID! Ignoring\n";
		$exclude_this_feat .= "Feature with no associated gene ID.  ";
	    }
	    if(exists( ${$gene_info{$name[0]}}[0] )){# compare to previously generated list of primary_tag=='gene' entries
		foreach(@{$gene_info{$name[0]}}){
		    if (($exclude_dubious eq 'yes') && (index($_, 'Dubious') == 0) && (index($_, 'verlap') != -1)) {  # Note: this may only work for S. cerevisiae 3 - *begins* with'Dubious'
		        #print "\n Warning: encountered dubious ORF feature in transcript: $name[0] overlapping another feature. Ignoring\n";
		        $exclude_this_feat .= "Dubious ORF feature overlapping another feature.  ";
		    }
		}
	    } else {
			print "\n Warning: encountered feature with no associated gene object: $name[0]\n";
			#$exclude_this_feat .= "Feature with no associated gene OBJECT!  "; # un-comment to exclude these
	    }
	    
        my @transcript;
	    if ($feat_ob -> has_tag('transcript_id')){
		@transcript = $feat_ob -> get_tag_values('transcript_id')  
	    } elsif ($feat_ob->has_tag('note')) {
                foreach ($feat_ob->get_tag_values('note')){
		    #print " NOTE: $_";
                    if ($_ =~ /transcript_id\=(.+)/){push @transcript, $1};  # sometimes the transcript name is hidden away as a "note" 
                }
            }
	    if ($transcript[0] ne $name[0]) {
		print "Warning: transcript name and gene name different. Transcript: $transcript[0]  Gene: $name[0] \n";  #<--- This is not normal for yeast genome
		#code
	    }
	    
	    
            $exclude_this_feat .= "Feature has no transcript field.  " unless (exists($transcript[0])); # double-check that feature has a 'transcript' field 
	    # $exclude_this_feat = "1" unless(exists($keywords{$name[0]}));  # SKIP EVERYTHING EXCEPT WHAT IS IN %KEYWORDS
	    
	    
	    $exclude_this_feat .= "Feature name found in externally provided blacklist.  " if(exists($blacklisted_genes{$transcript[0]}) or exists($blacklisted_genes{$name[0]})) ;
	    
	    if($exclude_this_feat ne ' '){
		#print "Excluding feature: $name[0]\n";
		if (exists($excluded_feats{$name[0]})){
		    $excluded_feats{$name[0]} = "$excluded_feats{$name[0]} \t $exclude_this_feat";
		} else {
		    $excluded_feats{$name[0]} = $exclude_this_feat; 
		}
                next;
            }
	    
	    if (exists($included_feats{$name[0]})){
		$included_feats{$name[0]} ++;
	    }else {
		$included_feats{$name[0]} = 1;
	    }
	    
	    $counter++;
	    my $strand = $feat_ob->location->strand;
	    
	    ##### get UTRs #####
	    my $UTR_5p = $fpdf;
            my $UTR_3p = $tpdf;
            if (exists($gene_UTRs{$name[0]})){ #override defaults if UTR lengths for GENE name specified in UTR table file
		$UTR_5p = $gene_UTRs{$name[0]}[0];
                $UTR_3p = $gene_UTRs{$name[0]}[1];
		$UTR_count ++;
            }
            my $UTRseq_5p="";
            my $UTRseq_3p="";
	    if ($UTR_5p > 0 or $UTR_3p >0){
		my $start = $feat_ob->location->start;       
		my $end = $feat_ob->location->end;
		if ($strand == -1){
		    $UTRseq_5p = rev_comp(substr ($contig_seq, $end, $UTR_5p));
		    $UTRseq_3p = rev_comp(substr ($contig_seq, $start-1-$UTR_3p, $UTR_3p));   
		} else {
		    $UTRseq_5p = substr ($contig_seq, $start-1-$UTR_5p, $UTR_5p);
		    $UTRseq_3p = substr ($contig_seq, $end, $UTR_3p);
		}
            
            }
	    	    
	    my (@exons) = get_exons($feat_ob, $contig_seq);
	    my $feat_seq = pop @exons; 
	    $feat_seq = $UTRseq_5p.$feat_seq.$UTRseq_3p;  	#add UTRs to $feat_seq
	    
            
	    if ($counter % 50 == 0){ # update display every n transcripts                                                                                                    ";
		#print "\nTranscripts analysed on chromosome: $counter  \t   of which $UTR_count had user-defined UTRs          "; # Gene: $name[0] \ttranscript: $transcript[0] \texons: ".scalar @{$exons[0]}. "  \t strand: $strand\t                          ";
	    }
	    my $orf_seqobj = Bio::Seq->new( -seq => $feat_seq,
                                 -id  => $transcript[0] );
    
	    $seq_out->write_seq($orf_seqobj);    
        }
	print "\r                                                                                                      ";
	print "\rTranscripts analysed on chromosome: $counter  \t   of which $UTR_count had user-defined UTRs         \n";
	print "Excluded feature under gene $_ : $excluded_feats{$_} \n" foreach(keys(%excluded_feats));
    } #next in file
} # next file

print "\n\n";

################################################################################
################################# END OF MAIN ##################################
################################################################################


########################### reverse complement ##############################

sub rev_comp {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}

########################### get exons ##############################

sub get_exons { # make a 2-element array in @exons for converting spliced coords to genomic coords,
                # $exons[0][$i]= for exon i: end in spliced coords and
                # $exons[1][$i]= for exon i: cumulative length of preceeding introns plus the start coordinate for + strand,
                # OR for - strand, the end coordinate (ie the start from the transcript's perspective)
                # MINUS the cumulative lengths of preceeding introns
    my $featob = shift;
    my $chromosome_seq = shift;
    my @exonobs;
    my @exons;
    my $spliced_seqstr;
    my $strand = $featob -> location -> strand;
    
    if($featob->location->isa('Bio::Location::SplitLocationI')) {
		@exonobs = $featob->location->sub_Location;
		my $prev_ex_right;
		for(my $i = 0; $i <= $#exonobs; $i++) {
			my $ex_left = $exonobs[$i]->start; # increment cumulative spliced exon end pos by current exon size
			my $ex_right = $exonobs[$i]->end;
			my $ex_len = ($ex_right-$ex_left)+1;
			my $exseq = substr $chromosome_seq, $ex_left-1, $ex_len;
			$exseq = rev_comp($exseq) if ($exonobs[$i]->strand eq '-1');
			
			if ($i==0){
				#$prev_ex_right = $exonobs[0]->end; # in this case $prev_ex_end is for the current exon
				#my $ex_left = $exonobs[0]->start;
				$prev_ex_right = $ex_right; # in the case of exon 1, $prev_ex_end is for the current exon
				$spliced_seqstr = $exseq;
				$exons[0][0] = ($ex_len-1); #exon 1 length for spliced coord of exon 1 end, not affected by strand
				if ($strand eq '1'){
					$exons[1][0] = $ex_left;
				} else{
					$exons[1][0] = $ex_right; # exon end is actually start of last exon from transcript's perspective
					$prev_ex_right=$ex_left;
				}
			}else{
				#my $ex_left = $exonobs[$i]->start; # increment cumulative spliced exon end pos by current exon size
				#my $ex_right = $exonobs[$i]->end; 
				
				#my $exseq = substr $chromosome_seq, $ex_left-1, ($ex_right-$ex_left)+1;
				#$exseq = rev_comp($exseq) if ($exonobs[$i]->strand eq '-1');
				
				$exons[0][$i] = $exons[0][$i-1] + ($ex_len-1);  # add new exon length to cumulative spliced length
				if ($strand eq '1'){
					$spliced_seqstr .= $exseq;
					$exons[1][$i] = $exons[1][$i-1] -1 + $ex_left - $prev_ex_right;
					$prev_ex_right = $ex_right;  # save value for next time to reduce computation load
				} else{ 
					$spliced_seqstr = $exseq.$spliced_seqstr;  # add it on at the START of string; this can make it consistent with Bioperl for mm_alt_Mm_Celera_chr*.gbk
					$exons[1][$i] = $exons[1][$i-1] +1 - abs($prev_ex_right - $ex_right); # $prev_exon_end is LEFT side of previous exon if strand = -1;  abs($prev_ex_end - $ex_right) is the length of the preceding intron from transcript's prespective 
					$prev_ex_right=$ex_left;
				}
	
			}
		}
    }
    else { # only 1 exon
	my $ex_left = $featob->location->start;
	my $ex_right = $featob->location->end;
	$spliced_seqstr = substr $chromosome_seq, $ex_left-1, ($ex_right-$ex_left)+1;
	$spliced_seqstr = rev_comp($spliced_seqstr) if ($strand eq '-1');
	$exons[0][0] = $ex_right - $ex_left; #exon 1 length for spliced coord of exon 1 end, not affected by strand
	if ($strand eq '1'){
		$exons[1][0] = $ex_left;
	} else{
		$exons[1][0] = $ex_right; # exon end is actually start of first exon from transcript's perspective
	}   
	    
    }
    push @exons, $spliced_seqstr;
    return @exons;
}

