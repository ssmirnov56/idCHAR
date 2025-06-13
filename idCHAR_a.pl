#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Text;
##use Data::Dumper;


########################################################

my ($help, $verbose, $debug) = (0,0,0);

my ($fasta_file, $charge, $aromatics, $prolines) = ('', 0,0,0);

my $rh_charge = {"E"=>-1, "D"=>-1, "R"=>1, "K"=>1};
my $rh_aromatics = {"W"=>1, "Y"=>1, "F"=>1};
my $rh_prolines = {"P"=>1};
my $min_segment = 1; 

my $args = {
	"h" => \$help,
	"v" => \$verbose,
	"d" => \$debug,
	"fasta:s"  => \$fasta_file,
	"charge"   => \$charge,
   	"aroma"    => \$aromatics,
	"pro"	   => \$prolines,
};

GetOptions(%$args);


my $msg2usr='';

$msg2usr .= "\nProvide fasta filename; -fasta" if not $fasta_file;
$msg2usr .= "\nProvide an analysis mode: -charge or -aroma or -pro" if not ($charge or $aromatics or $prolines);
$msg2usr .= "\nProvide a single analysis mode only: -charge or -aroma or -pro" if ( ($charge and $aromatics) or ($charge and $prolines) or ($aromatics and $prolines) );

if ($msg2usr or $help) {

warn "$msg2usr\n\n" if $msg2usr and not $help;

print "version D. Processing or separtion of prolines added.\n";
print "version C. Two reporting bugs fixed: for residue position and for single-digit % values.\n";
print "version B. Percent-score is added for the output: score/#aa\n";
print "version A. Processes partitioning by charge (KR vs ED) and by aromatics (FYW)\n";
print "\nPurpose:\n";
print "Finds the point in aa sequence which conveys the greatest partitioning (difference) in target feature of aa composition.\n";
print "\nUsage:\n";

print "\t$0 [-h] [-v] [-d] -fasta {fasta_filename} [-charge] [-aroma] [-pro]\n";
print "Params:\n";
print "\t-fasta:\t\tProtein sequence fasta filename;\n";
print "\t-charge:\tSequence partitioning by charge: Lys+Arg vs Asp+Glu content;\n";
print "\t-aroma:\tSequence partitioning by aromatics content: Trp+Tyr+Phe content vs. their absence;\n"; 
print "\t-pro:\tSequence partitioning by proline content: Pro vs its absence;\n";

print "(C):\n";
print "\tSerge Smirnov, Western Washington University, smirnos\@wwu.edu, 2020\n";
exit(0) if not $debug; 
}

report($args) if $verbose;


die "Fhheewww, made it here !" if (0);

#################################################
#
# M A I N    
#
#################################################


###########################################
##     Reading in the input files        ##
###########################################

my (@arr_fasta_seq, @arr_fasta_header, @arr_sites);



## Read in the one-letter aa seqeunce as $fasta_sequence. 

open (FASTA, $fasta_file) or die "could not open $fasta_file: $!";

my $curr_fasta_seq = -1;

while (<FASTA>) {
    chomp;


    if (/^\>/) {

    	push @arr_fasta_header, $_;

    	$curr_fasta_seq++;

    	$arr_fasta_seq[$curr_fasta_seq]='';

    	next
    } else {

        s/\s//g;

    	my $seq = $arr_fasta_seq[$curr_fasta_seq].$_;

	    $arr_fasta_seq[$curr_fasta_seq] = $seq;
    }

}
close (FASTA);


my $rh_property;

$rh_property = $rh_charge if ($charge);
$rh_property = $rh_aromatics if ($aromatics);
$rh_property = $rh_prolines if ($prolines);


## Run through ALL FASTA entries
my $fastaNum = scalar (@arr_fasta_header);


for my $seqN ( 0..$fastaNum-1 ) {

## Convert $fasta_sequence into @arr_sequence and get its size. 

	my ($fasta_sequence, $sequence_description) = ($arr_fasta_seq[$seqN], $arr_fasta_header[$seqN]); 

	my @arr_sequence = split("", $fasta_sequence); 

	my $sequence_size = scalar (@arr_sequence);

	next if ($sequence_size<2*$min_segment);
		
	my $rh_score;
## Run through the sequence and check for site pattern matches given the specified number of mutations

my ($max_score, $percent_score, $percentN, $percentC) = (-100, 0, 0, 0);


## Calc the score difference and determine the max_score value
for my $i ($min_segment .. $sequence_size-1-$min_segment) {

#Calc score for the N-terminus, positions 0 to $i-1

	my $scoreN = score ($fasta_sequence, 0, $i-1);

#Calc score for the C-terminus, positions $i+1 to $sequence_size-1

	my $scoreC = score ($fasta_sequence, $i+1, $sequence_size-1);

#Calc score difference and assign it to position $i


	my $diff = abs($scoreN-$scoreC);
	$rh_score->{$i} = $diff;

	if ($diff>$max_score) {
		$max_score = $diff;
		
		$percent_score = $max_score / $sequence_size * 100.0; 
		$percent_score =~ /(^\d\d*)/; $percent_score = $1;
		
		$percentN = abs ($scoreN / ($i+1) * 100.0);
		$percentN =~ /(^\d\d*)/; $percentN = $1;

		$percentC = abs ($scoreC / ($sequence_size - $i + 1) * 100.0);
		$percentC =~ /(^\d\d*)/; $percentC = $1;
	}
}
##        warn Data::Dumper->new([$rh_score])->Dump,"\n" if $debug;


## Report residues with greatest difference score = max_score value

	print "\n$arr_fasta_header[$seqN]\n";
	print "Max score difference: $max_score\;\t%score= $percent_score\t%scoreN= $percentN\t%scoreC=$percentC\t\tResidues:\t";

for my $i ($min_segment .. $sequence_size-1-$min_segment) {

	my $diff = $rh_score->{$i}; 
	if ( $diff == $max_score) {

		print $arr_sequence[$i],$i+1," ";
	}

}
	print "\n\n";
}


##print "hProt\n" if $verbose;

##print Data::Dumper->new([\%hProt])->Dump,"\n" if $verbose;

##warn Data::Dumper->new([$aliph->{'VAL'}])->Dump,"\n";


##print Data::Dumper->new([$rh])->Dump,"\n" if 1;

##warn Data::Dumper->new([\%best_pos])->Dump,"\n" if $verbose;
##warn Data::Dumper->new([\%best_misses])->Dump,"\n" if $verbose;

##warn Data::Dumper->new([\%BMRB_CS])->Dump,"\n";
##warn Data::Dumper->new([\%BMRB_SD])->Dump,"\n";


##print STDERR "Delta\tMisses\tOne\tTwo\tMany\tAll\n";
##print "$peak_delta\t$misses\t$ones\t$twos\t$many\t$total\n";





#################################################
#
# S U B s
#
#################################################


sub report { 
	my ($rh_args) = @_;
	
	foreach my $key (keys %$rh_args) {

		print "$key => ",${$rh_args->{$key}}, "\n";
	}
}

sub score {
	my ($sequence, $start, $stop)= @_;

	my @arr_sequence = split("", $sequence);

	my $score = 0;

	for my $j ($start .. $stop) {

		my $aa_score=0;
		
		$aa_score = $rh_property->{$arr_sequence[$j]} if ( defined $rh_property->{$arr_sequence[$j]} );

		$score+=$aa_score;
	}

	return $score
}
