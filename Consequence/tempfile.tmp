#!/usr/bin/perl
use strict;

#parameter checking
my $lenArgv = scalar @ARGV;
unless( $lenArgv >1 && $lenArgv<13 ){
	print "Error: Wrong number of parameters\n";
	&usage();
	exit 1;
}
	
if($lenArgv%2==1){
	print "Error: Wrong pair of options\n";
	&usage();
	exit 2;
}

open IN,"$ARGV[0]" or die "Can not find the file $ARGV[0]";
open OUT,">$ARGV[1]" or die "Can not open the output file $ARGV[1]";
my $unique = 1;
my $maxLength = 10000;
my $minLength = 1;
my $miscleavage = 0;
my $outputFormat = "fasta";

for (my $i=2;$i<$lenArgv;$i+=2){
	if(lc($ARGV[$i]) eq "-unique"){
		if(lc($ARGV[$i+1]) eq "yes"){
			$unique = 1;
		}elsif(lc($ARGV[$i+1]) eq "no"){
			$unique = 0;
		}else{
			print "unrecognized option value for -unique: $ARGV[$i+1]\n";
			&usage();
			exit 4;
		}
	}elsif (lc($ARGV[$i]) eq "-maxlength"){
		$maxLength = $ARGV[$i+1];
		checkIntegers($maxLength,"maxLength");
	}elsif (lc($ARGV[$i]) eq "-minlength"){
		$minLength = $ARGV[$i+1];
		checkIntegers($minLength,"minLength");
	}elsif (lc($ARGV[$i]) eq "-miscleavage"){
		$miscleavage = $ARGV[$i+1];
		checkIntegers($miscleavage,"miscleavage");
	}elsif(lc($ARGV[$i]) eq "-outputformat"){
		$outputFormat = lc($ARGV[$i+1]);
		unless ($outputFormat eq "fasta" || $outputFormat eq "tsv"){
			print "Error: Unrecognized value for parameter outputFormat which only allows fasta or tsv";
			exit 7;
		}
	}else{
		print "Error: Unrecognized option: $ARGV[$i]\n";
		&usage();
		exit 3;
	}
}

if($minLength > $maxLength){
	print "Error: The minimum length is greater than the maximum length\n";
	exit 5;
}

print OUT "Protein\tPeptide\n" if ($outputFormat eq "tsv");

#real business
my %peptides;
my %peptidesCount;
my $id;
my $seq="";
my $countProtein=0;

while(my $line=<IN>){
	chomp $line;
	if($line=~/^>/){
		unless(length $seq==0){
			&digest($id,$seq);
			$countProtein++;
			$seq="";
		}
		$id=substr($line,1);
	}else{
        	$seq .= $line; # add sequence
	}
}
unless(length $seq==0){
	&digest($id,$seq);
	$countProtein++;
}

exit if ($outputFormat eq "tsv");

my @uniquePep = sort {
	if((length $a)<(length $b)) {return -1;}
	elsif((length $a)>(length $b)) {return 1;}
	else {return $a cmp $b}
}  keys %peptides;

my %distri;
foreach my $pep(@uniquePep){
	my $len = length $pep;
	last if ($len > $maxLength);
	next if ($len < $minLength);
	$distri{$len}++;
	next if($peptidesCount{$pep}>1 && $unique==1);
	print OUT ">$peptides{$pep}\n$pep\n";
}
close OUT;
print "Peptides distribution between the limit $minLength-$maxLength:\n";
my @lens = sort {$a<=>$b} keys %distri;
foreach (@lens){
	print "$_\t$distri{$_}\n";
}

sub checkIntegers(){
	my ($value, $name) = @_;
	unless ($value=~/^\d+$/){
		print "The value $value given for the parameter $name is not valid, only integers allowed.\n\n";
		&usage();
		exit 6;
	}
}

sub digest(){
	my ($header,$seq) = @_;
	$header =~ s/\s.*//;#only keeps the part before the first white space
	$seq=~s/\s+//g;
	my @frags = split /(?<=[KR])(?=[^P])/,$seq;
	my $len = scalar @frags;
	for (my $mis = 0; $mis <= $miscleavage; $mis++){
		for(my $i=1;$i<=$len-$mis;$i++){
			my $pep = $frags[$i-1];
			if ($mis == 0 && $outputFormat eq "tsv"){
				my $len = length $pep;
				print OUT "$header\t$pep\n" if($minLength <= $len && $len <= $maxLength);
			}
			for (my $curr=1;$curr<=$mis;$curr++){
				$pep .= $frags[$i-1+$curr];
			}
			my $pepHeader = "$header-mis$mis-pep$i";
#			my $pepHeader = "$header-pep$i";
			if(exists $peptides{$pep}){ #not unique
				$peptidesCount{$pep}++;
				$peptides{$pep} .= ";$pepHeader";
			}else{
				$peptidesCount{$pep} = 1;
				$peptides{$pep} = "$pepHeader";
			}
		}
	}
}

sub usage(){
	print "This script is to theoritically digest proteins into tryptic peptides\n";
	print "Usage: perl digest.pl <protein fasta file> <output peptide file> [options]\n";
	print "Options include:\n";
	print "-unique [yes|no]  whether only output proteotypic peptides. The default value is no\n";
	print "-maxLength integer to limit the maximum length of peptide, default is 10000 which means include every peptide.\n";
	print "-minLength integer to limit the minimum length of peptide, default is 1 which means include every peptide.\n";
	print "-miscleavage integer to indicate the allowed maximum cleavage, default is 0 which means no miscleavage allowed.\n";
	print "-outputFormat determines what file format the ouput file will be, either fasta (default) or tsv.\n";
	print "Example 1: perl digest.pl uniprot_sprot_human.fasta\t\t this will digest the human swissprot data and output all sequences\n";
	print "Example 2: perl digest.pl uniprot_sprot_human.fasta -unique yes -minLength 4 -maxLength 24\t\t this will only output proteotypic peptides having AA between 4 and 24 inclusive, which is more likely to be observed on a MS machine\n";
	print "Example 3: perl digest.pl uniprot_sprot_human.fasta -miscleavage 1\t\t this will allow maximum one miscleavage.\n";
}