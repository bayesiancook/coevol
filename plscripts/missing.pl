
use strict;

my $infile = shift;

open (INFILE , $infile) or die "input error\n";

my $ntaxa;
my $nsite;

my $header = <INFILE>;
my @t = split(/\s/,$header);
my $dim = (@t);
if ($dim != 2)	{
	die "error with header : $header (@t)\n";
}
($ntaxa,$nsite) = (@t);
print "ntaxa : $ntaxa\n";
print "nsite : $nsite\n";
print "\n";

for (my $i=0; $i<$ntaxa; $i++)	{
	my $line = <INFILE>;
	chomp $line;
	my @t = split(/\s+/, $line);
	(my $tax, my $seq) = (@t);
	#my $tax = <INFILE>;
	#chomp $tax;
	#my $seq = <INFILE>;
	if (length($seq) != $nsite)	{
		die "error : non matching sequence length : " , length($seq), " instead of $nsite\n";
	}
	my $count = ($seq =~ tr/\?//);
	my $count += ($seq =~ tr/X//);
	my $count += ($seq =~ tr/N//);
	my $count += ($seq =~ tr/-//);
	my $mis = int (100 * $count / $nsite);
	print "$tax\t$mis\n";
}

