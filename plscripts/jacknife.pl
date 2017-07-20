
use strict;

my $inlist = shift;
my $nrep = shift;
my $nfold = shift;
my $basename = shift;

for (my $rep=0; $rep<$nrep; $rep++)	{
	open (INLIST, $inlist) or die "input error\n";

	my @genelist = (<INLIST>);
	open (OUTFILE, '>'.$basename.$rep) or die "output error\n";
	print OUTFILE "$nfold\n";
	for (my $i=0; $i<$nfold; $i++)	{
		my $size = (@genelist);
		my $choose = int(rand($size));
		if ($choose >= $size)	{
			die "error in random choice\n";
		}
		print OUTFILE $genelist[$choose];
		splice(@genelist,$choose,1);
	}
}


