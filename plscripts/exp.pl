use strict;

my $infile = shift;
my $outfile = shift;

open (INFILE, $infile) or die "input error: $infile\n";
open (OUTFILE, '>'.$outfile) or die "output error: $outfile\n";

my $header = <INFILE>;
print OUTFILE $header;

foreach my $line (<INFILE>)	{
	chomp $line;
	my @a = split(/\t/,$line);
	my $name = $a[0];
	my $n = (@a);
	print OUTFILE "$name";
	for (my $i=1; $i<$n; $i++)	{
		my $val = $a[$i];
		my $expval = exp($val);
		if ($val == -1)	{
			$expval = -1;
		}
		print OUTFILE "\t$expval";
	}
	print OUTFILE "\n";
		
	
	## (my $name, my $val) = split(/\s+/,$line);
	## my $expval = exp($val);
	## if ($val == -1)	{
	## 	$expval = -1;
	## }
	## print OUTFILE "$name\t$expval\n";
}

