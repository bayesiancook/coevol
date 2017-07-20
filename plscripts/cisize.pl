use strict;

my $infile = shift;

open (INFILE, $infile) or die "input error\n";

my $mean = 0;
my $count = 0;
foreach my $line (<INFILE>)	{

	chomp $line;
	my @a = split('\t',$line);
	my $tmp = $a[4] - $a[3];
	$mean += $tmp;
	$count++;
}

$mean /= $count;

print "mean CI size : $mean\n";

	
