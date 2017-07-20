use strict;

my $infile = shift;

open (INFILE, $infile) or die "input error\n";

my $count = 0;
my $mean = 0;
my $var = 0;
foreach my $line (<INFILE>)	{

	chomp $line;
	my @a = split('\t',$line);
	my $tmp1 = $a[2];
	my $tmp2 = $a[7];

	my $tmp = $tmp2 - $tmp1;

	$mean += $tmp;
	$var += $tmp * $tmp;
	$count++;

}

$mean /= $count;
$var /= $count;
$var -= $mean * $mean;

my $dist = sqrt($var);

print "mean quadratic distance : $dist\n";

	
	
