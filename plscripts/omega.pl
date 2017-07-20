use strict;

my $infile = shift;
my $burnin = shift;
my $cutoff = shift;

open (INFILE, $infile) or die "input error: $infile\n";

my $line;
for (my $i=0; $i<$burnin; $i++)	{
	$line = <INFILE>;
}

chomp $line;
my @a = split('\t',$line);
my $nsite = (@a);

my $count = 0;
my $possel = 0;
my @mean = ();
my @pp = ();
for (my $i=0; $i<$nsite; $i++)	{
	$mean[$i] = 0;
	$pp[$i] = 0;
}

foreach my $line (<INFILE>)	{

	chomp $line;
	my @a = split('\t',$line);
	for (my $i=0; $i<$nsite; $i++)	{
		my $tmp = $a[$i];
		$mean[$i] += $tmp;
		if ($tmp > 1)	{
			$pp[$i] ++;
			$possel++;
		}
	}
	$count++;
}

my $nsupp = 0;
for (my $i=0; $i<$nsite; $i++)	{
	$pp[$i] /= $count;
	$mean[$i] /= $count;
	if ($pp[$i] > $cutoff)	{
		$nsupp++;
	}
}

$possel /= $count * $nsite;
print "proportion above 1 : $possel\n";
print "number of sites significantly above 1 : $nsupp out of $nsite\n";

 
		

