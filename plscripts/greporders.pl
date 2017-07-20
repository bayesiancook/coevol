use strict;

my $infile = shift;
my $cladefile = shift;

open (INFILE, $infile) or die "input error\n";
open (CLADEFILE, $cladefile) or die "input error\n";

my %list;
foreach my $line (<CLADEFILE>)	{
	chomp $line;
	my @a = split('\t',$line);
	my $name = $a[0];
	my $left = $a[1];
	my $right = $a[2];

	my $leftright = $left . '\t' . $right;

	$list{$leftright} = $name;
}

foreach my $line (<INFILE>)	{

	chomp $line;
	my @a = split('\t',$line);
	my $left = $a[0];
	my $right = $a[1];
	my $leftright = $left . '\t' . $right;

	if (exists $list{$leftright})	{
		print "$list{$leftright}\t$a[2]\t$a[3]\t$a[4]\n";
	}
}

	
		
