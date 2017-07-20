
use strict;

my $file1 = shift;
my $file2 = shift;
my $outfile = shift;

my %tab1;

open (FILE1, $file1) or die "input error\n";
open (FILE2, $file2) or die "input error\n";

my $h1 = <FILE1>;
chomp $h1;
my @t1 = split(/\s/,$h1);
(my $ntax1, my $ncont1) = (@t1);

my $h2 = <FILE2>;
chomp $h2;
my @t2 = split(/\s/,$h2);
(my $ntax2, my $ncont2) = (@t2);
print "$ncont1\t$ncont2\n";


my $ncont = $ncont1 + $ncont2;

my @list1 = <FILE1>;

foreach my $line (@list1)	{
	chomp $line;
	my @t = split(/\s+/,$line);
	my $dim = (@t);
	if ($dim == $ncont1 + 1)	{
		my $tax = shift(@t);
		my $cont = join("\t",@t);
		$tab1{$tax} = $cont;
	}
}

my @list2 = <FILE2>;

my %inbothfiles;

foreach my $line (@list2)	{
	chomp $line;
	my @t = split(/\s+/,$line);
	my $dim = (@t);
	if ($dim == $ncont2 + 1)	{
		my $tax = shift(@t);
		if (exists $tab1{$tax})	{
			my $cont = join("\t",@t);
			$tab1{$tax} .= "\t".$cont;
			$inbothfiles{$tax} = 1;
		}
	}
}

my $ntax = (keys %inbothfiles);
	
open(OUTFILE,'>'.$outfile) or die "output error\n";

print OUTFILE "$ntax\t$ncont\n";
foreach my $tax (keys %inbothfiles)	{
	print OUTFILE "$tax\t$tab1{$tax}\n";
}


