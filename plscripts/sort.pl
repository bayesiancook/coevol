use strict;

my $infile = shift;
my $outfile = $infile . '.sorted';
my $index = shift;

open (INFILE, $infile) or die "input error:$infile\n";

open(OUTFILE , '>'.$outfile) or die "output error:$outfile\n";

my %hash;

foreach my $line (<INFILE>)	{

	chomp $line;
	my @a = split('\t', $line);
	$hash{$line} = $a[$index];
}

foreach my $line (sort {$hash{$a} <=> $hash{$b}} keys %hash)	{
	print OUTFILE "$line\n";
}

	
