
use strict;

my $infile = shift;
my $outfile = shift;

open (INFILE, $infile) or die "input error\n";
open (OUTFILE, '>',$outfile) or die "output error\n";

my $name ;
my $seq;

my %name2seq;

my $ntax = 0;
my $nsite = 0;

foreach my $line (<INFILE>)	{

	chomp $line;
	if ($line)	{
		if ($line =~ m/^>(.*)$/)	{
			# new sequence
			if ($name)	{
				$name2seq{$name} = $seq;
				if (! $nsite)	{
					$nsite = length($seq);
				}
				else	{
					if ($nsite != length($seq))	{
						die "error : " , len($seq) , " instead of $nsite\n";
					}
				}
				$seq = "";
				print "$1\n";
				$ntax++;
				
			}
			$name = $1;
		}
		else {
			$seq = $seq . $line;
		}
	}
}

$name2seq{$name} = $seq;
$ntax++;

print OUTFILE "$ntax\t$nsite\n";

foreach my $name (keys %name2seq)	{
	my $seq = $name2seq{$name};
	print OUTFILE "$name\t$seq\n";
}
