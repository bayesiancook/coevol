use strict;

my $listfile = shift;
my $taxfile = shift;

open (LISTFILE, $listfile) or die "input error\n";
open (TAXFILE, $taxfile.'.tax') or die "input error\n";

my %taxmap;
my $ntax = 0;
my $givenntax = <TAXFILE>;
foreach my $tax (<TAXFILE>)	{
	chomp $tax;
	$taxmap{$tax} = 1;
	$ntax++;
}
if ($givenntax != $ntax)	{
	die "error : non matching number of taxa in tax file\n";
}

my $ngene = <LISTFILE>;
foreach my $genefile (<LISTFILE>)	{
	chomp $genefile;
	open(ALIFILE, $genefile) or die "input error: $genefile\n";
	my $header = <ALIFILE>;
	my @a = split(/\s+/,$header);
	my $l = (@a);
	if ($l != 2)	{
		die "error when parsing header\n";
	}
	my $nsite = $a[1];
	open(OUTFILE, '>'.$taxfile.$genefile) or die "output error\n";
	print OUTFILE "$ntax\t$nsite\n";
	my $checkntax = 0;
	foreach my $line (<ALIFILE>)	{
		chomp $line;
		if ($line !~ /^\s*$/)	{
		my @a = split(/\s+/,$line);
		my $l = (@a);
		if ($l != 2)	{
			print "error when parsing sequence\n";
			print "$line\n";
			print "from $genefile\n";
			die;
		}
		my $tax = $a[0];
		my $seq = $a[1];
		if (exists $taxmap{$tax})	{
			print OUTFILE "$tax\t$seq\n";
			$checkntax++;
		}
		}
	}
	if ($ntax != $checkntax)	{
		die "error: non matching number of taxa: $checkntax instead of $ntax\n";
	}
}

	
