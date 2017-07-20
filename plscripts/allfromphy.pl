
use strict;

my $genefile = shift;
my $taxfile = shift;
my $base = shift;

open (GENEFILE, $genefile) or die "input error: $genefile\n";
open (TAXFILE, $taxfile) or die "intput error: $taxfile\n";

my $ntax = <TAXFILE>;
my @tax = <TAXFILE>;

my $ntaxa = (@tax);
if ($ntax != $ntaxa)	{
	die "error: non matching number of taxa\n";
}

my $ngene = <GENEFILE>;
chomp $ngene;
my @genelist = <GENEFILE>;

open (OUTLIST, '>'.$base.$genefile) or die "output error\n";
print OUTLIST "$ngene\n";

foreach my $gene (@genelist)	{

	chomp $gene;
	open (INFILE, $gene) or die "input error: $gene\n";
	my $header = <INFILE>;
	my @ali = <INFILE>;
	my %tax2seq;
	my $l = 0;
	foreach my $line (@ali)	{
		chomp $line;
		if ($line =~ /^(\w+)(\s+)([ACGT\-\?]+)$/gi)	{
			my $name = $1;
			my $seq = $3;
			if (! $l)	{
				$l = length($seq);
			}
			else	{
				if ($l != length($seq))	{
					my $ll = length($seq);
					die "error : non matching lengths : $l $ll\n";
				}
			}
			$tax2seq{$name} = $seq;
		}
		else	{
			if ($line !~ /^\s*$/)	{
				die "does not recognize\n$line\n";
			}
		}
	}

	print OUTLIST "$base$gene\n";

	open (OUTFILE, '>'.$base.$gene) or die "output error : $base$gene\n";
	print OUTFILE "$ntaxa\t$l\n";
	foreach my $species (@tax)	{
		chomp $species;
		print OUTFILE "$species\t";
		if (! exists $tax2seq{$species})	{
			for (my $i=0; $i<$l; $i++)	{
				print OUTFILE "?";
			}
			print OUTFILE "\n";
		}
		else	{
			print OUTFILE "$tax2seq{$species}\n";
		}
	}
}

			
