
use strict;

my $genefile = shift;
my $taxfile = shift;
my $base = shift;

open (GENEFILE, $genefile) or die "input error: $genefile\n";
open (TAXFILE, $taxfile) or die "intput error: $taxfile\n";

my @tax = <TAXFILE>;

my $ntaxa = (@tax);

my @genelist = <GENEFILE>;
foreach my $gene (@genelist)	{

	chomp $gene;
	open (INFILE, $gene) or die "input error: $gene\n";
	my $header = <INFILE>;
	my @ali = <INFILE>;
	my %tax2seq;
	my $l = 0;
	foreach my $line (@ali)	{
		chomp $line;
		my @t = split(/\t/,$line);
		my $n = (@t);
		if ($n != 2)	{
			print "$n\n";
			print "$line\n";
			print @t , "\n";
			die "error when splitting\n";
		}
		#(my $name,my $seq) = split($line,/\s*/);

		my $name = $t[0];
		my $seq = $t[1];
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

			
