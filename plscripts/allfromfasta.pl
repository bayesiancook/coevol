
use strict;
use Bio::SeqIO;

my $genefile = shift;
my $taxfile = shift;
my $base = shift;

open (GENEFILE, $genefile) or die "input error: $genefile\n";
open (TAXFILE, $taxfile) or die "intput error: $taxfile\n";

my @tax = <TAXFILE>;
my $ntaxa = (@tax);

foreach my $gene (<GENEFILE>)	{

	chomp $gene;
	my $seqio = Bio::SeqIO->new(-file => $gene);

	my %ali;
	my $l = 0;
	my %tax2seq;
	while (my $seq = $seqio->next_seq)	{
		if (! $l)	{
			$l = $seq->length;
		}
		else	{
			if ($l != $seq->length)	{
				my $ll = $seq->length;
				die "error : non matching lengths : $l $ll\n";
			}
		}
		$tax2seq{$seq->display_id} = $seq->seq;
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

			
