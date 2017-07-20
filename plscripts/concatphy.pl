use strict;


my $listfile = shift;
my $outfile = shift;

open(INFILE, $listfile) or die "error when opening input file\n";

my @list = <INFILE>;

my $count = 0;
my %seqhash;

my $currentlength = 0;

foreach my $filename (@list)	{

	open(INFILE, $filename) or die "input error: $filename\n";
	my $header = <INFILE>;
	foreach my $line (<INFILE>)	{
		chomp $line;
		my $name = "";
		my $seq = "";
		if ($line !~ /^\s*$/g)	{
			if ($line =~ /^(\S+)\s+(\S+)$/)	{
				$name = $1;
				$seq = $2;
			}
			else	{
				die "error: $line\n";
			}
			if (! $count)	{
				if (exists $seqhash{$name})	{
					die "error: found $name twice\n";
				}
				$seqhash{$name} = $seq;
			}
			else	{
				if (!exists $seqhash{$name})	{
					my $seq = "";
					for (my $i=0; $i<$currentlength; $i++)	{
						$seq .= '?';
					}
					$seqhash{$name} = $seq;
				}
				$seqhash{$name} .= $seq;
			}
		}
	}
	my $l = 0;
	foreach my $name (keys %seqhash)	{
		my $seq = $seqhash{$name};
		if ($l < length($seq))	{
			$l = length($seq);
		}
	}
	$currentlength = $l;
	foreach my $name (keys %seqhash)	{
		my $seq = $seqhash{$name};
		for (my $i=length($seq); $i<$l; $i++)	{
			$seq .= '?';
		}
		$seqhash{$name} = $seq;
	}
	$count++;
}

my $len = 0;

my $ntaxa = 0;
my $nex = 0;
foreach my $name (keys %seqhash)	{
	my $newlen =  length($seqhash{$name});
	if (! $len)	{
		$len =  $newlen;
	}
	else	{
		if ($len != $newlen)	{
			print "error with $name: non matching length : $len / $newlen\n";
			$nex++;
		}
	}
	$ntaxa ++;
}

print "total length of sequences : $len\n";

open(OUTFILE, '>'.$outfile) or die "error when opening output file\n";
$ntaxa -= $nex;
print OUTFILE "$ntaxa\t$len\n";
foreach my $name (keys %seqhash)	{
	my $newlen =  length($seqhash{$name});
	if ($len == $newlen)	{
		print OUTFILE "$name\t$seqhash{$name}\n";
	}
}

