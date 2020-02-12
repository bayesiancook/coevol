
#include <iostream>
#include <cstdlib>
using namespace std;

#include "CodonSequenceAlignment.h"
#include "Exception.h"
#include "Random.h"

CodonSequenceAlignment::CodonSequenceAlignment(SequenceAlignment* from, bool force_stops, GeneticCodeType type)	{

		try	{
			DNAsource = from;

			if (from->Nsite % 3 != 0)	{
				cerr << "not multiple of three\n";
				exit(1);
			}
			Nsite = from->Nsite/3;
			Ntaxa = from->Ntaxa;
			CodonStateSpace* tempstatespace = new CodonStateSpace(type);
			statespace = tempstatespace;
			DNAStateSpace* nucspace = tempstatespace->GetDNAStateSpace();

			taxset = DNAsource->GetTaxonSet();

			// make my own arrays
			// make translation
			Data = new int*[Ntaxa];
			for (int i=0; i<Ntaxa; i++)	{
				Data[i] = new int[Nsite];
				for (int j=0; j<Nsite; j++)	{
					try {
						Data[i][j] = GetCodonStateSpace()->GetCodonFromDNA(DNAsource->GetState(i, 3*j), DNAsource->GetState(i, 3*j+1), DNAsource->GetState(i, 3*j+2));
						if (Data[i][j] == -1)	{
							if ((DNAsource->GetState(i, 3*j) != -1) && (DNAsource->GetState(i, 3*j+1) != -1) && (DNAsource->GetState(i, 3*j+2) != -1))	{
								cerr << "stop codon: taxon " << taxset->GetTaxon(i) << " and codon " << j+1 << " (site " << 3*j+1 << ") :";
								cerr << nucspace->GetState(DNAsource->GetState(i, 3*j)) <<  nucspace->GetState(DNAsource->GetState(i, 3*j+1)) << nucspace->GetState(DNAsource->GetState(i, 3*j+2)) << '\n';
							}
						}
					}
					catch(...)	{
					// catch(Exception e)	{
						cerr << "stop codon: taxon " << taxset->GetTaxon(i) << " and codon " << j << " (site " << 3*j << ")\n";
						if (force_stops)	{
							// Data[i][j] = -2;
							Data[i][j] = -1;
						}
						else	{
							throw;
						}
					}
				}
			}

		}
		catch(Exception)	{
			cerr << "Codon Sequence Alignment: failed to read the datafile\n";
			exit(1);
		}
	}


void CodonSequenceAlignment::ToStreamRandomJackknife(ostream& os, double p)	{

	int* included = new int[Nsite];
	int ninc = 0;
	for (int i=0; i<Nsite; i++)	{
		if (Random::Uniform() < p)	{
			included[i] = 1;
			ninc ++;
		}
	}
	os << Ntaxa << '\t' << 3 * ninc << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}

	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			if (included[j])	{
				os << statespace->GetState(GetState(i,j));
			}
		}
		os << '\n';
	}
	os << '\n';
	delete[] included;
}

void CodonSequenceAlignment::ToStream(ostream& os)	{

	os << Ntaxa << '\t' << 3 * Nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}

	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			os << statespace->GetState(GetState(i,j));
		}
		os << '\n';
	}
	os << '\n';
}

void CodonSequenceAlignment::ToStream(ostream& os, int pos)	{

	int included = 0;
	for (int i=0; i<Nsite; i++)	{
		if (! AllMissingColumn(i))	{
			included++;
		}
	}
	os << Ntaxa << '\t' << included << '\n';
	// os << Ntaxa << '\t' << Nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}

	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			if (! AllMissingColumn(j))	{
				os << GetCodonStateSpace()->GetDNAStateSpace()->GetState(GetCodonStateSpace()->GetCodonPosition(pos,GetState(i,j)));
			}
		}
		os << '\n';
	}
	os << '\n';
}

void CodonSequenceAlignment::ToStreamFourFoldThird(ostream& os)	{

	int nsite = 0;
	int* included = new int[Nsite];
	for (int j=0; j<Nsite; j++)	{
		int i = 0;
		while ((i<Ntaxa) && (GetState(i,j) != -2) && ((GetState(i,j) == -1) || (GetCodonStateSpace()->GetDegeneracy(GetState(i,j)) == 4)))	{
			i++;
		}
		if ((i == Ntaxa) && (! AllMissingColumn(j)))	{
			included[j] = 1;
			nsite++;
		}
		else	{
			included[j] = 0;
		}
	}

	os << Ntaxa << '\t' << nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}

	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			if (included[j])	{
				os << GetCodonStateSpace()->GetDNAStateSpace()->GetState(GetCodonStateSpace()->GetCodonPosition(2,GetState(i,j)));
			}
		}
		os << '\n';
	}
	os << '\n';
}

void CodonSequenceAlignment::ToStreamFourFoldTriplet(ostream& os)	{

	int nsite = 0;
	int* included = new int[Nsite];
	for (int j=0; j<Nsite; j++)	{
		int i = 0;
		while ((i<Ntaxa) && (GetState(i,j) != -2) && ((GetState(i,j) == -1) || (GetCodonStateSpace()->GetDegeneracy(GetState(i,j)) == 4)))	{
			i++;
		}
		if (i == Ntaxa)	{
			included[j] = 1;
			nsite++;
		}
		else	{
			included[j] = 0;
		}
	}

	os << Ntaxa << '\t' << 3*nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}

	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			if (included[j])	{
				os << GetCodonStateSpace()->GetState(GetState(i,j));
			}
		}
		os << '\n';
	}
	os << '\n';
}

void CodonSequenceAlignment::ToStreamFourFoldThirdwoCpG(ostream& os)	{

	int nsite = 0;
	int* included = new int[Nsite];
	for (int j=0; j<Nsite-1; j++)	{
		bool keep = true;
		for (int i=0; i<Ntaxa; i++)	{
			if (GetState(i,j) != -1)	{
				if (GetState(i,j) == -2)	{
					keep = false;
				}
				else if (GetCodonStateSpace()->GetDegeneracy(GetState(i,j)) == 4)	{
						// if ( (GetState(i,j-1) != -2) && (GetState(i,j-1) != -1) && (GetCodonStateSpace()->GetCodonPosition(2,GetState(i,j-1)) != 1) 
						//  && (GetState(i,j+1) != -2) && (GetState(i,j+1) != -1) && (GetCodonStateSpace()->GetCodonPosition(0,GetState(i,j+1)) != 2) )  {
						// if (((GetState(i,j-1) == -1) || (GetCodonStateSpace()->GetCodonPosition(2,GetState(i,j-1)) != 1)) 
						//  && ((GetState(i,j+1) == -1) || (GetCodonStateSpace()->GetCodonPosition(0,GetState(i,j+1)) != 2) ) ) {
						if ( (GetCodonStateSpace()->GetCodonPosition(1,GetState(i,j)) == 1) || (GetCodonStateSpace()->GetCodonPosition(0,GetState(i,j+1)) == 2)) {
							keep = false;
						}
						/*
						else	{
							keep = false;
						}
						*/
				}
				else	{
					keep = false;
				}
			}
		}

		/*
		while (((i<Ntaxa) && (GetState(i,j) != -2) && ((GetState(i,j) == -1) || (GetCodonStateSpace()->GetDegeneracy(GetState(i,j)) == 4)))	&&
			// testing for absence of C in 3d position of previous codon and absence of G in 1st position of next codon
			( (GetState(i,j-1) != -2) && (GetState(i,j-1) != -1) && (GetCodonStateSpace()->GetCodonPosition(2,GetState(i,j-1)) != 1) 

			 && (GetState(i,j+1) != -2) && (GetState(i,j+1) != -1) && (GetCodonStateSpace()->GetCodonPosition(0,GetState(i,j+1)) != 2) ) ) {
			i++;
		}
		if (i == Ntaxa)	{
		*/
		if (keep && (! AllMissingColumn(j)))	{
			included[j] = 1;
			nsite++;
		}
		else	{
			included[j] = 0;
		}
	}

	os << Ntaxa << '\t' << nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}

	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite-1; j++)	{
			if (included[j])	{
				os << GetCodonStateSpace()->GetDNAStateSpace()->GetState(GetCodonStateSpace()->GetCodonPosition(2,GetState(i,j)));
			}
		}
		os << '\n';
	}
	os << '\n';
}
/*
void CodonSequenceAlignment::GetPairwiseDiff(int tax1, int tax2, double& tss, double& tsn, double& tvs, double& tvn)	{

	tss = 0;
	tsn = 0;
	tvs = 0;
	tvn = 0;

	int count = 0;
	for (int i=0; i<GetNsite(); i++)	{
		int c1 = Data(tax1,i);
		int c2 = Data(tax2,i);

		if ((c1 != -1) && (c2 != -1))	{
			count ++;




}

*/
