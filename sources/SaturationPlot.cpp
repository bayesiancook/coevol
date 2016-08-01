#include "SequenceAlignment.h"
#include "Tree.h"

class NucSequenceAlignment : public SequenceAlignment	{

	public:

	NucSequenceAlignment(SequenceAlignment* from)	{

		Ntaxa = from->Ntaxa;
		Nsite = from->Nsite;
		taxset = from->taxset;
		statespace = from->statespace;

		nucstatespace = dynamic_cast<DNAStateSpace*>(statespace);
		if (! nucstatespace)	{
			cerr << "error : not nucleotide\n";
			exit(1);
		}

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = from->Data[i][j];
			}
		}
	}

	double GetTransitionDistance(int tax1, int tax2)	{
		double total = 0;
		int count = 0;

		for (int i=0; i<Nsite; i++)	{
			int a1 = Data[tax1][i];
			int a2 = Data[tax2][i];
			if ((a1 != -1) && (a2 != -1))	{
				count++;
				if (((a1 == 0) && (a2 == 2)) || ((a1 == 2) && (a2 == 0)) || ((a1 == 1) && (a2 == 3)) || ((a1 == 3) && (a2 == 1)))	{
					total ++;
				}
			}
		}

		return total / count;
	}

	double GetTransversionDistance(int tax1, int tax2)	{
		double total = 0;
		int count = 0;

		for (int i=0; i<Nsite; i++)	{
			int a1 = Data[tax1][i];
			int a2 = Data[tax2][i];
			if ((a1 != -1) && (a2 != -1))	{
				count++;
				if (a1 != a2)	{
					if (!(((a1 == 0) && (a2 == 2)) || ((a1 == 2) && (a2 == 0)) || ((a1 == 1) && (a2 == 3)) || ((a1 == 3) && (a2 == 1))))	{
						total ++;
					}
				}
			}
		}

		return total / count;
	}

	protected:
	DNAStateSpace* nucstatespace;
};


class LTree : public Tree	{

	public:

	LTree(string treefile) : Tree(treefile) {}

	double GetLength(const Branch* branch)	{
		return atof(branch->GetName().c_str());
	}

	double GetLengthSinceLCA(const Link* down, const Link* up)	{
		const Link* link = down;
		double total = 0;
		while (link != up)	{
			total += GetLength(link->GetBranch());

			const Link* anc = link->Out();
			while (anc->Next() != link->Out())	{
				anc = anc->Next();
			}

			link = anc;
		}
		return total;
	}

	double GetLengthBetweenTaxonPair(string tax1, string tax2)	{
		/*
		const Link* leaf1 = GetLCA(tax1,tax1);
		const Link* leaf2 = GetLCA(tax2,tax2);
		const Link* anc = GetLCA(tax1,tax2);
		return GetLengthSinceLCA(leaf1,anc) + GetLengthSinceLCA(leaf2,anc);
		*/

		double totlength = 0;
		bool found1 = false;
		bool found2 = false;
		RecursiveGetLCAWithLength(GetRoot(),tax1,tax2,found1,found2,totlength);
		return totlength;
	}
	// returns 0 if not found
	// returns link if found (then found1 and found2 must

	const Link* RecursiveGetLCAWithLength(const Link* from, string tax1, string tax2, bool& found1, bool& found2, double& totlength)	{
		const Link* ret= 0;
		if (from->isLeaf())	{
			// found1 |= (from->GetNode()->GetName() == tax1);
			// found2 |= (from->GetNode()->GetName() == tax2);
			found1 |= (GetLeafNodeName(from) == tax1);
			found2 |= (GetLeafNodeName(from) == tax2);
			if (! ret)	{
				if (found1 && found2)	{
					ret = from;
				}
			}
		}
		else	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				bool tmp1 = false;
				bool tmp2 = false;
				const Link* ret2 = RecursiveGetLCAWithLength(link->Out(),tax1,tax2,tmp1,tmp2,totlength);
				found1 |= tmp1;
				found2 |= tmp2;
				if (ret2)	{
					if (ret)	{
						cerr << "error : found node twice\n";
						cerr << tax1 << '\t' << tax2 << '\n';
						ToStream(cerr,ret2->Out());
						cerr << '\n';
						ToStream(cerr,ret->Out());
						cerr << '\n';
						exit(1);
					}
					ret = ret2;
				}
			}
			if (! ret)	{
				if (found1 && found2)	{
					ret = from;
				}
			}
		}
		if ((found1 || found2) && ! (found1 && found2))	{
			totlength += GetLength(from->GetBranch());
		}
		return ret;
	}

	void ReadFromStream(istream& is);
};


int main(int argc, char* argv[])	{


	string datafile = argv[1];
	string treefile = argv[2];

	FileSequenceAlignment* data = new FileSequenceAlignment(datafile);
	NucSequenceAlignment* nucdata = new NucSequenceAlignment(data);
	TaxonSet* taxset = data->GetTaxonSet();
	LTree* tree = new LTree(treefile);

	tree->RegisterWith(taxset);


	int Ntaxa = taxset->GetNtaxa();
	int Nsite = data->GetNsite();

	for (int i=0; i<Ntaxa; i++)	{
		for (int j=i+1; j<Ntaxa; j++)	{
			string tax1 = taxset->GetTaxon(i);
			string tax2 = taxset->GetTaxon(j);
			cerr << taxset->GetTaxon(i) << '\t' << taxset->GetTaxon(j) << '\t' << tree->GetLengthBetweenTaxonPair(tax1,tax2) << '\t' << nucdata->GetTransitionDistance(i,j) << '\t' << nucdata->GetTransversionDistance(i,j) << '\n';
		}
	}
}

