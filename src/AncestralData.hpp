#ifndef ANCESTRALDATA_H
#define ANCESTRALDATA_H

#include "Tree.h"
#include "TaxonSet.h"
#include <fstream>
#include <cmath>
#include <map>
#include "SumConstrained.h"

class AncestralData	{

	public:

	AncestralData(Tree* intree) : tree(intree) {}

	// the list of taxa
	TaxonSet* GetTaxonSet()	{
		return taxset;
	}

	int GetNsite()	{
		return Nsite;
	}

	Tree* GetTree()	{
		return tree;
	}

	int GetNentry()	{
		return data.size();
	}

	bool Exists(const Link* link)	{
		return (data.find(link) != data.end());
	}

	void SetState(const Link* link, int site, double state)	{
		if (data.find(link) == data.end())	{
			cerr << "error in ancestral data: did not find link in hash\n";
			tree->ToStream(cerr,link);
			cerr << '\n';
			exit(1);
		}
		data[link][site] = state;
	}

	double GetState(const Link* link, int site)	{
		if (data.find(link) == data.end())	{
			cerr << "error in ancestral data: did not find link in hash\n";
			tree->ToStream(cerr,link);
			cerr << '\n';
			exit(1);
		}
		return data[link][site];
	}

	void SetState(string taxon1, string taxon2, int site, double state)	{
		const Link* link = tree->GetLCA(taxon1,taxon2);
		if (data.find(link) == data.end())	{
			cerr << "error in ancestral data: did not find link in hash\n";
			cerr << taxon1 << '\t' << taxon2 << '\n';
			exit(1);
		}
		data[link][site] = state;
	}

	double GetState(string taxon1, string taxon2, int site)	{
		const Link* link = tree->GetLCA(taxon1,taxon2);
		if (data.find(link) == data.end())	{
			cerr << "error in ancestral data: did not find link in hash\n";
			cerr << taxon1 << '\t' << taxon2 << '\n';
			exit(1);
		}
		return data[link][site];
	}

	bool isMissing(string taxon1, string taxon2, int site)	{
		return (GetState(taxon1,taxon2,site) == -1);
	}

	void ToStream(ostream& os)	{
		int count = 0;
		for (map<const Link*, double*>::iterator i = data.begin(); i!= data.end(); i++)	{
			count++;
		}
		os << count << '\t' << Nsite << '\n';
		RecursiveToStream(tree->GetRoot(),os);
	}

	void RecursiveToStream(const Link* from, ostream& os)	{

		os << tax1[from] << '\t' << tax2[from];
		double* tmp = data[from];
		for (int i=0; i<Nsite; i++)	{
			os << '\t' << tmp[i];
		}
		os << '\n';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveToStream(link->Out(),os);
		}
	}

	int Nsite;
	TaxonSet* taxset;
	Tree* tree;
	map<const Link*, double*> data;
	map<const Link*, string> tax1;
	map<const Link*, string> tax2;

};

class FileAncestralData : public AncestralData {

	public:
	FileAncestralData(Tree* intree, istream& is)	: AncestralData(intree) {
		ReadDataFromFile(is);
	}

	FileAncestralData(Tree* intree, string filename)	: AncestralData(intree)	{
		ifstream is(filename.c_str());
		ReadDataFromFile(is);
	}

	FileAncestralData(Tree* intree, string filename, int inNsite)	: AncestralData(intree)	{
		Nsite = inNsite;
		ifstream is(filename.c_str());
		ReadTemplateFromFile(is);
	}

	private:

	int  ReadDataFromFile(istream& is)	{

		map<string,int> namemap;
		int nentry;
		is >> nentry >> Nsite;
		string taxon1, taxon2;
		for (int i=0; i<nentry; i++)	{
			is >> taxon1 >> taxon2;
			namemap[taxon1] = 1;
			namemap[taxon2] = 1;
			const Link* link = tree->GetLCA(taxon1, taxon2);
			if (data.find(link) != data.end())	{
				cerr << "error in ancestral data: found node twice in hash table\n";
				cerr << taxon1 << '\t' << taxon2 << '\n';
				cerr << tax1[link] << '\t' << tax2[link] << '\n';
				cerr << link << '\n';

				/*
				cerr << '\n';
				tree->ToStreamSimplified(cerr,link);
				cerr << '\n';
				*/

				exit(1);
			}
			if (link)	{
				double* tmp = new double[Nsite];
				for (int j=0; j<Nsite; j++)	{
					is >> tmp[j];
				}
				data[link] = tmp;
				tax1[link] = taxon1;
				tax2[link] = taxon2;
			}
			else	{
				cerr << "did not find ancestor of: " << taxon1 << '\t' << taxon2 << '\n';
				double tmp;
				for (int j=0; j<Nsite; j++)	{
					is >> tmp;
				}
			}
		}
		int ntaxa = namemap.size();
		string* name = new string[ntaxa];
		int count = 0;
		for (map<string,int>::iterator i = namemap.begin(); i!= namemap.end(); i++)	{
			name[count++] = i->first;
		}
		if (count != ntaxa)	{
			cerr << "error in file ancestral data: non matching tax size\n";
			exit(1);
		}
		taxset = new TaxonSet(name,ntaxa);
		delete[] name;
		return 1;
	}

	int  ReadTemplateFromFile(istream& is)	{

		map<string,int> namemap;
		int nentry;
		is >> nentry;
		string taxon1, taxon2;
		for (int i=0; i<nentry; i++)	{
			is >> taxon1 >> taxon2;
			namemap[taxon1] = 1;
			namemap[taxon2] = 1;
			const Link* link = tree->GetLCA(taxon1, taxon2);
			if (data.find(link) != data.end())	{
				cerr << "error in ancestral data: found node twice in hash table\n";
				cerr << taxon1 << '\t' << taxon2 << '\n';
				cerr << tax1[link] << '\t' << tax2[link] << '\n';
				cerr << link << '\n';

				/*
				cerr << '\n';
				tree->ToStreamSimplified(cerr,link);
				cerr << '\n';
				*/

				exit(1);
			}
			if (link)	{
				double* tmp = new double[Nsite];
				for (int j=0; j<Nsite; j++)	{
					tmp[j] = 0;
				}
				data[link] = tmp;
				tax1[link] = taxon1;
				tax2[link] = taxon2;
			}
			else	{
				cerr << "did not find ancestor of : " << taxon1 << '\t' << taxon2 << '\n';
				exit(1);
			}
		}
		int ntaxa = namemap.size();
		string* name = new string[ntaxa];
		int count = 0;
		for (map<string,int>::iterator i = namemap.begin(); i!= namemap.end(); i++)	{
			name[count++] = i->first;
		}
		if (count != ntaxa)	{
			cerr << "error in file ancestral data: non matching tax size\n";
			exit(1);
		}
		taxset = new TaxonSet(name,ntaxa);
		delete[] name;
		return 1;
	}
};

class SumConstrainedAncestralData : public AncestralData	{

	SumConstrainedMapping* basis;

	public:

	SumConstrainedMapping* GetBasis() {return basis;}

	SumConstrainedAncestralData(AncestralData* from) : AncestralData(from->GetTree())	{

		Nsite = from->GetNsite() - 1;
		int fromnsite = from->GetNsite();
		taxset = from->GetTaxonSet();

		basis = new SumConstrainedMapping(fromnsite);

		double** b = basis->base;
		double* x = new double[fromnsite];
		double* y = new double[fromnsite];
		map<const Link*, double*>& fromdata = from->data;
		for (map<const Link*,double*>::iterator linki=fromdata.begin(); linki!=fromdata.end(); linki++)	{

			double total = 0;
			double* fromd = fromdata[linki->first];
			for (int i=0; i<fromnsite; i++)	{
				x[i] = log(fromd[i]);
				total += x[i];
			}
			total /= fromnsite;
			for (int i=0; i<fromnsite; i++)	{
				x[i] -= total;
			}
			for (int i=0; i<fromnsite; i++)	{
				double tmp = 0;
				for (int j=0; j<fromnsite; j++)	{
					tmp += b[i][j] * x[j];
				}
				y[i] = tmp;
				if (!i)	{
					if (fabs(y[i]) > 1e-6)	{
						cerr << "error : " << y[i] << '\n';
					}
				}
			}

			double* d = new double[Nsite];

			for (int i=1; i<fromnsite; i++)	{
				d[i-1] = y[i];
			}

			tax1[linki->first] = from->tax1[linki->first];
			tax2[linki->first] = from->tax2[linki->first];
			data[linki->first] = d;
		}
		delete[] x;
		delete[] y;
	}
};


#endif

