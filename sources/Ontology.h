
#ifndef ONTO_H
#define ONTO_H

#include "CovMatrix.h"

class Ontology	{

	public:

	Ontology(string filename)	{

		Lambda = 0;

		ifstream is(filename.c_str());
		is >> Ngene >> Nconcept;
		genename = new string[Ngene];
		conceptname = new string[Nconcept];
		for (int l=0; l<Nconcept; l++)	{
			is >> conceptname[l];
			// cerr << conceptname[l] << '\n';
			// conceptmap[conceptname[l]] = l;
		}
		z = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			z[gene] = new int[Nconcept];
			is >> genename[gene];
			genemap[genename[gene]] = gene;
			string tmp;
			is >> tmp;
			// cerr << genename[gene];
			// cerr << tmp << '\n';
			for (int l=0; l<Nconcept; l++)	{
				if (tmp[l] == '0')	{
					z[gene][l] = 0;
				}
				else if (tmp[l] == '1')	{
					z[gene][l] = 1;
				}
				else	{
					cerr << "error when parsing ontology\n";
					cerr << "gene : " << genename[gene] << '\n';
					cerr << tmp << '\n';
					exit(1);
				}
			}
		}
	}

	Ontology(Ontology* from, int inNgene, string* ingenename, int mingene, int minconcept, int* include)	{

		Lambda = 0;

		int* genecount = new int[from->GetNgene()];
		for (int gene=0; gene<from->GetNgene(); gene++)	{
			genecount[gene] = 0;
		}
		int* conceptcount = new int[from->GetNconcept()];
		for (int k=0; k<from->GetNconcept(); k++)	{
			conceptcount[k] = 0;
		}

		int* includegene = new int[from->GetNgene()];
		for (int gene=0; gene<from->GetNgene(); gene++)	{
			includegene[gene] = 0;
		}
		int m = 0;
		int missing = 0;
		for (int gene=0; gene<inNgene; gene++)	{
			int geneindex = from->GetGeneIndex(ingenename[gene]);
			if (geneindex == -1)	{
				// cerr << "in ontology reduction: did not find gene " << ingenename[gene] << '\n';
				// includegene[geneindex] = 0;
				missing ++;
			}
			else	{
				// cerr << ingenename[gene] << '\n';
				// cerr << geneindex << '\t' << from->GetGeneName(geneindex) << '\n';
				includegene[geneindex] = 1;
				m++;
			}
		}
		cerr << "number of missing genes     : " << missing << '\n';
		cerr << "initial number of genes kept: " << m << '\n';

		int* includeconcept = new int[from->GetNconcept()];
		for (int k=0; k<from->GetNconcept(); k++)	{
			includeconcept[k] = 1;
		}

		int cont = 1;
		while (cont)	{

			int tmpngene = 0;
			for (int gene=0; gene<from->GetNgene(); gene++)	{
				tmpngene += includegene[gene];
			}
			int tmpnconcept = 0;
			for (int k=0; k<from->GetNconcept(); k++)	{
				tmpnconcept += includeconcept[k];
			}

			for (int gene=0; gene<from->GetNgene(); gene++)	{
				genecount[gene] = 0;
			}
			for (int k=0; k<from->GetNconcept(); k++)	{
				conceptcount[k] = 0;
			}
			for (int gene=0; gene<from->GetNgene(); gene++)	{
				if (includegene[gene])	{
					for (int k=0; k<from->GetNconcept(); k++)	{
						if (includeconcept[k] && from->GetZ(gene,k))	{
							conceptcount[k]++;
							genecount[gene]++;
						}
					}
				}
			}
			cont = 0;
			for (int gene=0; gene<from->GetNgene(); gene++)	{
				if (includegene[gene] && (genecount[gene] < mingene))	{
					cont = 1;
					includegene[gene] = 0;
				}
			}
			for (int k=0; k<from->GetNconcept(); k++)	{
				if (includeconcept[k] && (conceptcount[k] < minconcept))	{
					cont = 1;
					includeconcept[k] = 0;
				}
			}
		}

		Ngene = 0;
		for (int gene=0; gene<from->GetNgene(); gene++)	{
			Ngene += includegene[gene];
		}
		Nconcept = 0;
		for (int k=0; k<from->GetNconcept(); k++)	{
			Nconcept += includeconcept[k];
		}

		genename = new string[Ngene];
		conceptname = new string[Nconcept];
		z = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			z[gene] = new int[Nconcept];
		}

		int concept = 0;
		int* conceptindexmap = new int[Nconcept];
		for (int k=0; k<from->GetNconcept(); k++)	{
			if (includeconcept[k])	{
				conceptname[concept] = from->GetConceptName(k);
				conceptindexmap[concept] = k;
				concept++;
			}
		}

		int gene = 0;
		for (int g=0; g<inNgene; g++)	{
			int geneindex = from->GetGeneIndex(ingenename[g]);
			if ((geneindex == -1) || (includegene[geneindex] == 0))	{
				include[g] = 0;
			}
			else	{
				include[g] = 1;
				genename[gene] = ingenename[g];
				genemap[genename[gene]] = gene;
				int concept = 0;
				for (int l=0; l<Nconcept; l++)	{
					z[gene][l] = from->GetZ(geneindex,conceptindexmap[l]);
				}
				gene++;
			}
		}

		delete[] conceptindexmap;
		delete[] conceptcount;
		delete[] genecount;
		delete[] includegene;
		delete[] includeconcept;
	}

	~Ontology()	{
		for (int gene=0; gene<Ngene; gene++)	{
			delete[] z[gene];
		}
		delete[] z;
		delete[] genename;
		delete[] conceptname;
	}

	int GetNgene() {return Ngene;}
	int GetNconcept() {return Nconcept;}

	int GetZ(int gene, int l) {return z[gene][l];}
	int* GetZ(int gene)	{return z[gene];}
	int** GetZ() {return z;}


	double** GetLambda()	{
		if (! Lambda)	{
			CreateAndComputeLambda();
		}
		return Lambda;
	}

	double** GetZmatrix()	{
		return Zmatrix;
	}

	double** GetTransZmatrix()	{
		return transZmatrix;
	}

	void CreateAndComputeLambda()	{
		Lambda = new double*[GetNconcept()];
		for (int k=0; k<GetNconcept(); k++)	{
			Lambda[k] = new double[GetNconcept()];
			for (int l=0; l<GetNconcept(); l++)	{
				Lambda[k][l] = 0;
				for (int i=0; i<GetNgene(); i++)	{
					Lambda[k][l] += GetZ(i,k) * GetZ(i,l);
				}
			}
		}
		Zmatrix = new double*[GetNconcept()];
		for (int k=0; k<GetNconcept(); k++)	{
			Zmatrix[k] = new double[GetNgene()];
			for (int i=0; i<GetNgene(); i++)	{
				Zmatrix[k][i] = GetZ(i,k);
			}
		}
		transZmatrix = new double*[GetNgene()];
		for (int i=0; i<GetNgene(); i++)	{
			transZmatrix[i] = new double[GetNconcept()];
			for (int k=0; k<GetNconcept(); k++)	{
				transZmatrix[i][k] = GetZ(i,k);
			}
		}
	}

	string GetGeneName(int gene) {return genename[gene];}
	string GetConceptName(int l) {return conceptname[l];}

	int GetGeneIndex(string name)	{
		if (genemap.find(name) == genemap.end())	{
			return -1;
		}
		return genemap[name];
	}

	void ToStream(ostream& os)	{
		os << GetNgene() << '\t' << GetNconcept() << '\n';
		for (int i=0; i<GetNconcept(); i++)	{
			os << GetConceptName(i) << '\t';
		}
		os << '\n';
		os << '\n';
		for (int i=0; i<GetNgene(); i++)	{
			os << GetGeneName(i) << '\t';
			for (int k=0; k<GetNconcept(); k++)	{
				os << GetZ(i,k);
			}
			os << '\n';
		}
	}

	protected:

	double** Lambda;
	double** Zmatrix;
	double** transZmatrix;

	int Ngene;
	int Nconcept;
	int** z;
	string* genename;
	string* conceptname;
	map<string,int> genemap;
	// map<string,int> conceptmap;
};

#endif
