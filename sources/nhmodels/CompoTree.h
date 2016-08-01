

#ifndef COMPOTREE_H
#define COMPOTREE_H

#include "BiologicalSequences.h"
#include "DrawTree.h"

class CompoTree : public DrawTree	{

	public:

	CompoTree(Tree* from, int innstate, char* inalphabet = 0) : DrawTree(from)	{
		alphabet = inalphabet;
		nstate = innstate;
		compo = new double[nstate];
		composcale = 0.8;
		barwidth = 0.05;
		rowsep = -1.4;
		colsep = -2.1;
		alphabetfontsize = 4;
		if (! alphabet)	{
			alphabet = new char[nstate];
			if (nstate == Nnuc)	{
				for (int i=0; i<Nnuc; i++)	{
					alphabet[i] = DNAset[i];
				} 
			}
			else if (nstate == Naa)	{
				for (int i=0; i<Naa; i++)	{
					alphabet[i] = AminoAcids[i];
				}
				composcale = 0.2;
				barwidth = 0.03;
				rowsep = -1.4;
				colsep = -2.5;
				alphabetfontsize = 3;
			}
			else	{
				cerr << "error in CompoTree : should specify an alphabet\n";
				exit(1);
			}
		}
	}

	void GetComposition(Link* link)	{
		string from = link->GetNode()->GetName();
		if (link->isLeaf())	{
			ostringstream s;
			unsigned int k = 0;
			while ((k < from.length()) && (from[k] != '@'))	{
				k++;
			}
			k++;
			while (k < from.length())	{
				s << from[k];
				k++;
			}
			cerr << s.str() << '\n';
			GetComposition(s.str());
		}
		else	{
			GetComposition(from);
		}
		maxplus = 0;
		maxminus = 0;
		for (int i=0; i<nstate; i++)	{
			/*
			if (fabs(compo[i]) < 1e-4)	{
				compo[i] = 0;
			}
			*/
			if (compo[i] > 0)	{
				if (maxplus < compo[i])	{
					maxplus = compo[i];
				}
			}
			if (compo[i] < 0)	{
				if (maxminus < -compo[i])	{
					maxminus = -compo[i];
				}
			}
		}
	}

	void GetComposition(string from)	{
		for (int i=0; i<nstate; i++)	{
			compo[i] = 0;
		}
		unsigned int k = 0;
		if (from[0] == 'X')	{
			k++;
		}
		for (int i=0; i<nstate; i++)	{
			ostringstream s;
			while ((k < from.length()) && (from[k] != '&'))	{
				s << from[k];
				k++;
			}
			compo[i] = atof(s.str().c_str());
			k++;
		}
	}

	void	ToTikz(string target, double sizeX=6, double sizeY=10, double fontsize = 4);
	double	DrawTikz(Link* from, ostream& os, double X, double Y, double xscale, double yscale);
	void DrawProfile(ostream& os);

	private:

	char* alphabet;
	int nstate;
	double composcale;
	double* compo;
	double maxplus;
	double maxminus;
	double barwidth;
	double alphabetfontsize;
	double colsep;
	double rowsep;

};

#endif



