
#include "Tree.h"
#include "SequenceAlignment.h"
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>

class Simulator : public NewickTree {

	public:

	Simulator(string datafile, string treefile, string paramfile)	{

		// get data from file

		ifstream is(datafile.c_str());
		is >> Ngene;
		cerr << Ngene << '\n';
		Nsite = new int[Ngene];
		nucdata = new FileSequenceAlignment*[Ngene];
		genename = new string[Ngene];
		// codondata = new CodonSequenceAlignment*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			string filename;
			is >> filename;
			cerr << filename << '\n';
			genename[gene] = filename;
			nucdata[gene] = new FileSequenceAlignment(filename);
			// codondata[gene] = new CodonSequenceAlignment(nucdata[gene]);
			int nsite = nucdata[gene]->GetNsite();
			Nsite[gene] = nsite;
			cerr << filename << '\t' << Nsite[gene] << '\n';
		}

		taxonset = nucdata[0]->GetTaxonSet();
		statespace = nucdata[0]->GetStateSpace();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);
		Ntaxa = taxonset->GetNtaxa();
		RecursiveSetBranchLengths(GetRoot());

		// get parameters from file
		ifstream prmis(paramfile.c_str());
		string tmp;
		prmis >> tmp;
		if (tmp == "Brownian")	{
			diffusion = 0;
			prmis >> sigma >> rootlogB;
		}
		else if (tmp == "OUP")	{
			diffusion = 1;
			prmis >> sigma >> phi >> meanlogB;
		}
		else	{
			cerr << "error: diffusion model not recognized\n";
			cerr << tmp << '\n';
			exit(1);
		}

		prmis >> alpha >> karyorate;

		prmis >> rho >> pi >> h;

		prmis >> pia >> pic >> pig >> pit;
		double total = pia + pic + pig + pit;
		if (fabs(total - 1)>1e-10)	{
			cerr << "error : nuc frequencies do not sum to 1 : " << total << '\n';
			exit(1);
		}
		q[0] = pia;
		q[1] = pia + pic;
		q[2] = pia + pic + pig;
		q[3] = 1;

		prmis >> mu;
		prmis >> at2cg >> at2gc >> at2ta >> cg2at >> cg2gc >> cg2ta;
		prmis >> cpgrate;
		mutrate[0][0] = 0;
		mutrate[0][1] = at2cg;
		mutrate[0][2] = at2gc;
		mutrate[0][3] = at2ta;
		mutrate[1][0] = cg2at;
		mutrate[1][1] = 0;
		mutrate[1][2] = cg2gc;
		mutrate[1][3] = cg2ta;
		mutrate[2][0] = cg2ta;
		mutrate[2][1] = cg2gc;
		mutrate[2][2] = 0;
		mutrate[2][3] = cg2at;
		mutrate[3][0] = at2ta;
		mutrate[3][1] = at2gc;
		mutrate[3][2] = at2cg;
		mutrate[3][3] = 0;

		strength[0] = 0;
		strength[1] = 1;
		strength[2] = 1;
		strength[3] = 0;

		prmis >> scale;
		prmis >> epsilon;
		prmis >> bmax;

		// aux array
		maxnsite = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			if (maxnsite < Nsite[gene])	{
				maxnsite = Nsite[gene];
			}
		}
		subrate = new double*[maxnsite];
		for (int i=0; i<maxnsite; i++)	{
			subrate[i] = new double[Nnuc+1];
		}
	}

	StateSpace* GetStateSpace()	{
		return statespace;
	}

	const Tree* GetTree()	const {
		return tree;
	}

	Tree* GetTree()	{
		return tree;
	}

	Link* GetRoot()	{
		return GetTree()->GetRoot();
	}

	const Link* GetRoot()	const {
		return GetTree()->GetRoot();
	}

	string GetNodeName(const Link* link)	const {
		ostringstream s;
		if (link->isLeaf())	{
			s << link->GetNode()->GetName() << "_";
		}
		s << exp(instantlogB.find(link->GetNode())->second);
		return s.str();
	}

	string GetBranchName(const Link* link)	const {
		if (link->isRoot())	{
			return "";
		}
		return GetTree()->GetBranchName(link);
	}

	void RecursiveSetBranchLengths(const Link* from)	{
		if (! from->isRoot())	{
			double l = atof(from->GetBranch()->GetName().c_str());
			bl[from->GetBranch()] = l;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetBranchLengths(link->Out());
		}
	}

	void Simulate()	{
		RecursiveSimulate(GetRoot());
	}

	void RecursiveSimulate(const Link* from)	{

		cerr << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\n';
		if (from->isRoot())	{
			RootSimulate();
		}
		else	{
			BranchSimulate(from);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSimulate(link->Out());
		}
	}

	void RootSimulate()	{

		// root value of B
		double logb = 0;
		if (diffusion == 0)	{
			logb = rootlogB;
		}
		else if (diffusion == 1)	{
			logb = meanlogB + Random::sNormal() * sqrt((sigma / (2 * phi)));

		}

		instantlogB[GetRoot()->GetNode()] = logb;
		double* gamma = new double[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			gamma[gene] = Random::Gamma(alpha,alpha);
		}
		genegamma[GetRoot()->GetNode()] = gamma;

		// coding sequences at root
		int** seq = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			seq[gene] = new int[Nsite[gene]];
			for (int i=0; i<Nsite[gene]; i++)	{
				seq[gene][i] = drawNucFromStat();
			}
		}
		geneseq[GetRoot()->GetNode()] = seq;

		// hot spot states
		int* hotspot = new int[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			double u = Random::Uniform();
			if (u < pi)	{
				hotspot[gene] = 1;
			}
			else	{
				hotspot[gene] = 0;
			}
		}
		genehotspotstate[GetRoot()->GetNode()] = hotspot;
	}

	void BranchSimulate(const Link* from)	{

		double l = bl[from->GetBranch()];

		// define current state as equal to state at ancestral node;
		// state includes
		// logb
		// genegammas
		// genehotspotstates
		// geneseqs
		double logb = instantlogB[from->Out()->GetNode()];
		double* ancgamma = genegamma[from->Out()->GetNode()];
		double* gamma = new double[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			gamma[gene] = ancgamma[gene];
		}

		int** ancseq = geneseq[from->Out()->GetNode()];
		int** seq = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			seq[gene] = new int[Nsite[gene]];
			for (int i=0; i<Nsite[gene]; i++)	{
				seq[gene][i] = ancseq[gene][i];
			}
		}

		int* anchotspot = genehotspotstate[from->Out()->GetNode()];
		int* hotspot = new int[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			hotspot[gene] = anchotspot[gene];
		}

		double meanb = 0;
		double* genemeanb = new double[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			genemeanb[gene] = 0;
		}

		double meanvarcoeff = 0;

		int count = 0;
		double t = 0;

		while (t < l)	{

			double b0 = exp(logb);
			double bcold = b0 / (pi*h + (1-pi));
			double bhot = h*b0 / (pi*h + (1-pi));

			double mean = 0;
			double var = 0;

			for (int gene=0; gene<Ngene; gene++)	{

				double b = 1;
				if (hotspot[gene])	{
					b = bhot * gamma[gene];
				}
				else	{
					b = bcold * gamma[gene];
				}

				mean += b;
				var += b*b;

				// draw a substitution event for current sequence
				int nsite = Nsite[gene];
				int* tmpseq = seq[gene];
				double totrate = 0;
				for (int i=0; i<nsite; i++)	{
					int prev = -1;
					int next = -1;
					int init = tmpseq[i];
					if (i > 0)	{
						prev = tmpseq[i-1];
					}
					if (i < nsite-1)	{
						next = tmpseq[i+1];
					}
					double totpos = 0;
					for (int final=0; final<Nnuc; final++)	{
						double tmp = epsilon * GetSubRate(init,final,prev,next,b);
						subrate[i][final] = tmp;
						totrate += tmp;
						totpos += tmp;
					}
					subrate[i][Nnuc] = totpos;
				}

				if (totrate > 0.1)	{
					cerr << "error, tot rate toolarge : " << totrate << '\n';
					exit(1);
				}

				double u = Random::Uniform();

				if (u < totrate)	{ // these is a substitution

					// choose position
					double tot = subrate[0][Nnuc];
					int pos = 0;
					while ((pos < nsite) && (u > tot))	{
						pos++;
						tot += subrate[pos][Nnuc];
					}
					if (pos == nsite)	{
						cerr << "error: overflow when choosing substituting position\n";
						exit(1);
					}

					// choose final nucleotide
					double v = subrate[pos][Nnuc] * Random::Uniform();
					tot = subrate[pos][0];
					int final = 0;
					while ((final < Nnuc) && (v > tot))	{
						final++;
						tot += subrate[pos][final];
					}
					if (final == Nnuc)	{
						cerr << "error: overflow when choosing final nucleotide\n";
						exit(1);
					}

					tmpseq[pos] = final;
				}

				genemeanb[gene] += b;

				// hotspot
				if (hotspot[gene])	{
					if (Random::Uniform() < epsilon * rho * (1-pi))	{
						hotspot[gene] = 0;
					}
				}
				else	{
					if (Random::Uniform() < epsilon * rho * pi)	{
						hotspot[gene] = 1;
					}
				}

				// karyotype
				if (Random::Uniform() < epsilon * karyorate)	{
					gamma[gene] = Random::Gamma(alpha,alpha);
				}

			}

			mean /= Ngene;
			var /= Ngene;
			var -= mean * mean;
			meanvarcoeff += sqrt(var / mean / mean);

			// global logb
			if (diffusion == 0)	{
				// brownian
				logb += sqrt(sigma * epsilon) * Random::sNormal();
			}
			else if (diffusion == 1)	{
				// oup
				logb += -phi*(logb - meanlogB) + sqrt(sigma * epsilon) * Random::sNormal();
			}
			else	{
				cerr << "unrecognized diffusion\n";
				exit(2);
			}

			// maintain averages

			meanb += b0;
			t += epsilon;
			count++;

		}

		// store new values
		instantlogB[from->GetNode()] = logb;
		geneseq[from->GetNode()] = seq;
		genegamma[from->GetNode()] = gamma;
		genehotspotstate[from->GetNode()] = hotspot;

		// averages
		meanb /= count;
		meanB[from->GetBranch()] = meanb;
		genemeanB[from->GetBranch()] = genemeanb;

		meanvarcoeff /= count;
		instantvarcoeff[from->GetBranch()] = meanvarcoeff;

		double mean = 0;
		double var = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			genemeanb[gene] /= count;
			mean += genemeanb[gene];
			var += genemeanb[gene] * genemeanb[gene];
		}
		mean /= Ngene;
		var /= Ngene;
		var -= mean * mean;
		varcoeff[from->GetBranch()] = sqrt(var / mean / mean);

	}

	double GetMutRate(int init, int final, int prev, int next)	{
		if (init == final)	{
			return 0;
		}
		if ((prev == 1) && (init == 2) && (final == 0))	{
			return mu * cpgrate;
		}
		if ((next == 3) && (init == 1) && (final == 3))	{
			return mu * cpgrate;
		}
		return mu * mutrate[init][final];
	}

	double GetFixProb(int init, int final, double B)	{
		double b = B * (strength[final] - strength[init]);
		if (b > bmax)	{
			b = bmax;
		}
		if (b < -bmax)	{
			b = -bmax;
		}
		double fixprob = 1;
		if (fabs(b) > 1e-6)	{
			fixprob = b / (1 - exp(-b));
		}
		return fixprob;
	}

	double GetSubRate(int init, int final, int prev, int next, double B)	{
		return GetMutRate(init,final,prev,next) * GetFixProb(init,final,B);
	}

	int drawNucFromStat()	{
		double u = Random::Uniform();
		int k = 0;
		while ((k<Nnuc) && (u>q[k]))	{
			k++;
		}
		if (k == Nnuc)	{
			cerr << "error in draw nuc: overflow\n";
			exit(1);
		}
		return k;
	}

	void WriteSimu(string basename)	{

		// dataset
		ofstream os((basename + ".list").c_str());
		os << Ngene << '\n';
		for (int gene=0; gene<Ngene; gene++)	{
			os << basename << genename[gene] << '\n';
		}

		for (int gene=0; gene<Ngene; gene++)	{
			ofstream os((basename + genename[gene]).c_str());
			os << Ntaxa << '\t' << Nsite[gene] << '\n';
			RecursiveWriteData(GetRoot(),gene,os);
		}

		// true value of b
		ofstream bos((basename + ".trueb").c_str());
		ToStream(bos);

		// global means
		ofstream sos((basename + ".simu").c_str());
		RecursiveWriteGlobalMeans(GetRoot(),sos);

		// gene means
		ofstream gos((basename + ".geneb").c_str());
		RecursiveWriteGeneMeans(GetRoot(),gos);

	}

	void RecursiveWriteData(const Link* from, int gene, ostream& os)	{
		if (from->isLeaf())	{
			os << from->GetNode()->GetName() << '\t';
			int* seq = geneseq[from->GetNode()][gene];
			for (int i=0; i<Nsite[gene]; i++)	{
				os << GetStateSpace()->GetState(seq[i]);
			}
			os << '\n';
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveWriteData(link->Out(), gene, os);
		}
	}

	void RecursiveWriteGlobalMeans(const Link* from, ostream& os)	{
		if (! from->isRoot())	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
			os << meanB[from->GetBranch()] << '\t' << varcoeff[from->GetBranch()] << '\t' << instantvarcoeff[from->GetBranch()] << '\n';
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveWriteGlobalMeans(link->Out(), os);
		}
	}

	void RecursiveWriteGeneMeans(const Link* from, ostream& os)	{
		if (! from->isRoot())	{
			double* mean = genemeanB[from->GetBranch()];
			for (int gene=0; gene<Ngene; gene++)	{
				os << mean[gene] << '\t';
			}
			os << '\n';
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveWriteGeneMeans(link->Out(), os);
		}
	}

	// brownian for logB
	// or OUP
	int diffusion;

	double sigma;
	double meanlogB;
	double rootlogB;
	double phi;

	// gene long-term karyotypic fluctuations
	double alpha;
	double karyorate;

	// hot spot process
	double rho;
	double pi;
	double h;

	// number of nucleotide positions
	int* Nsite;
	int Ntaxa;

	double pia, pic, pig, pit;
	double q[Nnuc];
	// absolute mutation rate;
	double mu;
	// all relative mutation rates
	double at2cg, at2gc, at2ta, cg2at, cg2gc, cg2ta;
	double mutrate[Nnuc][Nnuc];
	int strength[Nnuc];

	// including CpG effect
	double cpgrate;

	int Ngene;

	Tree* tree;
	// with branch lengths, measured in time
	map<const Branch*, double> bl;
	map<const Node*, double> instantlogB;

	map<const Node*, double*> genegamma;
	map<const Node*, int*> genehotspotstate;
	map<const Node*, int**> geneseq;

	// branch averages
	map<const Branch*, double> meanB;
	map<const Branch*, double*> genemeanB;
	map<const Branch*, double> varcoeff;
	map<const Branch*, double> instantvarcoeff;

	// a global time scale for easy translation
	double scale; // 100 Myr
	double epsilon; // 0.1 Myr ?
	double bmax;

	int maxnsite;
	double** subrate;

	string* genename;
	FileSequenceAlignment** nucdata;
	SequenceAlignment** data;
	TaxonSet* taxonset;
	StateSpace* statespace;

};

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string paramfile = argv[3];
	string basename = argv[4];

	cerr << "new simulation\n";

	Simulator* sim = new Simulator(datafile,treefile,paramfile);

	cerr << "simulate\n";
	sim->Simulate();

	cerr << "output to file\n";
	sim->WriteSimu(basename);

	cerr << "ok\n";


}

