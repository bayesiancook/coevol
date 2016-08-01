
#include "Tree.h"
#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>

class Simulator : public NewickTree {

	public:

	Simulator(string datafile, string treefile, string paramfile, double inrho, double inpi, double inh, double incpgrate)	{

		// get data from file
		withextb = false;

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
			if (nsite % 3)	{
				cerr << "error : " << nsite << " not multiple of 3\n";
				exit(1);
			}
			Nsite[gene] = nsite;
			cerr << filename << '\t' << Nsite[gene] << '\n';
		}

		taxonset = nucdata[0]->GetTaxonSet();
		statespace = nucdata[0]->GetStateSpace();
		codonstatespace = new CodonStateSpace(Universal);

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
		if (inrho != -1)	{
			rho = inrho;
		}
		if (inpi != -1)	{
			pi = inpi;
		}
		if (inh != -1)	{
			h = inh;
		}

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
		if (incpgrate != -1)	{
			cpgrate = incpgrate;
		}
		prmis >> omegamean;
		prmis >> omegagenealpha;
		prmis >> omegasitealpha;


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

	void LoadExtB(string filename)	{
		ifstream is(filename.c_str());
		int N;
		is >> N;
		for (int i=0; i<N; i++)	{
			string tmp1, tmp2;
			double temp1, temp2;
			is >> tmp1 >> tmp2 >> temp1 >> temp2;
			if ((taxonset->GetTaxonIndex(tmp1) != -1) && (taxonset->GetTaxonIndex(tmp2) != -1))	{
			const Link* link = tree->GetLCA(tmp1,tmp2);
			if (link->isRoot())	{
				cerr << tmp1 << '\t' << tmp2 << '\t' << tree->GetLeftMost(link) << '\t' << tree->GetRightMost(link) << '\n';
				cerr << temp2 << '\n';
			}
			extB[link->GetNode()] = temp2;
			}
		}
		withextb = true;
	}

	/*
	double GetNodeVal(const Node* node)	{
		string name = node->GetName();
		ostringstream s;
		unsigned int i = 0;
		while ((i < name.length()) && (name[i] != '_'))	{
			s << name[i];
			i++;
		}
		double ret = atof(s.str().c_str());
		// cerr << ret << '\n';
		return ret;
	}

	void RecursiveGetNodeVal(const Link* from)	{
		cerr << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\t' << GetNodeVal(from->GetNode()) << '\n';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetNodeVal(link->Out());
		}
	}
	*/

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
		omega = new double*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			omega[gene] = new double[Nsite[gene] / 3];
			double omegagene = omegamean * Random::Gamma(omegagenealpha,omegagenealpha);
			for (int i=0; i<Nsite[gene] / 3; i++)	{
				omega[gene][i] = omegagene * Random::Gamma(omegasitealpha,omegasitealpha);
			}
		}
		Link* link = new Link();
		Node* node = new Node();
		Branch* branch = new Branch();
		link->SetBranch(branch);
		GetRoot()->SetBranch(branch);
		link->SetNode(node);
		bl[branch] = 100;
		link->SetOut(GetRoot());
		GetRoot()->SetOut(link);
		extB[node] = extB[GetRoot()->GetNode()];
		RootSimulate(link);
		bksigma = sigma;
		if (diffusion == 0)	{
			sigma = 0;
		}
		RecursiveSimulate(GetRoot());
		GetRoot()->SetOut(GetRoot());
		GetRoot()->SetBranch(0);
		delete link;
		delete node;
		delete branch;
	}

	void RecursiveSimulate(const Link* from)	{

		cerr << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\n';
		/*
		if (from->isRoot())	{
			RootSimulate();
		}
		else	{
			BranchSimulate(from);
		}
		*/
		BranchSimulate(from);
		sigma = bksigma;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSimulate(link->Out());
		}
	}

	void RootSimulate(const Link* link)	{

		// root value of B
		double logb = 0;
		if (withextb)	{
			logb = log(extB[link->GetNode()]);
		}
		else	{
			if (diffusion == 0)	{
				logb = rootlogB;
			}
			else if (diffusion == 1)	{
				logb = meanlogB + Random::sNormal() * sqrt((sigma / (2 * phi)));
			}
		}

		instantlogB[link->GetNode()] = logb;
		// instantlogB[GetRoot()->GetNode()] = logb;
		double* gamma = new double[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			gamma[gene] = Random::Gamma(alpha,alpha);
		}
		genegamma[link->GetNode()] = gamma;
		// genegamma[GetRoot()->GetNode()] = gamma;

		// coding sequences at root
		int** seq = new int*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			seq[gene] = new int[Nsite[gene]];
			for (int i=0; i<Nsite[gene]/3; i++)	{
				int pos1,pos2,pos3;
				drawCodonFromStat(pos1,pos2,pos3);
				seq[gene][3*i+0] = pos1;
				seq[gene][3*i+1] = pos2;
				seq[gene][3*i+2] = pos3;
			}
			/*
			for (int i=0; i<Nsite[gene]; i++)	{
				seq[gene][i] = drawNucFromStat();
			}
			*/
		}
		geneseq[link->GetNode()] = seq;
		// geneseq[GetRoot()->GetNode()] = seq;

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
		genehotspotstate[link->GetNode()] = hotspot;
		// genehotspotstate[GetRoot()->GetNode()] = hotspot;
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
		double logb1 = 0;
		double logb2 = 0;
		if (withextb)	{
			if (!from->isRoot())	{
				logb1 = log(extB[from->Out()->GetNode()]);
				logb2 = log(extB[from->GetNode()]);
			}
		}
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
			if (withextb)	{
				b0 = exp(((l-t)*logb1 + t*logb2) / l);
			}
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

				double w2s = 1;
				double s2w = 1;
				if (b > bmax)	{
					w2s = b;
					s2w = 0;
				}
				else	{
					w2s = b / (1 - exp(-b));
					s2w = -b / (1 - exp(b));
				}

				// draw a substitution event for current sequence
				int nsite = Nsite[gene];
				int* tmpseq = seq[gene];
				double totrate = 0;
				for (int i=0; i<nsite; i++)	{
					int prev = -1;
					int next = -1;
					int prev2 = -1;
					int next2 = -1;
					int init = tmpseq[i];
					int jprev = i-1;
					if (jprev >=0)	{
						prev = tmpseq[jprev];
					}
					jprev--;
					if (jprev >=0)	{
						prev2 = tmpseq[jprev];
					}
					int jnext = i+1;
					if (jnext < nsite)	{
						next = tmpseq[jnext];
					}
					jnext++;
					if (jnext < nsite)	{
						next2 = tmpseq[jnext];
					}
					int phase = i % 3;
					double o = omega[gene][i/3];
					double totpos = 0;
					for (int final=0; final<Nnuc; final++)	{
						double tmp = epsilon * GetSubRate(init,final,prev,next,prev2,next2,phase,w2s,s2w,o);
						subrate[i][final] = tmp;
						totrate += tmp;
						totpos += tmp;
					}
					subrate[i][Nnuc] = totpos;
				}

				if (totrate > 0.5)	{
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
		if (withextb)	{
			instantlogB[from->GetNode()] = logb2;
		}
		else	{
			instantlogB[from->GetNode()] = logb;
		}
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

	double GetSubRate(int init, int final, int prev, int next, int prev2, int next2, int phase, double w2s, double s2w, double omega)	{
		double om = 1;
		if (phase == 0)	{
			if (codonstatespace->CheckStop(final,next,next2))	{
				om = 0;
			}
			int initcodon = codonstatespace->GetCodonFromDNA(init,next,next2);
			int finalcodon = codonstatespace->GetCodonFromDNA(final,next,next2);
			if (! codonstatespace->Synonymous(initcodon,finalcodon))	{
				om *= omega;
			}
		}
		else if (phase == 1)	{
			if (codonstatespace->CheckStop(prev,final,next))	{
				om = 0;
			}
			int initcodon = codonstatespace->GetCodonFromDNA(prev,init,next);
			int finalcodon = codonstatespace->GetCodonFromDNA(prev,final,next);
			if (! codonstatespace->Synonymous(initcodon,finalcodon))	{
				om *= omega;
			}
		}
		else	{
			if (codonstatespace->CheckStop(prev2,prev,final))	{
				om = 0;
			}
			int initcodon = codonstatespace->GetCodonFromDNA(prev2,prev,init);
			int finalcodon = codonstatespace->GetCodonFromDNA(prev2,prev,final);
			if (! codonstatespace->Synonymous(initcodon,finalcodon))	{
				om *= omega;
			}
		}
		double fixprob = 1;
		if (strength[final] - strength[init] == 1)	{
			fixprob = w2s;
		}
		if (strength[final] - strength[init] == -1)	{
			fixprob = s2w;
		}
		return om * GetMutRate(init,final,prev,next) * fixprob;
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

	void drawCodonFromStat(int& pos1, int& pos2, int& pos3)	{
		bool stop = true;
		while (stop)	{
			pos1 = drawNucFromStat();
			pos2 = drawNucFromStat();
			pos3 = drawNucFromStat();
			stop = codonstatespace->CheckStop(pos1,pos2,pos3);
		}
	}

	void WriteSimu(string basename)	{

		// dataset
		ofstream allos((basename + "_all.list").c_str());
		allos << Ngene << '\n';
		ofstream thirdos((basename + "_pos3.list").c_str());
		thirdos << Ngene << '\n';
		ofstream fourfoldos((basename + "_4fold.list").c_str());
		fourfoldos << Ngene << '\n';
		ofstream wocpgos((basename + "_4foldwocpg.list").c_str());
		wocpgos << Ngene << '\n';
		for (int gene=0; gene<Ngene; gene++)	{
			allos << basename << genename[gene] << "_all.ali\n";
			thirdos << basename << genename[gene] << "_pos3.ali\n";
			fourfoldos << basename << genename[gene] << "_4fold.ali\n";
			wocpgos << basename << genename[gene] << "_4foldwocpg.ali\n";
		}

		for (int gene=0; gene<Ngene; gene++)	{
			int** data = new int*[Ntaxa];
			string* names = new string[Ntaxa];
			int n = 0;
			RecursiveMakeData(GetRoot(),gene,data,names,n);
			if (n != Ntaxa)	{
				cerr << "error when making sequence alignment: wrong number of taxa : " << n << '\t' << " instead of " << Ntaxa << '\n';
				exit(1);
			}

			cerr << "make new ali\n";
			SequenceAlignment* simuali = new SequenceAlignment(data,names,Nsite[gene],statespace,taxonset);

			cerr << "make codon ali\n";
			CodonSequenceAlignment* codonsimuali = new CodonSequenceAlignment(simuali,true,Universal);

			cerr << "output\n";
			ofstream allos((basename + genename[gene] + "_all.ali").c_str());
			codonsimuali->ToStream(allos);

			ofstream thirdos((basename + genename[gene] + "_pos3.ali").c_str());
			codonsimuali->ToStream(thirdos,2);

			ofstream fourfoldos((basename + genename[gene] + "_4fold.ali").c_str());
			codonsimuali->ToStreamFourFoldThird(fourfoldos);

			ofstream wocpgos((basename + genename[gene] + "_4foldwocpg.ali").c_str());
			codonsimuali->ToStreamFourFoldThirdwoCpG(wocpgos);

			delete[] data;
			delete[] names;
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

	void RecursiveMakeData(const Link* from, int gene, int** data, string* names, int& n)	{
		if (from->isLeaf())	{
			names[n] = from->GetNode()->GetName();
			data[n] = geneseq[from->GetNode()][gene];
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveMakeData(link->Out(), gene, data, names, n);
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
	double bksigma;
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
	double** omega;
	double omegamean;
	double omegasitealpha;
	double omegagenealpha;

	map<const Node*, double> extB;
	bool withextb;

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
	CodonStateSpace* codonstatespace;

};

int main(int argc, char* argv[])	{

	/*
	CodonStateSpace* cod = new CodonStateSpace(Universal);
	for (int pos1=0; pos1<Nnuc; pos1++)	{
	for (int pos2=0; pos2<Nnuc; pos2++)	{
	for (int pos3=0; pos3<Nnuc; pos3++)	{
			cerr << cod->GetState(cod->GetCodonFromDNA(pos1,pos2,pos3)) << '\t' << cod->CheckStop(pos1,pos2,pos3) << '\n';
	}
	}
	}
	exit(1);
	*/
	string datafile = argv[1];
	string treefile = argv[2];
	string paramfile = argv[3];
	double rho = -1;
	double pi = -1;
	double h = -1;
	double cpgrate = -1;
	string basename = "";
	string extbfile = "";
	if (argc == 10)	{
		rho = atof(argv[4]);
		pi = atof(argv[5]);
		h = atof(argv[6]);
		cpgrate = atof(argv[7]);
		extbfile = argv[8];
		basename = argv[9];
	}
	else if (argc == 9)	{
		rho = atof(argv[4]);
		pi = atof(argv[5]);
		h = atof(argv[6]);
		cpgrate = atof(argv[7]);
		basename = argv[8];
	}
	else if (argc == 6)	{
		extbfile = argv[4];
		basename = argv[5];
	}
	else	{
		basename = argv[4];
	}

	cerr << "new sim\n";
	Simulator* sim = new Simulator(datafile,treefile,paramfile,rho,pi,h,cpgrate);

	if (extbfile != "")	{
		cerr << "load ext b\n";
		sim->LoadExtB(extbfile);
	}
	cerr << "simu\n";
	sim->Simulate();
	sim->WriteSimu(basename);

}

