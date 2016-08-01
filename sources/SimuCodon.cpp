
#include "Tree.h"
#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>

#include "StringStreamUtils.h"

class Simulator : public NewickTree {

	public:

	Simulator(string datafile, int inNsite, string treefile, string paramfile)	{


		nucdata = new FileSequenceAlignment(datafile);
		int fromNsite = nucdata->GetNsite();
		// 3* : nuc -> codon
		// +1 : final stop codon (removed at the end)
		if (inNsite == -1)	{
			Nsite = fromNsite;
		}
		else	{
			Nsite = 3*(inNsite+1);
		}

		codondata = new CodonSequenceAlignment(nucdata);
		protdata = new ProteinSequenceAlignment(codondata);

		taxonset = nucdata->GetTaxonSet();
		statespace = nucdata->GetStateSpace();
		codonstatespace = new CodonStateSpace(Universal);

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);
		Ntaxa = taxonset->GetNtaxa();

		Create();

		// get parameters from file
		ifstream prmis(paramfile.c_str());

		// read file
		string tmp;

		betaav = 0;
		betacm = 0;

		epibeta = 0;

		// alpha
		prmis >> tmp;

		if (tmp == "Structure")	{
			ReadStructure(prmis);
			SetStructure();
			prmis >> tmp;
		}

		if (tmp == "Epistasy")	{
			prmis >> epibeta >> epifrac >> episigma;
			if (epibeta)	{
				SetEpistasy();
			}
			prmis >> tmp;
		}

		if (tmp == "SiteFitness")	{
			prmis >> tmp;
			if (tmp == "Random")	{
				prmis >> alphasigma;
				MakeRandomAlpha();
			}
			else if (tmp == "File")	{
				string alphafile;
				prmis >> alphafile;
				ReadAlphaFromFile(alphafile);
			}
			else	{
				cerr << "error: site fitness model not recognized\n";
				cerr << tmp << '\n';
				exit(1);
			}
			prmis >> tmp;
		}

		if (tmp != "Fluctuations")	{
			cerr << "error: missing Fluctuations keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> flucrho >> flucstat >> flucsigma >> Khidden;
		SetupHiddenZ();
		MakeRandomFluctuations();

		prmis >> tmp;
		if (tmp != "MutRates")	{
			cerr << "error: missing MutRates keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
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

		for (int i=0; i<Nnuc; i++)	{
			double total = 0;
			for (int j=0; j<Nnuc; j++)	{
				total += mutrate[i][j];
			}
			mutrate[i][i] = -total;
		}

		CalculateNucStat();

		// calculate total mut rate
		totmutrate = 0;
		for (int i=0; i<Nnuc; i++)	{
			cerr << nucstat[i] << '\t';
			for (int j=0; j<Nnuc; j++)	{
				if (i != j)	{
					totmutrate += nucstat[i] * mutrate[i][j];
				}
			}
		}
		cerr << '\n';

		prmis >> tmp;
		if (tmp != "Ne")	{
			cerr << "error: missing Ne keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> Ne;

		prmis >> tmp;
		if (tmp != "Scale")	{
			cerr << "error: missing Scale keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> scale;
		prmis >> tmp;
		if (tmp != "CodonUsage")	{
			cerr << "error: missing CodonUSage keyword\n";
			cerr << tmp << '\n';
			exit(1);
		}
		prmis >> cusigma;
		MakeRandomCodonUsage();
		totlength = 0;
		RecursiveSetBranchLengths(GetRoot());

		cerr << "total subs length per codon: " << mu * totlength * totmutrate * 3 << '\n';
		cerr << "mutrate : " << totmutrate << '\n';
		cerr << "length  : " << totlength << '\n';
	}

	void ReadStructure(istream& is)	{

		is >> betaav;
		is >> betacm;

		int Ngene;
		is >> Ngene;
		TotNsite = 1;

		string seqname;
		is >> seqname;
		ifstream sis(seqname.c_str());
		int* nsite = new int[Ngene];
		string temp;
		for (int gene=0; gene<Ngene; gene++)	{
			sis >> temp >> nsite[gene] >> temp;
			TotNsite += nsite[gene];
		}
		TotNsite++;
		cerr << "structure: total number of sites: " << TotNsite << '\n';

		string avname;
		is >> avname;
		ifstream ais(avname.c_str());
		for (int i=0; i<5; i++)	{
			ReadLine(ais);
		}

		av = new int[TotNsite];
		// start and stop codons
		av[0] = 0;
		av[TotNsite-1] = 0;

		int index = 1;
		for (int gene=0; gene<Ngene; gene++)	{
			int checknsite;
			ais >> temp >> checknsite;
			if (checknsite != nsite[gene])	{
				cerr << "error when reading av file\n";
				exit(1);
			}
			for (int i=0; i<nsite[gene]; i++)	{
				ais >> av[index];
				index++;
			}
		}

		string cmname;
		is >> cmname;
		ifstream cis(cmname.c_str());

		ReadLine(cis);

		int MaxNconc = 20;
		ncont = new int[TotNsite];
		cm = new int*[TotNsite];
		for (int i=0; i<TotNsite; i++)	{
			ncont[i] = 0;
			cm[i] = new int[MaxNconc];
		}
		int totnsite = 1;
		for (int gene=0; gene<Ngene; gene++)	{
			string temp;
			int checknsite, ncontact;
			cis >> temp >> checknsite >> ncontact >> temp >> temp >> temp;
			if (checknsite != nsite[gene])	{
				cerr << "error when reading structure: non matching gene length\n";
				exit(1);
			}

			for (int j=0; j<ncontact; j++)	{
				int tmp1 = 0;
				int tmp2 = 0;
				cis >> tmp1 >> tmp2;
				// cerr << totnsite << '\t' << tmp1 << '\t' << tmp2 << '\t' << TotNsite << '\n';
				if (ncont[totnsite+tmp1] == MaxNconc)	{
					cerr << "contact overflow\n";
					exit(1);
				}
				if (ncont[totnsite+tmp2] == MaxNconc)	{
					cerr << "contact overflow\n";
					exit(1);
				}
				cm[totnsite+tmp1][ncont[totnsite+tmp1]] = totnsite+tmp2;
				ncont[totnsite+tmp1]++;
				cm[totnsite+tmp2][ncont[totnsite+tmp2]] = totnsite+tmp2;
				ncont[totnsite+tmp2]++;
			}

			totnsite += nsite[gene];
		}

		string avpotname;
		is >> avpotname;

		ifstream avpis(avpotname.c_str());
		int Nacc;
		avpis >> temp >> Nacc;

		avpot = new double*[Nacc];
		for (int acc=0; acc<Nacc; acc++)	{
			avpot[acc] = new double[Naa];
			for (int a=0; a<Naa; a++)	{
				avpis >> avpot[acc][a];
			}
		}
		
		string cmpotname;
		is >> cmpotname;

		ifstream cmpis(cmpotname.c_str());
		cmpis >> temp >> temp;
		cmpot = new double*[Naa];
		for (int i=0; i<Naa; i++)	{
			cmpot[i] = new double[Naa];
			for (int j=0; j<Naa; j++)	{
				cmpis >> cmpot[i][j];
			}
		}

		delete[] nsite;
	}

	void SetStructure()	{

		for (int codpos=1; codpos<Nsite/3-1; codpos++)	{
			for (int a=0; a<Naa; a++)	{
				alpha[codpos][a] = - betaav * avpot[av[codpos]][a];
			}
		}

	}

	void SetEpistasy()	{

		epicont = new int*[Nsite/3];
		for (int i=0; i<Nsite/3; i++)	{
			epicont[i] = new int[Nsite/3];
			for (int j=0; j<Nsite/3; j++)	{
				epicont[i][j] = 0;
			}
		}

		epipot = new double***[Nsite/3];
		for (int i=0; i<Nsite/3; i++)	{
			epipot[i] = new double**[Nsite/3];
			for (int j=0; j<Nsite/3; j++)	{
				epipot[i][j] = new double*[Naa];
				for (int a=0; a<Naa; a++)	{
					epipot[i][j][a] = new double[Naa];
				}
			}
		}
		
		for (int i=0; i<Nsite/3; i++)	{
			for (int j=0; j<Nsite/3; j++)	{
				for (int a=0; a<Naa; a++)	{
					for (int b=0; b<Naa; b++)	{
						epipot[i][j][a][b] = 0;
					}
				}
			}
		}

		for (int i=1; i<Nsite/3-1; i++)	{
			for (int j=i+1; j<Nsite/3-1; j++)	{
				if (Random::Uniform() < epifrac)	{
					epicont[i][j] = epicont[j][i] = 1;
					for (int a=0; a<Naa; a++)	{
						for (int b=0; b<Naa; b++)	{
							epipot[i][j][a][b] = episigma * Random::sNormal();
							epipot[j][i][b][a] = epipot[i][j][a][b];
						}
					}
				}
			}
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
		return GetTree()->GetNodeName(link);
		/*
		ostringstream s;
		if (link->isLeaf())	{
			s << link->GetNode()->GetName();
		}
		// s << exp(instantlogB.find(link->GetNode())->second);
		return s.str();
		*/
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
			if (! l)	{
				cerr << "null branch length : " << from->GetBranch()->GetName() << '\n';
				exit(1);
			}
			bl[from->GetBranch()] = l * scale;
			totlength += l*scale;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetBranchLengths(link->Out());
		}
	}

	void Simulate()	{

		count = 0;
		dscount = 0;
		dncount = 0;
		hiddencount = 0;

		// RecursiveSimulate(GetRoot());

		// add a stem before the root
		Link* link = new Link();
		Node* node = new Node();
		Branch* branch = new Branch();
		link->SetBranch(branch);
		GetRoot()->SetBranch(branch);
		link->SetNode(node);
		bl[branch] = 10000;
		if ((!flucsigma) && (! betacm) && (! epibeta))	{
			bl[branch] = 1;
		}
		link->SetOut(GetRoot());
		GetRoot()->SetOut(link);

		RootSimulate();
		StoreCurrentValues(link);

		RecursiveSimulate(GetRoot());
		GetRoot()->SetOut(GetRoot());
		GetRoot()->SetBranch(0);
		delete link;
		delete node;
		delete branch;

		NormalizeFitnessStats();
	}

	void RecursiveSimulate(const Link* from)	{

		// cerr << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\n';

		cerr << '.';

		/*
		if (from->isRoot())	{
			RootSimulate(from);
		}
		else	{
			BranchSimulate(from);
		}
		*/

		BranchSimulate(from);

		StoreCurrentValues(from);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSimulate(link->Out());
		}

	}

	void StoreCurrentValues(const Link* from)	{

		// store new values
		int* newseq = new int[Nsite];
		int* newhidden = new int[Nsite/3];
		for (int i=0; i<Nsite; i++)	{
			newseq[i] = currentseq[i];
		}
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			newhidden[codpos] = currenthidden[codpos];
		}
		nodeseq[from->GetNode()] = newseq;
		nodehidden[from->GetNode()] = newhidden;
	}

	void RootSimulate()	{

		// draw hidden states
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			if (hiddenz[codpos])	{
				currenthidden[codpos] = (int) (Khidden * Random::Uniform());
			}
			else	{
				currenthidden[codpos] = 0;
			}
		}

		// draw codon sequence
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			DrawCodonFromStat(codpos);
		}
	}

	void BranchSimulate(const Link* from)	{

		double l = bl[from->GetBranch()];

		// define current state as equal to state at ancestral node;

		int* ancseq = nodeseq[from->Out()->GetNode()];
		int* anchidden = nodehidden[from->Out()->GetNode()];

		for (int i=0; i<Nsite; i++)	{
			currentseq[i] = ancseq[i];
		}
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			currenthidden[codpos] = anchidden[codpos];
		}

		// setup array of rates
		UpdateAllSubRates();

		double t = 0;

		while (t < l)	{

			double totsubrate = 0;
			for (int i=0; i<Nsite; i++)	{
				totsubrate += subrate[i][Nnuc];
			}
			double flucrate = Nhidden * flucrho;
			double totrate = totsubrate + flucrate;

			double dt = Random::sExpo() / totrate;
			t += dt;

			IncrementTimeCounters(dt);

			if (t < l)	{

				double u = totrate * Random::Uniform();

				if (u > totsubrate)	{
					u-= totsubrate;
					// choose site under fluctuating selection

					int j = (int) (u / flucrho);
					if (j >= Nhidden)	{
						cerr << "error when choosing fluct. pos: " << j << '\t' << Nhidden << '\n';
						exit(1);
					}
					// change category at this site
					int pos = 0;
					if (hiddenz[pos])	{
						j--;
					}
					while ((pos < Nsite/3) && (j>0))	{
						pos++;
						if (pos == Nsite/3)	{
							cerr << "error when choosing fluctuating coding position\n";
							exit(1);
						}
						if (hiddenz[pos])	{
							j--;
						}
					}

					int newh = (int) ((Khidden-1) * Random::Uniform());
					if (newh >= currenthidden[pos])	{
						newh ++;
					}
					currenthidden[pos] = newh;

					// recalculate rates at all three sites of the same codon
					UpdateSubRate(3*pos);
					UpdateSubRate(3*pos + 1);
					UpdateSubRate(3*pos + 2);

					UpdateFitnessStats(pos);

					// UpdateAllSubRates();

					if (from != GetRoot())	{
						hiddencount++;
					}
				}
				else	{

					// choose position

					double prevtot = 0;
					double tot = subrate[0][Nnuc];
					int pos = 0;
					while ((pos < Nsite) && (u > tot))	{
						pos++;
						prevtot = tot;
						tot += subrate[pos][Nnuc];
					}
					if (pos == Nsite)	{
						cerr << "error: overflow when choosing substituting position\n";
						exit(1);
					}

					// choose final nuc state
					u -= prevtot;
					tot = subrate[pos][0];
					int final = 0;
					while ((final < Nnuc) && (u > tot))	{
						final++;
						tot += subrate[pos][final];
					}
					if (final == Nnuc)	{
						cerr << "error: overflow when choosing final nucleotide\n";
						exit(1);
					}

					int codpos = pos/3;

					int a1 = currentseq[3*codpos];
					int a2 = currentseq[3*codpos + 1];
					int a3 = currentseq[3*codpos + 2];
					int c = codonstatespace->GetCodonFromDNA(a1,a2,a3);
					currentseq[pos] = final;
					int b1 = currentseq[3*codpos];
					int b2 = currentseq[3*codpos + 1];
					int b3 = currentseq[3*codpos + 2];
					int d = codonstatespace->GetCodonFromDNA(b1,b2,b3);
					if (from != GetRoot())	{
						if (codonstatespace->Synonymous(c,d))	{
							dscount++;
						}
						else	{
							dncount++;
						}
					}

					// recalculate rates at all three sites of the same codon
					UpdateSubRate(3*codpos);
					UpdateSubRate(3*codpos + 1);
					UpdateSubRate(3*codpos + 2);

					UpdateFitnessStats(codpos);

					// if codon first site: recalculate rate at previous nuc site
					if (pos % 3 == 0)	{
						if (pos)	{
							UpdateSubRate(pos-1);
							UpdateFitnessStats(codpos-1);
						}
					}

					// if codon last site: recalculate rate at next nuc site
					if (pos % 3 == 2)	{
						if (pos < Nsite-1)	{
							UpdateSubRate(pos+1);
							UpdateFitnessStats(codpos+1);
						}
					}

					// if epistasy
					// recalculate rates at all codons in contact
					if (betacm)	{
						for (int j=0; j<ncont[codpos]; j++)	{
							/*
							if ((cm[codpos][j] < 0) || (cm[codpos][j] >= Nsite/3))	{
								cerr << "contact map : " << cm[codpos][j] << '\t' << Nsite << '\n';
								exit(1);
							}
							*/
							if (cm[codpos][j] < Nsite/3)	{
								int codpos2 = cm[codpos][j];
								UpdateSubRate(3*codpos2);
								UpdateSubRate(3*codpos2+1);
								UpdateSubRate(3*codpos2+2);

								UpdateFitnessStats(codpos2);
							}
						}
					}

					if (epibeta)	{
						for (int codpos2=0; codpos2<Nsite/3; codpos2++)	{
							if (epicont[codpos][codpos2])	{
								UpdateSubRate(3*codpos2);
								UpdateSubRate(3*codpos2+1);
								UpdateSubRate(3*codpos2+2);

								UpdateFitnessStats(codpos2);
							}
						}
					}

					if (from != GetRoot())	{
						count++;
					}
					// UpdateAllSubRates();
					/*
					if (! (count % 30))	{
						cout << GetTotalFitness() << '\n';
					}
					*/
				}

				// CheckSubRates();
			}
		}
		// exit(1);
	}

	double GetTotalFitness()	{

		double total = 0;
		for (int i=0; i<Nsite/3; i++)	{
			total += GetFitness(i,0);
			total += GetFitness(i,1);
		}
		return 0.5 * total;
	}
		
	void UpdateAllSubRates()	{
		for (int i=0; i<Nsite; i++)	{
			UpdateSubRate(i);
		}
	}

	double GetMutRate(int init, int final, int prev, int next)	{
		if ((prev == 1) && (init == 2) && (final == 0))	{
			return mu * mutrate[2][0] * cpgrate;
		}
		if ((next == 3) && (init == 1) && (final == 3))	{
			return mu * mutrate[1][3] * cpgrate;
		}
		return mu * mutrate[init][final];
	}

	void CheckSubRates()	{
		for (int i=0; i<Nsite; i++)	{
			CheckSubRate(i);
		}
	}

	void CheckSubRate(int site)	{

		double bk[Nnuc];
		for (int n=0; n<Nnuc; n++)	{
			bk[n] = subrate[site][n];
		}
		UpdateSubRate(site);
		for (int n=0; n<Nnuc; n++)	{
			if (fabs(bk[n] - subrate[site][n]) > 1e-6)	{
				cerr << "error: subrates not updated : " << bk[n] << '\t' << subrate[site][n] << '\n';
				exit(1);
			}
		}
	}

	void UpdateSubRate(int site)	{

		// first and last codon position stay at ATG and TGA resp.
		if ((site < 3) || (site >= Nsite-3))	{
			for (int n=0; n<Nnuc; n++)	{
				subrate[site][n] = 0;
			}
			subrate[site][Nnuc] = 0;
		}
		else	{
			int codpos = site / 3;
			int current = currentseq[site];
			subrate[site][current] = 0;
			double tot = 0;
			for (int n=0; n<Nnuc; n++)	{
				if (n != current)	{

					double tmp = GetMutRate(current,n,currentseq[site-1],currentseq[site+1]);
					double initfitness = GetFitness(codpos);
					currentseq[site] = n;
					double finalfitness = GetFitness(codpos);
					currentseq[site] = current;
					double fixprob = GetFixProb(Ne*(finalfitness-initfitness));
					tmp *= fixprob;
					subrate[site][n] = tmp;
					tot += tmp;
				}
			}
			subrate[site][Nnuc] = tot;
		}
	}

	void IncrementTimeCounters(double dt)	{

		tottime += dt;
		for (int i=0; i<Nsite/3; i++)	{
			currenttime[i] += dt;
		}
	}

	void CreateFitnessStats()	{

		currenttime = new double[Nsite/3];
		currentdist = new double[Nsite/3];
		for (int i=0; i<Nsite/3; i++)	{
			currenttime[i] = 0;
			currentdist[i] = 0;
		}

		currentfitness = new double*[Nsite/3];
		meanfitness = new double*[Nsite/3];
		varfitness = new double*[Nsite/3];
		for (int i=0; i<Nsite/3; i++)	{
			currentfitness[i] = new double[Naa];
			meanfitness[i] = new double[Naa];
			varfitness[i] = new double[Naa];
			for (int a=0; a<Naa; a++)	{
				currentfitness[i][a] = 0;
				meanfitness[i][a] = 0;
				varfitness[i][a] = 0;
			}
		}
	}

	void UpdateAllFitnessStats()	{
		for (int i=0; i<Nsite/3; i++)	{
			UpdateFitnessStats(i);
		}
	}

	void UpdateFitnessStats(int codpos, double beta=1)	{

		double* oldfit = new double[Naa];
		for (int a=0; a<Naa; a++)	{
			oldfit[a] = currentfitness[codpos][a];
		}
		double* newfit = currentfitness[codpos];

		double t = currenttime[codpos];

		currenttime[codpos] = 0;

		for (int a=0; a<Naa; a++)	{
			newfit[a] = GetAAFitness(codpos,a, beta);
		}

		// mean
		for (int a=0; a<Naa; a++)	{
			meanfitness[codpos][a] += t * oldfit[a];
			varfitness[codpos][a] += t * oldfit[a] * oldfit[a];
		}
		double dist = 0;
		for (int a=0; a<Naa; a++)	{
			double tmp = newfit[a] - oldfit[a];
			dist += tmp * tmp;
		}
		currentdist[codpos] += sqrt(dist);

		delete[] oldfit;

	}

	void NormalizeFitnessStats()	{

		globalmeanfitness = 0;
		globalvarfitness = 0;
		globaldist = 0;
		for (int i=0; i<Nsite/3; i++)	{
			for (int a=0; a<Naa; a++)	{
				meanfitness[i][a] /= tottime;
				varfitness[i][a] /= tottime;
				varfitness[i][a] -= meanfitness[i][a] * meanfitness[i][a];

				globalmeanfitness += meanfitness[i][a];
				globalvarfitness += varfitness[i][a];
			}

			currentdist[i] /= tottime;
			globaldist += currentdist[i];
		}

		globaldist /= Nsite/3/Naa;
		globalmeanfitness /= Nsite/3/Naa;
		globalvarfitness /= Nsite/3/Naa;

	}

	double GetFitness(int codpos, double beta = 1)	{

		int nuc1 = currentseq[3*codpos];
		int nuc2 = currentseq[3*codpos + 1];
		int nuc3 = currentseq[3*codpos + 2];
		if (codonstatespace->CheckStop(nuc1,nuc2,nuc3))	{
			return -200;
		}
		int codon = codonstatespace->GetCodonFromDNA(nuc1,nuc2,nuc3);
		int aa = codonstatespace->Translation(codon);

		double fitness = GetAAFitness(codpos,aa,beta);
		fitness += cu[codon];

		return fitness;
	}

	double GetAAFitness(int codpos, int aa, double beta = 1)	{

		double fitness = alpha[codpos][aa];
		if (hiddenz[codpos])	{
			fitness += dalpha[codpos][currenthidden[codpos]][aa];
		}

		if (betacm && beta)	{
			for (int j=0; j<ncont[codpos]; j++)	{
				int codpos2 = cm[codpos][j];
				if ((codpos2 < 0) || (codpos2 >= TotNsite))	{
					cerr << "contact map : " << codpos2 << '\t' << Nsite << '\n';
					exit(1);
				}
				if (codpos2 < Nsite/3)	{
					int nuc1 = currentseq[3*codpos2];
					int nuc2 = currentseq[3*codpos2 + 1];
					int nuc3 = currentseq[3*codpos2 + 2];
					if (codonstatespace->CheckStop(nuc1,nuc2,nuc3))	{
						cerr << "error: position in contact has stop codon\n";
						exit(1);
					}
					int codon = codonstatespace->GetCodonFromDNA(nuc1,nuc2,nuc3);
					int bb = codonstatespace->Translation(codon);
					fitness -= betacm * cmpot[aa][bb];
				}
			}
		}

		if (epibeta)	{
			for (int codpos2=0; codpos2<Nsite/3; codpos2++)	{
				if (epicont[codpos][codpos2])	{
					int nuc1 = currentseq[3*codpos2];
					int nuc2 = currentseq[3*codpos2 + 1];
					int nuc3 = currentseq[3*codpos2 + 2];
					if (codonstatespace->CheckStop(nuc1,nuc2,nuc3))	{
						cerr << "error: position in contact has stop codon\n";
						exit(1);
					}
					int codon = codonstatespace->GetCodonFromDNA(nuc1,nuc2,nuc3);
					int bb = codonstatespace->Translation(codon);
					fitness += epibeta * epipot[codpos][codpos2][aa][bb];
				}
			}
		}


		return fitness;
	}

	double GetFixProb(double s)	{

		double smax = 100;
		if (s > smax)	{
			s = smax;
		}
		if (s < -smax)	{
			s = -smax;
		}
		double fixprob = 1;
		if (fabs(s) > 1e-6)	{
			fixprob = s / (1 - exp(-s));
		}
		return fixprob;
	}

	void DrawCodonFromStat(int codpos)	{

		// assumes no epistasy
		// ideally, if epistasy is present, should run simulator forward for some time to equilibrate sequence

		if (codpos == 0)	{

			currentseq[0] = 0;
			currentseq[1] = 3;
			currentseq[2] = 2;
		}
		else if (codpos == Nsite/3 -1)	{

			currentseq[3*codpos] = 3;
			currentseq[3*codpos+1] = 2;
			currentseq[3*codpos+2] = 0;
		}
		else	{

			double* cumulprob = new double[64];
			double tot = 0;
			int i = 0;
			for (int n1=0; n1<Nnuc; n1++)	{
				for (int n2=0; n2<Nnuc; n2++)	{
					for (int n3=0; n3<Nnuc; n3++)	{
						currentseq[3*codpos] = n1;
						currentseq[3*codpos + 1] = n2;
						currentseq[3*codpos + 2] = n3;
						double f = GetFitness(codpos, 0);
						double tmp = nucstat[n1] * nucstat[n2] * nucstat[n3] * exp(f);
						tot += tmp;
						cumulprob[i] = tot;
						i++;
					}
				}
			}
			
			double u = tot * Random::Uniform();
			i = 0;
			while ((i<64) && (u>cumulprob[i]))	{
				i++;
			}
			if (i == 64)	{
				cerr << "error in draw codon from stat: overflow\n";
				exit(1);
			}
			int nuc3 = i % 4;
			i /= 4;
			int nuc2 = i % 4;
			i /= 4;
			int nuc1 = i;
			if (codonstatespace->CheckStop(nuc1,nuc2,nuc3))	{
				cerr << "error in draw codon from stat: stop codon\n";
				exit(1);
			}

			currentseq[3*codpos] = nuc1;
			currentseq[3*codpos + 1] = nuc2;
			currentseq[3*codpos + 2] = nuc3;
		}
	}

	void WriteSimu(string basename)	{

		// dataset
		int** data = new int*[Ntaxa];
		string* names = new string[Ntaxa];
		int n = 0;
		Nsite -= 3;
		RecursiveMakeData(GetRoot(),data,names,n);
		if (n != Ntaxa)	{
			cerr << "error when making sequence alignment: wrong number of taxa : " << n << '\t' << " instead of " << Ntaxa << '\n';
			exit(1);
		}

		cerr << "make new ali\n";
		SequenceAlignment* simuali = new SequenceAlignment(data,names,Nsite,statespace,taxonset);

		/*
		ofstream allos((basename + ".ali").c_str());
		simuali->ToStream(allos);
		*/

		CodonSequenceAlignment* codonali = new CodonSequenceAlignment(simuali);
		CodonSequenceAlignment* codontemplate = new CodonSequenceAlignment(nucdata);
		// codonali->GetMissingCellsFromTemplate(codontemplate);

		codonali->Mask(codondata);

		ofstream callos((basename + ".ali").c_str());
		codonali->ToStream(callos);

	
		ProteinSequenceAlignment* protali = new ProteinSequenceAlignment(codonali);
		ofstream protos((basename + "_prot.ali").c_str());
		protali->ToStream(protos);

		double** empfreq = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			empfreq[i] = new double[Naa];
		}
		protali->GetSiteEmpiricalFreq(empfreq);
		ofstream fos((basename + ".freq").c_str());
		fos << Nsite << '\t' << Naa << '\n';
		for (int i=0; i<Nsite; i++)	{
			for (int a=0; a<Naa; a++)	{
				fos << empfreq[i][a] << '\t';
			}
			fos << '\n';
		}
		fos.close();

		ofstream sos((basename + ".summary").c_str());
		
		sos << "mean number of subs per site : " << ((double) count) / Nsite * 3 << '\n';
		sos << "mean number of syns per site : " << ((double) dscount) / Nsite * 3 << '\n';
		sos << "mean number of reps per site : " << ((double) dncount) / Nsite * 3<< '\n';
		sos << "mean number of fluctuations per site : " << ((double) hiddencount) / Nsite * 3 << '\n';
		sos << "mean diversity : " << protali->GetMeanDiversity() << '\n';
		sos << "ref  diversity : " << protdata->GetMeanDiversity() << '\n';
		sos << "global stdev fitness : " << sqrt(globalvarfitness) << '\n';
		sos << "rate of change      : " << globaldist / sqrt(globalvarfitness) / tottime / 3.0 / mu << '\n';

		cerr << "mean number of subs per site : " << ((double) count) / Nsite * 3 << '\n';
		cerr << "mean number of syns per site : " << ((double) dscount) / Nsite * 3 << '\n';
		cerr << "mean number of reps per site : " << ((double) dncount) / Nsite * 3 << '\n';
		cerr << "mean number of fluctuations per site : " << ((double) hiddencount) / Nsite * 3 << '\n';
		cerr << "mean diversity : " << protali->GetMeanDiversity() << '\n';
		cerr << "ref  diversity : " << protdata->GetMeanDiversity() << '\n';
		cerr << "global stdev fitness : " << sqrt(globalvarfitness)  / Naa  << '\n';
		cerr << "rate of change      : " << globaldist / sqrt(globalvarfitness) / tottime  / 3.0 / mu << '\n';

		delete[] data;
		delete[] names;
	}

	void RecursiveMakeData(const Link* from, int** data, string* names, int& n)	{
		if (from->isLeaf())	{
			names[n] = from->GetNode()->GetName();
			data[n] = nodeseq[from->GetNode()];
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveMakeData(link->Out(), data, names, n);
		}
	}

	void SetupHiddenZ()	{

		Nhidden = 0;
		hiddenz = new int[Nsite/3];
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			if (Random::Uniform() < flucstat)	{
				hiddenz[codpos] = 1;
				Nhidden++;
			}
			else	{
				hiddenz[codpos] = 0;
			}
		}
		cerr << "total number of fluctuating sites: " << Nhidden << '\n';

	}

	void MakeRandomCodonUsage()	{

		cu = new double[codonstatespace->GetNstate()];
		for (int c=0; c<codonstatespace->GetNstate(); c++)	{
			cu[c] = cusigma * Random::sNormal();
		}
	}

	void MakeRandomFluctuations()	{

		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			for (int k=0; k<Khidden; k++)	{
				if (hiddenz[codpos])	{
					for (int a=0; a<Naa; a++)	{
						dalpha[codpos][k][a] = flucsigma * Random::sNormal();
					}
				}
				else	{
					for (int a=0; a<Naa; a++)	{
						dalpha[codpos][k][a] = 0;
					}
				}
			}
		}
	}

	void MakeRandomAlpha()	{

		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			for (int a=0; a<Naa; a++)	{
				alpha[codpos][a] = alphasigma * Random::sNormal();
			}
		}
	}

	void ReadAlphaFromFile(string filename)	{

		ifstream is(filename.c_str());
		int nsite, nstate;
		is >> nsite >> nstate;
		if (nsite < Nsite/3)	{
			cerr << "error: too small\n";
			exit(1);
		}
		if (nstate != Naa)	{
			cerr << "error when reading " << filename << '\n';
			exit(1);
		}
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			double tot = 0;
			for (int a=0; a<Naa; a++)	{
				double tmp;
				is >> tmp;
				tot += tmp;
				if (tmp <= 0)	{
					cerr << "error: negative fitness\n";
					exit(1);
				}
				alpha[codpos][a] = log(tmp);
			}
			if (fabs(tot - 1) > 1e-4)	{
				cerr << "error: profile does not sum to one : " << tot << '\n';
				exit(1);
			}
		}
	}

	void CalculateNucStat()	{

		double r[Nnuc][Nnuc];
		double r2[Nnuc][Nnuc];
		double max = 0;
		for (int i=0; i<Nnuc; i++)	{
			if (max < fabs(mutrate[i][i]))	{
				max = fabs(mutrate[i][i]);
			}
		}
		cerr << "max : " << max << '\n';
		for (int i=0; i<Nnuc; i++)	{
			for (int j=0; j<Nnuc; j++)	{
				r[i][j] = mutrate[i][j] / max;
			}
		}
		for (int i=0; i<Nnuc; i++)	{
			r[i][i] += 1;
		}

		double diff = 1;
		while (diff > 1e-6)	{
			for (int i=0; i<Nnuc; i++)	{
				for (int j=0; j<Nnuc; j++)	{
					double tmp = 0;
					for (int k=0; k<Nnuc; k++)	{
						tmp += r[i][k] * r[k][j];
					}
					r2[i][j] = tmp;
				}
			}
			for (int i=0; i<Nnuc; i++)	{
				for (int j=0; j<Nnuc; j++)	{
					double tmp = 0;
					for (int k=0; k<Nnuc; k++)	{
						tmp += r2[i][k] * r2[k][j];
					}
					r[i][j] = tmp;
				}
			}
			diff = 0;
			for (int i=0; i<Nnuc; i++)	{
				for (int j=0; j<Nnuc; j++)	{
					for (int k=j+1; k<Nnuc; k++)	{
						double tmp = fabs(r[j][i] - r[k][i]);
						if (diff < tmp)	{
							diff = tmp;
						}
					}
				}
			}
		}
			
		for (int i=0; i<Nnuc; i++)	{
			nucstat[i] = r[0][i];
		}
	}

	void Create()	{

		currentseq = new int[Nsite];
		currenthidden = new int[Nsite/3];
		hiddenz = new int[Nsite/3];

		subrate = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			subrate[i] = new double[Nnuc+1];
		}

		alpha = new double*[Nsite/3];
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			alpha[codpos] = new double[Naa];
			for (int a=0; a<Naa; a++)	{
				alpha[codpos][a] = 0;
			}
		}

		int Kmax = 10;
		dalpha = new double**[Nsite/3];
		for (int codpos=0; codpos<Nsite/3; codpos++)	{
			dalpha[codpos] = new double*[Kmax];
			for (int k=0; k<Kmax; k++)	{
				dalpha[codpos][k] = new double[Naa];
			}
		}

		CreateFitnessStats();
	}

	// statistical potentials and related structural desciptions
	int TotNsite;
	int* av;
	int* ncont;
	int** cm;
	double** cmpot;
	double** avpot;

	double betaav;
	double betacm;

	double epibeta;
	double epifrac;
	double episigma;
	int** epicont;
	double**** epipot;

	double** subrate;
	int count;
	int dscount;
	int dncount;
	int hiddencount;
	double Ne;
	
	double flucrho;
	double flucstat;
	double flucsigma;

	double alphasigma;

	double** alpha;
	double*** dalpha;

	int Nhidden;
	int Khidden;
	int* currenthidden;
	int* hiddenz;

	int Nsite;
	int* currentseq;

	double* currenttime;
	double** currentfitness;
	double** meanfitness;
	double** varfitness;
	double* currentdist;

	double tottime;

	double globalmeanfitness;
	double globalvarfitness;
	double globaldist;

	int Ntaxa;

	double nucstat[Nnuc];
	// absolute mutation rate;
	double mu;
	// all relative mutation rates
	double at2cg, at2gc, at2ta, cg2at, cg2gc, cg2ta;
	double mutrate[Nnuc][Nnuc];

	// including CpG effect
	double cpgrate;

	double cusigma;
	double* cu;

	Tree* tree;

	// with branch lengths, measured in time
	map<const Branch*, double> bl;
	double totlength;
	double totmutrate;

	// map<const Node*, double> instantNe;

	map<const Node*, int*> nodehidden;
	map<const Node*, int*> nodeseq;

	// a global time scale for easy translation
	double scale; // 100 Myr

	FileSequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	ProteinSequenceAlignment* protdata;
	TaxonSet* taxonset;
	StateSpace* statespace;
	CodonStateSpace* codonstatespace;

};

int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cerr << "simucodon datafile treefile paramfile Nsite basename\n";
		exit(1);
	}

	string datafile = argv[1];
	string treefile = argv[2];
	string paramfile = argv[3];
	int Nsite = atoi(argv[4]);
	string basename = argv[5];

	cerr << "new sim\n";
	Simulator* sim = new Simulator(datafile,Nsite,treefile,paramfile);

	cerr << "simu\n";
	sim->Simulate();

	cerr << "write simu\n";
	sim->WriteSimu(basename);

	cerr << "simu written in " << basename << '\n';
	cerr << '\n';

}

