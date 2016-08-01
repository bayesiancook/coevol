
#include "Random.h"
#include <map>
#include <iostream>
#include <sstream>

#include "BiologicalSequences.h"
#include "StringStreamUtils.h"


class StructStat	{

	public:

	int Nsite;
	int* av;
	int* ncont;
	int** cm;
	double** cmpot;
	double** avpot;

	double betaav;
	double betacm;

	double** freq;
	int* seq;

	double** alpha;

	StructStat(string filename)	{

		ifstream is(filename.c_str());
		ReadStructure(is);
		SetStructure();
	}

	void ReadStructure(istream& is)	{

		cerr << "Read\n";
		string temp;
		is >> temp;
		is >> betaav;
		is >> betacm;

		int Ngene;
		is >> Ngene;
		cerr << "Ngene : " << Ngene << '\n';
		Nsite = 0;

		string seqname;
		is >> seqname;
		ifstream sis(seqname.c_str());
		int* nsite = new int[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			sis >> temp >> nsite[gene] >> temp;
			Nsite += nsite[gene];
		}
		cerr << "structure: total number of sites: " << Nsite << '\n';

		string avname;
		is >> avname;
		ifstream ais(avname.c_str());
		for (int i=0; i<5; i++)	{
			ReadLine(ais);
		}

		av = new int[Nsite];

		int index = 0;
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
		ncont = new int[Nsite];
		cm = new int*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			ncont[i] = 0;
			cm[i] = new int[MaxNconc];
		}
		int totnsite = 0;
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
				// cerr << totnsite << '\t' << tmp1 << '\t' << tmp2 << '\t' << Nsite << '\n';
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

		alpha = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			alpha[i] = new double[Naa];
		}
		is >> temp;
		if (temp == "SiteFitness")	{
			cerr << "site fitnesses\n";
			is >> temp;
			if (temp == "Random")	{
				double sigma;
				is >> sigma;
				MakeRandomAlpha(sigma);
			}
			else if (temp == "File")	{
				string alphafile;
				is >> alphafile;
				ReadAlphaFromFile(alphafile);
			}
			else	{
				cerr << "error: site fitness model not recognized\n";
				cerr << temp << '\n';
				exit(1);
			}
		}
		else	{
			cerr << "using solvent accessibility potential for site fitnesses\n";
			SetAlphaFromPot();
		}
		
	}

	void SetStructure()	{

		
		freq = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			freq[i] = new double[Naa];
			for (int a=0; a<Naa; a++)	{
				freq[i][a] = 0;
			}
		}

		seq = new int[Nsite];
		for (int i=0; i<Nsite; i++)	{
			double tot = 0;
			double cumul[Naa];
			for (int a=0; a<Naa; a++)	{
				double tmp = alpha[i][a];
				tot += exp(tmp);
				cumul[a] = tot;
			}
			double u = tot * Random::Uniform();
			int b = 0;
			while ((b<Naa) && (cumul[b] < u))	{
				b++;
			}
			seq[i] = b;
		}

	}

	void MakeRandomAlpha(double sigma)	{

		for (int i=0; i<Nsite; i++)	{
			for (int a=0; a<Naa; a++)	{
				alpha[i][a] = sigma * Random::sNormal();
			}
		}
	}

	void SetAlphaFromPot()	{

		for (int i=0; i<Nsite; i++)	{
			for (int a=0; a<Naa; a++)	{
				alpha[i][a] = - betaav * avpot[av[i]][a];
			}
		}
	}

	void ReadAlphaFromFile(string filename)	{

		ifstream is(filename.c_str());
		int nsite, nstate;
		is >> nsite >> nstate;
		if (nsite < Nsite)	{
			cerr << "error: not enough sites in " << filename << '\n';
			exit(1);
		}
		if (nstate != Naa)	{
			cerr << "error when reading " << filename << '\n';
			exit(1);
		}
		for (int i=0; i<Nsite; i++)	{
			for (int a=0; a<Naa; a++)	{
				double tmp;
				is >> tmp;
				if (tmp <= 0)	{
					cerr << "error: negative fitness\n";
					exit(1);
				}
				alpha[i][a] = log(tmp);
			}
		}
	}

	void Run(int burnin, int nrep)	{

		for (int rep=0; rep<nrep; rep++)	{

			for (int ii=0; ii<Nsite; ii++)	{
			int i = (Nsite * (Random::Uniform()));

			double tot = 0;
			double cumul[Naa];
			for (int a=0; a<Naa; a++)	{
				double tmp = alpha[i][a];
				// double tmp = betaav * avpot[av[i]][a];
				for (int k=0; k<ncont[i]; k++)	{
					int j = cm[i][k];
					tmp -= betacm * cmpot[a][seq[j]];
				}

				tot += exp(tmp);
				if (rep > burnin)	{
					freq[i][a] += exp(tmp);
				}
				cumul[a] = tot;
			}
			double u = tot * Random::Uniform();
			int b = 0;
			while ((b<Naa) && (cumul[b] < u))	{
				b++;
			}
			seq[i] = b;
			/*
			if (rep > burnin)	{
				freq[i][b] ++;
			}
			*/
			}

			cout << Energy() << '\n';
			cout.flush();
		}
	}

	double Energy()	{

		double tot = 0;
		for (int i=0; i<Nsite; i++)	{
			tot -= alpha[i][seq[i]];
			// tot += betaav * avpot[av[i]][seq[i]];
			for (int k=0; k<ncont[i]; k++)	{
				int j = cm[i][k];
				tot += 0.5 * betacm * cmpot[seq[i]][seq[j]];
			}
		}
		return tot;
	}

	void WriteStat(string outfile)	{

		ofstream os(outfile.c_str());

		os << Nsite << '\t' << Naa << '\n';
		for (int i=0; i<Nsite; i++)	{
			double tot = 0;
			for (int a=0; a<Naa; a++)	{
				tot += freq[i][a];
			}
			for (int a=0; a<Naa; a++)	{
				freq[i][a] /= tot;
				os << freq[i][a] << '\t';
			}
			os << '\n';
		}

	}
};

int main(int argc, char* argv[])	{

	int burnin = atoi(argv[1]);
	int nrep = atoi(argv[2]);
	string filename = argv[3];
	string outfile = argv[4];

	StructStat* s = new StructStat(filename);

	s->Run(burnin,nrep);
	s->WriteStat(outfile);

}
