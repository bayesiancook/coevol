
#include "Sample.h"
#include "UniformMixOmegaModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"

class MixOmegaSample : public Sample	{

	private:

	int burnindone;
	int myid;
	int nprocs;

	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	int clampdiag;
	int clamptree;
	int meanexp;

	GeneticCodeType type;

	int conjpath;

	bool normalise;
	int nrep;

	int df;
	int clampreg;

	double mappingfreq;
	double mappingevery;
	string omegafile;

	public:

	string GetModelType() {return modeltype;}

	MixOmegaModel* GetModel() {return (MixOmegaModel*) model;}

	MixOmegaSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag;
		is >> conjpath;
		is >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> type;
		is >> clamptree;
		is >> suffstatfile;
		is >> rootfile;
		is >> clampreg;
		is >> mappingfreq;
		is >> mappingevery;
		is >> omegafile;

		nprocs = 1;
		myid = 0;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "UNIMIXOMEGA")	{
			model = new MixOmegaModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,df,clampdiag,conjpath,clamptree,meanexp,normalise,nrep,suffstatfile,rootfile,clampreg,mappingfreq,mappingevery,omegafile,type,myid,nprocs,false);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);

		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		cerr << "read\n";

		// ofstream dos((GetName() + ".reg").c_str());

		int Ngene = GetModel()->GetNgene();
		int Ncont = GetModel()->GetNcont();

		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}


		double meansynreg[Ngene][Ncont];
		double varsynreg[Ngene][Ncont];
		double ppsynreg[Ngene][Ncont];
		double meanomreg[Ngene][Ncont];
		double varomreg[Ngene][Ncont];
		double ppomreg[Ngene][Ncont];

		for (int g=0; g<Ngene; g++)	{
			for (int c=0; c<Ncont; c++)	{
				meansynreg[g][c] = 0;
				varsynreg[g][c] = 0;
				ppsynreg[g][c] = 0;
				varomreg[g][c] = 0;
				meanomreg[g][c] = 0;
				ppomreg[g][c] = 0;
			}
		}

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			for (int g=0; g<Ngene; g++)	{
				for (int c=0; c<Ncont; c++)	{
					double tmp1 = GetModel()->GetSynRegCoeff(g,c);
					meansynreg[g][c] += tmp1;
					varsynreg[g][c] += tmp1 * tmp1;
					if (tmp1 > 0)	{
						ppsynreg[g][c] ++;
					}
					double tmp2 = GetModel()->GetRegCoeff(g,c);
					// dos << tmp2 << '\t';
					meanomreg[g][c] += tmp2;
					varomreg[g][c] += tmp2 * tmp2;
					if (tmp2 > 0)	{
						ppomreg[g][c] ++;
					}
				}
				// dos << '\n';
			}

			meansynrate->Add(GetModel()->GetProcess(),GetModel()->GetChronogram(),Ncont);
			meanomega->Add(GetModel()->GetProcess(),GetModel()->GetChronogram(),Ncont+1);
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetProcess(),GetModel()->GetChronogram(),k);
			}
		}

		cerr << '\n';
		cerr << "normalise\n";

		for (int g=0; g<Ngene; g++)	{
			for (int c=0; c<Ncont; c++)	{
				meansynreg[g][c] /= size;
				varsynreg[g][c] /= size;
				varsynreg[g][c] -= meansynreg[g][c] * meansynreg[g][c];
				ppsynreg[g][c] /= size;
				meanomreg[g][c] /= size;
				varomreg[g][c] /= size;
				varomreg[g][c] -= meanomreg[g][c] * meanomreg[g][c];
				ppomreg[g][c] /= size;
			}
		}

		ofstream os((GetName() + ".postmeanreg").c_str());
		for (int g=0; g<Ngene; g++)	{
			for (int c=0; c<Ncont; c++)	{
				os << GetModel()->GetGeneName(g) << '\t' << meansynreg[g][c] << '\t' << sqrt(varsynreg[g][c]) << '\t' << ppsynreg[g][c] << '\t' << meanomreg[g][c] << '\t' << sqrt(varomreg[g][c]) << '\t' << ppomreg[g][c] << '\n';
			}
		}

		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
		/*
		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);
		cerr << "reconstructed variation in Ks in " << name << ".postmeansynrate.tre\n";
		*/

		meanomega->Normalise();
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);
		/*
		cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
		meanomega->SetWithLeaf(true);
		meanomega->SetWithInternal(true);
		ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
		meanomega->Tabulate(ooos);
		ooos.close();
		*/

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
		}
		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	// bool nans = false;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if (s == "+stdev")	{
				printstdev = true;
			}
			else if (s == "-stdev")	{
				printstdev = false;
			}
			else if (s == "+ci")	{
				printci = true;
			}
			else if (s == "-ci")	{
				printci = false;
			}
			else if (s == "+leaf")	{
				withleaf = true;
			}
			else if (s == "-leaf")	{
				withleaf = false;
			}
			else if (s == "+internal")	{
				withinternal = true;
			}
			else if (s == "-internal")	{
				withinternal = false;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readcoevol [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	MixOmegaSample sample(name,burnin,every,until);
	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}



