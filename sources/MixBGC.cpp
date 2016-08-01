
#include "Chain.h"
#include "MixBGCModel.h"
#include "StringStreamUtils.h"

// int BGCMutSelSubMatrix::bgccount = 0;

class BGCChain : public Chain	{

	private:
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

	bool clampdiag;
	bool autoregressive;
	bool clamptree;
	bool meanexp;

	double priorsigma;
	double fixalpha;

	int conjpath;
	int contdatatype;

	bool normalise;
	int nrep;

	int nsplit;

	int df;

	int clampbgcoffset;
	int flexrho;

	double lambda;
	double at2cg;
	double at2gc;
	double gc2cg;

	double cg2at;
	double cg2gc;
	double cg2ta;
	int discn;
	int discgam;
	int triplet;
	int clampreg;

	public:

	string GetModelType() {return modeltype;}

	BGCChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inpriorsigma, int indf, string filename, bool inclampdiag, bool inautoregressive, int inconjpath, int incontdatatype, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int innsplit, string insuffstatfile, string inrootfile, double infixalpha, int inclampbgcoffset, int inflexrho, double inlambda, double inat2cg, double inat2gc, double ingc2cg, double incg2at, double incg2gc, double incg2ta, int indiscgam, int indiscn, int intriplet, int inclampreg, int force = 1)	{
		modeltype = "BGC";
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		calibfile = incalibfile;
		suffstatfile = insuffstatfile;
		rootfile = inrootfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		name = filename;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		priorsigma = inpriorsigma;
		fixalpha = infixalpha;
		df = indf;
		conjpath = inconjpath;
		contdatatype = incontdatatype;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		normalise = innormalise;
		nrep = innrep;
		nsplit = innsplit;
		clampbgcoffset = inclampbgcoffset;
		flexrho = inflexrho;
		lambda = inlambda;
		at2cg = inat2cg;
		at2gc = inat2gc;
		gc2cg = ingc2cg;
		cg2at = incg2at;
		cg2gc = incg2gc;
		cg2ta = incg2ta;
		discgam = indiscgam;
		discn = indiscn;
		triplet = intriplet;
		clampreg = inclampreg;
		New(force);
	}

	BGCChain(string filename)	{
		conjpath = 1;
		contdatatype = 0;
		name = filename;
		priorsigma = 1;
		fixalpha = 0;
		df = 2;
		clamptree = false;
		meanexp = false;
		autoregressive = false;
		nsplit = 1;
		cerr << "open\n";
		suffstatfile = "None";
		rootfile = "None";
		clampbgcoffset = 0;
		flexrho = 0;
		lambda = 2;
		at2cg = -1;
		at2gc = -1;
		cg2at = -1;
		cg2gc = -1;
		cg2ta = -1;
		discgam = 0;
		discn = 0;
		triplet = 0;
		clampreg = 0;
		Open();
	}

	void New(int force)	{
		if (modeltype == "BGC")	{
			model = new BGCModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,clampbgcoffset,flexrho,lambda,at2cg,at2gc,gc2cg,cg2at,cg2gc,cg2ta,discgam,discn,triplet,clampreg,true);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		Reset(force);
		// model->Update();
		cerr << "ln prob = " << model->GetLogProb() << "\n";
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag;
		is >> conjpath;
		is >> contdatatype;
		is >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> nsplit;
		is >> clamptree;
		is >> suffstatfile;
		is >> autoregressive;
		is >> rootfile;

		int check;
		is >> check;
		if (check)	{
			is >> fixalpha;
			is >> check;
			if (check)	{
				is >> clampbgcoffset;
				is >> flexrho;
				is >> check;
				if (check)	{
					is >> lambda;
					is >> at2cg;
					is >> at2gc;
					is >> gc2cg;
					is >> check;
					if (check)	{
						is >> discgam;
						is >> cg2at;
						is >> cg2gc;
						is >> cg2ta;
						is >> discn;
						is >> check;
						if (check)	{
							is >> triplet;
							is >> check;
							if (check)	{
								is >> clampreg;
								is >> check;
								if (check)	{
									cerr << "error when reading model\n";
									exit(1);
								}
							}
						}
					}
				}
			}
		}

		is >> every >> until >> size;

		if (modeltype == "BGC")	{
			model = new BGCModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,clampbgcoffset,flexrho,lambda,at2cg,at2gc,gc2cg,cg2at,cg2gc,cg2ta,discgam,discn,triplet,clampreg,true);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		cerr << "update\n";
		model->Update();
		cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << datafile << '\t' << treefile << '\t' << contdatafile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << chronoprior << '\t' << meanchi << '\t' << meanchi2 << '\n';
		param_os << clampdiag << '\n';
		param_os << conjpath << '\n';
		param_os << contdatatype << '\n';
		param_os << meanexp << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << df << '\n';
		param_os << priorsigma << '\n';
		param_os << nsplit << '\n';
		param_os << clamptree << '\n';
		param_os << suffstatfile << '\n';
		param_os << autoregressive << '\n';
		param_os << rootfile << '\n';
		param_os << 1 << '\n';
		param_os << fixalpha << '\n';
		param_os << 1 << '\n';
		param_os << clampbgcoffset << '\n';
		param_os << flexrho << '\n';
		param_os << 1 << '\n';
		param_os << lambda << '\n';
		param_os << at2cg << '\n';
		param_os << at2gc << '\n';
		param_os << gc2cg << '\n';
		param_os << 1 << '\n';
		param_os << discgam << '\n';
		param_os << cg2at << '\n';
		param_os << cg2gc << '\n';
		param_os << cg2ta << '\n';
		param_os << discn << '\n';
		param_os << 1 << '\n';
		param_os << triplet << '\n';
		param_os << 1 << '\n';
		param_os << clampreg << '\n';
		param_os << 0 << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';
		model->ToStream(param_os);
	}

	void Move()	{
		for (int i=0; i<every; i++)	{
			model->Move(1);
		}
		SavePoint();
		Save();
		Monitor();
	}
};

int main(int argc, char* argv[])	{

	if (argc == 2)	{
		string name = argv[1];
		BGCChain* chain = new BGCChain(name);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
	else	{

		string datafile = "";
		string treefile = "";
		string contdatafile = "None";
		string calibfile = "None";
		string suffstatfile = "None";
		string rootfile = "None";
		double rootage = 0;
		double rootstdev = 0;

		int chronoprior = 0;
		double meanchi = 1e-3;
		double meanchi2 = 1e-3;

		string name = "";
		bool clampdiag = false;
		bool autoregressive = false;
		bool clamptree = false;
		bool meanexp = false;

		double priorsigma = -1;
		double fixalpha = 0;
		int df = 0;

		int conjpath = -1;
		int contdatatype = 0;

		int every = 1;
		int until = -1;

		int nrep = 0;
		bool normalise = false;

		int geneoffset = 1;
		int clampbgc = 0;
		int flexrho = 0;

		int nsplit = 1;

		int force = 0;

		double lambda = 2;
		double at2cg = -1;
		double at2gc = -1;
		double gc2cg = -1;

		double cg2at = -1;
		double cg2gc = -1;
		double cg2ta = -1;
		int discn = 0;
		int discgam = 0;
		int triplet = 0;
		int clampreg = 0;

		try	{

			if (argc == 1)	{
				throw(0);
			}

			int i = 1;
			while (i < argc)	{
				string s = argv[i];

				if (s == "-d")	{
					i++;
					datafile = argv[i];
				}
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				else if (s == "-c")	{
					i++;
					contdatafile = argv[i];
				}
				else if (s == "-suffstat")	{
					i++;
					suffstatfile = argv[i];
				}
				else if (s == "-cpg")	{
					triplet = 1;
				}
				else if (s == "-clampreg")	{
					clampreg = 1;
				}
				else if (s == "-root")	{
					i++;
					rootfile = argv[i];
				}
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-linear")	{
					contdatatype = 2;
				}
				else if (s == "-logit")	{
					contdatatype = 1;
				}
				else if (s == "-split")	{
					i++;
					s = argv[i];
					i++;
					nsplit = atoi(argv[i]);
					if (s == "length")	{
						nsplit = -nsplit;
					}
				}
				else if (s == "-cal")	{
					i++;
					calibfile = argv[i];
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootage = atof(argv[i]);
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootstdev = atof(argv[i]);
				}
				else if (s == "-bd")	{
					chronoprior = 1;
				}
				else if (s == "-cbd")	{
					chronoprior = 2;
				}
				else if (s == "-rbd")	{
					chronoprior = 3;
				}
				else if (s == "-bdhyperprior")	{
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -bd <p1> <p2>\n";
						cerr << '\n';
						exit(1);
					}
					meanchi = atof(argv[i]);
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -bd <p1> <p2>\n";
						cerr << '\n';
						exit(1);
					}
					meanchi2 = atof(argv[i]);
				}
				else if ((s == "-priorsigma") || (s == "-priornu"))	{
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -priorsigma <priorsigma>\n";
						cerr << '\n';
						exit(1);
					}
					priorsigma = atof(argv[i]);
				}
				else if (s == "-df")	{
					i++;
					df = atoi(argv[i]);
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-prior")	{
					conjpath = 2;
				}
				else if (s == "-norm")	{
					normalise = true;
				}
				else if (s == "-nonnorm")	{
					normalise = false;
				}
				else if (s == "-nrep")	{
					i++;
					nrep = atoi(argv[i]);
				}
				else if (s == "-fixbl")	{
					clamptree = true;
				}
				else if (s == "-meanexp")	{
					meanexp = true;
				}
				else if (s == "-arith")	{
					meanexp = false;
				}
				else if (s == "-geod")	{
					meanexp = true;
				}
				else if (s == "-diag")	{
					clampdiag = true;
				}
				else if (s == "-oup")	{
					autoregressive = true;
				}
				else if (s == "-freealpha")	{
					fixalpha = -1;
				}
				else if (s == "-freealphabeta")	{
					fixalpha = -2;
				}
				else if (s == "-hotspot")	{
					fixalpha = -3;
				}
				else if (s == "-clamphotspot")	{
					fixalpha = -4;
				}
				else if (s == "-fixalpha")	{
					i++;
					fixalpha = atof(argv[i]);
				}
				else if (s == "-geneoffset")	{
					geneoffset = 0;
				}
				else if (s == "-clampbgc")	{
					clampbgc = 1;
				}
				else if (s == "-flexrho")	{
					flexrho = 1;
				}
				else if (s == "-gammagene")	{
					flexrho = -1;
				}
				else if (s == "-wngene")	{
					flexrho = -2;
				}
				else if (s == "-lngene")	{
					flexrho = -3;
				}
				else if (s == "-mixedgene")	{
					flexrho = -4;
				}
				else if (s == "-lambda")	{
					i++;
					lambda = atof(argv[i]);
				}
				else if (s == "-discgam")	{
					i++;
					discgam = atoi(argv[i]);
				}
				else if (s == "-nonrev")	{
					discn = 10;
				}
				else if (s == "-discn")	{
					i++;
					discn = atoi(argv[i]);
				}
				else if (s == "-mut")	{
					if (!discn)	{
						i++;
						at2cg = atof(argv[i]);
						i++;
						at2gc = atof(argv[i]);
						i++;
						gc2cg = atof(argv[i]);
					}
					else	{
						i++;
						at2cg = atof(argv[i]);
						i++;
						at2gc = atof(argv[i]);
						i++;
						cg2at = atof(argv[i]);
						i++;
						cg2gc = atof(argv[i]);
						i++;
						cg2ta = atof(argv[i]);
					}
				}
				else if ( (s == "-x") || (s == "-extract") )	{
					i++;
					if (i == argc) throw(0);
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					until = atoi(argv[i]);
				}
				else	{
					if (i != (argc -1))	{
						throw(0);
					}
					name = argv[i];
				}
				i++;
			}
			if ((datafile == "") || (treefile == "") || (name == ""))	{
				throw(0);
			}
			if (contdatafile == "")	{
				contdatafile = "None";
			}
		}
		catch(...)	{
			cerr << "coevol -d <alignment> -t <tree> [-c <phenodata>] [-diag -dNdS -dSomega] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		int clampbgcoffset = 2 * clampbgc + geneoffset;

		BGCChain* chain = new BGCChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,name,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,clampbgcoffset,flexrho,lambda,at2cg,at2gc,gc2cg,cg2at,cg2gc,cg2ta,discgam,discn,triplet,clampreg,force);
		chain->SetEvery(every);
		chain->SetUntil(until);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
}

