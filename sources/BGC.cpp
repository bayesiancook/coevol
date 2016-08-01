
#include "Chain.h"
#include "BGCModel.h"
// #include "BGCModelwoLinReg.h"
#include "StringStreamUtils.h"


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

	double lambda;
	double at2cg;
	double at2gc;
	double gc2cg;

	public:

	string GetModelType() {return modeltype;}

	BGCChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inpriorsigma, int indf, string filename, bool inclampdiag, bool inautoregressive, int inconjpath, int incontdatatype, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int innsplit, string insuffstatfile, string inrootfile, double infixalpha, double inlambda, double inat2cg, double inat2gc, double ingc2cg, int force = 1)	{
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
		lambda = inlambda;
		at2cg = inat2cg;
		at2gc = inat2gc;
		gc2cg = ingc2cg;
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
		lambda = 2;
		at2cg = -1;
		at2gc = -1;
		gc2cg = -1;
		Open();
	}

	void New(int force)	{
		if (modeltype == "BGC")	{
			model = new BGCModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,lambda,at2cg,at2gc,gc2cg,true);
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
				is >> lambda;
				is >> at2cg;
				is >> at2gc;
				is >> gc2cg;
				is >> check;
				if (check)	{
					cerr << "error when reading model\n";
					exit(1);
				}
			}
		}

		is >> every >> until >> size;

		if (modeltype == "BGC")	{
			model = new BGCModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,lambda,at2cg,at2gc,gc2cg,true);
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
		param_os << lambda << '\n';
		param_os << at2cg << '\n';
		param_os << at2gc << '\n';
		param_os << gc2cg << '\n';
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

		int nsplit = 1;

		int force = 0;

		double lambda = 2;
		double at2cg = -1;
		double at2gc = -1;
		double gc2cg = -1;

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
				else if (s == "-fixalpha")	{
					i++;
					fixalpha = atof(argv[i]);
				}
				else if (s == "-lambda")	{
					i++;
					lambda = atof(argv[i]);
				}
				else if (s == "-mut")	{
					i++;
					at2cg = atof(argv[i]);
					i++;
					at2gc = atof(argv[i]);
					i++;
					gc2cg = atof(argv[i]);
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

		BGCChain* chain = new BGCChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,name,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,lambda,at2cg,at2gc,gc2cg,force);
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

