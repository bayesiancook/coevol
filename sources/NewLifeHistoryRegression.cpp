
#include "Chain.h"
#include "NewLifeHistoryRegressionModel.h"
#include "StringStreamUtils.h"

class LifeHistoryRegressionChain : public Chain	{

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

	double N0;
	double N1;
	double N2;
	double T0;

	bool clampdiag;
	bool autoregressive;
	bool clamproot;
	bool clamptree;
	int withdrift;
	int withpulse;
	double kt;

	GeneticCodeType type;

	double priorsigma;
	int conjpath;
	int contdatatype;

	int omegaratiotree;

	bool normalise;
	int nrep;

	int nsplit;

	int df;

	public:

	string GetModelType() {return modeltype;}

	LifeHistoryRegressionChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inN0, double inN1, double inN2, double inT0, double inpriorsigma, int indf, string filename, bool inclampdiag, bool inautoregressive, GeneticCodeType intype, int inconjpath, int incontdatatype, int inomegaratiotree, bool inclamproot, bool inclamptree, bool innormalise, int innrep, int innsplit, int inwithdrift, int inwithpulse, double inkt, string insuffstatfile, string inrootfile, int force = 1)	{
		modeltype = "LIFEHISTORYREGRESSION2";
		type = intype;
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
		N0 = inN0;
		N1 = inN1;
		N2 = inN2;
		T0 = inT0;
		name = filename;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		priorsigma = inpriorsigma;
		df = indf;
		conjpath = inconjpath;
		contdatatype = incontdatatype;
		omegaratiotree = inomegaratiotree;
		clamproot = inclamproot;
		clamptree = inclamptree;
		normalise = innormalise;
		nrep = innrep;
		nsplit = innsplit;
		withdrift = inwithdrift;
		withpulse = inwithpulse;
		kt = inkt;
		New(force);
	}

	LifeHistoryRegressionChain(string filename)	{
		conjpath = 1;
		contdatatype = 0;
		name = filename;
		priorsigma = 1;
		df = 2;
		omegaratiotree = 0;
		clamproot = false;
		clamptree = false;
		autoregressive = false;
		nsplit = 1;
		cerr << "open\n";
		withdrift = 0;
		withpulse = 0;
		kt = 0;
		suffstatfile = "None";
		rootfile = "None";
		Open();
	}

	void New(int force)	{
		if (modeltype == "LIFEHISTORYREGRESSION2")	{
			model = new LifeHistoryRegressionModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,N0,N1,N2,T0,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,clamptree,normalise,nrep,nsplit,withdrift,withpulse,kt,suffstatfile,rootfile,true,type);
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
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		if (chronoprior == 4)	{
			is >> N0 >> N1 >> N2 >> T0;
		}
		is >> clampdiag;
		is >> conjpath;
		is >> contdatatype;
		is >> omegaratiotree;
		is >> clamproot;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> nsplit;
		is >> withdrift;
		is >> withpulse;
		is >> kt;
		is >> clamptree;
		is >> suffstatfile;
		is >> autoregressive;
		is >> rootfile;

		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "LIFEHISTORYREGRESSION2")	{
			model = new LifeHistoryRegressionModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,N0,N1,N2,T0,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,clamptree,normalise,nrep,nsplit,withdrift,withpulse,kt,suffstatfile,rootfile,true,type);
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
		param_os << type << '\n';
		param_os << datafile << '\t' << treefile << '\t' << contdatafile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << chronoprior << '\t' << meanchi << '\t' << meanchi2 << '\n';
		if (chronoprior == 4)	{
			param_os << N0 << '\t' << N1 << '\t' << N2 << '\t' << T0 << '\n';
		}
		param_os << clampdiag << '\n';
		param_os << conjpath << '\n';
		param_os << contdatatype << '\n';
		param_os << omegaratiotree << '\n';
		param_os << clamproot << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << df << '\n';
		param_os << priorsigma << '\n';
		param_os << nsplit << '\n';
		param_os << withdrift << '\n';
		param_os << withpulse << '\n';
		param_os << kt << '\n';
		param_os << clamptree << '\n';
		param_os << suffstatfile << '\n';
		param_os << autoregressive << '\n';
		param_os << rootfile << '\n';
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
		LifeHistoryRegressionChain* chain = new LifeHistoryRegressionChain(name);
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

		double N0;
		double N1;
		double N2;
		double T0;

		string name = "";
		bool clampdiag = false;
		bool autoregressive = false;
		bool clamproot = false;
		bool clamptree = false;
		GeneticCodeType type = Universal;

		double priorsigma = -1;
		int df = 0;

		int conjpath = -1;
		int contdatatype = 0;

		int every = 1;
		int until = -1;

		int omegaratiotree = 0;

		int nrep = 0;
		bool normalise = true;

		int nsplit = 1;
		int withdrift = 0;
		int withpulse = 0;
		double kt = 0;

		int force = 0;
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
				else if (s == "-pulse")	{
					withpulse = 1;
				}
				else if (s == "-kt")	{
					i++;
					kt = atof(argv[i]);
				}
				else if (s == "-drift")	{
					withdrift = 1;
				}
				else if (s == "-expdrift")	{
					withdrift = 2;
				}
				else if (s == "-doubleexpdrift")	{
					withdrift = 3;
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
				else if (s == "-coal")	{
					chronoprior = 4;
					i++;
					meanchi = atof(argv[i]);
					i++;
					meanchi2 = atof(argv[i]);
					i++;
					N0 = atof(argv[i]);
					i++;
					N1 = atof(argv[i]);
					i++;
					N2 = atof(argv[i]);
					i++;
					T0 = atof(argv[i]);
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
				else if ((s == "-dSdN") || (s == "-dsdn"))	{
					omegaratiotree = 2;
				}
				else if ((s == "-dSomega") || (s == "-dsomega") || (s == "-dsom"))	{
					omegaratiotree = 1;
				}
				else if (s == "-clamproot")	{
					clamproot = true;
				}
				else if (s == "-fixbl")	{
					clamptree = true;
				}
				else if (s == "-diag")	{
					clampdiag = true;
				}
				else if (s == "-oup")	{
					autoregressive = true;
				}
				else if (s == "-mtmam")	{
					type = MtMam;
				}
				else if (s == "-mtinv")	{
					type = MtInv;
				}
				else if (s == "-uni")	{
					type = Universal;
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

		LifeHistoryRegressionChain* chain = new LifeHistoryRegressionChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,N0,N1,N2,T0,priorsigma,df,name,clampdiag,autoregressive,type,conjpath,contdatatype,omegaratiotree,clamproot,clamptree,normalise,nrep,nsplit,withdrift,withpulse,kt,suffstatfile,rootfile,force);
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

