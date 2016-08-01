
#include "Chain.h"
#include "ConjugateTamuraModel.h"

class TamuraChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool clampdiag;
	bool autoregressive;
	bool clamproot;
	bool meanexp;
	GeneticCodeType type;

	double priorsigma;

	int conjpath;
	int contdatatype;

	int omegaratiotree;

	bool normalise;
	int nrep;

	public:

	string GetModelType() {return modeltype;}

	TamuraChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inpriorsigma, string filename, bool inclampdiag, bool inautoregressive, GeneticCodeType intype, bool conjugate, int inconjpath, int incontdatatype, int inomegaratiotree, bool inclamproot, bool inmeanexp, bool innormalise, int innrep, int force = 1)	{
		if (conjugate)	{
			modeltype = "CONJUGATETAMURA";
		}
		else	{
			modeltype = "TAMURA";
		}
		type = intype;
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		calibfile = incalibfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		name = filename;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		priorsigma = inpriorsigma;
		conjpath = inconjpath;
		contdatatype = incontdatatype;
		omegaratiotree = inomegaratiotree;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		normalise = innormalise;
		nrep = innrep;
		New(force);
	}

	TamuraChain(string filename)	{
		conjpath = true;
		contdatatype = 0;
		name = filename;
		priorsigma = 1;
		omegaratiotree = 0;
		clamproot = false;
		meanexp = false;
		cerr << "open\n";
		Open();
	}

	void New(int force)	{
		if (modeltype == "TAMURA")	{
			model = new TamuraModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,clampdiag,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep,true,type);
		}
		else if (modeltype == "CONJUGATETAMURA")	{
			model = new ConjugateTamuraModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep,true,type);
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
		is >> clampdiag >> autoregressive;
		is >> conjpath;
		is >> contdatatype;
		is >> omegaratiotree;
		is >> clamproot >> meanexp;
		is >> normalise >> nrep;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "TAMURA")	{
			model = new TamuraModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,clampdiag,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep,true,type);
		}
		else if (modeltype == "CONJUGATETAMURA")	{
			model = new ConjugateTamuraModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep,true,type);
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
		param_os << clampdiag << '\t' << autoregressive << '\n';
		param_os << conjpath << '\n';
		param_os << contdatatype << '\n';
		param_os << omegaratiotree << '\n';
		param_os << clamproot << '\t' << meanexp << '\n';
		param_os << normalise << '\t' << nrep << '\n';
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
		TamuraChain* chain = new TamuraChain(name);
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
		double rootage = 0;
		double rootstdev = 0;

		int chronoprior = 0;
		double meanchi = 1e-3;
		double meanchi2 = 1e-3;

		string name = "";
		bool autoregressive = false;
		bool clampdiag = false;
		bool clamproot = false;
		bool meanexp = false;
		bool conjugate = true;
		GeneticCodeType type = Universal;

		double priorsigma = 1;

		int conjpath = -1;
		int contdatatype = 0;

		int every = 1;
		int until = -1;

		int omegaratiotree = 0;

		int nrep = 0;
		bool normalise = false;

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
				else if (s == "-linear")	{
					contdatatype = 2;
				}
				else if (s == "-logit")	{
					contdatatype = 1;
				}
				else if (s == "-cal")	{
					i++;
					calibfile = argv[i];
					i++;
					rootage = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
				}
				else if (s == "-bd")	{
					chronoprior = 1;
				}
				else if (s == "-bdhyperprior")	{
					i++;
					meanchi = atof(argv[i]);
					i++;
					meanchi2 = atof(argv[i]);
				}
				else if (s == "-priornu")	{
					i++;
					priorsigma = atof(argv[i]);
				}
				else if (s == "-conj")	{
					conjugate = true;
				}
				else if (s == "-nonconj")	{
					conjugate = false;
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-norm")	{
					normalise = true;
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
				else if (s == "-oup")	{
					autoregressive = true;
				}
				else if (s == "-bp")	{
					autoregressive = false;
				}
				else if (s == "-clamproot")	{
					clamproot = true;
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
					conjugate = false;
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
			if (clampdiag && conjugate)	{
				cerr << "error : conjugate sampling and diagonal model are not compatible\n";
				exit(1);
			}
		}
		catch(...)	{
			cerr << "coevol -d <alignment> -t <tree> [-c <phenodata>] [-diag -dNdS -dSomega] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		TamuraChain* chain = new TamuraChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,name,clampdiag,autoregressive,type,conjugate,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep);
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

