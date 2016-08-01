
#include "Chain.h"
#include "ConjugateRateMultivariateModel.h"

class RateMultivariateChain : public Chain	{

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

	double priorsigma;
	bool gc;

	int conjpath;
	int contdatatype;

	public:

	string GetModelType() {return modeltype;}

	RateMultivariateChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inpriorsigma, bool ingc, string filename, bool inclampdiag, bool inautoregressive, bool conjugate, int inconjpath, int incontdatatype, bool inclamproot, bool inmeanexp, int force = 1)	{
		if (conjugate)	{
			modeltype = "CONJUGATERATEMULTIVARIATE";
		}
		else	{
			modeltype = "RATEMULTIVARIATE";
		}
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
		gc = ingc;
		conjpath = inconjpath;
		contdatatype = incontdatatype;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		New(force);
	}

	RateMultivariateChain(string filename)	{
		conjpath = true;
		contdatatype = 0;
		name = filename;
		priorsigma = 1;
		clamproot = false;
		meanexp = false;
		cerr << "open\n";
		Open();
	}

	void New(int force)	{
		if (modeltype == "RATEMULTIVARIATE")	{
			model = new RateMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,clampdiag,autoregressive,conjpath,contdatatype,clamproot,meanexp,true);
		}
		else if (modeltype == "CONJUGATERATEMULTIVARIATE")	{
			model = new ConjugateRateMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,autoregressive,conjpath,contdatatype,clamproot,meanexp,true);
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
		is >> clampdiag >> autoregressive >> gc;
		is >> conjpath >> contdatatype >> clamproot >> meanexp;
		is >> every >> until >> size;

		if (modeltype == "RATEMULTIVARIATE")	{
			model = new RateMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,clampdiag,autoregressive,conjpath,contdatatype,clamproot,meanexp,true);
		}
		else if (modeltype == "CONJUGATERATEMULTIVARIATE")	{
			model = new ConjugateRateMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,autoregressive,conjpath,contdatatype,clamproot,meanexp,true);
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
		param_os << clampdiag << '\t' << autoregressive << '\t' << gc << '\n';
		param_os << conjpath << '\t' << contdatatype << '\t' << clamproot << '\t' << meanexp << '\n';
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
		RateMultivariateChain* chain = new RateMultivariateChain(name);
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

		double priorsigma = 1;
		bool gc = false;

		int conjpath = -1;
		int contdatatype = 0;

		int every = 1;
		int until = -1;

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
				else if (s == "-gc")	{
					gc = true;
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
				else if (s == "-geod")	{
					meanexp = true;
				}
				else if (s == "-diag")	{
					clampdiag = true;
					conjugate = false;
				}
				else if (s == "-diag")	{
					clampdiag = true;
					conjugate = false;
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
			cerr << "ratecorrel -d <alignment> -t <tree> [-c <phenodata>] [-conj -nonconj -bp -oup -diag] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		RateMultivariateChain* chain = new RateMultivariateChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,name,clampdiag,autoregressive,conjugate,conjpath,contdatatype,clamproot,meanexp);
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

