
#include "Chain.h"
#include "ConjugateGammaTreeOmegaMultivariateModel.h"

class GammaTreeOmegaMultivariateChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	bool clampdiag;
	bool autoregressive;
	GeneticCodeType type;

	double priorsigma;
	bool gc;

	public:

	string GetModelType() {return modeltype;}

	GammaTreeOmegaMultivariateChain(string indata, string intree, string incontdata, double inpriorsigma, bool ingc, string filename, bool inclampdiag, bool inautoregressive, GeneticCodeType intype, bool conjugate, int force = 1)	{
		if (conjugate)	{
			modeltype = "CONJUGATEGAMMATREEOMEGAMULTIVARIATE";
		}
		else	{
			modeltype = "GAMMATREEOMEGAMULTIVARIATE";
		}
		type = intype;
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		name = filename;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		priorsigma = inpriorsigma;
		gc = ingc;
		New(force);
	}

	GammaTreeOmegaMultivariateChain(string filename)	{
		name = filename;
		priorsigma = 1;
		cerr << "open\n";
		Open();
	}

	void New(int force)	{
		if (modeltype == "GAMMATREEOMEGAMULTIVARIATE")	{
			model = new GammaTreeOmegaMultivariateModel(datafile,treefile,contdatafile,priorsigma,gc,clampdiag,autoregressive,true,type);
		}
		else if (modeltype == "CONJUGATEGAMMATREEOMEGAMULTIVARIATE")	{
			model = new ConjugateGammaTreeOmegaMultivariateModel(datafile,treefile,contdatafile,priorsigma,gc,autoregressive,true,type);
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
		is >> clampdiag >> autoregressive >> gc;

		bool check;
		is >> check;
		if (check)	{
			cerr << "error when reading file\n";
			exit(1);
		}
		is >> every >> until >> size;

		if (modeltype == "GAMMATREEOMEGAMULTIVARIATE")	{
			model = new GammaTreeOmegaMultivariateModel(datafile,treefile,contdatafile,priorsigma,gc,clampdiag,autoregressive,true,type);
		}
		else if (modeltype == "CONJUGATEGAMMATREEOMEGAMULTIVARIATE")	{
			model = new ConjugateGammaTreeOmegaMultivariateModel(datafile,treefile,contdatafile,priorsigma,gc,autoregressive,true,type);
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
		param_os << clampdiag << '\t' << autoregressive << '\t' << gc << '\n';
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
		GammaTreeOmegaMultivariateChain* chain = new GammaTreeOmegaMultivariateChain(name);
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
		string name = "";
		bool autoregressive = false;
		bool clampdiag = false;
		bool conjugate = true;
		GeneticCodeType type = Universal;

		double priorsigma = 1;
		bool gc = false;

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
				else if (s == "-priornu")	{
					i++;
					priorsigma = atof(argv[i]);
				}
				else if (s == "-conj")	{
					conjugate = true;
				}
				else if (s == "-gc")	{
					gc = true;
				}
				else if (s == "-nonconj")	{
					conjugate = false;
				}
				else if (s == "-oup")	{
					autoregressive = true;
				}
				else if (s == "-bp")	{
					autoregressive = false;
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
			cerr << "mulomega -d <alignment> -t <tree> [-c <phenodata>] [-conj -nonconj -bp -oup -diag] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		GammaTreeOmegaMultivariateChain* chain = new GammaTreeOmegaMultivariateChain(datafile,treefile,contdatafile,priorsigma,gc,name,clampdiag,autoregressive,type,conjugate);
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

