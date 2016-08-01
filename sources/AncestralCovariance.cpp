
#include "Chain.h"
#include "AncestralCovarianceModel.h"
#include "StringStreamUtils.h"

class AncestralCovarianceChain : public Chain	{

	private:
	string modeltype;
	string treefile;
	string contdatafile;
	string datafile;
	string ancdatafile;

	bool clampdiag;
	int autoregressive;
	bool clamproot;
	bool meanexp;
	bool withdrift;
	int kalman;
	double priorsigma;

	int contdatatype;
	int ancdatatype;
	int df;

	public:

	string GetModelType() {return modeltype;}

	AncestralCovarianceChain(string intree, string incontdata, string inancdata, double inpriorsigma, int indf, string filename, bool inclampdiag, int inautoregressive, int incontdatatype, int inancdatatype, bool inclamproot, bool inmeanexp, bool inwithdrift, int inkalman, int force = 1)	{
		modeltype = "ANCESTRALCOV";
		treefile = intree;
		contdatafile = incontdata;
		ancdatafile = inancdata;
		name = filename;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		priorsigma = inpriorsigma;
		df = indf;
		contdatatype = incontdatatype;
		ancdatatype = inancdatatype;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		withdrift = inwithdrift;
		kalman = inkalman;
		New(force);
	}

	AncestralCovarianceChain(string filename)	{
		contdatatype = 0;
		ancdatatype = 0;
		name = filename;
		priorsigma = 1;
		df = 2;
		clamproot = false;
		meanexp = false;
		withdrift = false;
		kalman = 0;
		cerr << "open\n";
		Open();
	}

	void New(int force)	{
		if (modeltype == "ANCESTRALCOV")	{
			model = new AncestralCovarianceModel(treefile,contdatafile,ancdatafile,priorsigma,df,clampdiag,autoregressive,contdatatype,ancdatatype,clamproot,meanexp,withdrift,kalman,true);
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
		is >> treefile >> contdatafile >> ancdatafile ;
		is >> clampdiag >> autoregressive;
		is >> contdatatype;
		is >> ancdatatype;
		is >> clamproot >> meanexp;
		is >> df;
		is >> priorsigma;
		is >> withdrift;
		is >> kalman;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "ANCESTRALCOV")	{
			model = new AncestralCovarianceModel(treefile,contdatafile,ancdatafile,priorsigma,df,clampdiag,autoregressive,contdatatype,ancdatatype,clamproot,meanexp,withdrift,kalman,true);
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
		param_os << treefile << '\t' << contdatafile << '\t' << ancdatafile << '\n';
		param_os << clampdiag << '\t' << autoregressive << '\n';
		param_os << contdatatype << '\n';
		param_os << ancdatatype << '\n';
		param_os << clamproot << '\t' << meanexp << '\n';
		param_os << df << '\n';
		param_os << priorsigma << '\n';
		param_os << withdrift << '\n';
		param_os << kalman << '\n';
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
		AncestralCovarianceChain* chain = new AncestralCovarianceChain(name);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
	else	{

		string treefile = "";
		string contdatafile = "None";
		string ancdatafile = "None";

		string name = "";
		int autoregressive = 0;
		bool clampdiag = false;
		bool clamproot = false;
		bool meanexp = false;
		bool withdrift = false;
		bool kalman = 0;

		double priorsigma = -1;
		int df = 0;

		int contdatatype = 0;
		int ancdatatype = 0;

		int every = 1;
		int until = -1;

		int force = 0;

		bool conjugate = true;

		try	{

			if (argc == 1)	{
				throw(0);
			}

			int i = 1;
			while (i < argc)	{
				string s = argv[i];

				if (s == "-anc")	{
					i++;
					ancdatafile = argv[i];
				}
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				else if (s == "-c")	{
					i++;
					contdatafile = argv[i];
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
				else if (s == "-log")	{
					contdatatype = 0;
				}
				else if (s == "-anclinear")	{
					ancdatatype = 2;
				}
				else if (s == "-anclogit")	{
					ancdatatype = 1;
				}
				else if (s == "-anclog")	{
					ancdatatype = 0;
				}
				else if (s == "-ancfreq")	{
					ancdatatype = -1;
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
				else if (s == "-conj")	{
					conjugate = true;
				}
				else if (s == "-nonconj")	{
					conjugate = false;
				}
				else if (s == "-oup")	{
					autoregressive = 1;
				}
				else if (s == "-boup")	{
					autoregressive = 2;
					clampdiag = true;
					conjugate = false;
				}
				else if (s == "-bp")	{
					autoregressive = 1;
				}
				else if (s == "-clamproot")	{
					clamproot = true;
				}
				else if (s == "-drift")	{
					withdrift = true;
				}
				else if (s == "-k")	{
					kalman = 1;
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
			if ((contdatafile == "") || (treefile == "") || (name == ""))	{
				throw(0);
			}
		}
		catch(...)	{
			cerr << "ancov -t <tree> -c <contdata> [-anc <ancestral_data>] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		AncestralCovarianceChain* chain = new AncestralCovarianceChain(treefile,contdatafile,ancdatafile,priorsigma,df,name,clampdiag,autoregressive,contdatatype,ancdatatype,clamproot,meanexp,withdrift,kalman,force);
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

