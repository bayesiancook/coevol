
#include "Chain.h"
#include "CompoModel.h"
#include "StringStreamUtils.h"

class CompoChain : public Chain	{

	private:

	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double softa;

	double meanchi;
	double meanchi2;

	int df;

	int clampdiag;
	int contdatatype;
	string rrtype;
	int clamptree;
	int meanexp;

	int normalise;
	int nrep;
	string mix;
	string rootfile;
	int separatesyn;

	public:

	string GetModelType() {return modeltype;}

	CompoChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double insofta, double inmeanchi, double inmeanchi2, int indf, int inclampdiag, int incontdatatype, string inrrtype, int inclamptree, int inmeanexp, int innormalise, int innrep, string inmix, string inrootfile, int inseparatesyn, string filename, int force = 1)	{
		modeltype = "COMPO";

		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		calibfile = incalibfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		softa = insofta;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		df = indf;
		clampdiag = inclampdiag;
		contdatatype = incontdatatype;
		rrtype = inrrtype;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		normalise = innormalise;
		nrep = innrep;
		mix = inmix;
		rootfile = inrootfile;
		separatesyn =  inseparatesyn;

		name = filename;
		New(force);
	}

	CompoChain(string filename)	{
		contdatatype = 0;
		rrtype = "lg";
		name = filename;
		df = 2;
		clamptree = false;
		meanexp = false;
		Open();
	}

	void New(int force)	{
		if (modeltype == "COMPO")	{
			model = new CompoModel(datafile, treefile, contdatafile, calibfile, rootage, rootstdev, chronoprior, softa,  meanchi, meanchi2, df, clampdiag, contdatatype, rrtype, clamptree, meanexp, normalise, nrep, mix, rootfile, separatesyn, true);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		Reset(force);
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
		is >> chronoprior >> softa >> meanchi >> meanchi2;
		is >> df;
		is >> clampdiag;
		is >> contdatatype;
		is >> rrtype;
		is >> clamptree;
		is >> meanexp;
		is >> normalise;
		is >> nrep;
		is >> mix;
		is >> rootfile;
		is >> separatesyn;

		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "COMPO")	{
			model = new CompoModel(datafile, treefile, contdatafile, calibfile, rootage, rootstdev, chronoprior, softa,  meanchi, meanchi2, df, clampdiag, contdatatype, rrtype, clamptree, meanexp, normalise, nrep, mix, rootfile, separatesyn, true);
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
		param_os << chronoprior << '\t' << softa << '\t' << meanchi << '\t' << meanchi2 << '\n';
		param_os << df << '\n';
		param_os << clampdiag << '\n';
		param_os << contdatatype << '\n';
		param_os << rrtype << '\n';
		param_os << clamptree << '\n';
		param_os << meanexp << '\n';
		param_os << normalise << '\n';
		param_os << nrep << '\n';
		param_os << mix << '\n';
		param_os << rootfile << '\n';
		param_os << separatesyn << '\n';
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
		CompoChain* chain = new CompoChain(name);
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
		double softa = 0.025;
		double meanchi = 1e-3;
		double meanchi2 = 1e-3;

		string name = "";
		int clampdiag = 0;
		int clamptree = 0;
		int meanexp = 0;
		int normalise = 0;

		int separatesyn = 0;
		int df = 0;

		int contdatatype = 0;
		string rrtype = "lg";

		int nrep = 0;

		string mix = "None";
		string rootfile = "None";

		int every = 1;
		int until = -1;

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
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-linear")	{
					contdatatype = 2;
				}
				else if (s == "-logit")	{
					contdatatype = 1;
				}
				else if (s == "-rr")	{
					i++;
					rrtype = argv[i];
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
				else if (s == "-unif")	{
					chronoprior = 0;
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
				else if (s == "-sbd")	{
					chronoprior = 4;
					i++;
					softa = atof(argv[i]);
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
				else if (s == "-uncons")	{
					calibfile = "Unconstrained";
				}
				else if (s == "-df")	{
					i++;
					df = atoi(argv[i]);
				}
				else if (s == "-norm")	{
					normalise = 0;
				}
				else if (s == "-nonnorm")	{
					normalise = 0;
				}
				else if (s == "-nrep")	{
					i++;
					nrep = atoi(argv[i]);
				}
				else if ((s == "-fixtimes") || (s == "-fixbl"))	{
					clamptree = 1;
				}
				else if (s == "-meanexp")	{
					meanexp = 1;
				}
				else if (s == "-arith")	{
					meanexp = 0;
				}
				else if (s == "-geod")	{
					meanexp = 1;
				}
				else if (s == "-diag")	{
					clampdiag = 1;
				}
				else if (s == "-sepsyn")	{
					separatesyn = 1;
				}
				else if (s == "-mix")	{
					i++;
					mix = argv[i];
				}
				else if (s == "-root")	{
					i++;
					rootfile = argv[i];
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
			if (contdatafile == "")	{
				contdatafile = "None";
			}
		}
		catch(...)	{
			cerr << "coevol -d <alignment> -t <tree> [-c <phenodata>] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		CompoChain* chain = new CompoChain(datafile, treefile, contdatafile, calibfile, rootage, rootstdev, chronoprior, softa, meanchi, meanchi2, df, clampdiag, contdatatype, rrtype, clamptree, meanexp, normalise, nrep, mix, rootfile, separatesyn, name, force);
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

