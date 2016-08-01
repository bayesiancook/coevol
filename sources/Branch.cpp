
#include "Chain.h"
#include "ConjugateBranchModel.h"
// #include "BranchModel.h"
#include "StringStreamUtils.h"

class BranchChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string calibfile;

	string suffstatfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool clamptree;
	GeneticCodeType type;

	int gc;
	int codonmodel;

	int conjpath;
	int fullconj;
	double mappingfreq;
	int contdatatype;

	bool normalise;
	int nrep;

	public:

	string GetModelType() {return modeltype;}

	BranchModel* GetBranchModel()	{
		BranchModel* tmp = dynamic_cast<BranchModel*>(model);
		if (! tmp)	{
			cerr << "error in BranchChain: null model pointer\n";
			exit(1);
		}
		return tmp;
	}

	BranchChain(string indata, string intree, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, int ingc, int incodonmodel, string filename, GeneticCodeType intype, int inconjpath, int infullconj, double inmappingfreq, bool inclamptree, bool innormalise, int innrep, string insuffstatfile, int force = 1)	{
		modeltype = "BRANCHWISE";

		type = intype;
		datafile = indata;
		treefile = intree;
		calibfile = incalibfile;
		suffstatfile = insuffstatfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		name = filename;
		gc = ingc;
		codonmodel = incodonmodel;
		conjpath = inconjpath;
		fullconj = infullconj;
		mappingfreq = inmappingfreq;
		clamptree = inclamptree;
		normalise = innormalise;
		nrep = innrep;
		New(force);
	}

	BranchChain(string filename)	{
		conjpath = true;
		fullconj = true;
		mappingfreq = 1;
		name = filename;
		codonmodel = 0;
		gc = 0;
		clamptree = false;
		suffstatfile = "None";
		Open();
	}

	void New(int force)	{
		if (modeltype == "BRANCHWISE")	{
			model = new BranchModel(datafile,treefile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,gc,codonmodel,conjpath,fullconj,mappingfreq,clamptree,normalise,nrep,suffstatfile,true,type);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		cerr << "reset\n";
		Reset(force);
		ofstream os((name + ".suffstat").c_str());
		// ofstream osat((name + ".nansgc").c_str());
		// ofstream osgc((name + ".nansat").c_str());
		// model->Update();
		cerr << "log prob\n";
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
		is >> datafile >> treefile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> gc >> codonmodel;
		is >> conjpath;
		is >> normalise >> nrep;
		is >> clamptree;
		int check;
		is >> check;
		if (check)	{
			is >> suffstatfile;
			is >> check;
			if (check)	{
				is >> mappingfreq;
				is >> check;
				if (check)	{
					is >> fullconj;
					is >> check;
					if (check)	{
						cerr << "error when reading model\n";
						exit(1);
					}
				}
			}
		}

		is >> every >> until >> size;

		if (modeltype == "BRANCHWISE")	{
			model = new BranchModel(datafile,treefile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,gc,codonmodel,conjpath,fullconj,mappingfreq,clamptree,normalise,nrep,suffstatfile,true,type);
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
		param_os << datafile << '\t' << treefile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << chronoprior << '\t' << meanchi << '\t' << meanchi2 << '\n';
		param_os << gc << '\t' << codonmodel << '\n';
		param_os << conjpath << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << clamptree << '\n';
		param_os << 1 << '\n';
		param_os << suffstatfile << '\n';
		param_os << 1 << '\n';
		param_os << mappingfreq << '\n';
		param_os << 1 << '\n';
		param_os << fullconj << '\n';
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
		ofstream os((name + ".suffstat").c_str(), ios_base::app);
		GetBranchModel()->PrintSuffStat(os);
	}
};

int main(int argc, char* argv[])	{

	if (argc == 2)	{
		string name = argv[1];
		BranchChain* chain = new BranchChain(name);
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
		string calibfile = "None";
		string suffstatfile = "None";
		double rootage = 0;
		double rootstdev = 0;

		int chronoprior = 0;
		double meanchi = 1e-3;
		double meanchi2 = 1e-3;

		string name = "";
		bool clamptree = false;
		GeneticCodeType type = Universal;

		int gc = 0;
		int codonmodel = 0;

		int conjpath = -1;
		double mappingfreq = -1;

		int fullconj = 1;

		int every = 1;
		int until = -1;

		int nrep = 0;
		bool normalise = true;

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
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-suffstat")	{
					i++;
					suffstatfile = argv[i];
				}
				else if (s == "-mapfreq")	{
					i++;
					mappingfreq = atof(argv[i]);
				}
				else if (s == "-uncons")	{
					calibfile = "Unconstrained";
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
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-conjom")	{
					fullconj = 1;
				}
				else if (s == "-nonconjom")	{
					fullconj = 0;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
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
				else if ((s == "-dSomega") || (s == "-dsomega") || (s == "-dsom"))	{
					codonmodel = 1;
				}
				else if (s == "-dsom3")	{
					codonmodel = 3;
				}
				else if (s == "-gc")	{
					gc = 1;
				}
				else if (s == "-gc3")	{
					gc = 3;
				}
				else if (s == "-fixbl")	{
					clamptree = true;
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
		}
		catch(...)	{
			cerr << "branchmodel -d <alignment> -t <tree> [-dsom -gc -gc3] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		BranchChain* chain = new BranchChain(datafile,treefile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,gc,codonmodel,name,type,conjpath,fullconj,mappingfreq,clamptree,normalise,nrep,suffstatfile,force);
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

