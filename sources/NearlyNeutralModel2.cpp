
#include "Chain.h"
#include "NearlyNeutralModel2.h"
#include "StringStreamUtils.h"

class BranchOmegaMultivariateChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	double rootage;
	double rootstdev;

	bool clamptree;
	bool meanexp;
	bool sameseq;
	bool noadapt;

	GeneticCodeType type;

	double priorsigma;

	int contdatatype;
	int nrep;

	int df;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaMultivariateChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, double inpriorsigma, int indf, GeneticCodeType intype, int incontdatatype, bool insameseq, bool innoadapt, bool inclamptree, bool inmeanexp, int innrep, string filename, int force = 1)	{
		modeltype = "CONJUGATEBRANCHOMEGAMULTIVARIATE";

		type = intype;
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		calibfile = incalibfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		priorsigma = inpriorsigma;
		df = indf;
		contdatatype = incontdatatype;
		sameseq = insameseq;
		noadapt = innoadapt;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		nrep = innrep;
		name = filename;
		New(force);
	}

	BranchOmegaMultivariateChain(string filename)	{
		contdatatype = 0;
		name = filename;
		priorsigma = 1;
		df = 2;
		clamptree = false;
		meanexp = false;
		Open();
	}

	void New(int force)	{
		if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,priorsigma,df,contdatatype,sameseq,noadapt,clamptree,meanexp,nrep,true,type);
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
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> contdatatype;
		is >> sameseq;
		is >> noadapt;
		is >> meanexp;
		is >> nrep;
		is >> df;
		is >> priorsigma;
		is >> clamptree;

		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,priorsigma,df,contdatatype,sameseq,noadapt,clamptree,meanexp,nrep,true,type);
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
		param_os << contdatatype << '\n';
		param_os << sameseq << '\n';
		param_os << noadapt << '\n';
		param_os << meanexp << '\n';
		param_os << nrep << '\n';
		param_os << df << '\n';
		param_os << priorsigma << '\n';
		param_os << clamptree << '\n';
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
		BranchOmegaMultivariateChain* chain = new BranchOmegaMultivariateChain(name);
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

		string name = "";
		bool sameseq = false;
		bool noadapt = false;
		bool clamptree = false;
		bool meanexp = false;
		GeneticCodeType type = Universal;

		double priorsigma = -1;
		int df = 0;

		int contdatatype = 0;

		int every = 1;
		int until = -1;

		int nrep = 0;

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
				else if (s == "-calspe")	{
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootage = atof(argv[i]);
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
				else if (s == "-nrep")	{
					i++;
					nrep = atoi(argv[i]);
				}
				else if (s == "-fixbl")	{
					clamptree = true;
				}
				else if (s == "-rootage")	{
					i++;
					rootage = atof(argv[i]);
				}
				else if (s == "-sameseq")	{
					sameseq = true;
				}
				else if (s == "-noadapt")	{
					noadapt = true;
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
			if (contdatafile == "")	{
				contdatafile = "None";
			}
		}
		catch(...)	{
			cerr << "coevol -d <alignment> -t <tree> [-c <phenodata>] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		BranchOmegaMultivariateChain* chain = new BranchOmegaMultivariateChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,priorsigma,df,type,contdatatype,sameseq,noadapt,clamptree,meanexp,nrep,name,force);
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

