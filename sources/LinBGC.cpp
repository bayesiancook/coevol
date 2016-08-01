
#include "Chain.h"
#include "LinBGCModel.h"
#include "StringStreamUtils.h"

// int BGCMutSelSubMatrix::bgccount = 0;

class BGCChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string suffstatfile;
	string mix;
	int withbranch;
	int withgene;
	int withbranchgene;
	double priorsigma;
	int conjpath;
	bool normalise;
	int nrep;
	int df;
	int discn;
	GeneticCodeType type;

	LinBGCModel* linbgcmodel;

	public:

	string GetModelType() {return modeltype;}

	BGCChain(string indata, string intree, double inpriorsigma, int indf, string filename, int inconjpath, bool innormalise, int innrep, string insuffstatfile, string inmix, int inwithbranch, int inwithgene, int inwithbranchgene, int indiscn, GeneticCodeType intype, int force = 1)	{
		modeltype = "LINBGC";
		datafile = indata;
		treefile = intree;
		suffstatfile = insuffstatfile;
		mix = inmix;
		withbranch = inwithbranch;
		withgene = inwithgene;
		withbranchgene = inwithbranchgene;
		name = filename;
		priorsigma = inpriorsigma;
		df = indf;
		conjpath = inconjpath;
		normalise = innormalise;
		nrep = innrep;
		discn = indiscn;
		type = intype;
		New(force);
	}

	BGCChain(string filename)	{
		name = filename;
		/*
		conjpath = 1;
		priorsigma = 1;
		df = 2;
		suffstatfile = "None";
		mix = "branch";
		withbranch = 1;
		withgene = 1;
		withbranchgene = 1;
		discn = 0;
		*/
		Open();
	}

	void New(int force)	{
		if (modeltype == "LINBGC")	{
			linbgcmodel = new LinBGCModel(datafile,treefile,priorsigma,df,conjpath,normalise,nrep,suffstatfile,mix,withbranch,withgene,withbranchgene,discn,type,true);
			model = linbgcmodel;
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		Reset(force);
		// model->Update();
		ofstream trace_os((name + ".sigma").c_str());
		linbgcmodel->TraceSigmaHeader(trace_os);
		cerr << "ln prob = " << model->GetLogProb() << "\n";
	}

	void Monitor()	{
		ofstream trace_os((name + ".sigma").c_str(), ios_base::app);
		linbgcmodel->TraceSigma(trace_os);
		Chain::Monitor();
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> datafile >> treefile;
		is >> conjpath;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> suffstatfile;
		is >> mix;
		is >> withbranch;
		is >> withgene;
		is >> withbranchgene;
		is >> type;
		is >> discn;

		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "LINBGC")	{
			linbgcmodel = new LinBGCModel(datafile,treefile,priorsigma,df,conjpath,normalise,nrep,suffstatfile,mix,withbranch,withgene,withbranchgene,discn,type,true);
			model = linbgcmodel;
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
		param_os << datafile << '\t' << treefile << '\n';
		param_os << conjpath << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << df << '\n';
		param_os << priorsigma << '\n';
		param_os << suffstatfile << '\n';
		param_os << mix << '\n';
		param_os << withbranch << '\n';
		param_os << withgene << '\n';
		param_os << withbranchgene << '\n';
		param_os << type << '\n';
		param_os << discn << '\n';
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
		cerr << name << '\n';
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
		string suffstatfile = "None";
		string mix = "branch";

		int withbranch = 1;
		int withgene = 1;
		int withbranchgene = 1;

		string name = "";

		double priorsigma = -1;
		int df = 0;

		int conjpath = -1;

		int every = 1;
		int until = -1;

		int nrep = 0;
		bool normalise = false;

		int force = 0;

		int discn = 10;
		GeneticCodeType type = Universal;

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
				else if (s == "-suffstat")	{
					i++;
					suffstatfile = argv[i];
				}
				else if (s == "-mix")	{
					i++;
					mix = argv[i];
				}
				else if (s == "-fixbranch")	{
					withbranch = 0;
					withbranchgene = 0;
				}
				else if (s == "-fixgene")	{
					withgene = 0;
					withbranchgene = 0;
				}
				else if (s == "-fixbranchgene")	{
					withbranchgene = 0;
				}
				else if (s == "-f")	{
					force = 1;
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
				else if (s == "-discn")	{
					i++;
					discn = atoi(argv[i]);
				}
				else if (s == "-mtvert")	{
					type = MtMam;
				}
				else if (s == "-mtmam")	{
					type = MtMam;
				}
				else if (s == "-mtinv")	{
					type = MtInv;
				}
				else if (s == "-univ")	{
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
			cerr << "error in command\n";
			cerr << '\n';
			exit(1);
		}

		BGCChain* chain = new BGCChain(datafile,treefile,priorsigma,df,name,conjpath,normalise,nrep,suffstatfile,mix,withbranch,withgene,withbranchgene,discn,type,force);
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

