
#include "Chain.h"
#include "SiteOmegaModel.h"
#include "StringStreamUtils.h"

class SiteOmegaChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	GeneticCodeType type;

	bool clamptree;

	int conjpath;
	int conjom;
	double mappingfreq;
	bool normalise;
	int nrep;

	public:

	string GetModelType() {return modeltype;}

	SiteOmegaModel* GetSiteOmegaModel()	{
		SiteOmegaModel* tmp = dynamic_cast<SiteOmegaModel*>(model);
		if (! tmp)	{
			cerr << "error in SiteOmegaChain: null model pointer\n";
			exit(1);
		}
		return tmp;
	}

	SiteOmegaChain(string indata, string intree, string filename, GeneticCodeType intype, int inconjpath, int inconjom, double inmappingfreq, bool inclamptree, bool innormalise, int innrep, int force = 1)	{

		modeltype = "SITEWISE";

		type = intype;
		datafile = indata;
		treefile = intree;
		name = filename;
		conjpath = inconjpath;
		conjom = inconjom;
		mappingfreq = inmappingfreq;
		clamptree = inclamptree;
		normalise = innormalise;
		nrep = innrep;
		New(force);
	}

	SiteOmegaChain(string filename)	{

		conjpath = true;
		conjom = true;
		mappingfreq = 1;
		name = filename;
		clamptree = false;
		Open();
	}

	void New(int force)	{
		if (modeltype == "SITEWISE")	{
			model = new SiteOmegaModel(datafile,treefile,conjpath,conjom,mappingfreq,clamptree,normalise,nrep,true,type);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		cerr << "reset\n";
		Reset(force);
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
		is >> conjpath;
		is >> conjom;
		is >> normalise >> nrep;
		is >> clamptree;
		is >> mappingfreq;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "SITEWISE")	{
			model = new SiteOmegaModel(datafile,treefile,conjpath,conjom,mappingfreq,clamptree,normalise,nrep,true,type);
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
		param_os << conjpath << '\n';
		param_os << conjom << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << clamptree << '\n';
		param_os << mappingfreq << '\n';
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
		SiteOmegaChain* chain = new SiteOmegaChain(name);
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

		string name = "";
		bool clamptree = false;
		GeneticCodeType type = Universal;

		int conjpath = -1;
		int conjom = 1;
		double mappingfreq = -1;

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
				else if (s == "-mapfreq")	{
					i++;
					mappingfreq = atof(argv[i]);
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-conjom")	{
					conjom = 1;
				}
				else if (s == "-nonconjom")	{
					conjom = 0;
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
			cerr << "siteom -d <alignment> -t <tree> [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		SiteOmegaChain* chain = new SiteOmegaChain(datafile,treefile,name,type,conjpath,conjom,mappingfreq,clamptree,normalise,nrep,force);
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

