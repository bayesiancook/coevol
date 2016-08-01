
#include "Chain.h"
#include "TimeLineClockModel.h"
#include "StringStreamUtils.h"

class TimeLineClockChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string ancdatafile;
	string calibfile;
	string timelinefile;
	int flextimeline;

	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool clampdiag;
	bool clamproot;
	bool clamptree;
	bool meanexp;

	double priorsigma;

	int conjpath;
	int contdatatype;
	int ancdatatype;

	bool normalise;
	int nrep;

	int nsplit;

	int df;

	string mix;

	public:

	string GetModelType() {return modeltype;}

	TimeLineClockChain(string indata, string intree, string incontdata, string inancdata, string intimelinefile, int inflextimeline, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inpriorsigma, int indf, string filename, bool inclampdiag, int inconjpath, int incontdatatype, int inancdatatype, bool inclamproot, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, string inmix, int innsplit, string insuffstatfile, string inrootfile, int force = 1)	{
		modeltype = "TIMELINE";
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		ancdatafile = inancdata;
		calibfile = incalibfile;
		timelinefile = intimelinefile;
		flextimeline = inflextimeline;
		suffstatfile = insuffstatfile;
		rootfile = inrootfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		name = filename;
		clampdiag = inclampdiag;
		priorsigma = inpriorsigma;
		df = indf;
		conjpath = inconjpath;
		contdatatype = incontdatatype;
		ancdatatype = inancdatatype;
		clamproot = inclamproot;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		normalise = innormalise;
		nrep = innrep;
		mix = inmix;
		nsplit = innsplit;
		New(force);
	}

	TimeLineClockChain(string filename)	{
		cerr << "open\n";
		Open();
	}

	void New(int force)	{
		if (modeltype == "TIMELINE")	{
			model = new TimeLineClockModel(datafile,treefile,contdatafile,ancdatafile,timelinefile,flextimeline,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,conjpath,contdatatype,ancdatatype,clamproot,clamptree,meanexp,normalise,nrep,mix,nsplit,suffstatfile,rootfile,true);
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
		is >> datafile >> treefile >> contdatafile >> ancdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> timelinefile >> flextimeline;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag;
		is >> conjpath;
		is >> contdatatype;
		is >> ancdatatype;
		is >> clamproot >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> mix;
		is >> nsplit;
		is >> clamptree;
		is >> suffstatfile;
		is >> rootfile;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;

		if (modeltype == "TIMELINE")	{
			model = new TimeLineClockModel(datafile,treefile,contdatafile,ancdatafile,timelinefile,flextimeline,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,conjpath,contdatatype,ancdatatype,clamproot,clamptree,meanexp,normalise,nrep,mix,nsplit,suffstatfile,rootfile,true);
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
		param_os << datafile << '\t' << treefile << '\t' << contdatafile << '\t' << ancdatafile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << timelinefile << '\t' << flextimeline << '\n';
		param_os << chronoprior << '\t' << meanchi << '\t' << meanchi2 << '\n';
		param_os << clampdiag << '\n';
		param_os << conjpath << '\n';
		param_os << contdatatype << '\n';
		param_os << ancdatatype << '\n';
		param_os << clamproot << '\t' << meanexp << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << df << '\n';
		param_os << priorsigma << '\n';
		param_os << mix << '\n';
		param_os << nsplit << '\n';
		param_os << clamptree << '\n';
		param_os << suffstatfile << '\n';
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
		TimeLineClockChain* chain = new TimeLineClockChain(name);
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
		string ancdatafile = "None";
		string calibfile = "None";
		string timelinefile = "None";
		int flextimeline = 0;
		string suffstatfile = "None";
		string rootfile = "None";
		double rootage = 0;
		double rootstdev = 0;

		int chronoprior = 0;
		double meanchi = 1e-3;
		double meanchi2 = 1e-3;

		string name = "";
		bool clampdiag = false;
		bool clamproot = false;
		bool clamptree = false;
		bool meanexp = false;

		double priorsigma = -1;
		int df = 0;

		int conjpath = -1;
		int contdatatype = 0;
		int ancdatatype = 0;

		int every = 1;
		int until = -1;

		int nrep = 0;
		bool normalise = true;

		string mix = "None";

		int nsplit = 1;

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
				else if (s == "-a")	{
					i++;
					ancdatafile = argv[i];
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
				else if (s == "-clinear")	{
					contdatatype = 2;
				}
				else if (s == "-clogit")	{
					contdatatype = 1;
				}
				else if (s == "-alinear")	{
					ancdatatype = 2;
				}
				else if (s == "-alogit")	{
					ancdatatype = 1;
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
				else if (s == "-mix")	{
					i++;
					mix = argv[i];
				}
				else if (s == "-timeline")	{
					i++;
					timelinefile = argv[i];
				}
				else if (s == "-flextimeline")	{
					i++;
					timelinefile = argv[i];
					flextimeline = 1;
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
				else if (s == "-clamproot")	{
					clamproot = true;
				}
				else if (s == "-fixbl")	{
					clamptree = true;
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
			cerr << "coevol -d <alignment> -t <tree> [-c <phenodata>] [-diag -dNdS -dSomega] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		TimeLineClockChain* chain = new TimeLineClockChain(datafile,treefile,contdatafile,ancdatafile,timelinefile,flextimeline,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,name,clampdiag,conjpath,contdatatype,ancdatatype,clamproot,clamptree,meanexp,normalise,nrep,mix,nsplit,suffstatfile,rootfile,force);
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

