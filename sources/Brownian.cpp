
#include "Chain.h"
#include "BrownianModel.h"
#include "StringStreamUtils.h"




class BrownianChain : public Chain	{

	private:
	string modeltype;
	string treefile;
	string datafile;
	string contdatafile;
	int conjpath;
	int conjsigma;
	bool mapSegment;
	int nSegments;
	Segmentation segm;
	bool commut;
	int gc;
	bool omega;
	bool fixtimes;
	bool paral;

	public:

	string GetModelType() {return modeltype;}

	BrownianModel* GetModel() {return (BrownianModel*) model;}

	BrownianChain(string intree, string indata, string incontdatafile, bool inconjpath, bool inmapSegment, int innSegments, Segmentation insegm, int incommut, int ingc, int inomega, int inconjsigma, bool infixtimes, bool inparal, string filename, int force = 1)	{
		modeltype = "BrownianNonCommutative";
		treefile = intree;
		datafile = indata;
		contdatafile = incontdatafile;
		conjpath = inconjpath;
		mapSegment = inmapSegment;
		nSegments = innSegments;
		segm = insegm;
		commut = incommut;
		gc = ingc;
		omega = inomega;
		conjsigma = inconjsigma;
		fixtimes = infixtimes;
		paral = inparal;
		name = filename;
		New(force);
	}

	BrownianChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "BrownianNonCommutative")	{
			model = new BrownianModel(datafile,treefile,contdatafile,conjpath,mapSegment,nSegments,segm,commut,gc,omega,conjsigma,fixtimes, paral,true);
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
		is >> treefile >> datafile >> contdatafile;
		is >> conjpath;
				is >> mapSegment;
		is >> nSegments;
				bool b;
				is >> b;
				segm = b ? SEGM_REGULAR : SEGM_ABSOLUTE;
				is >> commut;
				is >> gc;
				is >> omega;
				is >> conjsigma;
				is >> fixtimes;
				is >> paral;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> every >> until >> size;
		if (modeltype == "BrownianNonCommutative")	{
			model = new BrownianModel(datafile,treefile, contdatafile, conjpath, mapSegment, nSegments, segm, commut, gc, omega, conjsigma,fixtimes,paral,true);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		model->test();
		cerr << "update\n";
		model->Update();
		cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << treefile << '\t' << datafile << '\t' << contdatafile << '\n';

		param_os << conjpath << '\n';
				param_os << mapSegment << '\n';
		param_os << nSegments << '\n';
				param_os << (segm == SEGM_REGULAR) << '\n';
				param_os << commut << '\n';
				param_os << gc << '\n';
				param_os << omega << '\n';
				param_os << conjsigma << '\n';
				param_os << fixtimes << '\n';
				param_os << paral << '\n';
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
		BrownianChain* chain = new BrownianChain(name);
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
		string datafile = "";
		string contdatafile = "None";
		string name = "";
		bool conjpath = true;
		bool mapSegment = false;
		int nSegments = 100;
		Segmentation segm = SEGM_ABSOLUTE;
		// Segmentation segm = SEGM_REGULAR;
		bool commut = false;
		int gc = 0;
		bool omega = false;
		int conjsigma = 1;
		bool fixtimes = false;
		bool paral = false;

		string simufile = "None";
		string simuname = "None";
		int nrep = 1;

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
				else if (s == "-commut")	{
					commut = true;
				}
				else if (s == "-gc")	{
					gc = 1;
				}
				else if (s == "-compo")	{
					gc = 2;
				}
				else if (s == "-omega")	{
					omega = true;
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-fullconjsigma")	{
					conjsigma = 2;
				}
				else if (s == "-conjsigma")	{
					conjsigma = 1;
				}
				else if (s == "-nonconjsigma")	{
					conjsigma = 0;
				}
				else if (s == "-mapsegment")	{
					mapSegment = true;
					segm = SEGM_REGULAR;
				}
				else if (s == "-mapbranch")	{
					mapSegment = false;
				}
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-ns")	{
					i++;
					nSegments = atoi(argv[i]);
				}
				else if ((s == "-sAbs") || (s == "-sabs"))	{
					segm = SEGM_ABSOLUTE;
				}
				else if ((s == "-sReg") || (s == "-srel"))	{
					segm = SEGM_REGULAR;
				}
				else if (s == "-fixtimes")	{
					fixtimes = true;
				}
				else if (s == "-paral")	{
					paral = true;
				}
				else if (s == "-simu")	{
					i++;
					simufile = argv[i];
					i++;
					nrep = atoi(argv[i]);
				}
				else if (s == "-x")	{
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
			cerr << "brownianPath -d <datafile> -t <tree> [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		BrownianChain* chain = new BrownianChain(treefile,datafile,contdatafile,conjpath,mapSegment, nSegments, segm, commut, gc, omega, conjsigma, fixtimes, paral, name,force);
		if (simufile != "None")	{
			chain->GetModel()->Load(simufile);
			for (int rep=0; rep<nrep; rep++)	{
				ostringstream s;
				s << name << "_" << rep;
				cerr << '.';
				chain->GetModel()->PostPredSimu(s.str(),true,true,-1);
			}
			cerr << '\n';
			exit(1);
		}
		chain->SetEvery(every);
		chain->SetUntil(until);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
				//chain->GetModel()->test();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
}
