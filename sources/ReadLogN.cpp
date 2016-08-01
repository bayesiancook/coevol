
#include "Sample.h"
#include "LogNModel.h"
#include "MeanValTree.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"

class LogNSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string datafile;

	string calibfile;
	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	double N0;
	double N1;
	double N2;
	double T0;
	double softa;

	bool conjpath;

	public:

	string GetModelType() {return modeltype;}

	GTRLogNormalModel* GetModel() {return (GTRLogNormalModel*) model;}

	LogNSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> treefile >> datafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		if (chronoprior == 4)	{
			is >> softa;
		}
		if (chronoprior == 5)	{
			is >> N0 >> N1 >> N2 >> T0;
		}
		is >> conjpath;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		// make a new model depending on the type obtained from the file
		if (modeltype == "LOGN")	{
			model = new GTRLogNormalModel(datafile,treefile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,N0,N1,N2,T0,softa,conjpath,false);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		cerr << "read\n";

		MeanExpNormTree* meanratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			meanratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetLogNormalProcess(), GetModel()->GetChronogram());
			meanchrono->Add(GetModel()->GetChronogram());

		}
		cerr << '\n';
		cerr << "normalise\n";

		meanratetree->Normalise();
		ofstream ros((name + ".rates").c_str());
		meanratetree->ToStream(ros);
		cerr << "rates in " << name << ".rates\n";

		meanchrono->Normalise();
		ofstream dos((name + ".dates").c_str());
		meanchrono->ToStream(dos);
		cerr << "dates in " << name << ".dates\n";

		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if (s == "+stdev")	{
				printstdev = true;
			}
			else if (s == "-stdev")	{
				printstdev = false;
			}
			else if (s == "+ci")	{
				printci = true;
			}
			else if (s == "-ci")	{
				printci = false;
			}
			else if (s == "+leaf")	{
				withleaf = true;
			}
			else if (s == "-leaf")	{
				withleaf = false;
			}
			else if (s == "+internal")	{
				withinternal = true;
			}
			else if (s == "-internal")	{
				withinternal = false;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readlogn [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	LogNSample sample(name,burnin,every,until);

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}



