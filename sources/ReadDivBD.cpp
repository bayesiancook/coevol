
#include "Sample.h"
#include "DivBDModel.h"
#include "MeanValTree.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"

class DivBDSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string datafile;
	string calibfile;
	int divmodel;
	double initlambda;
	double finallambda;
	double initmu;
	double finalmu;

	int conjpath;

	int chainburnin;

	public:

	string GetModelType() {return modeltype;}

	DivBDModel* GetModel() {return (DivBDModel*) model;}

	DivBDSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	bool IsAIS()	{
		return chainburnin;
	}

	void Open()	{

		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> treefile >> datafile;
		is >> calibfile;
		is >> divmodel;
		is >> initlambda >> finallambda >> initmu >> finalmu;
		is >> conjpath;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainburnin >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "DIVBD")	{
			model = new DivBDModel(datafile,treefile,calibfile,divmodel,initlambda,finallambda,initmu,finalmu,conjpath,true);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void GetSuffStat()	{

		ofstream os((name + ".suffstat").c_str());

		os << size << '\t' << GetModel()->GetNtaxa() << '\n';
		for (int i=0; i<size; i++)	{
			cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();
			GetModel()->GetChronogram()->specialUpdate();
			GetModel()->WriteSuffStat(os);
		}
		cerr << '\n';
		cerr << "suff stats in " << name << ".suffstat\n";
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
			GetModel()->GetChronogram()->specialUpdate();

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

	bool ss = false;

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
			else if (s == "-ss")	{
				ss = true;
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
		cerr << "readdivbd [-ss -x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	DivBDSample sample(name,burnin,every,until);

	if (sample.IsAIS())	{
		cerr << "sample is from AIS\n";
	}
	else if (ss)	{
		sample.GetSuffStat();
	}
	else	{
		sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);
	}

}



