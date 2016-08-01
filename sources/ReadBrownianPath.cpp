
#include "Sample.h"
#include "BrownianModel.h"
#include "MeanValTree.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanDiscretization.h"


class BrownianPathSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string datafile;

	int nSegments;
		int typemodel;
	bool conjpath;

	public:

	string GetModelType() {return modeltype;}

	BrownianModel* GetModel() {return (BrownianModel*) model;}

	BrownianPathSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
 
		is >> conjpath;
		is >> nSegments;
				is >> typemodel;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		// make a new model depending on the type obtained from the file
		if (modeltype == "BrownianPath")	{
			model = new BrownianModel(datafile,treefile,conjpath,nSegments,typemodel,false);
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
				MeanDiscretization* meandiscret = new MeanDiscretization(GetModel()->GetTree());
	  
		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			GetModel()->GetBrownianProcess()->specialUpdate();
			meanratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetBrownianProcess(), GetModel()->GetChronogram());
			meanchrono->Add(GetModel()->GetChronogram());
						meandiscret->Add(GetModel()->GetBrownianProcess()->GetPureBrownianProcess());


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

				meandiscret->Normalise();
		ofstream dscrtos((name + ".discret").c_str());
		meandiscret->ToStream(dscrtos);
		cerr << "discretization in " << name << ".discret\n";

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
		cerr << "readbrownianPath [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	BrownianPathSample sample(name,burnin,every,until);

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}



