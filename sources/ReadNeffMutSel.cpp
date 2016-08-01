
#include "Sample.h"
#include "NeffAAProfileMutSelModel.h"
#include "MeanValTree.h"

class NeffAAProfileMutSelSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	string profilefile;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	NeffAAProfileMutSelModel* GetModel() {return (NeffAAProfileMutSelModel*) model;}

	NeffAAProfileMutSelSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> type;
		is >> datafile >> treefile;
		is >> P;
		is >> profilefile;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "NEFFAAPROFILEMUTSEL")	{
			model = new NeffAAProfileMutSelModel(datafile,treefile,P,true,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		cerr << "UPDATE\n";
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read()	{

		// MeanExpNormTree* meannefftree = new MeanExpNormTree(GetModel()->GetTree());
		// MeanExpNormTree* meannefftree = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meannefftree = new MeanExpNormTree(GetModel()->GetTree(),false,false,false,true,false,true,true);
		meannefftree->Reset();

		// cycle over the sample
		cout << "reading chain...\n";
		for (int i=0; i<size; i++)	{
			GetNextPoint();
			meannefftree->Add((NodeVarTree<Real>*) GetModel()->GetNeffTree(),GetModel()->GetChronogram());
			//cout << *GetModel()->AAProfileConcentration << '\n';
			//cout << *GetModel()->
		}
		meannefftree->Normalise();
		meannefftree->ToStream(cout);

		cout << '\n';
	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int postpred = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				until = atoi(argv[i]);
			}
			else if (s == "-pp")	{
				postpred = 1;
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
		cerr << "readpb [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	NeffAAProfileMutSelSample sample(name,burnin,every,until);
	if (postpred)	{
		//sample.PostPredCompositionalHeterogeneity();
		cout << "No posterior predictive yet...\n";
		cout.flush();
	}
	else	{
		sample.Read();
	}

}



