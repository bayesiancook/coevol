
#include "Sample.h"
#include "NeffAAProfileInfMixMutSelModel.h"
#include "MeanValTree.h"

class NeffAAProfileInfMixMutSelSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	NeffAAProfileInfMixMutSelModel* GetModel() {return (NeffAAProfileInfMixMutSelModel*) model;}

	NeffAAProfileInfMixMutSelSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "NEFFAAPROFILEINFMIXMUTSEL")	{
			model = new NeffAAProfileInfMixMutSelModel(datafile,treefile,P,true,type);
			//model = new NeffAAProfileInfMixMutSelModel(datafile,treefile,P,false,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		//cerr << "UPDATE\n";
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read()	{

		// create array in which to load profiles for every site
		int Nsite = GetModel()->aaprofilemutselmix->GetSize();
		double** meanAAProfile = new double*[Nsite];
		double* tempAAProfile = new double[Naa];
		for (int i=0; i<Nsite; i++)	{
			meanAAProfile[i] = new double[Naa];
			for (int j=0; j<Naa; j++)	{
				meanAAProfile[i][j] = 0;
			}
		}
		// cycle over the sample
		for (int i=0; i<size; i++)	{
			GetNextPoint();
			//cout << GetModel()->aaprofilemutselmix->GetComponentNumber() << "\n";

			for (int site=0; site<Nsite; site++)	{
				for (int aa=0; aa<Naa; aa++)	{
					tempAAProfile[aa] = (*GetModel()->aaprofilemutselmix->GetRandomVariable(site))[aa];
					meanAAProfile[site][aa] += tempAAProfile[aa];
				}
			}
			//cout << (*GetModel()->aaprofilemutselmix->GetRandomVariable(0))[0] << "\n"; // site 0, just to check
			//tempAAProfile = GetModel()->aaprofilemutselmix->GetRandomVariable(0) << "\n"; // site 0, just to check

			//cout << *GetModel()->AAProfileConcentration << '\n';
			//cout << *GetModel()->
		}

		ofstream os((name + ".ssprofiles").c_str());
		os << Nsite << '\n';
		for (int site=0; site<Nsite; site++)	{
			for (int aa=0; aa<Naa; aa++)	{
				meanAAProfile[site][aa] /= size;
				os << meanAAProfile[site][aa] << '\t';
			}
			os << '\n';
		}
		cout << '\n';

		// eventually delete meanAAProfile and tempAAProfile

	}
	void ReadNeffTree()	{

		// MeanExpNormTree* meannefftree = new MeanExpNormTree(GetModel()->GetTree());
		// MeanExpNormTree* meannefftree = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meannefftree = new MeanExpNormTree(GetModel()->GetTree(),false,false,false,true,false,true,true);
		meannefftree->Reset();

		// cycle over the sample
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

	void PostPredNonsynSubs()	{

		int obs;
		int pred;
		int countPredGreaterThanObs = 0;
		ofstream os((name + ".nsobspred").c_str());
		os << "#obs\tpred\n";

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			GetNextPoint();
			GetModel()->Update();
			obs = GetModel()->ObservedNonsynSubCount();
			os << obs << '\t';
			pred =  GetModel()->PostPredNonsynSubCount();
			os << pred << '\n';
			if (pred >= obs)	{
				countPredGreaterThanObs++;
			}
		}
		double pvalue = (double)(countPredGreaterThanObs)/size;
		os << "#pvalue: " << pvalue << "\n";

		cout << "MyStatistic p-value: " << pvalue << '\n';
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

	NeffAAProfileInfMixMutSelSample sample(name,burnin,every,until);
	if (postpred)	{
		sample.PostPredNonsynSubs();
	}
	else	{
		//sample.Read();
		sample.ReadNeffTree();
	}

}



