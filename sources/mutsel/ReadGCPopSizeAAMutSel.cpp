
#include "Sample.h"
#include "GCPopSizeAAMutSelMatMixModel.h"
#include "MeanValTree.h"

class GCPopSizeAAMutSelSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	GCPopSizeAAMutSelModel* GetModel() {return (GCPopSizeAAMutSelModel*) model;}

	GCPopSizeAAMutSelSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		if (modeltype == "GCPOPSIZEAAMUTSEL")	{
			model = new GCPopSizeAAMutSelModel(datafile,treefile,P,true,type);
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

		/*
		MeanExpNormTree* meanpopsize = new MeanExpNormTree(GetModel()->GetTree());
		meanpopsize->Reset();
		*/

		// true means : logit and not exponential
		MeanExpNormTree* meangc = new MeanExpNormTree(GetModel()->GetTree(),true);
		meangc->Reset();

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			// get total length
			// note that "model", although a ProbModel*, points in fact to a GTRModel
			// GetModel() makes en explicit cast, and return a pointer of the right type (GTRModel*)

			// meanpopsize->Add(GetModel()->GetPopSizeTree(), GetModel()->GetChronogram());
			meangc->Add(GetModel()->GetGCTree(), GetModel()->GetChronogram());
		}
		cerr << '\n';

		/*
		meanpopsize->Normalise();
		ofstream oos((GetName() + ".postmeanpopsize.tre").c_str());
		meanpopsize->ToStream(oos);
		*/

		meangc->Normalise();
		ofstream gcos((GetName() + ".postmeangc.tre").c_str());
		meangc->ToStream(gcos);

	}

	// A framework for posterior and posterior predictive stochastic mappings
	void PostPredCompositionalHeterogeneity()	{

		// prepare files
		/*
		ofstream post_os((name + ".post").c_str());
		ofstream postpred_os((name + ".postpred").c_str());
		*/

		ofstream compos((name + ".comp").c_str());
		compos << size << '\t' << GetModel()->GetNtaxa() << '\t' << Naa << '\n';
		compos << '\n';

		double obs = GetModel()->ObservedCompositionalHeterogeneity(compos);
		cerr << "obs : " << obs << '\n';

		double mean = 0;
		double var = 0;
		double pp = 0;
		// cycle over the sample
		cerr << "total number of points : " << size << '\n';
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> from now on, accessible through GetModel()
			GetNextPoint();
			
			// update point (to be sure that stochastic mapping will work based on correct model specs)
			GetModel()->Update();

			// sample a mapping for site <site>, matching the data at the leaves (posterior)
			// and write it into <name>.post
			double tmp = GetModel()->PostPredCompositionalHeterogeneity(compos);

			mean += tmp;
			var += tmp * tmp;
			if (obs < tmp)	{
				pp++;
			}
		}
		mean /= size;
		var /= size;
		var -= mean * mean;
		pp /= size;

		cout << '\n';
		cout << "observed  : " << obs << '\n';
		cout << "predictive: " << mean << '\t' << sqrt(var) << '\n';
		cout << "pp        : " << pp << '\n';
		cout << "z-score   : " << (obs - mean) / sqrt(var) << '\n';
		cout << '\n';
		/*
		cerr << '\n';
		cerr << "posterior mappings in " << name << ".post\n";
		cerr << "post pred mappings in " << name << ".postpred\n";
		cerr << '\n';
		*/
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

	GCPopSizeAAMutSelSample sample(name,burnin,every,until);
	if (postpred)	{
		sample.PostPredCompositionalHeterogeneity();
	}
	else	{
		sample.Read();
	}

}

