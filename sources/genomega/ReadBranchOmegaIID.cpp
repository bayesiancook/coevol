
#include "Sample.h"
#include "BranchOmegaIIDModel.h"
#include "MeanValTree.h"

class BranchOmegaIIDSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaIIDModel* GetModel() {return (BranchOmegaIIDModel*) model;}

	BranchOmegaIIDSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> datafile >> treefile;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "BRANCHOMEGAIID")	{
			model = new BranchOmegaIIDModel(datafile,treefile,false);
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

		MeanBranchTree* meanomega = new MeanBranchTree(GetModel()->GetTree());
		meanomega->Reset();

		cerr << '\n';
		// prepare the mean and variance
		double meanlength = 0;
		double varlength = 0;

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			// get total length
			double tmp = GetModel()->GetLength();
			// note that "model", although a ProbModel*, points in fact to a GTRModel
			// GetModel() makes en explicit cast, and return a pointer of the right type (GTRModel*)
	
			// increment mean and var
			meanlength += tmp;
			varlength += tmp * tmp;

			meanomega->Add(GetModel()->GetOmegaTree());
		}
		cerr << '\n';

		// normalize mean and var
		meanlength /= size;
		varlength /= size;
		varlength -= meanlength * meanlength;

		meanomega->Normalise();
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);

		// can also make a brief output of essential results onto the standard output
		cout << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
		cout << '\n';
	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int ppred = 0;
	int site = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-pp")	{
				ppred = 1;
				i++;
				site = atoi(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
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

	BranchOmegaIIDSample sample(name,burnin,every,until);
	sample.Read();

}

