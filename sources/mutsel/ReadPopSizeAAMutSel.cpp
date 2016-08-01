
#include "Sample.h"
#include "PopSizeAAMutSelModel.h"
#include "MeanValTree.h"

class PopSizeAAMutSelSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	PopSizeAAMutSelModel* GetModel() {return (PopSizeAAMutSelModel*) model;}

	PopSizeAAMutSelSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		if (modeltype == "POPSIZEAAMUTSEL")	{
			model = new PopSizeAAMutSelModel(datafile,treefile,P,false,type);
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

		MeanExpNormTree* meanpopsize = new MeanExpNormTree(GetModel()->GetTree());
		meanpopsize->Reset();

		double** tmps = new double*[Naa];
		double** means = new double*[Naa];
		double** vars = new double*[Naa];
		for (int i=0; i<Naa; i++)	{
			tmps[i] = new double[Naa];
			means[i] = new double[Naa];
			vars[i] = new double[Naa];
			for (int j=0; j<Naa; j++)	{
				tmps[i][j] = 0;
				means[i][j] = 0;
				vars[i][j] = 0;
			}
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			// get total length
			// note that "model", although a ProbModel*, points in fact to a GTRModel
			// GetModel() makes en explicit cast, and return a pointer of the right type (GTRModel*)

			GetModel()->GetAASelectionCoefficients(tmps);
			for (int i=0; i<Naa; i++)	{
				for (int j=0; j<Naa; j++)	{
					means[i][j] += tmps[i][j];
					vars[i][j] += tmps[i][j] * tmps[i][j];
				}
			}

			meanpopsize->Add(GetModel()->GetPopSizeTree(), GetModel()->GetChronogram());
		}
		cerr << '\n';

		ofstream os((GetName() + ".aaselcoeff").c_str());
		for (int i=0; i<Naa; i++)	{
			for (int j=i+1; j<Naa; j++)	{
				means[i][j] /= size;
				vars[i][j] /= size;
				vars[i][j] -= means[i][j] * means[i][j];
				os << AminoAcids[i] << ' ' << AminoAcids[j] << '\t' << means[i][j] << '\t' << sqrt(vars[i][j]) << '\n';
			}
		}

		meanpopsize->Normalise();
		ofstream oos((GetName() + ".postmeanpopsize.tre").c_str());
		meanpopsize->ToStream(oos);
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

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

	PopSizeAAMutSelSample sample(name,burnin,every,until);
	sample.Read();

}

