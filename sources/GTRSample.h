
#include "GTRModel.h"
#include "Sample.h"

class GTRSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;

	public:

	string GetModelType() {return modeltype;}

	GTRModel* GetModel() {return (GTRModel*) model;}

	GTRSample(string name, int burnin, int every, int until) : Sample(name,burnin,every,until)	{
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
		if (modeltype == "GTR")	{
			model = new GTRModel(datafile,treefile);
		}
		else if (modeltype == "CAT")	{
			// to be implemeted
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		cerr << "UPDATE\n";
		model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read()	{

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
		}
		cerr << '\n';

		// normalize mean and var
		meanlength /= size;
		varlength /= size;
		varlength -= meanlength * meanlength;

		// output result into a file, using the usual file name and extension standards
		ofstream os((GetName() + ".postmeanlength").c_str());
		os << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';

		// can also make a brief output of essential results onto the standard output
		cout << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
		cout << '\n';
	}

	// A framework for posterior and posterior predictive stochastic mappings
	void PostPredMapping(int site, bool redundant = true)	{

		// prepare files
		ofstream post_os((name + ".post").c_str());
		ofstream postpred_os((name + ".postpred").c_str());

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> from now on, accessible through GetModel()
			GetNextPoint();

			// update point (to be sure that stochastic mapping will work based on correct model specs)
			GetModel()->Update();

			// sample a mapping for site <site>, matching the data at the leaves (posterior)
			// and write it into <name>.post
			GetModel()->PostMapping(site,post_os, redundant);

			// sample a mapping for site <site>, not trying to match the data at the leaves (posterior predictive)
			// and write it into <name>.postpred
			GetModel()->PostPredMapping(site,postpred_os,redundant);
		}
		cerr << '\n';
		cerr << "posterior mappings in " << name << ".post\n";
		cerr << "post pred mappings in " << name << ".postpred\n";
		cerr << '\n';
	}
};
