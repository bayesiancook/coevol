
#include "Sample.h"
#include "BranchGenOmegaModel.h"
#include "MeanValTree.h"

class BranchGenOmegaSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	BranchGenOmegaModel* GetModel() {return (BranchGenOmegaModel*) model;}

	BranchGenOmegaSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> chainevery >> chainuntil >> chainsize;
		string contdatafile = "nuc42.maturity";
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "BRANCHGENOMEGA")	{
			model = new BranchGenOmegaModel(datafile,treefile,contdatafile,false);
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

		MeanRateBranchTree* meanomega = new MeanRateBranchTree(GetModel()->GetTree());
		meanomega->Reset();

		MeanRateBranchTree* meangentime = new MeanRateBranchTree(GetModel()->GetTree());
		meangentime->Reset();

		MeanRateBranchTree* meanrate = new MeanRateBranchTree(GetModel()->GetTree(), MEAN);
		meanrate->Reset();

		MeanRateBranchTree* meanratepergen = new MeanRateBranchTree(GetModel()->GetTree(), MEAN);
		meanratepergen->Reset();

		MeanPopSizeTree* meanpopsize = new MeanPopSizeTree(GetModel()->GetTree(), 0.16);
		meanpopsize->Reset();

		cerr << '\n';
		// prepare the mean and variance
		double meanlength = 0;
		double varlength = 0;

		double meansigma = 0;
		double meantheta = 0;
		double meantau = 0;

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

			meansigma += *(GetModel()->sigma);
			meantheta += *(GetModel()->theta);
			meantau += *(GetModel()->tau);

			meanomega->Add(GetModel()->GetOmegaTree(), GetModel()->GetChronogram());
			meangentime->Add(GetModel()->GetGenTimeTree(), GetModel()->GetChronogram());
			GetModel()->GetLogSynRateTree()->specialUpdate();
			GetModel()->GetSynRateTree()->specialUpdate();
			meanrate->Add(GetModel()->GetSynRateTree(), GetModel()->GetChronogram());
			meanratepergen->Add(GetModel()->GetSynRatePerGenTree(), GetModel()->GetChronogram());
			meanpopsize->Add(GetModel()->GetOmegaTree(), GetModel()->GetChronogram());
		}
		cerr << '\n';

		// normalize mean and var
		meanlength /= size;
		varlength /= size;
		varlength -= meanlength * meanlength;

		meansigma /= size;
		meantheta /= size;
		meantau /= size;

		meanomega->Normalise();
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);

		meangentime->Normalise();
		ofstream gos((GetName() + ".postmeangentime.tre").c_str());
		meangentime->ToStream(gos);

		meanrate->Normalise();
		ofstream sos((GetName() + ".postmeands.tre").c_str());
		meanrate->ToStream(sos);

		meanratepergen->Normalise();
		ofstream tos((GetName() + ".postmeandspergen.tre").c_str());
		meanratepergen->ToStream(tos);

		meanpopsize->Normalise();
		ofstream ppos((GetName() + ".postmeanpopsize.tre").c_str());
		meanpopsize->ToStream(ppos);

		// output result into a file, using the usual file name and extension standards
		ofstream os((GetName() + ".postmeanlength").c_str());
		os << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
		
		// can also make a brief output of essential results onto the standard output
		cout << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
		cout << '\n';

		cout << "sigma : " << meansigma << '\n';
		cout << "theta : " << meantheta << '\n';
		cout << "tau   : " << meantau << '\n';
	}

	/*
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
	*/
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

	BranchGenOmegaSample sample(name,burnin,every,until);
	sample.Read();

}

