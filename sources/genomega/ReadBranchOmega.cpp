
#include "Sample.h"
#include "BranchOmegaModel.h"
#include "MeanValTree.h"

class BranchOmegaSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	GeneticCodeType type;

	double alpha;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaModel* GetModel() {return (BranchOmegaModel*) model;}

	BranchOmegaSample(string filename, int inburnin, int inevery, int inuntil, double inalpha = 0.1) : Sample(filename,inburnin,inevery,inuntil)	{
		alpha = inalpha;
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
		is >> datafile >> treefile >> contdatafile;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "BRANCHOMEGA")	{
			model = new BranchOmegaModel(datafile,treefile,contdatafile,false,type);
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

		MeanExpNormTree* meanomega = new MeanExpNormTree(GetModel()->GetTree());
		meanomega->Reset();

		MeanRateBranchTree* meanbranchomega = new MeanRateBranchTree(GetModel()->GetTree());
		meanbranchomega->Reset();

		MeanRateBranchTree* meanrate = new MeanRateBranchTree(GetModel()->GetTree(), MEAN);
		meanrate->Reset();

		MeanRateBranchTree* meangentime = new MeanRateBranchTree(GetModel()->GetTree());
		meangentime->Reset();

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

			meanomega->Add(GetModel()->GetOmegaTree(), GetModel()->GetChronogram());
			meanbranchomega->Add(GetModel()->GetOmegaTree(), GetModel()->GetChronogram());
			meanrate->Add(GetModel()->GetLogNormalTree(), GetModel()->GetChronogram());
			meangentime->Add(GetModel()->GetGenTimeTree(), GetModel()->GetChronogram());
		}
		cerr << '\n';

		// normalize mean and var
		meanlength /= size;
		varlength /= size;
		varlength -= meanlength * meanlength;

		meanomega->Normalise();
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);

		meanbranchomega->Normalise();
		ofstream boos((GetName() + ".postmeanbranchomega.tre").c_str());
		meanbranchomega->ToStream(boos);

		meanrate->Normalise();
		ofstream sos((GetName() + ".postmeands.tre").c_str());
		meanrate->ToStream(sos);

		meangentime->Normalise();
		ofstream ppos((GetName() + ".postmeangentime.tre").c_str());
		meangentime->ToStream(ppos);

		// output result into a file, using the usual file name and extension standards
		ofstream os((GetName() + ".postmeanlength").c_str());
		os << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
		
		// can also make a brief output of essential results onto the standard output
		cout << "mean tree length : " << meanlength << " +/- " << sqrt(varlength) << '\n';
		cout << '\n';

		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meanrate->Tabulate(ssos,false);
		ofstream lssos((GetName() + ".leafpostmeansynrate.tab").c_str());
		meanrate->Tabulate(lssos,true);

		ofstream ggos((GetName() + ".postmeangentime.tab").c_str());
		meangentime->Tabulate(ggos,false);
		ofstream lggos((GetName() + ".leafpostmeangentime.tab").c_str());
		meangentime->Tabulate(lggos,true);

		ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
		meanomega->Tabulate(ooos,false);
		ofstream looos((GetName() + ".leafpostmeanomega.tab").c_str());
		meanomega->Tabulate(looos,true);

		ssos.close();
		lssos.close();
		ooos.close();
		looos.close();
		ggos.close();
		lggos.close();
		string com = "paste " + GetName() + ".postmeansynrate.tab " + GetName() + ".postmeanomega.tab " + GetName() + ".postmeangentime.tab > " + GetName() + ".tab";
		string lcom = "paste " + GetName() + ".leafpostmeansynrate.tab " + GetName() + ".leafpostmeanomega.tab " + GetName() + ".leafpostmeangentime.tab > " + GetName() + ".leaftab";
		system(com.c_str());
		system(lcom.c_str());


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

	double alpha = 0.1;

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
			else if (s == "-a")	{
				i++;
				alpha = atof(argv[i]);
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

	BranchOmegaSample sample(name,burnin,every,until,alpha);
	sample.Read();

}

