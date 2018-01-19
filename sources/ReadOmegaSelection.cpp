
#include "Sample.h"
#include "OmegaSelectionModel.h"
#include "GTRSubMatrix.h"

class OmegaSelectionSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string contdatafile;
	string treefile;
	int category;
	double alpha;

	public:

	string GetModelType() {return modeltype;}

	OmegaSelectionModel* GetModel() {return (OmegaSelectionModel*) model;}

	OmegaSelectionSample(string filename, int inburnin, int inevery, int inuntil, double inalpha = 0.1) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> datafile >> contdatafile >> treefile >> category;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size


		// make a new model depending on the type obtained from the file
		if (modeltype == "OMEGADIFF")	{
			model = new OmegaSelectionModel(datafile,contdatafile,treefile,category,false);
		}
		else	{
			cerr << "error when opening file "  << name << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		cerr << "from stream\n";
		model->FromStream(is);

		/*
		cerr << "UPDATE\n";
		model->Update();
		*/

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read(double cutoff)	{

		// prepare the mean and variance
		
		int Nsite = GetModel()->GetSite();

		int K = GetModel()->GetCategory();

		double omegaarray[K][Nsite];
		double omegapp[K][Nsite];


		for(int k=0;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				omegaarray[k][i] = 0;
				omegapp[k][i] = 0;
			}
		}

  
		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '*';

			for(int k=0;k<K;k++)	{
				for(int i=0;i<Nsite;i++)	{
					double tmp = GetModel()->GetOmega(k,i);
					if (isnan(tmp))	{
						cerr << "omega is nan\n";
						exit(1);
					}
					if (isinf(tmp))	{
						cerr << "omega is inf\n";
						exit(1);
					}
					if (tmp < 0)	{
						cerr << "negative omega\n";
						exit(1);
					}
					if (tmp > 100)	{
						cerr << "large omega: " << tmp << '\n';
						exit(1);
					}
					omegaarray[k][i] += tmp;
					if (tmp > 1)	{
						omegapp[k][i]++;
					}
				}
			}

		}
		cerr << '\n';

		
		for(int k=0;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				omegaarray[k][i] /= size;
				omegapp[k][i] /= size;
			}
		}


		ofstream os((name + ".postmeanom").c_str());
		for(int i=0;i<Nsite;i++)	{
			os << i << '\t';
			for(int k=0;k<K;k++)	{
				os << omegaarray[k][i] << '\t' << omegapp[k][i] << '\t';
				if (omegapp[k][i] > cutoff)	{
					cerr << k << '\t' << i << '\t' << omegaarray[k][i] << '\t' << omegapp[k][i] << '\n';
				}
			}
			os << '\n';
		}
		/*
		for(int k=0;k<K;k++)	{
			ostringstream s1;
			s1 << GetName() + "_" << k << ".omega";
			string name1 = s1.str();
			ofstream os(name1.c_str());


			ostringstream s2;
			s2 << GetName() + "_" << k << ".omegapp";
			string name2 = s2.str();
			ofstream ppos(name2.c_str());

			os << Nsite << '\n';
			for(int i=0;i<Nsite;i++)	{
				os << omegaarray[k][i] << '\n';
				ppos << omegapp[k][i] << '\n';
				if (omegapp[k][i] > 0.9)	{
					cerr << k << '\t' << i << '\t' << omegaarray[k][i] << '\t' << omegapp[k][i] << '\n';
				}
			}
		}
		*/
	}

	void PostPred(string basename)	{

		for (int i=0; i<size; i++)	{
			cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			GetModel()->Update();

			ostringstream s;
			s << basename << "_" << i;

			GetModel()->PostPredAli(s.str());
		}

		cerr << "simulated data in " << basename << "_[rep].ali\n";
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	double alpha = 0.1;
	double cutoff = 0.7;

	int ppred = 0;
	string basename = "";

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-a")	{
				i++;
				alpha = atof(argv[i]);
			}
			else if (s == "-c")	{
				i++;
				cutoff = atof(argv[i]);
			}
			else if (s == "-ppred")	{
				ppred = 1;
				i++;
				basename = argv[i];
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

	OmegaSelectionSample sample(name,burnin,every,until,alpha);
	if (ppred)	{
		cerr << "posterior predictive simulation\n";
		sample.PostPred(basename);
	}
	else	{
		cerr << "read\n";
		sample.Read(cutoff);
	}

}


