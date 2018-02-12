
#include "Sample.h"
#include "DirichletCodonUsageSelectionModelMS.h"
#include "GTRSubMatrix.h"

class DirichletCodonUsageSelectionMSSample : public Sample	{

	private:
	string modeltype;
	string datafile;
    string contdatafile;
	string treefile;
    int ncond;
    int fixglob;
    int codonmodel;

    double alpha;

	public:

	string GetModelType() {return modeltype;}

	DirichletCodonUsageSelectionModelMS* GetModel() {return (DirichletCodonUsageSelectionModelMS*) model;}

	DirichletCodonUsageSelectionMSSample(string filename, int inburnin, int inevery, int inuntil, double inalpha = 0.1) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> datafile >> contdatafile >> treefile;
        is >> ncond;
        is >> fixglob >> codonmodel;
		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "DIFFSELDIR")	{
			cerr << "CREATE\n";
			model = new DirichletCodonUsageSelectionModelMS(datafile,contdatafile,treefile,ncond,fixglob,codonmodel,false);
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

	void ReadFalsePositives(int ncat)	{

		int Nsite = GetModel()->GetSite();
		int Naa = GetModel()->GetaaState();
		int K = GetModel()->GetCategory();

		double*** meansel = new double**[K];
		double*** ppsel = new double**[K];

		int kmin = 1;

		for (int k=kmin; k<K; k++)	{
			meansel[k] = new double*[Nsite];
			ppsel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				meansel[k][i] = new double[Naa];
				ppsel[k][i] = new double[Naa];
				for (int j=0; j<Naa; j++)	{
					meansel[k][i][j] = 0;
					ppsel[k][i][j] = 0;
				}
			}
		}

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '.';

			//get selection profile
			for(int k=kmin;k<K;k++)	{
				for(int i=0;i<Nsite;i++)	{
					double logsum = 0; 
					for(int j=0;j<Naa;j++)	{ 
						double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j));
						logsum += tmp;
					}
					for(int j=0;j<Naa;j++)	{  
						double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j)) - (1.0/20) * (logsum);
						meansel[k][i][j] += tmp;
						if(tmp > 0)	{
							ppsel[k][i][j]++;
						}
					}
				}
			}
		}
		cerr << '\n';

		double** falsepositives = new double*[K];
		for(int k=kmin;k<K;k++)	{
			falsepositives[k] = new double[ncat];
			for (int cat=0; cat<ncat; cat++)	{
				falsepositives[k][cat] = 0;
			}
		}
		for(int k=kmin;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{ 
					meansel[k][i][j] /= size;
					ppsel[k][i][j] /= size;
					int cat = (int) (10*(ppsel[k][i][j] - 0.5));
					if (cat < 0)	{
						cat = -cat;
					}
					if (cat == ncat)	{
						cat--;
					}
					falsepositives[k][cat]++;
				}
			}
		}

		for(int k=kmin;k<K;k++)	{
			for (int cat=0; cat<ncat; cat++)	{
				for (int cat2=cat+1; cat2<ncat; cat2++)	{
					falsepositives[k][cat] +=falsepositives[k][cat2];
				}
			}
		}

		for(int k=kmin;k<K;k++)	{
			for (int cat=0; cat<ncat; cat++)	{
				falsepositives[k][cat] /= Nsite*Naa;
			}
		}


		ofstream cout((GetName() + ".fp").c_str());
		cout << '\n';
		cout << "false positive rates:\n";
		cout << '\n';
		for(int k=kmin;k<K;k++)	{
			cout << "category : " << k << '\n';
			for (int cat=0; cat<ncat; cat++)	{
				cout << 0.5 + ((double) cat) / 2 / ncat << '\t' << 0.5 + ((double) (cat+1)) / 2 / ncat << '\t';
				cout << 100 * falsepositives[k][cat] << '\t' << Nsite*Naa * falsepositives[k][cat] << '\n';
			}
			cout << '\n';
		}
	}

	void Read(double cutoff)	{

		// prepare the mean and variance
		
		int Nsite = GetModel()->GetSite();
		int K = GetModel()->GetCategory();

		double*** meansel = new double**[K];
		for (int k=1; k<K; k++)	{
			meansel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				meansel[k][i] = new double[Naa];
				for (int a=0; a<Naa; a++)	{
					meansel[k][i][a] = 0;
				}
			}
		}

		double*** ppsel = new double**[K];
		for (int k=1; k<K; k++)	{
			ppsel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				ppsel[k][i] = new double[Naa];
				for (int a=0; a<Naa; a++)	{
					ppsel[k][i][a] = 0;
				}
			}
		}

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '.';

			for (int k=1; k<K; k++)	{
				for (int i=0; i<Nsite; i++)	{

					double delta[Naa];
					for(int j=0;j<Naa;j++)	{ 
						delta[j] = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j));
					}

					double mean = 0;
					for (int a=0; a<Naa; a++)	{
						mean += delta[a];
					}
					mean /= Naa;

					for (int a=0; a<Naa; a++)	{
						double tmp = delta[a] - mean;
						meansel[k][i][a] += tmp;
						if (tmp > 0)	{
							ppsel[k][i][a] ++;
						}
					}
				}
			}
		}
		cerr << '\n';

		for (int k=1; k<K; k++)	{

			ostringstream s1,s2;
			s1 << GetName() << "_" << k << ".meandiffsel";
			s2 << GetName() << "_" << k << ".signdiffsel";
			ofstream os1(s1.str().c_str());
			ofstream os2(s2.str().c_str());

			for (int i=0; i<Nsite; i++)	{

				os1 << i;
				for (int a=0; a<Naa; a++)	{
					ppsel[k][i][a] /= size;
					os1 << '\t' << ppsel[k][i][a];
					
				}
				for (int a=0; a<Naa; a++)	{
					meansel[k][i][a] /= size;
					os1 << '\t' << meansel[k][i][a];
				}
				os1 << '\n';

				for (int a=0; a<Naa; a++)	{
                    if ((ppsel[k][i][a] > cutoff) || (ppsel[k][i][a] < 1-cutoff)) {
						os2 << i << '\t' << AminoAcids[a] << '\t' << ppsel[k][i][a] << '\t' << meansel[k][i][a] << '\n';
					}
				}
			}
		}

		ofstream os((GetName() + ".postmeansel").c_str());
		for (int i=0; i<Nsite; i++)	{
			os << i << '\t';
			int k = K-1;
			// for (int k=1; k<K; k++)	{
			for (int a=0; a<Naa; a++)	{
				os << meansel[k][i][a] << '\t' << ppsel[k][i][a] << '\t';
			}
			// }
			os << '\n';
		}
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
	double cutoff = 0.9;

	int ppred = 0;
	string basename = "";

	int fp = 0;
	int ncat = 10;

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
			else if (s == "-fp")	{
				fp = 1;
			}
			else if (s == "-ncat")	{
				i++;
				ncat = atoi(argv[i]);
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

	DirichletCodonUsageSelectionMSSample sample(name,burnin,every,until,alpha);
	if (ppred)	{
		cerr << "posterior predictive simulation\n";
		sample.PostPred(basename);
	}
	else if (fp)	{
		sample.ReadFalsePositives(ncat);
	}
	else    {
		sample.Read(cutoff);
	}
}


