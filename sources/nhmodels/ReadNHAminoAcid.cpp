
#include "Sample.h"
#include "NHLGAminoAcidModel.h"
#include "MeanValTree.h"
#include "BranchMeanTree.h"

class NHAminoAcidSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	bool withsep;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	// NHAminoAcidModel* GetModel() {return (NHAminoAcidModel*) model;}
	NHAminoAcidModel* GetModel() {return (NHAminoAcidModel*) model;}

	NHAminoAcidSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> withsep;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "NHAMINOACID")	{
			model = new NHAminoAcidModel(datafile,treefile,contdatafile,withsep,true,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read()	{

		int Ncont = GetModel()->Ncont;

		MeanExpNormTree* meangc = new MeanExpNormTree(GetModel()->GetTree(),true);
		meangc->Reset();

		BranchMeanTree<Profile>* meanaa = new BranchMeanTree<Profile>(GetModel()->GetTree(),Naa);
		meanaa->Reset();

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree());
			tree[k]->Reset();
		}

		double** mean = new double*[Ncont+1];
		double** var = new double*[Ncont+1];
		double** ppp = new double*[Ncont+1];
		for (int k=0; k<=Ncont; k++)	{
			mean[k] = new double[Naa];
			var[k] = new double[Naa];
			ppp[k] = new double[Naa];
			for (int i=0; i<Naa; i++)	{
				mean[k][i] = 0;
				var[k][i] = 0;
				ppp[k][i] = 0;
			}
		}

		// prepare the mean and variance


		CovMatrix* mat = new CovMatrix(GetModel()->GetCovMatrix()->GetDim());
		mat->Reset();
		int dim = mat->GetDim();
		double** pp = new double*[dim];
		for (int k=0; k<dim; k++)	{
			pp[k] = new double[dim];
			for (int l=0; l<dim; l++)	{
				pp[k][l] = 0;
			}
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			meangc->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetGamTree(), GetModel()->Kgc);
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetGamTree(), GetModel()->Kcont + k);
			}

			GetModel()->GetStatTree()->specialUpdate();
			GetModel()->GetBranchStatTree()->specialUpdate();
			meanaa->Add(GetModel()->GetBranchStatTree(), GetModel()->GetGamTree());

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			(*mat) += m;
			for (int i=0; i<dim; i++)	{
				for (int j=0; j<dim; j++)	{
					if (m[i][j] > 0)	{
						pp[i][j] ++;
					}
				}
			}

			for (int k=0; k<=Ncont; k++)	{
				for (int i=0; i<Naa; i++)	{
					double tmp = m[Naa+k][i];
					mean[k][i] += tmp;
					var[k][i] += tmp * tmp;
					if (tmp > 0)	{
						ppp[k][i] ++;
					}
				}
			}
			
		}
		cerr << '\n';

		mat->ScalarMultiplication(1.0/size);
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				pp[i][j] /= size;
			}
		}

		ofstream cout((GetName() + ".general").c_str());
		cout << "mean posterior covariances: \n";
		// cout << *mat << '\n';
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				if (i == j)	{
					// cout << '\t' << '-';
					cout << '\t' << ((double) ((int) (100 * (*mat)[i][j]))) / 100;
				}
				else	{
					cout << '\t' << ((double) ((int) (100 * (*mat)[i][j]))) / 100;
				}
			}
			cout << '\n';
		}

		cout << '\n';

		cout << "correlation coefficients: \n";
		// mat->PrintCorrelationCoefficients(cout);
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				if (i == j)	{
					cout << '\t' << '-';
				}
				else	{
					cout << '\t' << ((double) ((int) (100 * ( (*mat)[i][j] / sqrt((*mat)[i][i]) * (*mat)[j][j]) ))) / 100;
				}
			}
			cout << '\n';
		}
		cout << '\n';

		cout << "posterior probabilities of being >0: \n";
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				if (i == j)	{
					cout << '\t' << '-';
				}
				else	{
					cout << '\t' << ((double) ((int) (100 * pp[i][j]))) / 100;
				}
			}
			cout << '\n';
		}
		cout << '\n';

		/*
		for (int i=0; i<dim; i++)	{
			if (i == GetModel()->Kgc)	{
				cout << (*mat)[GetModel()->Kgc][i] << '\n';
				// cout << '\t' << '-';
			}
			else	{
				cout << (*mat)[GetModel()->Kgc][i] << '\n';
				}
		}
		for (int k=0; k<GetModel()->Ncont; k++)	{
			for (int i=0; i<dim; i++)	{
				if (i == GetModel()->Kcont + k)	{
					cout << (*mat)[GetModel()->Kcont + k][i] << '\n';
					// cout << '\t' << '-';
				}
				else	{
					cout << (*mat)[GetModel()->Kcont + k][i] << '\n';
				}
			}
		}
		*/

		ofstream aaos((GetName() + ".aaprofile").c_str());
		for (int i=0; i<Naa; i++)	{
			for (int k=0; k<=Ncont; k++)	{
				mean[k][i] /= size;
				var[k][i] /= size;
				ppp[k][i] /= size;
				aaos << mean[k][i] << '\t' << sqrt(var[k][i]) << '\t' << ppp[k][i] << '\t';
			}
			aaos << '\n';
		}

		/*
		mat->Diagonalise();
		cout << "eigenvectors :\n";
		mat->PrintEigenVectors(cout);
		cout << '\n';
		*/

		meanaa->Normalise();
		ofstream aatreeos((GetName() + ".postaa.tre").c_str());
		meanaa->ToStream(aatreeos);

		meangc->Normalise();
		ofstream gcos((GetName() + ".postmeangc.tre").c_str());
		meangc->ToStream(gcos);

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
		}

		cerr << "done\n";
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

	NHAminoAcidSample sample(name,burnin,every,until);
	sample.Read();

}

