
#include "Sample.h"
#include "ConjugateGammaTreeOmegaMultivariateModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"

class GammaTreeOmegaMultivariateSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	bool gc;
	bool clampdiag;
	bool autoregressive;
	GeneticCodeType type;


	public:

	string GetModelType() {return modeltype;}

	GammaTreeOmegaMultivariateModel* GetModel() {return (GammaTreeOmegaMultivariateModel*) model;}

	GammaTreeOmegaMultivariateSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> clampdiag >> autoregressive >> gc;

		bool check;
		is >> check;
		if (check)	{
			cerr << "error when reading file\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		double priorsigma = 1;
		// make a new model depending on the type obtained from the file
		if (modeltype == "GAMMATREEOMEGAMULTIVARIATE")	{
			model = new GammaTreeOmegaMultivariateModel(datafile,treefile,contdatafile,priorsigma,gc,clampdiag,autoregressive,false,type);
		}
		else if (modeltype == "CONJUGATEGAMMATREEOMEGAMULTIVARIATE")	{
			model = new ConjugateGammaTreeOmegaMultivariateModel(datafile,treefile,contdatafile,priorsigma,gc,autoregressive,false,type);
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

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void ReadIC()	{

		MeanICTree* ictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		ictree->Reset();
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ictree->Add();
		}
		cerr << '\n';
		ictree->Normalise();
		ofstream os((name + ".postmeanic").c_str());
		ictree->Tabulate(os,false);
		ofstream los((name + ".postmeanleafic").c_str());
		ictree->Tabulate(los,true);
		cerr << "mean independent contrasts in " << name  << ".postmeanic\n";
		cerr << "for terminal branches only in " << name  << ".postmeanleafic\n";
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		int Ncont = GetModel()->Ncont;

		MeanBranchTree* meantree = new MeanBranchTree(GetModel()->GetTree(),false);

		MeanExpNormTree* meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
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

		double*** proj = new double**[dim];
		double*** pproj = new double**[dim];
		double*** tempproj = new double**[dim];
		for (int i=0; i<dim; i++)	{
			proj[i] = new double*[dim-1];
			pproj[i] = new double*[dim-1];
			tempproj[i] = new double*[dim-1];
			for (int k=0; k<dim-1; k++)	{
				proj[i][k] = new double[dim-1];
				pproj[i][k] = new double[dim-1];
				tempproj[i][k] = new double[dim-1];
				for (int l=0; l<dim-1; l++)	{
					proj[i][k][l] = 0;
					pproj[i][k][l] = 0;
				}
			}
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			// cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			meantree->Add(GetModel()->GetGammaTree());

			meanomega->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetGammaTree(), 0);
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetGammaTree(), GetModel()->GetL()+k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			(*mat) += m;
			for (int i=0; i<dim; i++)	{
				for (int j=0; j<dim; j++)	{
					if (m[i][j] > 0)	{
						pp[i][j] ++;
					}
				}
			}

			for (int k=0; k<dim; k++)	{
				m.Project(k,tempproj[k]);

				for (int i=0; i<dim-1; i++)	{
					for (int j=0; j<dim-1; j++)	{
						proj[k][i][j] += tempproj[k][i][j];
						if (tempproj[k][i][j] > 0 )	{
							pproj[k][i][j] ++;
						}
					}
				}
			}
		}
		// cerr << '\n';

		mat->ScalarMultiplication(1.0/size);
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				pp[i][j] /= size;
			}
		}

		ofstream cout((GetName() + ".general").c_str());
		cout << "mean posterior covariances: \n";
		cout << *mat << '\n';
		cout << '\n';

		cout << "correlation coefficients: \n";
		mat->PrintCorrelationCoefficients(cout);
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

		mat->Diagonalise();
		cout << "eigenvectors :\n";
		mat->PrintEigenVectors(cout);
		cout << '\n';

		for (int k=0; k<dim; k++)	{
			cout << '\n';
			cout << "projection orthogonal to component " << k << "\n";
			for (int i=0; i<dim-1; i++)	{
				for (int j=0; j<dim-1; j++)	{
					proj[k][i][j]/= size;
					pproj[k][i][j]/= size;
				}
			}
			for (int i=0; i<dim-1; i++)	{
				for (int j=0; j<dim-1; j++)	{
					cout << proj[k][i][j] << '\t';
				}
				cout << '\t';
				for (int j=0; j<dim-1; j++)	{
					if (i == j)	{
						cout << '\t' << '-';
					}
					else	{
						cout << '\t' << ((double) ((int) (100 * pproj[k][i][j]))) / 100;
					}
				}
				cout << '\n';
			}
			cout << '\n';
		}

		cerr << "covariance matrix in " << name << ".general\n";
		cerr << '\n';

		meantree->Normalise();
		ofstream tos((GetName() + ".postmean.tre").c_str());
		meantree->ToStream(tos);

		meanomega->Normalise();
		/*
		ofstream popos((GetName() + ".popsize.tre").c_str());
		meanomega->PrintPopSizeTree(popos,0.18,10000,"HUMAN");
		ofstream popos2((GetName() + ".popsize2.tre").c_str());
		meanomega->PrintPopSizeTree2(popos2,10000,100000,"HUMAN","MOUSE");
		*/
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);
		cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k << " in "  << name << ".postmean" << k << ".tre\n";
		}

		ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
		meanomega->Tabulate(ooos);
		ooos.close();

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}

		/*
		ostringstream s;
		s << "paste " << GetName() << ".postmeanomega.tab ";
		for (int k=0; k<Ncont; k++)	{
			s << GetName() << ".postmean" << k << ".tab ";
		}
		s << " > " << GetName() << ".tab";
		// cerr << s.str() << '\n';
		system(s.str().c_str());
		*/

		cerr << '\n';
	}

	void Check(string truefile)	{

		int Ncheck = GetModel()->GetNcheck();
		cerr << Ncheck << '\n';
		cerr << "size : " << size << '\n';
		double** tab = new double*[size];
		for (int i=0; i<size; i++)	{
			tab[i] = new double[Ncheck];
		}

		double* trueval = new double[Ncheck];
		double* mean = new double[Ncheck];
		for (int k=0; k<Ncheck; k++)	{
			mean[k] = 0;
		}
		double* var = new double[Ncheck];
		for (int k=0; k<Ncheck; k++)	{
			var[k] = 0;
		}

		// read true value

		ifstream is(truefile.c_str());

		GetModel()->FromStream(is);
		GetModel()->GetCheck(trueval);

		for (int i=0; i<size; i++)	{
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();
			GetModel()->GetCheck(tab[i]);
			for (int k=0; k<Ncheck; k++)	{
				mean[k] += tab[i][k];
				var[k] += tab[i][k] * tab[i][k];
			}
		}
		ofstream os((name + ".check").c_str());
		ofstream vos((name + ".vcheck").c_str());
		for (int k=0; k<Ncheck; k++)	{
			mean[k] /= size;
			var[k] /= size;
			var[k] -= mean[k] * mean[k];
			int n = 0;
			for (int i=0; i<size; i++)	{
				if (trueval[k] > tab[i][k])	{
					n++;
				}
			}
			os << ((double) n) / size << '\t';
			vos << trueval[k] << '\t' << mean[k] << '\t' << sqrt(var[k]) << '\t' << ((double) n) / size << '\n';
		}
		os << '\n';
		vos << '\n';

		for (int k=0; k<Ncheck; k++)	{
			delete[] tab[k];
		}
		delete[] tab;
		delete[] trueval;
	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	bool check = false;
	string truefile = "";
	bool ic = false;

	int ppred = 0;
	int site = 0;

	bool printlog = true;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

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
			else if (s == "-ic")	{
				ic = true;
			}
			else if (s == "-ch")	{
				check = true;
				i++;
				truefile = argv[i];
			}
			else if (s == "+printlog")	{
				printlog = true;
			}
			else if (s == "-printlog")	{
				printlog = false;
			}
			else if (s == "+printmean")	{
				printmean = true;
			}
			else if (s == "-printmean")	{
				printmean = false;
			}
			else if (s == "+printstdev")	{
				printstdev = true;
			}
			else if (s == "-printstdev")	{
				printstdev = false;
			}
			else if (s == "+printci")	{
				printci = true;
			}
			else if (s == "-printci")	{
				printci = false;
			}
			else if (s == "+withleaf")	{
				withleaf = true;
			}
			else if (s == "-printwithleaf")	{
				withleaf = false;
			}
			else if (s == "+printinternal")	{
				withinternal = true;
			}
			else if (s == "-printinternal")	{
				withinternal = false;
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

	GammaTreeOmegaMultivariateSample sample(name,burnin,every,until);

	if (check)	{
		sample.Check(truefile);
		exit(1);
	}
	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}

