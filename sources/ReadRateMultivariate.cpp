
#include "Sample.h"
#include "ConjugateRateMultivariateModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"

class RateMultivariateSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool gc;
	bool clampdiag;
	bool autoregressive;
	bool clamproot;
	bool meanexp;

	int conjpath;
	int contdatatype;

	public:

	string GetModelType() {return modeltype;}

	RateMultivariateModel* GetModel() {return (RateMultivariateModel*) model;}

	RateMultivariateSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag >> autoregressive >> gc;
		is >> conjpath >> contdatatype >> clamproot >> meanexp;
		is >> chainevery >> chainuntil >> chainsize;

		double priorsigma = 1;
		// make a new model depending on the type obtained from the file
		if (modeltype == "BRANCHOMEGAMULTIVARIATE")	{
			model = new RateMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,clampdiag,autoregressive,conjpath,contdatatype,clamproot,meanexp,false);
		}
		else if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new ConjugateRateMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,gc,autoregressive,conjpath,contdatatype,clamproot,meanexp,false);
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

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());

		MeanExpNormTree* meanrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		// prepare the mean and variance


		CovMatrix* mat = new CovMatrix(GetModel()->GetCovMatrix()->GetDim());
		mat->Reset();
		CovMatrix* matom = new CovMatrix(GetModel()->GetCovMatrix()->GetDim());
		matom->Reset();
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

		/*
		double meanslope = 0;
		list<double> slopedist;
		*/

		// cycle over the sample
		for (int i=0; i<size; i++)	{

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			meanchrono->Add(GetModel()->GetChronogram());

			GetModel()->GetRateTree()->specialUpdate();
			meantree->Add(GetModel()->GetRateTree());

			meanrate->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 0);
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), GetModel()->GetL()+k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());

			/*
			double tmpslope = m[3][4] / m[4][4];
			meanslope += tmpslope;
			slopedist.push_front(tmpslope);
			*/

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
			(*mat) += m;
		}
		// cerr << '\n';


		/*
		meanslope /= size;
		cout << "meanslope : " << meanslope << '\n';
		slopedist.sort();
		int ss = slopedist.size();
		int t = ((int) (((double) ss) / 100 * 2.5));
		list<double>::iterator i = slopedist.begin();
		for (int j=0; j<t; j++)	{
			i++;
		}
		cout << *i << '\t';
		i = slopedist.begin();
		for (int j=0; j<ss-t; j++)	{
			i++;
		}
		cout << *i << '\n';
		*/

		mat->ScalarMultiplication(1.0/size);
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				pp[i][j] /= size;
			}
		}

		ofstream cout((GetName() + ".cov").c_str());
		cout << "entries are in the following order:\n";
		GetModel()->PrintEntries(cout);
		cout << '\n';
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

		/*
		mat->Diagonalise();
		cout << "eigenvectors :\n";
		mat->PrintEigenVectors(cout);
		cout << '\n';
		*/

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

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		meantree->Normalise();
		ofstream tos((GetName() + ".postmean.tre").c_str());
		meantree->ToStream(tos);

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);

		meanrate->Normalise();
		ofstream sos((GetName() + ".postmeanrate.tre").c_str());
		meanrate->ToStream(sos);
		cerr << "reconstructed variation in rate in " << name << ".postmeanrate.tre\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k +1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k + 1 << " in "  << name << ".postmean" << k + 1 << ".tre\n";
		}

		ofstream ssos((GetName() + ".postmeanrate.tab").c_str());
		meanrate->Tabulate(ssos);
		ssos.close();

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}

		cerr << '\n';
	}

	void PostPred()	{
		for (int i=0; i<size; i++)	{
			// cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			GetModel()->Update();

			ostringstream s;
			s << GetName() << "_" << i;
			// GetModel()->PostPredSimu(s.str());
		}
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	bool checkroot = false;
	bool check = false;
	string truefile = "";
	bool ic = false;

	int ppred = 0;

	bool printlog = false;
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
			}
			else if (s == "-cr")	{
				checkroot = true;
			}
			else if (s == "-ic")	{
				ic = true;
			}
			else if (s == "-ch")	{
				check = true;
				i++;
				truefile = argv[i];
			}
			else if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if (s == "+stdev")	{
				printstdev = true;
			}
			else if (s == "-stdev")	{
				printstdev = false;
			}
			else if (s == "+ci")	{
				printci = true;
			}
			else if (s == "-ci")	{
				printci = false;
			}
			else if (s == "+leaf")	{
				withleaf = true;
			}
			else if (s == "-leaf")	{
				withleaf = false;
			}
			else if (s == "+internal")	{
				withinternal = true;
			}
			else if (s == "-internal")	{
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
		cerr << "readratecorrel [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	RateMultivariateSample sample(name,burnin,every,until);

	if (ppred)	{
		sample.PostPred();
		exit(1);
	}
	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}

