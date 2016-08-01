
#include "Sample.h"
#include "ConjugateTamuraModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"

class TamuraSample : public Sample	{

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

	bool clampdiag;
	bool autoregressive;
	bool clamproot;
	bool meanexp;
	GeneticCodeType type;

	int conjpath;
	int contdatatype;

	int omegaratiotree;

	bool normalise;
	int nrep;

	public:

	string GetModelType() {return modeltype;}

	TamuraModel* GetModel() {return (TamuraModel*) model;}

	TamuraSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag >> autoregressive;
		is >> conjpath;
		is >> contdatatype;
		is >> omegaratiotree;
		is >> clamproot >> meanexp;
		is >> normalise >> nrep;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		double priorsigma = 1;
		// make a new model depending on the type obtained from the file
		if (modeltype == "TAMURA")	{
			model = new TamuraModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,clampdiag,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep,false,type);
		}
		else if (modeltype == "CONJUGATETAMURA")	{
			model = new ConjugateTamuraModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,nrep,false,type);
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

		MeanExpNormTree* meangcsynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanatsynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meantstv = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree* meanomega = 0;
		if (omegaratiotree)	{
			meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

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

		double** ppom = new double*[dim];
		for (int k=0; k<dim; k++)	{
			ppom[k] = new double[dim];
			for (int l=0; l<dim; l++)	{
				ppom[k][l] = 0;
			}
		}

		double*** proj = new double**[dim];
		double*** pproj = new double**[dim];
		double*** projom = new double**[dim];
		double*** pprojom = new double**[dim];
		double*** tempproj = new double**[dim];
		for (int i=0; i<dim; i++)	{
			proj[i] = new double*[dim-1];
			pproj[i] = new double*[dim-1];
			projom[i] = new double*[dim-1];
			pprojom[i] = new double*[dim-1];
			tempproj[i] = new double*[dim-1];
			for (int k=0; k<dim-1; k++)	{
				proj[i][k] = new double[dim-1];
				pproj[i][k] = new double[dim-1];
				projom[i][k] = new double[dim-1];
				pprojom[i][k] = new double[dim-1];
				tempproj[i][k] = new double[dim-1];
				for (int l=0; l<dim-1; l++)	{
					proj[i][k][l] = 0;
					pproj[i][k][l] = 0;
					projom[i][k][l] = 0;
					pprojom[i][k][l] = 0;
				}
			}
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{

			GetNextPoint();

			meanchrono->Add(GetModel()->GetChronogram());

			GetModel()->GetSynRateTree()->specialUpdate();
			meantree->Add(GetModel()->GetSynRateTree());

			meangcsynrate->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 0);
			meanatsynrate->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 1);
			meantstv->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 2);
			if (omegaratiotree)	{
				meanomega->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 3);
			}
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), GetModel()->GetL()+k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());

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

		if (omegaratiotree)	{
			meanomega->Normalise();
			ofstream oos((GetName() + ".postmeanomega.tre").c_str());
			meanomega->ToStream(oos);
			cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
		}

		meangcsynrate->Normalise();
		ofstream sos((GetName() + ".postmeangcsynrate.tre").c_str());
		meangcsynrate->ToStream(sos);
		cerr << "reconstructed variation in Ks in " << name << ".postmeangcsynrate.tre\n";

		meanatsynrate->Normalise();
		ofstream satos((GetName() + ".postmeanatsynrate.tre").c_str());
		meanatsynrate->ToStream(satos);
		cerr << "reconstructed variation in Ks in " << name << ".postmeanatsynrate.tre\n";

		meantstv->Normalise();
		ofstream ststvos((GetName() + ".postmeantstv.tre").c_str());
		meantstv->ToStream(ststvos);
		cerr << "reconstructed variation in Ks in " << name << ".postmeantstv.tre\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
		}

		ofstream ssos((GetName() + ".postmeangcsynrate.tab").c_str());
		meangcsynrate->Tabulate(ssos);
		ssos.close();

		ofstream ssatos((GetName() + ".postmeanatsynrate.tab").c_str());
		meanatsynrate->Tabulate(ssatos);
		ssatos.close();

		ofstream sststvos((GetName() + ".postmeantstv.tab").c_str());
		meantstv->Tabulate(sststvos);
		sststvos.close();

		if (omegaratiotree)	{
			ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
			meanomega->Tabulate(ooos);
			ooos.close();
		}

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
			// tree[k]->TabulateDistribution(os);
			// tree[k]->TabulateDistribution(os,false);
			/*
			ostringstream s2;
			s2 << GetName() << ".postmean" << k << ".loo";
			ofstream os2(s2.str().c_str());
			tree[k]->LooTabulate(os2);
			*/
		}
	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	string truefile = "";
	bool ic = false;

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

			if (s == "-ic")	{
				ic = true;
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
		cerr << "readcoevol [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	TamuraSample sample(name,burnin,every,until);

	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}



