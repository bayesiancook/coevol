
#include "Sample.h"
#include "AncestralCovarianceModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"


void confint(list<double>& l, double& mean, double& median, double& min, double& max, double c = 0.975)	{

	mean = 0;

	for (list<double>::iterator i=l.begin(); i!=l.end(); i++)	{
		mean += *i;
	}
	mean /= l.size();

	l.sort();
	int n = ((int) (((double) l.size()) * (1-c)));
	list<double>::const_iterator i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	min = *i;
	n = ((int) (((double) l.size()) * c));
	i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	max = *i;

	n = ((int) (((double) l.size()) * 0.5));
	i = l.begin();
	for (int j=0; j<n; j++)	{
		i++;
	}
	median = *i;
}

class AncestralCovarianceSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string contdatafile;
	string ancdatafile;
	bool clampdiag;
	int autoregressive;
	bool clamproot;
	bool meanexp;
	bool withdrift;
	int kalman;
	int contdatatype;
	int ancdatatype;
	int df;

	string bounds;
	string mix;

	public:

	string GetModelType() {return modeltype;}

	AncestralCovarianceModel* GetModel() {return (AncestralCovarianceModel*) model;}

	AncestralCovarianceSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		double priorsigma = 1;

		// read model type, and other standard fields
		is >> modeltype;
		is >> treefile >> contdatafile >> ancdatafile ;
		is >> clampdiag >> autoregressive;
		is >> contdatatype;
		is >> ancdatatype;
		is >> clamproot >> meanexp;
		is >> df;
		is >> priorsigma;
		is >> withdrift;
		is >> kalman;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		// make a new model depending on the type obtained from the file
		if (modeltype == "ANCESTRALCOV")	{
			model = new AncestralCovarianceModel(treefile,contdatafile,ancdatafile,priorsigma,df,clampdiag,autoregressive,contdatatype,ancdatatype,clamproot,meanexp,withdrift,kalman,true);
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

	void ReadSim()	{
	// void ReadSim(string refname)	{

		/*
		ifstream is(refname.c_str());
		double x0, x1, x2;
		is >> x0 >> x1 >> x2;
		if (x0 != 0)	{
			cerr << "error in readsim\n";
			exit(1);
		}
		*/

		list<double> y0;
		list<double> y1;
		list<double> y2;
		list<double> c;

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			double tmp0, tmp1, tmp2;
			GetModel()->GetRootValues(tmp0,tmp1,tmp2);
			y0.push_back(tmp0);
			// y1.push_back(tmp1 - tmp0);
			// y2.push_back(tmp2 - tmp0);
			y1.push_back(tmp1);
			y2.push_back(tmp2);
			double tmp = GetModel()->GetCorrelCoeff();
			c.push_back(tmp);
		}
		cerr << '\n';

		double mean, mean0, mean1, mean2;
		double med, med0, med1, med2;
		double min, min0, min1, min2;
		double max, max0, max1, max2;

		confint(y0,mean0,med0,min0,max0);
		confint(y1,mean1,med1,min1,max1);
		confint(y2,mean2,med2,min2,max2);

		confint(c,mean,med,min,max);

		ofstream os((name + ".postmeanroot").c_str());
		os << mean0 << '\t' << med0 << '\t' << min0 << '\t' << max0 << '\n';
		os << mean1 << '\t' << med1 << '\t' << min1 << '\t' << max1 << '\n';
		os << mean2 << '\t' << med2 << '\t' << min2 << '\t' << max2 << '\n';
		os << mean << '\t' << med << '\t' << min << '\t' << max << '\n';
		cerr << "stats in " << name << ".postmeanroot\n";
		cerr << '\n';
	}

	void Read(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal, bool withanc, string mulreg)	{

		cerr << "read\n";
		int Ncont = GetModel()->Ncont;
		int Nanc = GetModel()->Nanc;

		int N = withanc ? Ncont + Nanc : Ncont;

		MeanExpNormTree** tree = new MeanExpNormTree*[N];
		for (int k=0; k<N; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		int dim = GetModel()->GetCovMatrix()->GetDim();
		// MeanCovMatrix*  mat = new MeanCovMatrix(dim);

		SumConstrainedMeanCovMatrix* expandedmat = 0;
		MeanCovMatrix* mat = 0;
		if (GetModel()->GetBasis())	{
			int* kcomp = new int[1];
			int ncomp = 1;
			kcomp[0] = Ncont;
			SumConstrainedMapping** map = new SumConstrainedMapping*[1];
			map[0] = GetModel()->GetBasis();
			expandedmat = new SumConstrainedMeanCovMatrix(dim,ncomp,kcomp,map);
			mat = expandedmat;
		}
		else	{
			mat = new MeanCovMatrix(dim);
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), k);
			}
			if (withanc)	{
				for (int k=0; k<Nanc; k++)	{
					tree[Ncont + k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), Ncont + k);
				}
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			mat->Add(&m);

			/*
			if (expandedmat)	{
				expandedmat->Add(&m);
			}
			*/
		}
		cerr << '\n';
		cerr << "normalise\n";

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());
		cout << "entries are in the following order:\n";
		GetModel()->PrintEntries(cout);
		cout << '\n';
		cout << *mat;
		cout << '\n';
		mat->PrintPropVariances(cout);

		if (mulreg != "")	{
			ofstream cout((GetName() + ".marginalcov").c_str());
			if (mulreg.size() != ((unsigned int) mat->GetDim()))	{
				cerr << "error when specifying multiple regression : " << mulreg << '\n';
				exit(1);
			}

			int* indexarray = new int[mat->GetDim()];
			int dim = 0;
			for (int l=0; l<mat->GetDim(); l++)	{
				if (mulreg[l] == '1')	{
					indexarray[l] = 1;
					dim++;
				}
				else if (mulreg[l] == '0')	{
					indexarray[l] = 0;
				}
				else	{
					cerr << "error when specifying multiple regression : " << mulreg << '\n';
					exit(1);
				}
			}

			cout << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout, indexarray);

			cout << '\n';
			PartialMeanCovMatrix* partmat = new PartialMeanCovMatrix(mat,indexarray,dim);

			ReducedMeanCovMatrix* redmat = new ReducedMeanCovMatrix(mat,indexarray,dim);
			cout << '\n';
			cout << "partial correl\n";
			cout << '\n';
			cout << *partmat;

			ofstream cout2((GetName() + ".controlcov").c_str());
			cout2 << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout2, indexarray);
			cout2 << '\n';
			cout2 << "reduced correl\n";
			cout2 << '\n';
			redmat->PrintCovariances(cout2);
			redmat->PrintCorrel(cout2);
			redmat->PrintPosteriorProbs(cout2);
			// cout2 << *redmat;
			cout2 << '\n';
			delete partmat;
			delete redmat;
		}


		/*
		if (mulreg != "")	{
			if (mulreg.size() != ((unsigned int) mat->GetDim()))	{
				cerr << "error when specifying multiple regression : " << mulreg << '\n';
				exit(1);
			}

			int* indexarray = new int[mat->GetDim()];
			int dim = 0;
			for (int l=0; l<mat->GetDim(); l++)	{
				if (mulreg[l] == '1')	{
					indexarray[l] = 1;
					dim++;
				}
				else if (mulreg[l] == '0')	{
					indexarray[l] = 0;
				}
				else	{
					cerr << "error when specifying multiple regression : " << mulreg << '\n';
					exit(1);
				}
			}

			cout << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout, indexarray);

			cout << '\n';
			cerr << dim << '\t' << mat->GetDim() << '\n';
			PartialMeanCovMatrix* partmat = new PartialMeanCovMatrix(mat,indexarray,dim);
			// ReducedMeanCovMatrix* redmat = new ReducedMeanCovMatrix(mat,indexarray,dim);
			cout << '\n';
			// cout << *redmat;
			cout << *partmat;
			delete partmat;
			// delete redmat;
		}
		else	{
			// cout << "entries are in the following order:\n";
			// GetModel()->PrintEntries(cout);

			cout << '\n';
			cout << *mat;
			cout << '\n';
			mat->PrintPropVariances(cout);

			if (expandedmat)	{
				expandedmat->Normalize();
				ofstream cout((GetName() + ".expandedcov").c_str());
				cout << '\n';
				cout << *expandedmat;
				cout << '\n';
				expandedmat->PrintPropVariances(cout);
			}
		}
		*/


		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
		}
		if (withanc)	{
			for (int k=0; k<Nanc; k++)	{
				tree[Ncont + k]->Normalise();
				ostringstream s;
				s << GetName() << ".postmeananc" << k+1 << ".tre";
				ofstream os(s.str().c_str());
				tree[Ncont + k]->ToStream(os);
				cerr << "reconstructed variation in ancestral character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			}
		}

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->SetWithDepth(GetModel()->GetLengthTree());
			tree[k]->Tabulate(os);
		}

		if (withanc)	{
			for (int k=0; k<Nanc; k++)	{
				ostringstream s;
				s << GetName() << ".postmeananc" << k+1 << ".tab";
				ofstream os(s.str().c_str());
				tree[Ncont + k]->Tabulate(os);

			}
		}

		cerr << '\n';
	}

	void ReadWoCov(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal, bool withanc, string mulreg)	{

		cerr << "read\n";
		int Ncont = GetModel()->Ncont;
		int Nanc = GetModel()->Nanc;

		int N = withanc ? Ncont + Nanc : Ncont;

		MeanExpNormTree** tree = new MeanExpNormTree*[N];
		for (int k=0; k<N; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), k);
			}
			if (withanc)	{
				for (int k=0; k<Nanc; k++)	{
					tree[Ncont + k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), Ncont + k);
				}
			}

		}
		cerr << '\n';
		cerr << "normalise\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
		}
		if (withanc)	{
			for (int k=0; k<Nanc; k++)	{
				tree[Ncont + k]->Normalise();
				ostringstream s;
				s << GetName() << ".postmeananc" << k+1 << ".tre";
				ofstream os(s.str().c_str());
				tree[Ncont + k]->ToStream(os);
				cerr << "reconstructed variation in ancestral character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			}
		}

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->SetWithDepth(GetModel()->GetLengthTree());
			tree[k]->Tabulate(os);
		}

		if (withanc)	{
			for (int k=0; k<Nanc; k++)	{
				ostringstream s;
				s << GetName() << ".postmeananc" << k+1 << ".tab";
				ofstream os(s.str().c_str());
				tree[Ncont + k]->Tabulate(os);

			}
		}

		cerr << '\n';
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
    bool printmed = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;
	bool withanc = false;

	string mulreg = "";

	bool sim = false;
	bool wocov = false;

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
			else if (s == "-wocov")	{
				wocov = true;
			}
			else if (s == "+anc")	{
				withanc = true;
			}
			else if (s == "-anc")	{
				withanc = false;
			}
			else if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
                printmed = false;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if ((s == "+med") || (s == "+median"))	{
				printmed = true;
                printmean = false;
			}
			else if ((s == "-med") || (s == "-median"))	{
				printmed = false;
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
			else if (s == "-mulreg")	{
				i++;
				mulreg = argv[i];
			}
			else if (s == "-partial")	{
				i++;
				mulreg = argv[i];
			}
			else if (s == "-sim")	{
				sim = true;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
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

	AncestralCovarianceSample sample(name,burnin,every,until);

	if (wocov)	{
		sample.ReadWoCov(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,withanc,mulreg);
		exit(1);
	}

	if (sim)	{
		sample.ReadSim();
		exit(1);
	}

	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,withanc,mulreg);

}



