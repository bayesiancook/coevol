
#include "Sample.h"
#include "TimeLineAncestralCovarianceModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"

class TimeLineAncestralCovarianceSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string contdatafile;
	string ancdatafile;
	bool clampdiag;
	bool clamproot;
	bool meanexp;
	int contdatatype;
	int ancdatatype;
	int df;

	string bounds;
	string mix;

	public:

	string GetModelType() {return modeltype;}

	TimeLineAncestralCovarianceModel* GetModel() {return (TimeLineAncestralCovarianceModel*) model;}

	TimeLineAncestralCovarianceSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> clampdiag;
		is >> contdatatype;
		is >> ancdatatype;
		is >> clamproot >> meanexp;
		is >> df;
		is >> priorsigma;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		// make a new model depending on the type obtained from the file
		if (modeltype == "TIMELINEANCESTRALCOV")	{
			model = new TimeLineAncestralCovarianceModel(treefile,contdatafile,ancdatafile,priorsigma,df,clampdiag,contdatatype,ancdatatype,clamproot,meanexp,true);
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

	void ReadTimeLine()	{

		int timesize = GetModel()->GetTimeLine()->GetSize();
		int K = GetModel()->GetCovMatrix()->GetDim();
		double* date = new double[timesize];
		for (int i=0; i<timesize; i++)	{
			date[i] = 1 - GetModel()->GetTimeIntervals()->GetDate(i);
		}
		double** mean = new double*[timesize];
		double** var = new double*[timesize];
		for (int i=0; i<timesize; i++)	{
			mean[i] = new double[K];
			var[i] = new double[K];
			for (int k=0; k<K; k++)	{
				mean[i][k] = 0;
				var[i][k] = 0;
			}
		}

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			for (int i=0; i<timesize; i++)	{
				for (int k=0; k<K; k++)	{
					double tmp = (*GetModel()->GetTimeLine()->GetVal(i))[k];
					double tmp2 = (*GetModel()->GetTimeLine()->GetVal(2))[k];
					tmp-=tmp2;
					mean[i][k] += tmp;
					var[i][k] += tmp * tmp;
				}
			}
		}
		cerr << '\n';

		ofstream os((name + ".timeline").c_str());
		for (int i=0; i<timesize; i++)	{
			os << date[i];
			for (int k=0; k<K; k++)	{
				mean[i][k] /= size;
				var[i][k] /= size;
				var[i][k] -= mean[i][k] * mean[i][k];
				os << '\t' << mean[i][k] << '\t' << sqrt(var[i][k]);
			}
			os << '\n';
		}
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, bool withanc, string mulreg)	{

		cerr << "read\n";
		int Ncont = GetModel()->Ncont;
		int Nanc = GetModel()->Nanc;

		int N = withanc ? Ncont + Nanc : Ncont;

		MeanExpNormTree** tree = new MeanExpNormTree*[N];
		for (int k=0; k<N; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		int dim = GetModel()->GetCovMatrix()->GetDim();
		MeanCovMatrix*  mat = new MeanCovMatrix(dim);

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
		}
		cerr << '\n';
		cerr << "normalise\n";

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());

		cerr << "matrix ok\n";

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

			// cout << "entries are in the following order:\n";
			// GetModel()->PrintEntries(cout, indexarray);

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
		}

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
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;
	bool withanc = false;

	bool timeline = false;

	string mulreg = "";

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
			else if (s == "-tl")	{
				timeline = true;
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
			else if (s == "-mulreg")	{
				i++;
				mulreg = argv[i];
			}
			else if (s == "-partial")	{
				i++;
				mulreg = argv[i];
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

	TimeLineAncestralCovarianceSample sample(name,burnin,every,until);

	if (ic)	{
		sample.ReadIC();
		exit(1);
	}
	if (timeline)	{
		sample.ReadTimeLine();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,withanc,mulreg);

}



