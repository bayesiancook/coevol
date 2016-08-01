
#include "Sample.h"
// #include "BGCModelwoLinReg.h"
#include "BGCModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"


class BGCSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;
	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool clampdiag;
	bool autoregressive;
	bool clamptree;
	bool meanexp;

	double priorsigma;
	double fixalpha;

	int conjpath;
	int contdatatype;

	bool normalise;
	int nrep;

	int df;

	int nsplit;

	double lambda;
	double at2cg;
	double at2gc;
	double gc2cg;

	public:

	string GetModelType() {return modeltype;}

	BGCModel* GetModel() {return (BGCModel*) model;}

	BGCSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		suffstatfile = "None";
		rootfile = "None";
		autoregressive = false;
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

		priorsigma = 1;
		fixalpha = 0;

		nsplit = 1;

		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag;
		is >> conjpath;
		is >> contdatatype;
		is >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> nsplit;
		is >> clamptree;
		is >> suffstatfile;
		is >> autoregressive;
		is >> rootfile;

		int check;
		is >> check;
		if (check)	{
			is >> fixalpha;
			is >> check;
			if (check)	{
				is >> lambda;
				is >> at2cg;
				is >> at2gc;
				is >> gc2cg;
				is >> check;
				if (check)	{
					cerr << "error when reading model\n";
					exit(1);
				}
			}
		}

		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "BGC")	{
			model = new BGCModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,lambda,at2cg,at2gc,gc2cg,false);
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

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext)	{

		cerr << "df : " << GetModel()->df << '\n';
		int Ncont = GetModel()->GetContMatrix()->GetDim();
		MeanCovMatrix*  mat = new MeanCovMatrix(Ncont);
		MeanExpNormTree* meanbgc = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal,0,0);
		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		double* meanreg = new double[Ncont];
		double* varreg = new double[Ncont];
		double* ppreg = new double[Ncont];
		for (int i=0; i<Ncont; i++)	{
			meanreg[i] = 0;
			varreg[i] = 0;
			ppreg[i] = 0;
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			GetModel()->GetBGCProcess()->specialUpdate();
			meanbgc->Add(GetModel()->GetBGCProcess(), GetModel()->GetLengthTree(), false);

			CovMatrix& m = *(GetModel()->GetContMatrix());
			mat->Add(&m);

			for (int k=0; k<Ncont; k++)	{
				double tmp = GetModel()->GetRegCoef(k);
				meanreg[k] += tmp;
				varreg[k] += tmp * tmp;
				if (tmp > 0)	{
					ppreg[k]++;
				}
			}

		}
		cerr << '\n';
		cerr << "normalise\n";

		meanbgc->Normalise();
		ofstream bgcos((GetName() + ".postmeanbgc.tre").c_str());
		meanbgc->ToStream(bgcos);

		for (int k=0; k<Ncont; k++)	{
			meanreg[k] /= size;
			varreg[k] /= size;
			ppreg[k] /= size;
			varreg[k] -= meanreg[k] * meanreg[k];
			cerr << '\n';
			cerr << "mean reg: " << k << '\t' << meanreg[k] << " +/- " << sqrt(varreg[k]) << '\t' << ppreg[k] << '\n';
		}

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());
		cout << *mat;
		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	string truefile = "";

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	bool tex = false;
	double nodescale = 5;
	double nodepower = 1;
	double barwidth = 0.04;
	double x = 6;
	double y = 8;
	bool withnodebubbles = true;
	bool withtimeintervals = true;
	int fontsize = 4;

	bool bubbletext = false;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-tex")	{
				tex = true;
			}
			else if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-np")	{
				i++;
				nodepower = atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
			}
			else if (s == "+bt")	{
				bubbletext = true;
			}
			else if (s == "-bt")	{
				bubbletext = false;
			}
			else if (s == "-fs")	{
				i++;
				fontsize = atoi(argv[i]);
			}
			else if (s == "-xs")	{
				i++;
				x = atof(argv[i]);
			}
			else if (s == "-ys")	{
				i++;
				y = atof(argv[i]);
			}
			else if (s == "+nb")	{
				withnodebubbles = true;
			}
			else if (s == "-nb")	{
				withnodebubbles = false;
			}
			else if (s == "+ti")	{
				withtimeintervals = true;
			}
			else if (s == "-ti")	{
				withtimeintervals = false;
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

	BGCSample sample(name,burnin,every,until);

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext);

}



