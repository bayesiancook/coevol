
#include "Sample.h"
#include "TimeLineClockModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"


class TimeLineClockSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string ancdatafile;
	string calibfile;
	string timelinefile;
	int flextimeline;

	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool clampdiag;
	bool clamproot;
	bool clamptree;
	bool meanexp;

	double priorsigma;

	int conjpath;
	int contdatatype;
	int ancdatatype;

	bool normalise;
	int nrep;

	int nsplit;

	int df;

	string mix;

	public:

	string GetModelType() {return modeltype;}

	TimeLineClockModel* GetModel() {return (TimeLineClockModel*) model;}

	TimeLineClockSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> datafile >> treefile >> contdatafile >> ancdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> timelinefile >> flextimeline;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag;
		is >> conjpath;
		is >> contdatatype;
		is >> ancdatatype;
		is >> clamproot >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> mix;
		is >> nsplit;
		is >> clamptree;
		is >> suffstatfile;
		is >> rootfile;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "TIMELINE")	{
			model = new TimeLineClockModel(datafile,treefile,contdatafile,ancdatafile,timelinefile,flextimeline,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,conjpath,contdatatype,ancdatatype,clamproot,clamptree,meanexp,normalise,nrep,mix,nsplit,suffstatfile,rootfile,false);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		OpenChainFile();
		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext)	{

		int Ncont = GetModel()->GetNcont();
		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());
		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		int dim = GetModel()->GetCovMatrix()->GetDim();
		MeanCovMatrix*  mat = new MeanCovMatrix(dim);

		double meanalpha = 0;
		double varalpha = 0;
		double ppalpha = 0;

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			double tmp = GetModel()->GetAlpha();
			meanalpha += tmp;
			varalpha += tmp * tmp;
			if (tmp > 0) ppalpha++;

			GetModel()->GetChronogram()->specialUpdate();
			meanchrono->Add(GetModel()->GetChronogram());

			if (GetModel()->Split())	{
				GetModel()->GetSplitLengthTree()->specialUpdate();
			}

			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetProcess(), GetModel()->GetLengthTree(), k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			mat->Add(&m);
		}
		cerr << '\n';
		cerr << "normalise\n";

		meanalpha /= size;
		varalpha /= size;
		varalpha -= meanalpha * meanalpha;
		cerr << "global/local scaling factor : " << meanalpha << " +/- " << sqrt(varalpha) << '\n';
		cerr << "pp : " << ppalpha << '\n';
		cerr << '\n';

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());

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
			// cout << '\n';
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
			// cout << '\n';
			cout << *mat;
		}

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
			if (tex)	{
				MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,tree[k],xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext);
				textree->Draw((s.str() + ".tex").c_str());
			}
		}

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
			// tree[k]->TabulateDistribution(os);
			// tree[k]->TabulateDistribution(os,false);
			ostringstream s2;
			s2 << GetName() << ".postmean" << k << ".loo";
			ofstream os2(s2.str().c_str());
			tree[k]->LooTabulate(os2, GetModel()->GetContinuousData());
		}

		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

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

	string mulreg;

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

	TimeLineClockSample sample(name,burnin,every,until);
	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,mulreg,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext);

}



