
#include "Sample.h"
#include "BrownianModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"

class BrownianSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string datafile;
		string contdatafile;
	int conjpath;
		bool conjsigma;
	bool mapSegment;
		int nSegments;
		Segmentation segm;
		bool commut;
		bool gc;
		bool omega;
		bool fixtimes;
		bool paral;



	public:

	string GetModelType() {return modeltype;}

	BrownianModel* GetModel() {return (BrownianModel*) model;}

	BrownianSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> treefile >> datafile >> contdatafile;
		is >> conjpath;
				is >> mapSegment;
		is >> nSegments;
				bool b;
				is >> b;
				segm = b ? SEGM_REGULAR : SEGM_ABSOLUTE;
				is >> commut;
				is >> gc;
				is >> omega;
				is >> conjsigma;
				is >> fixtimes;
				is >> paral;
		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

				is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "BrownianNonCommutative")	{
			model = new BrownianModel(datafile,treefile, contdatafile, conjpath, mapSegment, nSegments, segm, commut, gc, omega, conjsigma, fixtimes, paral, false);
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

	void PostPred(int nrep, int resamplecov, int resampleprocess, int nsegments, string newname)	{

		if (newname == "")	{
			newname = name;
		}
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();

			for (int rep=0; rep<nrep; rep++)	{
				ostringstream s;
				if (nrep > 1)	{
					s << newname << "_" << i << "_" << rep;
				}
				else	{
					s << newname << "_" << i ;
				}
				GetModel()->PostPredSimu(s.str(),resamplecov,resampleprocess,nsegments);
			}
		}
		cerr << '\n';
	}

	void Read(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext, double meanreg, double stdevreg)	{

		int Ncont = GetModel()->GetNcont();

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());
		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,meanreg,stdevreg);
		MeanExpNormTree* meangc = 0;
				if(gc)
					meangc = new MeanExpNormTree(GetModel()->GetTree(),true,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		int dim = GetModel()->GetSigma()->GetDim();

		MeanCovMatrix*  mat = new MeanCovMatrix(dim);

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			//GetModel()->GetSynRateTree()->specialUpdate();
			GetModel()->GetChronogram()->specialUpdate();

			meanchrono->Add(GetModel()->GetChronogram());

			meansynrate->Add(GetModel()->GetBrownianProcess()->GetInstantProcess(), GetModel()->GetChronogram(), 0);
			if(gc)
							meangc->Add(GetModel()->GetBrownianProcess()->GetInstantProcess(), GetModel()->GetChronogram(), 1);


			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetBrownianProcess()->GetInstantProcess(), GetModel()->GetChronogram(), 2+k);
			}

			CovMatrix& m = *(GetModel()->GetSigma());
			mat->Add(&m);
		}
		cerr << '\n';
		cerr << "normalise\n";

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

			cout << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout);

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
			cout << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout);

			cout << '\n';
			cout << *mat;
		}

		cout.close();


		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);


		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
				if(gc) {
					meangc->Normalise();
					ofstream gcos((GetName() + ".postmeangc.tre").c_str());
					meangc->ToStream(gcos);
				}
		cerr << "reconstructed variations of Ks in " << name << ".postmeansynrate.tre\n";
		cerr << "pp of mean leaf values > root value : " << meansynrate->GetPPLeafRoot() << '\n';

		if (tex)	{
			ostringstream s;
			s << GetName() << ".postmeansynrate.tre";
			MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,meansynrate,xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext);
			textree->Draw((s.str() + ".tex").c_str());
		}

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variations of continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
			if (tex)	{
				MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,tree[k],xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext);
				textree->Draw((s.str() + ".tex").c_str());
			}
		}

		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);
		ssos.close();

		if(gc) {
			ofstream gcsos((GetName() + ".postmeangc.tab").c_str());
			meangc->Tabulate(gcsos);
			gcsos.close();
		}


		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}

		cerr << '\n';
	}

	void ReadShort(string outfile)	{

		GetNextPoint();
		ofstream os(outfile.c_str());
		GetModel()->ToStreamShort(os);
		os.close();
	}

	void CheckCov(string truefile)	{

		int dim = GetModel()->GetSigma()->GetDim();

		double trueval[dim][dim];
		double mean[dim][dim];
		double var[dim][dim];
		double meancorrel[dim][dim];
		double varcorrel[dim][dim];
		double truecorrel[dim][dim];
		double totval[dim][dim][size];
		double totcorrel[dim][dim][size];

		for (int j=0; j<dim; j++)	{
			for (int k=0; k<dim; k++)	{
				mean[j][k] = 0;
				var[j][k] = 0;
				meancorrel[j][k] = 0;
				varcorrel[j][k] = 0;
			}
		}

		cerr << "TRUE  " << truefile << '\n';
		ifstream is(truefile.c_str());
		// ifstream is((truefile + "cov").c_str());
		if (!is)	{
			cerr << "error: did not find file : " << truefile << '\n';
			exit(1);
		}

		/*
		int tmpdim;
		is >> tmpdim;
		if (tmpdim != dim)	{
			cerr << "error when reading true file\n";
			cerr << tmpdim << '\t' << dim << '\n';
			exit(1);
		}
		*/
		for (int j=0; j<dim; j++)	{
			for (int k=0; k<dim; k++)	{
				is >> trueval[j][k];
			}
		}
		for (int j=0; j<dim; j++)	{
			for (int k=j+1; k<dim; k++)	{
				truecorrel[j][k] = trueval[j][k] / sqrt(trueval[j][j] * trueval[k][k]);
			}
		}
		for (int i=0; i<size; i++)	{
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();
			for (int j=0; j<dim; j++)	{
				for (int k=j; k<dim; k++)	{
					double tmp = (*GetModel()->GetSigma())[j][k];
					mean[j][k] += tmp;
					var[j][k] += tmp * tmp;
					totval[j][k][i] = tmp;
				}
			}
			for (int j=0; j<dim; j++)	{
				for (int k=j+1; k<dim; k++)	{
					double tmp = (*GetModel()->GetSigma())[j][k] / sqrt( (*GetModel()->GetSigma())[j][j] * (*GetModel()->GetSigma())[k][k]);
					meancorrel[j][k] += tmp;
					totcorrel[j][k][i] = tmp;
				}
			}
		}

		ofstream vos((name + ".varcheck").c_str());
		for (int j=0; j<dim; j++)	{
			// for (int k=j; k<dim; k++)	{
				int k = j;
				mean[j][k] /= size;
				var[j][k] /= size;
				var[j][k] -= mean[j][k] * mean[j][k];
				int n = 0;
				for (int i=0; i<size; i++)	{
					if (trueval[j][k] > totval[j][k][i])	{
						n++;
					}
				}

				for (int i=0; i<size; i++)	{
					for (int l=size-1; l>i; l--)	{
						if (totval[j][k][i] > totval[j][k][l])	{
							double tmp = totval[j][k][i];
							totval[j][k][i] = totval[j][k][l];
							totval[j][k][l] = tmp;
						}
					}
				}
				int min = ((int) (((double) (size)) / 100 * 2.5));
				vos << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t' << var[j][k] << '\t' << ((double) n)/size << '\n';
				// os << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t';
			// }
		}
		// os << '\n';

		ofstream cos((name + ".covcheck").c_str());
		for (int j=0; j<dim; j++)	{
			for (int k=j+1; k<dim; k++)	{
				mean[j][k] /= size;
				var[j][k] /= size;
				var[j][k] -= mean[j][k] * mean[j][k];

				meancorrel[j][k] /= size;
				varcorrel[j][k] /= size;
				varcorrel[j][k] -= meancorrel[j][k] * meancorrel[j][k];

				int n = 0;
				for (int i=0; i<size; i++)	{
					if (trueval[j][k] > totval[j][k][i])	{
						n++;
					}
				}

				for (int i=0; i<size; i++)	{
					for (int l=size-1; l>i; l--)	{
						if (totval[j][k][i] > totval[j][k][l])	{
							double tmp = totval[j][k][i];
							totval[j][k][i] = totval[j][k][l];
							totval[j][k][l] = tmp;
						}
					}
				}

				for (int i=0; i<size; i++)	{
					for (int l=size-1; l>i; l--)	{
						if (totcorrel[j][k][i] > totcorrel[j][k][l])	{
							double tmp = totcorrel[j][k][i];
							totcorrel[j][k][i] = totcorrel[j][k][l];
							totcorrel[j][k][l] = tmp;
						}
					}
				}
				int min = ((int) (((double) (size)) / 100 * 2.5));
				cos << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t' << sqrt(var[j][k]) << '\t' << ((double) n)/size << '\t' << truecorrel[j][k] << '\t' << meancorrel[j][k] << '\t' << totcorrel[j][k][min] << '\t' << totcorrel[j][k][size-min-1] << '\t' << sqrt(varcorrel[j][k]) << '\n';
				// cos << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t' << var[j][k] << '\t' << ((double) n)/size << '\t' << truecorrel[j][k] << '\t' << meancorrel[j][k] << '\t' << totcorrel[j][k][min] << '\t' << totcorrel[j][k][size-min-1] << '\t' << varcorrel[j][k] << '\n';
				// os << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t';
			}
		}
		// os << '\n';

	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	bool ic = false;

	bool printlog = false;
	bool printmean = false;
    bool printmed = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	string mulreg = "";

	bool check = false;
	string truefile = "";

	bool tex = false;
	double nodescale = 5;
	double nodepower = 1;
	double barwidth = 0.04;
	double x = 6;
	double y = 8;
	// bool withnodebubbles = true;
	// bool withtimeintervals = true;
	int fontsize = 4;

	bool bubbletext = false;

	double meanreg = 0;
	double stdevreg = 0;

	int pp = 0;
	int ppcov = 0;
	int ppbrown = 0;
	int nrep = 1;
	int nsegments = -1;
	string ppname = "";

	string shortfile = "None";

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
				nodepower= atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
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
			else if (s == "+bt")	{
				bubbletext = true;
			}
			else if (s == "-bt")	{
				bubbletext = false;
			}
			else if (s == "-pointwise")	{
				i++;
				meanreg = atof(argv[i]);
				i++;
				stdevreg = atof(argv[i]);
			}
			/*
			else if (s == "+nb")	{
				withnodebubbles = true;
			}
			else if (s == "-nb")	{
				withnodebubbles = false;
			}
			*/
			else if (s == "-ic")	{
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
			else if ((s == "-pp") || (s == "-ppred"))	{
				pp = 1;
			}
			else if (s == "+ppcov")	{
				pp = 1;
				ppcov = 1;
			}
			else if (s == "-ppcov")	{
				pp = 1;
				ppcov = 0;
			}
			else if (s == "+ppbrown")	{
				pp = 1;
				ppbrown = 1;
			}
			else if (s == "-ppbrown")	{
				pp = 1;
				ppbrown = 0;
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-ns")	{
				i++;
				nsegments = atoi(argv[i]);
			}
			else if (s == "-ppname")	{
				i++;
				ppname = argv[i];
			}
			else if (s == "-ch")	{
				check = true;
				i++;
				truefile = argv[i];
			}
			else if (s == "-short")	{
				i++;
				shortfile = argv[i];
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
		cerr << "readbrownian [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	BrownianSample sample(name,burnin,every,until);

	if (pp)	{
		sample.PostPred(nrep,ppcov,ppbrown,nsegments,ppname);
	}
	else if (shortfile != "None")	{
		sample.ReadShort(shortfile);
	}
	else if (check)	{
		sample.CheckCov(truefile);
	}
	else	{
		sample.Read(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,mulreg,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext,meanreg,stdevreg);

	}

}




