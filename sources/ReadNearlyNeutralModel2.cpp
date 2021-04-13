
#include "Sample.h"
#include "NearlyNeutralModel2.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"

class BranchOmegaMultivariateSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	
	bool sameseq;
	bool noadapt;
	bool clamptree;
	bool meanexp;
    bool shiftages;
    double beta0;
	GeneticCodeType type;

	int contdatatype;

	int nrep;

	int df;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaMultivariateModel* GetModel() {return (BranchOmegaMultivariateModel*) model;}

	BranchOmegaMultivariateSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil) {	
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
        chronoprior = 0;

		is >> modeltype;
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> contdatatype;
		is >> sameseq;
		is >> noadapt;
		is >> meanexp;
		is >> nrep;
		is >> df;
		is >> priorsigma;
		is >> clamptree;

		int check;
		is >> check;
		if (check)	{
			is >> chronoprior;
			is >> check;
            shiftages = 0;
			if (check)	{
                is >> shiftages;
                is >> check;
                if (check)  {
                    is >> beta0;
                    is >> check;
                    if (check)  {
                        cerr << "error when reading model\n";
                        exit(1);
                    }
                }
			}
		}

		is >> chainevery >> chainuntil >> chainsize;

		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,priorsigma,df,contdatatype,sameseq,noadapt,shiftages,beta0,clamptree,meanexp,nrep,false,type);
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

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext, double meanreg, double stdevreg)	{

		int Ncont = GetModel()->Ncont;
		int dim = GetModel()->GetCovMatrix()->GetDim();

        vector<double> betadist;
        double meanbeta = 0;

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());
		MeanExpNormTree* meanadaptativeomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal,meanreg,stdevreg);
		MeanExpNormTree* meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanu = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanNe = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);	
		MeanExpNormTree* meanneutralomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		
        meanNe->SetLogScale(10.0);
        meanu->SetLogScale(10.0);
        meansynrate->SetLogScale(10.0);

		double alpha[dim];
		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		MeanCovMatrix*  mat = new MeanCovMatrix(dim);

        int nerr = 0;

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

            GetModel()->UpdateLengthTree();
            nerr += GetModel()->GetShiftNerr();

			GetModel()->GetSynRateTree()->specialUpdate();
			GetModel()->GetNeutralOmegaNodeTree()->specialUpdate();
			GetModel()->GetSynrateNodeTree()->specialUpdate();
			if (!noadapt) {GetModel()->GetOmegaNodeTree()->specialUpdate();}
			GetModel()->GetNeNodeTree()->specialUpdate();
			GetModel()->GetUNodeTree()->specialUpdate();
			
			double rootage = GetModel()->GetRootAge();

            betadist.push_back(GetModel()->GetBeta());
            meanbeta += GetModel()->GetBeta();

			meanchrono->Add(GetModel()->GetChronogram());

			if (!noadapt) {meanadaptativeomega->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(),0);}

            // dS: mutation rate per tree depth
            // tree depth in years: rootage * 365.10^6
            // mutation rate per year: r = dS / rootage / 10^6
            double synrate_offset = -log(rootage * 1000000);
			meansynrate->Add(GetModel()->GetSynrateNodeTree(), GetModel()->GetLengthTree(), synrate_offset);

			if (!noadapt) {meanomega->Add(GetModel()->GetOmegaNodeTree(), GetModel()->GetLengthTree());}
			meanu->Add(GetModel()->GetUNodeTree(), GetModel()->GetLengthTree());
			meanNe->Add(GetModel()->GetNeNodeTree(), GetModel()->GetLengthTree());
			meanneutralomega->Add(GetModel()->GetNeutralOmegaNodeTree(), GetModel()->GetLengthTree());

			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetL()+k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			mat->Add(&m);
		}
						
		cerr << '\n';
        cerr << "number of shifting errors: " << nerr << '\t' << ((double) nerr)  / size << '\n';
		cerr << '\n';

		cerr << "normalise\n";

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());

		cout << *mat;

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';
		
        meanbeta /= size;
        sort(betadist.begin(), betadist.end());
        ofstream bos((GetName() + ".beta").c_str());
        int l = (int) (0.025 * size);
        int m = (int) (0.5 * size);
        bos << "median\tCI95min\tCI95max\n";
        bos << betadist[m] << '\t' << betadist[l] << '\t' << betadist[size-1-l] << '\n';
		cerr << "shape parameter beta in " << name << ".beta\n";
		cerr << '\n';

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);

		if (!noadapt) {
			meanadaptativeomega->Normalise();
			ofstream aoos((GetName() + ".postmeanadaptativeomega.tre").c_str());
			meanadaptativeomega->ToStream(aoos);
			cerr << "reconstructed variations of adaptative omega in " << name << ".postmeanadaptativeomega.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanadaptativeomega->GetPPLeafRoot() << '\n';
		}
		
		meanu->Normalise();
		ofstream mos((GetName() + ".postmeanu.tre").c_str());
		meanu->ToStream(mos);
		cerr << "reconstructed variations of u in " << name << ".postmeanu.tre\n";
		cerr << "pp of mean leaf values > root value : " << meanu->GetPPLeafRoot() << '\n';
			

		if (!noadapt) {
			meanomega->Normalise();
			ofstream oos((GetName() + ".postmeanomega.tre").c_str());
			meanomega->ToStream(oos);
			cerr << "reconstructed variations of omega in " << name << ".postmeanomega.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanomega->GetPPLeafRoot() << '\n';
		}

		
		meanNe->Normalise();
		ofstream Neos((GetName() + ".postmeanNe.tre").c_str());
		meanNe->ToStream(Neos);
		cerr << "reconstructed variations of Ne in " << name << ".postmeanNe.tre\n";
		cerr << "pp of mean leaf values > root value : " << meanNe->GetPPLeafRoot() << '\n';
		
		
		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
		cerr << "reconstructed variations of Ks in " << name << ".postmeansynrate.tre\n";
		cerr << "pp of mean leaf values > root value : " << meansynrate->GetPPLeafRoot() << '\n';


		meanneutralomega->Normalise();
		ofstream aos((GetName() + ".postmeanneutral_omega.tre").c_str());
		meanneutralomega->ToStream(aos);
		cerr << "reconstructed variations of neutral_omega in " << name << ".postmeanneutral_omega.tre\n";
		cerr << "pp of mean leaf values > root value : " << meanneutralomega->GetPPLeafRoot() << '\n';

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variations of continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
		}

		if (!noadapt) {
			ofstream aoaoos((GetName() + ".postmeanadaptativeomega.tab").c_str());
			meanadaptativeomega->Tabulate(aoaoos);
			aoaoos.close();
		}
	
		ofstream mmos((GetName() + ".postmeanu.tab").c_str());
		meanu->Tabulate(mmos);
		mmos.close();
		
		
		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);
		ssos.close();

		if (!noadapt) {
			ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
			meanomega->Tabulate(ooos);
			ooos.close();
		}
		
		ofstream Neoos((GetName() + ".postmeanNe.tab").c_str());
		meanNe->Tabulate(Neoos);
		Neoos.close();
			

		ofstream aaos((GetName() + ".postmeanneutral_omega.tab").c_str());
		meanneutralomega->Tabulate(aaos);
		aaos.close();
		
		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}

		cerr << '\n';
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
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	string mulreg = "";

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
			else if (s == "+med")	{
				printmean = true;
			}
			else if (s == "-med")	{
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

	BranchOmegaMultivariateSample sample(name,burnin,every,until);

	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,mulreg,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext,meanreg,stdevreg);

}



