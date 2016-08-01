
#ifndef TIMELINEANCESTRAL
#define TIMELINEANCESTRAL

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "BranchProcess.h"
#include "ContinuousData.h"
#include "AncestralData.h"
#include "MeanExpTree.h"
#include "Normal.h"

#include "Jeffreys.h"
#include "Partition.h"
#include "MultiVarNormal.h"


#include "TimeLineMultiVariateTreeProcess.h"

class TimeLineAncestralCovarianceModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	ContinuousData* contdata;
	AncestralData* ancdata;
	TaxonSet* taxonset;

	int Ncont;
	int Nanc;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	Chronogram* fixtree;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;

	JeffreysIIDArray* TimeLineDiagArray;
	SigmaZero* TimeLineSigmaZero;
	Rvar<CovMatrix>* timelinesigma;

	TimeIntervals* timeintervals;
	TimeLine* timeline;

	MultiVariateTreeProcess* process;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool clamproot;

	bool meanexp;

	int df;


	public:

	ConjugateInverseWishart* GetConjugateInverseWishart() {
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(sigma);
		if (! tmp)	{
			cerr << "error : dynamic castt of conjugate inverse wishart : " << sigma << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ConjugateMultiVariateTreeProcess* GetConjugateMultiVariateTreeProcess() {
		ConjugateMultiVariateTreeProcess* tmp = dynamic_cast<ConjugateMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : in dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	TimeLineMultiVariateTreeProcess* GetTimeLineMultiVariateTreeProcess() {
		TimeLineMultiVariateTreeProcess* tmp = dynamic_cast<TimeLineMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : in dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	TimeLineAncestralCovarianceModel() {
	}

	TimeLineAncestralCovarianceModel(string treefile, string contdatafile, string ancdatafile, double priorsigma, int indf, bool inclampdiag, int contdatatype, int ancdatatype, bool inclamproot, bool inmeanexp, bool sample=true)	{

		df = indf;
		clampdiag = inclampdiag;
		clamproot = inclamproot;
		meanexp = inmeanexp;

		// get continuous data from file
		contdata = new FileContinuousData(contdatafile);
		Ncont = contdata->GetNsite();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// get ancestral data from file
		if (ancdatafile != "None")	{
			cerr << ancdatafile << "\n";
			ancdata = new FileAncestralData(tree,ancdatafile);
			Nanc = ancdata->GetNsite();
			cerr << "ok\n";
		}
		else	{
			ancdata = 0;
			Nanc = 0;
		}
		// taxonset = ancdata->GetTaxonSet();

		// check whether tree and data fit together
		// tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";
		cerr << '\n';

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);

		One = new Const<PosReal>(1);

		fixtree = new Chronogram(tree,One);

		double mindiag = 0.001;
		double maxdiag = 1000;

		TimeLineDiagArray = new JeffreysIIDArray(Ncont+Nanc,mindiag,maxdiag,Zero);
		TimeLineDiagArray->ClampAt(1.0);
		TimeLineSigmaZero = new SigmaZero(TimeLineDiagArray);
		timelinesigma = new DiagonalCovMatrix(TimeLineSigmaZero, Ncont+Nanc+df);
		timelinesigma->SetIdentity();

		timeintervals = new TimeIntervals(fixtree);
		timeline = new TimeLine(fixtree,timeintervals, timelinesigma);

		DiagArray = new JeffreysIIDArray(Ncont+Nanc,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		sigmaZero = new SigmaZero(DiagArray);
		if (clampdiag)	{
			sigma = new DiagonalCovMatrix(sigmaZero, Ncont+Nanc+df);
		}
		else	{
			sigma = new ConjugateInverseWishart(sigmaZero, Ncont+Nanc+df);
		}
		sigma->SetIdentity();

		if (clampdiag)	{
			process = new TimeLineMultiVariateTreeProcess(sigma,timeline,fixtree);
		}
		else	{
			cerr << "make conjugate process\n";
			process = new ConjugateTimeLineMultiVariateTreeProcess(GetConjugateInverseWishart(),timeline,fixtree);
		}

		if (contdata)	{
			cerr << "set and clamp continuous data\n";
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,i,i,contdatatype);
			}
			cerr << "ok\n";
		}
		// process->GetMultiNormal(GetTree()->GetRoot())->ClampAt(0.52,1);
		if (ancdata)	{
			cerr << "set and clamp ancestral data\n";
			process->SetAndClamp(ancdata,Ncont,ancdatatype);
			// process->SetAndClamp(ancdata,0,contdatatype);
			cerr << "ok\n";
		}

		if (clamproot)	{
			process->ClampRoot();
		}

		RootRegister(Zero);
		RootRegister(One);
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	~TimeLineAncestralCovarianceModel() {}

	Tree* GetTree() {return tree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	LengthTree* GetLengthTree() {return fixtree;}

	ContinuousData* GetContinuousData() {return contdata;}
	AncestralData* GetAncestralData() {if (! ancdata) {cerr << "get ancestral data : null\n"; exit(1);} return ancdata;}

	TimeLine* GetTimeLine() {return timeline;}
	TimeIntervals* GetTimeIntervals() {return timeintervals;}

	CovMatrix* GetCovMatrix() {return sigma;}

	double GetLogProb()	{
		double total = 0;

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += timelinesigma->GetLogProb();
		total += timeline->GetLogProb();
		total += process->GetLogProb();

		return total;
	}

	virtual void MakeScheduler()	{

		/*
		scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
		scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
		*/

		scheduler.Register(new SimpleMove(process,10),10,"multinormal");
		scheduler.Register(new SimpleMove(process,1),10,"multinormal");
		scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
		scheduler.Register(new SimpleMove(process,0.01),10,"multinormal");

		if (! clampdiag)	{
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,0),1,"conjugate sigma - process");
		}
		else	{
			scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
		}

		scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
		scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
		scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

		scheduler.Register(new SimpleMove(timelinesigma,10),100,"timelinesigma");
		scheduler.Register(new SimpleMove(timelinesigma,1),100,"timelinesigma");
		scheduler.Register(new SimpleMove(timelinesigma,0.1),100,"timelinesigma");
		scheduler.Register(new SimpleMove(timelinesigma,0.01),100,"timelinesigma");

		scheduler.Register(new SimpleMove(timeline,10),100,"timeline");
		scheduler.Register(new SimpleMove(timeline,1),100,"timeline");
		scheduler.Register(new SimpleMove(timeline,0.1),100,"timeline");
		scheduler.Register(new SimpleMove(timeline,0.01),100,"timeline");
	}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";
		DiagArray->Sample();

		sigma->Sample();
		timelinesigma->Sample();
		timeline->Sample();
		process->Sample();

		cerr << "ok\n";
	}

	double GetMeanTimeLine(int k)	{
		double mean = 0;
		for (int i=0; i<timeline->GetSize(); i++)	{
			mean += (*(timeline->GetVal(i)))[k];
		}
		mean /= timeline->GetSize();
		return mean;
	}

	double GetTimeLine(int i, int k)	{
		return (*(timeline->GetVal(i)))[k];
	}

	double GetVarTimeLine(int k)	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<timeline->GetSize(); i++)	{
			double tmp = (*(timeline->GetVal(i)))[k];
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= timeline->GetSize();
		var /= timeline->GetSize();
		var -= mean * mean;
		return var;
	}

	void TraceHeader(ostream& os)	{
		os << "#logprob";
		for (int k=0; k<Ncont+Nanc; k++)	{
			for (int l=k+1; l<Ncont+Nanc; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << "tmean_" << k;
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << "tvar_" << k;
		}
		os << "\tdim";
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << "root_" << k;
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			for (int l=k+1; l<Ncont+Nanc; l++)	{
				os << '\t' << "tsigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << "tsigma_" << k << '_' << k;
		}
		if (! DiagArray->GetVal(0)->isClamped())	{
			for (int k=0; k<Ncont+Nanc; k++)	{
				os << "\tdiag" << k;
			}
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogProb();
		for (int k=0; k<Ncont+Nanc; k++)	{
			for (int i=0; i<timeline->GetSize(); i++)	{
				os << '\t' << GetTimeLine(i,k);
			}
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			for (int l=k+1; l<Ncont+Nanc; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << GetMeanTimeLine(k);
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << GetVarTimeLine(k);
		}
		os << '\t' << GetMultiVariateProcess()->GetMultiNormal(GetTree()->GetRoot())->val();
		for (int k=0; k<Ncont+Nanc; k++)	{
			for (int l=k+1; l<Ncont+Nanc; l++)	{
				os << '\t' << (*timelinesigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+Nanc; k++)	{
			os << '\t' << (*timelinesigma)[k][k];
		}
		if (! DiagArray->GetVal(0)->isClamped())	{
			for (int k=0; k<Ncont+Nanc; k++)	{
				os << '\t' << DiagArray->GetVal(k)->val();
			}
		}
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *timelinesigma << '\n';
		os << *timeline << '\n';
		os << *process << '\n';
	}

	void FromStream(istream& is)	{
		is >> *DiagArray;
		is >> *sigma;
		is >> *timelinesigma;
		is >> *timeline;
		is >> *process;
	}
};

#endif
