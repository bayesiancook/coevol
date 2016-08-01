
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "AutoRegressiveMultiVariateTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "BranchProcess.h"
#include "ContinuousData.h"
#include "AncestralData.h"
#include "MeanExpTree.h"
#include "Normal.h"

#include "Jeffreys.h"
#include "Partition.h"
#include "MultiVarNormal.h"

#include "FixedLengthTree.h"


class AncestralCovarianceModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	ContinuousData* contdata;
	AncestralData** ancdata;
	TaxonSet* taxonset;

	int Ncont;
	int Nancdata;
	int* Nanc;
	int TotNanc;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	FixedLengthTree* fixtree;
	// Chronogram* fixtree;


	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;
	Gamma* phi;
	IIDUniform* mean;
	MultiVarNormal* drift;

	MultiVariateTreeProcess* process;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool autoregressive;

	bool clamproot;

	bool meanexp;

	int df;

	bool withdrift;

	bool kalman;

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
			cerr << "error : dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	AncestralCovarianceModel() {
	}

	AncestralCovarianceModel(string treefile, string contdatafile, string ancdatafile, double priorsigma, int indf, bool inclampdiag, int inautoregressive, int contdatatype, int ancdatatype, bool inclamproot, bool inmeanexp, bool inwithdrift, int inkalman, bool sample=true)	{

		df = indf;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		withdrift = inwithdrift;
		kalman = inkalman;

		// get continuous data from file
		contdata = new FileContinuousData(contdatafile);
		Ncont = contdata->GetNsite();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// get ancestral data from file
		if (ancdatafile != "None")	{
			cerr << ancdatafile << "\n";
			ifstream is(ancdatafile.c_str());
			if (! is)	{
				cerr << "error: cannot open ancdatafile\n";
				exit(1);
			}
			is >> Nancdata;
			Nanc = new int[Nancdata];
			ancdata = new AncestralData*[Nancdata];
			for (int a=0; a<Nancdata; a++)	{
				string name;
				is >> name;
				ancdata[a] = new FileAncestralData(tree,name);
				Nanc[a] = ancdata[a]->GetNsite();
				TotNanc += Nanc[a];
			}
			cerr << "ok\n";
		}
		else	{
			ancdata = 0;
			Nancdata = 0;
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

		// fixtree = new Chronogram(tree,One);
		fixtree = new FixedLengthTree(tree,One);


		double mindiag = 0.00001;
		double maxdiag = 100000;
		DiagArray = new JeffreysIIDArray(1 + Nancdata,mindiag,maxdiag,Zero);
		// DiagArray = new JeffreysIIDArray(Ncont+Nanc,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		Var<PosReal>** diagarray = new Var<PosReal>*[Ncont + TotNanc];
		int k = 0;
		for (int i=0; i<Ncont; i++)	{
			diagarray[k] = DiagArray->GetVal(0);
			k++;
		}
		for (int a=0; a<Nancdata; a++)	{
			for (int i=0; i<Nanc[a]; i++)	{
				diagarray[k] = DiagArray->GetVal(a+1);
				k++;
			}
		}
		if (k != TotNanc + Ncont)	{
			cerr << "error in diag array: non matching count\n";
			exit(1);
		}
		// sigmaZero = new SigmaZero(DiagArray);
		sigmaZero = new SigmaZero(diagarray, Ncont + TotNanc);
		if (clampdiag)	{
			sigma = new DiagonalCovMatrix(sigmaZero, Ncont+TotNanc+df);
		}
		else	{
			sigma = new ConjugateInverseWishart(sigmaZero, Ncont+TotNanc+df);
		}
		sigma->SetIdentity();

		// drift = new MultiVarNormal(Zero,priorOnSigmaZero);
		drift = new MultiVarNormal(Ncont+TotNanc,Zero,One);
		if (! withdrift)	{
			drift->ClampAtZero();
		}

		if (clampdiag)	{
			if (autoregressive)	{
				phi = new Gamma(One,One);
				mean = new IIDUniform(Zero,Ncont+TotNanc,100);
				process = new AutoRegressiveMultiVariateTreeProcess(sigma,mean,phi,fixtree);
			}
			else	{
				phi = 0;
				mean = 0;
				process = new MultiVariateTreeProcess(sigma,fixtree,0,drift);
			}
		}
		else	{
			if (autoregressive)	{
				phi = new Gamma(One,One);
				mean = new IIDUniform(Zero,Ncont+TotNanc,100);
				process = new ConjugateAutoRegressiveMultiVariateTreeProcess(GetConjugateInverseWishart(),mean,phi,fixtree);
			}
			else	{
				phi = 0;
				mean = 0;
				process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),fixtree,0,drift);
			}
		}

		if (contdata)	{
			cerr << "set and clamp continuous data\n";
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,i,i,contdatatype);
			}
			cerr << "ok\n";
		}
		// process->GetMultiNormal(GetTree()->GetRoot())->ClampAt(0.52,1);
		k = Ncont;
		for (int a=0; a<Nancdata; a++)	{
			cerr << "set and clamp ancestral data\n";
			process->SetAndClamp(ancdata[a],k,ancdatatype);
			k += Nanc[a];
			// process->SetAndClamp(ancdata,0,contdatatype);
			cerr << "ok\n";
		}

		if (clamproot)	{
			process->ClampRoot();
		}

		RootRegister(Zero);
		RootRegister(One);
		if (autoregressive)	{
			RootRegister(mean);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	~AncestralCovarianceModel() {}

	Tree* GetTree() {return tree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	LengthTree* GetLengthTree() {return fixtree;}

	ContinuousData* GetContinuousData() {return contdata;}
	AncestralData* GetAncestralData(int a) {
		return ancdata[a];
	}

	CovMatrix* GetCovMatrix() {return sigma;}

	void GetRootValues(double& x0, double& x1, double& x2)	{
		x0 = (*process->GetNodeVal(process->GetRoot()->GetNode()))[0];
		x1 = (*process->GetNodeVal(process->GetRoot()->Next()->Out()->GetNode()))[0];
		x2 = (*process->GetNodeVal(process->GetRoot()->Next()->Next()->Out()->GetNode()))[0];
	}

	double GetDrift(int i)	{
		return (*drift)[i];
	}

	double GetLogProb()	{
		double total = 0;

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += drift->GetLogProb();
		if (autoregressive)	{
			total += phi->GetLogProb();
			total += mean->GetLogProb();
		}
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

		if (kalman)	{
			scheduler.Register(new MultiVariateKalmanMove(process,Ncont,TotNanc),1,"kalman");
		}
		else	{
			/*
			scheduler.Register(new SimpleMove(process,10),10,"multinormal");
			scheduler.Register(new SimpleMove(process,1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),10,"multinormal");
			*/

			scheduler.Register(new MultiVariateSegmentMove(process,100,0,Ncont),10,"multinormal");
			scheduler.Register(new MultiVariateSegmentMove(process,10,0,Ncont),10,"multinormal");
			scheduler.Register(new MultiVariateSegmentMove(process,1,0,Ncont),10,"multinormal");
			scheduler.Register(new MultiVariateSegmentMove(process,0.1,0,Ncont),10,"multinormal");
			scheduler.Register(new MultiVariateSegmentMove(process,0.01,0,Ncont),10,"multinormal");
		}


		if (! clampdiag)	{
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,0),1,"conjugate sigma - process");
		}
		else	{
			scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
		}

		if (withdrift)	{
			scheduler.Register(new SimpleMove(drift,10),10,"drift");
			scheduler.Register(new SimpleMove(drift,1),10,"drift");
			scheduler.Register(new SimpleMove(drift,0.1),10,"drift");
			scheduler.Register(new SimpleMove(drift,0.01),10,"drift");
		}

		scheduler.Register(new SimpleMove(DiagArray,10),2,"theta");
		scheduler.Register(new SimpleMove(DiagArray,1),2,"theta");
		scheduler.Register(new SimpleMove(DiagArray,0.1),2,"theta");

		if (autoregressive)	{
			scheduler.Register(new SimpleMove(phi,10),100,"phi");
			scheduler.Register(new SimpleMove(phi,1),100,"phi");
			scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
			scheduler.Register(new SimpleMove(phi,0.01),100,"phi");

			scheduler.Register(new SimpleMove(mean,1),100,"mean");
			scheduler.Register(new SimpleMove(mean,0.1),100,"mean");
			scheduler.Register(new SimpleMove(mean,0.01),100,"mean");
		}
	}

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";
		DiagArray->Sample();

		sigma->Sample();
		drift->Sample();
		if (autoregressive)	{
			phi->Sample();
			mean->Sample();
		}

		process->Sample();

		cerr << "ok\n";
	}

	Var<RealVector>* GetMeanVector()	{
		return mean;
	}

	Var<RealVector>* GetDrift()	{
		return drift;
	}

	double GetCorrelCoeff()	{
		return (*sigma)[0][1] / sqrt((*sigma)[0][0] * (*sigma)[1][1]);
	}

	void TraceHeader(ostream& os)	{
		os << "#logprob";
		if (autoregressive)	{
			os << '\t' << "theta";
			os << '\t' << "dim";
			for (int k=0; k<Ncont+TotNanc; k++)	{
				os << '\t' << "mean" << k;
			}
		}
		else	{
			os << "\tdim";
			for (int k=0; k<Ncont+TotNanc; k++)	{
				os << '\t' << "root_" << k;
			}
		}
		for (int k=0; k<Ncont; k++)	{
			for (int l=k+1; l<Ncont+TotNanc; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		if (withdrift)	{
			os << "\tdim";
			for (int k=0; k<Ncont+TotNanc; k++)	{
				os << '\t' << "drift_" << k;
			}
		}
		if (! DiagArray->GetVal(0)->isClamped())	{
			for (int k=0; k<Nancdata+1; k++)	{
				os << "\tdiag" << k;
			}
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogProb();
		if (autoregressive)	{
			os << '\t' << phi->val();
			os << '\t' << *mean;
		}
		else	{
			os << '\t' << GetMultiVariateProcess()->GetMultiNormal(GetTree()->GetRoot())->val();
		}
		for (int k=0; k<Ncont; k++)	{
			for (int l=k+1; l<Ncont+TotNanc; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		if (withdrift)	{
			os << '\t' << *drift;
		}
		if (! DiagArray->GetVal(0)->isClamped())	{
			for (int k=0; k<Nancdata; k++)	{
				os << '\t' << DiagArray->GetVal(k)->val();
			}
		}
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *drift << '\n';
		if (autoregressive)	{
			os << *phi << '\n';
			os << *mean << '\n';
		}
		os << '\n';
		os << *process << '\n';
	}

	void FromStream(istream& is)	{
		is >> *DiagArray;
		is >> *sigma;
		is >> *drift;
		if (autoregressive)	{
			is >> *phi;
			is >> *mean;
		}
		is >> *process;
	}
};

#endif
