
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "GeneralConjugatePath.h"
#include "Jeffreys.h"
#include "PrecisionNormalTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "MeanChronogram.h"
#include "SumConstrained.h"

class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(Var<Profile>* inrelrate, BranchVarTree<Profile>* instattree, Var<Profile>* inrootstat) {
		SetWithRoot(true);
		stattree = instattree;
		relrate = inrelrate;
		rootstat = inrootstat;
		RecursiveCreate(GetRoot());
	}

	~MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new GTRRandomSubMatrixWithNormRates(relrate,rootstat,false);
		}
		return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()),true);
		// return new GTRRandomSubMatrixWithNormRates(relrate,stattree->GetBranchVal(link->GetBranch()),false);
	}

	Tree* GetTree() {return stattree->GetTree();}

	private:

	BranchVarTree<Profile>* stattree;
	Var<Profile>* relrate;
	Var<Profile>* rootstat;

};

class CombiCompoModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment** data;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	int Ndata;
	// number of columns
	int* Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int* Nstate;
	int Ncont;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	int chronoprior;
	double softa;

	double meanchi;
	double meanchi2;
	// 0 : uniform;
	// 1 : bd;
	// 2 : bd with cauchy proper lower bounds

	// chronogram
	Chronogram* chronogram;
	LengthTree** synratetree;
	MeanExpTreeFromMultiVariate** meanexpsynratetree;
	Jeffreys** synsigma;
	LogNormalTreeProcess** lognormalsyntree;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* sigma;
	ConjugateMultiVariateTreeProcess* process;

	JeffreysIIDArray* contDiagArray;
	SigmaZero* contSigmaZero;
	ConjugateInverseWishart* contsigma;
	ConjugateMultiVariateTreeProcess* contprocess;

	/*
	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;
	*/

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet** relrate;
	Dirichlet** rootstationary;

	SumConstrainedMapping** changeofbasis;
	SumConstrainedStatTree** stattree;

	MatrixTree** matrixtree;

	// phylo process
	PathConjugateTree** pathconjtree;
	PhyloProcess** phyloprocess;

	// if true: covariances are all set equal to 0
	int clampdiag;
	int clamptree;
	int meanexp;

	// total number of substitution parameters modelled as non homogeneous
	int nrep;

	int L;
	int M;
	int df;
	int jeffdim;

	int* Ksyn;
	int* Kcomp;
	int Kcont;

	bool iscalib;

	public:

	CombiCompoModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double insofta,  double inmeanchi, double inmeanchi2, int indf, int inclampdiag, int contdatatype, string rrtype, int inclamptree, int inmeanexp, int innormalise, int innrep, string inmix, string rootfile, int inseparatesyn, bool sample=true)	{

		contprocess = 0;

		clamptree = inclamptree;
		meanexp = inmeanexp;
		clampdiag = inclampdiag;

		chronoprior = inchronoprior;
		softa = insofta;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		cerr << "ali\n";
		// get data from file
		ifstream is(datafile.c_str());
		is >> Ndata;
		Nsite = new int[Ndata];
		Nstate = new int[Ndata];

		data = new SequenceAlignment*[Ndata];

		for (int d=0; d<Ndata; d++)	{
			string tmp;
			is >> tmp;
			data[d] = new FileSequenceAlignment(tmp);
			if (!d)	{
				taxonset = GetData(0)->GetTaxonSet();
			}
			else	{
				GetData(d)->RegisterWith(taxonset);
			}
			Nsite[d] = GetData(d)->GetNsite();
			Nstate[d] = GetData(d)->GetNstate();
		}

		nrep = innrep;
		if (nrep == 0)	{
			nrep = 30;
		}

		cerr << "tree\n";
		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		cerr << "cont data\n";
		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

		cerr << "tree and data ok\n";
		cerr << '\n';

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		double mindiag = 0.001;
		double maxdiag = 1000;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		iscalib = false;

		lognormalsyntree = 0;

		if (calibfile != "None")	{
			iscalib = true;

			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			if (rootage == -1)	{
				a = b = -1;
			}
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
			}
			else {
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,a,b,calibset,chronoprior,softa);
			}
		}
		else	{
			chronogram = new Chronogram(tree,One);
		}
		if (clamptree)	{
			chronogram->Clamp();
		}

		synratetree = new LengthTree*[Ndata];
		synsigma = 0;

		if (inseparatesyn)	{
			synsigma = new Jeffreys*[Ndata];
			lognormalsyntree = new LogNormalTreeProcess*[Ndata];
			for (int d=0; d<Ndata; d++)	{
				synsigma[d] = new Jeffreys(mindiag,maxdiag,Zero);
				synsigma[d]->setval(1);
				lognormalsyntree[d] = new LogNormalTreeProcess(chronogram,synsigma[d],INTEGRAL);
				synratetree[d] = lognormalsyntree[d];
			}
		}

		Ksyn = new int[Ndata];
		Kcomp = new int[Ndata];

		if (! inseparatesyn)	{
			int index = 0;
			for (int d=0; d<Ndata; d++)	{
				Ksyn[d] = index;
				index++;
				Kcomp[d] = index;
				index += Nstate[d] - 1;
			}
			L = index;
			jeffdim = 2 * Ndata;
		}
		else	{
			int index = 0;
			for (int d=0; d<Ndata; d++)	{
				Ksyn[d] = -1;
				Kcomp[d] = index;
				index += Nstate[d] - 1;
			}
			L = index;
			jeffdim = Ndata;
		}

		if (clampdiag)	{
			M = L;
			Kcont = -1;
		}
		else	{
			M = L + Ncont;
			jeffdim += Ncont;
			Kcont = L;
		}

		df = M + indf;

		// create an array of positive variables kappa_i, i=1..Ncont + L
		DiagArray = new JeffreysIIDArray(jeffdim,mindiag,maxdiag,Zero);
		for (int i=0; i<jeffdim; i++)	{
			// DiagArray->GetVal(i)->ClampAt(1.0);
			DiagArray->GetVal(i)->setval(1.0 + 0.01 * (Random::Uniform() - 0.5));
		}
		Var<PosReal>** diagarray = new Var<PosReal>*[M];
		int index = 0;
		int k = 0;
		cerr << "===\n";
		for (int d=0; d<Ndata; d++)	{
			if (Ksyn[d] != -1)	{
				cerr << index << '\t' << k << '\n';
				diagarray[index] = DiagArray->GetVal(k);
				index++;
				k++;
			}
			for (int i=0; i<Nstate[d]-1; i++)	{
				cerr << index << '\t' << k << '\n';
				diagarray[index] = DiagArray->GetVal(k);
				index++;
			}
			k++;
		}
		if (! clampdiag)	{
			for (int i=0; i<Ncont; i++)	{
				cerr << index << '\t' << k << '\n';
				diagarray[index] = DiagArray->GetVal(k);
				index++;
				k++;
			}
		}
		cerr << "===\n";
		if (index != M)	{
			cerr << "error when constructing diag array: index\n";
			exit(1);
		}
		if (k != jeffdim)	{
			cerr << "error when constructing diag array: k\n";
			exit(1);
		}

		cerr << Ncont << '\t' << L << '\t' << M << '\n';
		cerr << jeffdim << '\n';
		cerr << Kcont << '\n';
		for (int d=0; d<Ndata; d++)	{
			cerr << Ksyn[d] << '\t' << Kcomp[d] << '\n';
		}

		sigmaZero = new SigmaZero(diagarray, M);
		sigma = new ConjugateInverseWishart(sigmaZero, df);
		// sigma->SetIdentity();
		cerr << "process\n";
		process = new ConjugateMultiVariateTreeProcess(sigma,chronogram);
		// process->Reset();

		if (clampdiag)	{
			contDiagArray = new JeffreysIIDArray(Ncont,mindiag,maxdiag,Zero);
			for (int i=0; i<Ncont; i++)	{
				contDiagArray->GetVal(i)->setval(1.0 + 0.01 * (Random::Uniform() - 0.5));
			}
			contSigmaZero = new SigmaZero(contDiagArray);
			contsigma = new ConjugateInverseWishart(contSigmaZero, Ncont + indf);
			contprocess = new ConjugateMultiVariateTreeProcess(contsigma,chronogram);

			// condition the multivariate process
			// on the matrix of quantitative traits.
			for (int i=0; i<Ncont; i++)	{
				contprocess->SetAndClamp(contdata,i,i,contdatatype);
				(*contprocess->GetRootVal())[i] = contdata->GetMeanLog(i);
			}
		}
		else	{
			contDiagArray = 0;
			contSigmaZero = 0;
			contsigma = 0;
			contprocess = 0;

			// condition the multivariate process
			// on the matrix of quantitative traits.
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,Kcont+i,i,contdatatype);
				(*process->GetRootVal())[Kcont+i] = contdata->GetMeanLog(i);
			}
		}


		// just for numerical stability of the starting point
		for (int i=0; i<L; i++)	{
			process->CutOff(1,i);
		}

		if (! inseparatesyn)	{
			meanexpsynratetree = new MeanExpTreeFromMultiVariate*[Ndata];
			for (int d=0; d<Ndata; d++)	{
				cerr << "syn : " << Ksyn[d] << '\n';
				meanexpsynratetree[d] = new MeanExpTreeFromMultiVariate(process,Ksyn[d],INTEGRAL,false,meanexp);
				synratetree[d] = meanexpsynratetree[d];
			}
		}

		cerr << "matrix\n";
		relrate = new Dirichlet*[Ndata];
		rootstationary = new Dirichlet*[Ndata];
		changeofbasis = new SumConstrainedMapping*[Ndata];
		matrixtree = new MatrixTree*[Ndata];
		stattree = new SumConstrainedStatTree*[Ndata];

		pathconjtree = new PathConjugateTree*[Ndata];
		phyloprocess = new PhyloProcess*[Ndata];

		for (int d=0; d<Ndata; d++)	{
			relrate[d] = new Dirichlet(Nstate[d]*(Nstate[d]-1)/2);
			if ((Nstate[d] == Naa) && (rrtype == "lg"))	{
				double total = 0;
				for (int i=0; i<Nstate[d]*(Nstate[d]-1)/2; i++)	{
					(*relrate[d])[i] = LG_RR[i];
					total += (*relrate[d])[i];
				}
				for (int i=0; i<Nstate[d]*(Nstate[d]-1)/2; i++)	{
					(*relrate[d])[i] /= total;
				}
				relrate[d]->Clamp();
			}

			rootstationary[d] = new Dirichlet(Nstate[d]);
			changeofbasis[d] = new SumConstrainedMapping(Nstate[d]);

			cerr << "statree : " << d << '\t' << Kcomp[d] << '\t' << Nstate[d] << '\n';
			stattree[d] = new SumConstrainedStatTree(process,Kcomp[d],changeofbasis[d]);
			matrixtree[d] = new MatrixTree(relrate[d],stattree[d],rootstationary[d]);

			pathconjtree[d] = new BranchMatrixPathConjugateTree(synratetree[d], matrixtree[d], data[d]);
			phyloprocess[d] = new PathConjugatePhyloProcess(pathconjtree[d]);

			cerr << "unfold\n";
			phyloprocess[d]->Unfold();
			cerr << "sample\n";
			if (sample)	{
				if (phyloprocess[d])	{
					phyloprocess[d]->Sample();
				}
			}

		}
		cerr << "register\n";
		// register model
		RootRegister(Zero);
		RootRegister(One);
		if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		for (int d=0; d<Ndata; d++)	{
			RootRegister(relrate[d]);
			RootRegister(rootstationary[d]);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~CombiCompoModel() {}

	bool SeparateSyn()	{
		return lognormalsyntree;
	}

	bool SeparateCont()	{
		return clampdiag;
	}

	bool isCalibrated()	{
		return iscalib;
	}

	// accessors
	Tree* GetTree() {return tree;}

	SequenceAlignment* GetData(int d)	{
		return data[d];
	}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}

	MeanExpTreeFromMultiVariate* GetSynRateTree(int d) {return meanexpsynratetree[d];}

	Chronogram* GetChronogram() {
		return chronogram;
	}
	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix() {return sigma;}

	// probability computation

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		if (chronoprior)	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();

		if (SeparateSyn())	{
			for (int d=0; d<Ndata; d++)	{
				total += synsigma[d]->GetLogProb();
				total += lognormalsyntree[d]->GetLogProb();
			}
		}

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		if (clampdiag)	{
			total += contDiagArray->GetLogProb();
			total += contsigma->GetLogProb();
			total += contprocess->GetLogProb();
		}

		for (int d=0; d<Ndata; d++)	{
			total += relrate[d]->GetLogProb();
			total += rootstationary[d]->GetLogProb();
		}

		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		for (int d=0; d<Ndata; d++)	{
			ret += phyloprocess[d]->GetLogProb();
		}
		// double ret = pathconjtree->GetLogProb();
		return ret;
	}

	// MCMC schedule
	virtual void MakeScheduler()	{

		for (int d=0; d<Ndata; d++)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess[d],pathconjtree[d]),1,"mapping + sufficient stat");
		}

		for (int i=0; i<nrep; i++)	{

			if (! clamptree)	{
				if (chronoprior)	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.001),10,"chrono");
			}

			if (isCalibrated())	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");
			}

			if (SeparateSyn())	{
				for (int d=0; d<Ndata; d++)	{
					scheduler.Register(new SimpleMove(synsigma[d],10),100,"syn sigma");
					scheduler.Register(new SimpleMove(synsigma[d],1),100,"syn sigma");
					scheduler.Register(new SimpleMove(synsigma[d],0.1),100,"syn sigma");
					scheduler.Register(new SimpleMove(synsigma[d],0.01),100,"syn sigma");
					scheduler.Register(new SimpleMove(lognormalsyntree[d],10),10,"log normal syn");
					scheduler.Register(new SimpleMove(lognormalsyntree[d],1),10,"log normal syn");
					scheduler.Register(new SimpleMove(lognormalsyntree[d],0.1),10,"log normal syn");
					scheduler.Register(new SimpleMove(lognormalsyntree[d],0.01),10,"log normal syn");
				}
			}

			scheduler.Register(new SimpleMove(process,1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.3),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.03),10,"multinormal");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,1,1),1,"conjugate sigma - process");

			if (clampdiag)	{
				scheduler.Register(new SimpleMove(contprocess,1),10,"cont multinormal");
				scheduler.Register(new SimpleMove(contprocess,0.3),10,"cont multinormal");
				scheduler.Register(new SimpleMove(contprocess,0.1),10,"cont multinormal");
				scheduler.Register(new SimpleMove(contprocess,0.03),10,"cont multinormal");
				scheduler.Register(new ConjugateMultiVariateMove(contsigma,contprocess,1,1),1,"cont conjugate sigma - process");
				for (int c=0; c<Ncont; c++)	{
					scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(contprocess,10,c,1),10,"cont conjugate sigma - process - piecewise translation");
					scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(contprocess,1,c,1),10,"cont conjugate sigma - process - piecewise translation");
					scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(contprocess,0.1,c,1),10,"cont conjugate sigma - process - piecewise translation");
				}
			}
			else	{
				for (int c=0; c<Ncont; c++)	{
					scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,10,Kcont+c,1),10,"conjugate sigma - process - piecewise translation");
					scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,1,Kcont+c,1),10,"conjugate sigma - process - piecewise translation");
					scheduler.Register(new MultiVariateWholeTreePiecewiseTranslationMove(process,0.1,Kcont+c,1),10,"conjugate sigma - process - piecewise translation");
				}
			}

			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

			if (clampdiag)	{
				scheduler.Register(new SimpleMove(contDiagArray,10),10,"theta");
				scheduler.Register(new SimpleMove(contDiagArray,1),10,"theta");
				scheduler.Register(new SimpleMove(contDiagArray,0.1),10,"theta");
			}

			for (int d=0; d<Ndata; d++)	{
				scheduler.Register(new ProfileMove(relrate[d],0.1,1),10,"relrates");
				scheduler.Register(new ProfileMove(relrate[d],0.03,2),10,"relrates");
				scheduler.Register(new SimpleMove(relrate[d],0.01),10,"relrates");

				scheduler.Register(new ProfileMove(rootstationary[d],0.1,2),10,"stat");
				scheduler.Register(new ProfileMove(rootstationary[d],0.03,2),10,"stat");
				scheduler.Register(new ProfileMove(rootstationary[d],0.03,5),10,"stat");
				scheduler.Register(new SimpleMove(rootstationary[d],0.01),10,"stat");
			}
		}
	}

	void drawSample()	{
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	// summary statistics

	double GetMeanSynRate(int d)	{
		if (SeparateSyn())	{
			return lognormalsyntree[d]->GetMeanRate();
		}
		return GetSynRateTree(d)->GetTotal();
	}

	Var<PosReal>* GetScale()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	double GetMeanEntropy(int d)	{
		return stattree[d]->GetMeanEntropy();
	}

	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		for (int d=0; d<Ndata; d++)	{
			os << "\trate" << d << "\tstatent" << d;
		}
		for (int k=0; k<Ncont; k++)	{
			os << "\trootcont" << k;
		}
		if (isCalibrated())	{
			os << "\trootage";
		}
		os << "\tlogdet";
		if (clampdiag)	{
			os << "\tcontlogdet";
		}
		/*
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		*/
		for (int d=0; d<Ndata; d++)	{
			os << "\tstatent" << d;
			os << "\trrent" << d;
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		for (int d=0; d<Ndata; d++)	{
			os << '\t' << GetMeanSynRate(d);
			os << '\t' << GetMeanEntropy(d);
		}
		if (clampdiag)	{
			for (int k=0; k<Ncont; k++)	{
				os << '\t' << (*contprocess->GetRootVal())[k];
			}
		}
		else	{
			for (int k=0; k<Ncont; k++)	{
				os << '\t' << (*process->GetRootVal())[Kcont + k];
			}
		}
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		os << '\t' << sigma->GetLogDeterminant();
		if (clampdiag)	{
			os << '\t' << contsigma->GetLogDeterminant();
		}
		/*
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		*/
		for (int d=0; d<Ndata; d++)	{
			os << '\t' << rootstationary[d]->val().GetEntropy();
			os << '\t' << relrate[d]->val().GetEntropy();
		}
		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		if (chronoprior)	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		if (clampdiag)	{
			os << *contDiagArray << '\n';
			os << *contsigma << '\n';
			os << *contprocess << '\n';
		}
		for (int d=0; d<Ndata; d++)	{
			if (SeparateSyn())	{
				os << *synsigma[d] << '\n';
				os << *lognormalsyntree[d] << '\n';
			}
			changeofbasis[d]->ToStream(os);
			os << *relrate[d] << '\n';
			os << *rootstationary[d] << '\n';
		}
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		if (chronoprior)	{
			is >> *Chi >> *Chi2;
		}
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		if (clampdiag)	{
			is >> *contDiagArray;
			is >> *contsigma;
			is >> *contprocess;
		}
		for (int d=0; d<Ndata; d++)	{
			if (SeparateSyn())	{
				is >> *synsigma[d];
				is >> *lognormalsyntree[d];
			}
			changeofbasis[d]->FromStream(is);
			is >> *relrate[d];
			is >> *rootstationary[d];
		}
	}
};

#endif
