
#include "ConjugateMultiVariateTreeProcess.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"

#include "RateMultivariateModel.h"

class ConjugateRateMultivariateModel : public RateMultivariateModel {

	public:

	ConjugateRateMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, bool ingc, bool inautoregressive, int inconjpath, int contdatatype, bool inomegaratiotree, bool inclamproot, bool inmeanexp, bool sample=true)	{

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		autoregressive = inautoregressive;
		clampdiag = false;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		gc = ingc;
		L = gc ? 2 : 1;

		if (inconjpath == -1)	{
			/*
			if (Nsite > 5000)	{
				conjpath = true;
			}
			else	{
				conjpath = false;
			}
			*/
			conjpath = true;
		}
		else	{
			conjpath = inconjpath;
		}

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);
		Nsite = nucdata->GetNsite();	// # columns
		Nstate = nucdata->GetNstate();	// # states (20 for amino acids)

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

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

		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu);
		mu->ClampAt(1);

		RootAlpha = 0;
		RootBeta = 0;

		if (calibfile != "None")	{
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			cerr << "gamma params : " << a << '\t' << b << '\n';
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,mu,RootAlpha,RootBeta,calibset);
			}
			else if (chronoprior == 1)	{
				cerr << "bd\n";
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,RootAlpha,RootBeta,calibset);
			}
			else	{
				cerr << "error : chronoprior\n";
				exit(1);
			}
		}
		else	{
			chronogram = new Chronogram(tree,mu);
		}

		GammaDiag = new Gamma*[Ncont+L];
		diag = new Rvar<PosReal>*[Ncont+L];
		for (int k=0; k<Ncont+L; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(priorsigma);
		}
		cerr << "diag : " << priorsigma << '\n';
		priorOnSigmaZero = new RvarVec(diag,Ncont+L);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		sigma = new ConjugateInverseWishart(sigmaZero,Ncont+L+1);
		sigma->SetIdentity();

		if (autoregressive)	{
			phi = new Gamma(One,One);
			mean = new IIDUniform(Ncont+L,100);
			process = new ConjugateAutoRegressiveMultiVariateTreeProcess(GetConjugateInverseWishart(),mean,phi,chronogram);
		}
		else	{
			phi = 0;
			mean = 0;
			process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram);
		}

		// process->GetMultiNormal(tree->GetRoot())->ClampAt(log(100.0),3);

		if (clamproot)	{
			process->ClampRoot();
		}

		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,L+i,i,contdatatype);
			}
		}

		// cut off to avoid numerical errors
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		ratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);

		// make rate matrices
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);

		if (gc)	{
			gctree = new MeanLogitTreeFromMultiVariate(process,1,MEAN,false);
			rootgc = new Beta(One,One);
			stattree = new GCStatTree(gctree,rootgc);
			nucmatrixtree = new NucMatrixTree(relrate,stattree);
			nucmatrix = 0;
			stationary = 0;
		}
		else	{
			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);
			gctree = 0;
			rootgc = 0;
			stattree = 0;
			nucmatrixtree = 0;
		}

		// make substitution mappings
		if (conjpath)	{
			if (gc)	{
				pathconjtree = new BranchMatrixPathConjugateTree(ratetree, nucmatrixtree, nucdata);
			}
			else	{
				pathconjtree = new OneMatrixPathConjugateTree(ratetree,nucmatrix,nucdata);
			}
			phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
		}
		else	{
			pathconjtree = 0;
			if (gc)	{
				phyloprocess = new BranchMatrixPhyloProcess(ratetree, nucmatrixtree, nucdata);
			}
			else	{
				phyloprocess = new OneMatrixPhyloProcess(ratetree, nucmatrix, nucdata);
			}
		}

		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		RootRegister(PriorMu);
		RootRegister(One);
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		if (chronoprior == 1)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (autoregressive)	{
			RootRegister(Zero);
		}
		RootRegister(relrate);
		if (!gc)	{
			RootRegister(stationary);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();

			if (RootAlpha)	{
				cerr << "starting chrono : " << GetCalibratedChronogram()->GetLogProb() << '\n';
				cerr << "scale progeny : " << GetCalibratedChronogram()->GetScale()->down.size() << '\n';
			}
		}
	}

	~ConjugateRateMultivariateModel() {}

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

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	/*
	void MakeScheduler()	{

		scheduler.Register(new DSemiConjugateMove<Chronogram,PathConjugateTree>(chronogram,pathconjtree,1,10),1,"chrono");
		scheduler.Register(new DSemiConjugateMove<Chronogram,PathConjugateTree>(chronogram,pathconjtree,0.1,10),1,"chrono");
		scheduler.Register(new DSemiConjugateMove<Chronogram,PathConjugateTree>(chronogram,pathconjtree,0.01,10),1,"chrono");

		if (RootAlpha)	{
			scheduler.Register(new DSemiConjugateMove(GetCalibratedChronogram()->GetScale(),pathconjtree,1,10),1,"root age");
			scheduler.Register(new DSemiConjugateMove(GetCalibratedChronogram()->GetScale(),pathconjtree,0.1,10),1,"root age");
			scheduler.Register(new DSemiConjugateMove(GetCalibratedChronogram()->GetScale(),pathconjtree,0.01,10),1,"root age");
		}

		scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
		scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");

		if (autoregressive)	{
			scheduler.Register(new SimpleMove(phi,10),100,"phi");
			scheduler.Register(new SimpleMove(phi,1),100,"phi");
			scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
			scheduler.Register(new SimpleMove(phi,0.01),100,"phi");

			scheduler.Register(new SimpleMove(mean,1),100,"mean");
			scheduler.Register(new SimpleMove(mean,0.1),100,"mean");
			scheduler.Register(new SimpleMove(mean,0.01),100,"mean");
		}

		scheduler.Register(new DSemiConjugateMove(relrate,pathconjtree,0.1,10),1,"relrates");
		scheduler.Register(new DSemiConjugateMove(relrate,pathconjtree,0.01,10),1,"relrates");

		if (gc)	{
			scheduler.Register(new DSemiConjugateMove(rootgc,pathconjtree,1,10),1,"root gc");
			scheduler.Register(new DSemiConjugateMove(rootgc,pathconjtree,0.1,10),1,"root gc");
			scheduler.Register(new DSemiConjugateMove(rootgc,pathconjtree,0.01,10),1,"root gc");
		}
		else	{
			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new DSemiConjugateMove(stationary,pathconjtree,0.001),10,"stat");
		}


		scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,1,0),5,"root compensation");
		scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.1,0),5,"root compensation");
		scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.01,0),5,"root compensation");

		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}
	*/

	void MakeScheduler()	{

		if (conjpath)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		int nrep = conjpath ? 30 : 1;

		for (int i=0; i<nrep; i++)	{
			if (chronoprior == 1)	{
				scheduler.Register(new SimpleMove(Chi,1),10,"chi");
				scheduler.Register(new SimpleMove(Chi,0.1),10,"chi");
				scheduler.Register(new SimpleMove(Chi2,1),10,"chi2");
				scheduler.Register(new SimpleMove(Chi2,0.1),10,"chi2");
			}
			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			if (RootAlpha)	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
			}

			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");

			if (autoregressive)	{
				scheduler.Register(new SimpleMove(phi,10),100,"phi");
				scheduler.Register(new SimpleMove(phi,1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.01),100,"phi");

				scheduler.Register(new SimpleMove(mean,1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.01),100,"mean");
			}

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			if (gc)	{
				scheduler.Register(new SimpleMove(rootgc,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.01),10,"root gc");
			}
			else	{
				scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
				scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
				scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
				scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
			}

			// scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,1,0),5,"root compensation");
			// scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.1,0),5,"root compensation");
			// scheduler.Register(new MultiVariateRootCompensatoryMove(chronogram,process,phyloprocess,0.01,0),5,"root compensation");
		}
	}

};

