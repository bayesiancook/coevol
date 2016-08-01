
#include "ConjugateMultiVariateTreeProcess.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "Chronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"

#include "GammaTreeOmegaMultivariateModel.h"

class ConjugateGammaTreeOmegaMultivariateModel : public GammaTreeOmegaMultivariateModel {

	public:

	ConjugateGammaTreeOmegaMultivariateModel(string datafile, string treefile, string contdatafile, double priorsigma, bool ingc = false, bool inautoregressive = false, bool sample=true, GeneticCodeType type=Universal)	{

		autoregressive = inautoregressive;
		clampdiag = false;
		bool meanexp = false;
		gc = ingc;
		L = gc ? 2 : 1;

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(nucdata, true, type);
		Nsite = codondata->GetNsite();	// # columns
		Nstate = codondata->GetNstate();	// # states (20 for amino acids)

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

		PriorLambda = new Const<PosReal>(10);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);
		gamtree = new GammaTree(tree,One,lambda);

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
			mean = new IIDNormal(Ncont+L,Zero,One);
			process = new ConjugateAutoRegressiveMultiVariateTreeProcess(GetConjugateInverseWishart(),mean,phi,gamtree);
		}
		else	{
			phi = 0;
			mean = 0;
			process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),gamtree);
		}

		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,L+i,i,0);
			}
		}

		// cut off to avoid numerical errors
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		omegatree = new MeanExpTreeFromMultiVariate(process,0,MEAN,false,meanexp);

		// make rate matrices
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);

		if (gc)	{
			gctree = new MeanLogitTreeFromMultiVariate(process,1,MEAN,false);
			rootgc = new Beta(One,One);
			stattree = new GCStatTree(gctree,rootgc);
			nucmatrixtree = new NucMatrixTree(relrate,stattree);
			matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);

			nucmatrix = 0;
			stationary = 0;
		}
		else	{
			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);
			matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

			gctree = 0;
			rootgc = 0;
			stattree = 0;
			nucmatrixtree = 0;
		}

		// make substitution mappings
		// phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, codondata);
		pathconjtree = new BranchMatrixPathConjugateTree(gamtree, matrixtree, codondata);
		phyloprocess = new PathConjugatePhyloProcess(pathconjtree);

		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		RootRegister(PriorLambda);
		RootRegister(One);
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
		}
	}

	~ConjugateGammaTreeOmegaMultivariateModel() {}

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

	void MakeScheduler()	{

		// scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");

		for (int i=0; i<30; i++)	{
			scheduler.Register(new SimpleMove(lambda,1),100,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.1),100,"lambda");
			scheduler.Register(new SimpleMove(lambda,0.01),100,"lambda");
			scheduler.Register(new SimpleMove(gamtree,1),10,"bl");
			scheduler.Register(new SimpleMove(gamtree,0.1),10,"bl");
			scheduler.Register(new SimpleMove(gamtree,0.01),10,"bl");

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
		}
	}

};

