
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include <iostream>

#include <string> 
#include <cmath>

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "CalibratedChronogram.h"
#include "BDChronogram.h"
#include "BDCalibratedChronogram.h"
#include "ShiftedChronogram.h"
#include "BranchTimeTree.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "GeneralConjugatePath.h"
#include "CodonConjugatePath.h"
#include "MultiVariatePropagateMove.h"
#include "Jeffreys.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "LinearCombinationNodeTree2.h"
#include "MeanChronogram.h"
#include "AuxCoevol.h"

#include "ConstrainedLogKappa2.h"
#include "BetaKappaMove.h"

class BranchOmegaMultivariateModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	Tree* splittree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	ContinuousData* contdata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ncont;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;
	
    // prior on div times
	double meanchi;
	double meanchi2;
	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Rvar<PosReal>* Chi;
	Rvar<PosReal>* Chi2;

    // chronogram
	Chronogram* chronogram;

    // absolute root age
	Var<PosReal>* absrootage;
	
    ShiftedChronogram* shiftedchronogram;

    BranchTimeTree* branchtimetree;
    BranchTimeTree* shiftedbranchtimetree;

    // prior on scaling factors for Brownian process
    // (different scaling factor for each trait)
	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;

    // covariance matrix
	ConjugateInverseWishart* sigma;

    // prior for initial value of process at the root
	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

    // process of trait evolution (log-Brownian)
	ConjugateMultiVariateTreeProcess* process;
	
	// structural parameters of the DFE
    // (determine scaling of piN/piS and dN/dS w.r.t. Ne)
	Normal* beta;
	Normal* logkappa1;
	Normal* normal_logkappa2;	
    ConstrainedLogKappa2* constrained_logkappa2;
    Var<Real>* logkappa2;
	
    // indices of piS, piN/piS and generation times in data file
    int idxpiS;
    int idxpiNpiS;
    int idxgentime;

    // node-log-values of key molecular evolutionary quantities
    // omega0, omega, u, Ne and dS
    // all are linear combinations of components of the process (piS, piN/piS and generation time)
	NeutralOmegaLinearCombinationNodeTree* nodeneutralomegatree;
	OmegaLinearCombinationNodeTree* nodeomegatree;
	ULinearCombinationNodeTree* nodeutree;
	NeLinearCombinationNodeTree* nodeNetree;
	SynrateLinearCombinationNodeTree* nodesynratetree;
	
    // branch-specific values of dS and dN/dS
	MeanExpTree* synratetree;
	MeanExpTree* omegatree;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	BranchValPtrTree<RandomSubMatrix>* matrixtree;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

    // divergence times are fixed
	bool clamptree;

    // arithmetic / geodesic approximation for branch-averages of log-Brownian rates (see Poujol et al)
	bool meanexp;
	
    // same coding sequences are used in multiple sequence alignment and polymorphism data
	bool sameseq;
    // no adaptive component in omega (i.e omega = omega_0 = kappa1 N_e^-beta)
	bool noadapt;

	// number of entries of log-Brownian process
    // in addition to those given in data file
    // L == 1 if adaptive component included, 0 otherwise
	int L;

    // number of degrees of freedom of the inverse wishart
	int df;

    // whether tree is calibrated
	bool iscalib;

    // prior on div times
    // 0 : uniform prior
    // 1 : birth death prior
	int chronoprior;

    // number of repetitions for the whole MCMC cyle before saving new point
	int nrep;

    bool shiftages;

	public:

	BranchOmegaMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double priorsigma, int indf, int contdatatype, bool insameseq, bool innoadapt, bool inshiftages, bool inclamptree, bool inmeanexp, int innrep, bool sample=true, GeneticCodeType type=Universal)	{

		clamptree = inclamptree;
		chronoprior = inchronoprior;

		meanexp = inmeanexp;
		
		sameseq = insameseq;
		noadapt = innoadapt;

        shiftages = inshiftages;

		// here L = 1: adaptative omega
		if (!noadapt) {L = 1;}
		else {L = 0;}

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);

		// make codon alignment 
		codondata = new CodonSequenceAlignment(nucdata, true, type);

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		// number of MCMC cycles per point
		nrep = innrep;
		if (nrep == 0)	{
			nrep = 150;
		}

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

        cerr << "chronogram\n";

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;
		meanchi = -1;
		meanchi2 = -1;

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
				if (meanchi != -1)	{
					MeanChi = new Const<PosReal>(meanchi);
					MeanChi2 = new Const<PosReal>(meanchi2);
					Chi = new Exponential(MeanChi,Exponential::MEAN);
					Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				}
				else	{
					double min = 1e-6;
					double max = 1e6;
					Chi = new Jeffreys(min,max,Zero);
					Chi2 = new Jeffreys(min,max,Zero);
				}
				chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,a,b,calibset,chronoprior);
			}
		}
		else	{
			iscalib = false;
			chronogram = new Chronogram(tree,One);
		}
		if (clamptree)	{
			chronogram->Clamp();
		}

		if (iscalib) {
			absrootage = GetCalibratedChronogram()->GetScale();
		}
		else {
			absrootage = new Const<PosReal>(rootage);	
		}
		
        branchtimetree = new BranchTimeTree(chronogram, One);

		cerr << "sigma\n";

		// Ncont : number of quantitative traits
		// L : 1 if omega_a included, 0 otherwise
		// create an array of scaling factor k_i, i=1..Ncont + L
		double mindiag = 0.01;
		double maxdiag = 100;
		DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		// create a diagonal matrix, with the k_i's along the diagonal
		sigmaZero = new SigmaZero(DiagArray);

		// create covariance matrix
		// from an inverse wishart of parameter sigma0
		// and df degrees of freedom
		df = Ncont + L + indf;
		sigma = new ConjugateInverseWishart(sigmaZero, df);

		// create a multivariate brownian process (of dimension Ncont + L)
		// along the chronogram, and with covariance matrix sigma

		cerr << "log-Brownian process\n";
		process = new ConjugateMultiVariateTreeProcess(sigma,branchtimetree);
		process->Reset();
		
        idxpiS = -1;
        idxpiNpiS = -1;
        idxgentime = -1;

		for (int i=0; i<Ncont; i++) {
			if (GetContinuousData()->GetCharacterName(i)=="piS") {
                idxpiS = i;
			}	
			else if (GetContinuousData()->GetCharacterName(i)=="piNpiS") {
				idxpiNpiS = i;
			}	
			else if (GetContinuousData()->GetCharacterName(i)=="generation_time") {
				idxgentime = i;
			}	
		}	

        if (idxpiS == -1)   {
            cerr << "error: did not find piS in cont data file\n";
            exit(1);
        }
        if (idxpiNpiS == -1)    {
            cerr << "error: did not find piNpiS in cont data file\n";
            exit(1);
        }
        if (idxgentime == -1)   {
            cerr << "error: did not find generation_time in cont data file\n";
            exit(1);
        }
		
        // center process on mean observed values
        // necessary for gentle start (avoiding numerical errors)
		process->PiecewiseTranslation(GetContinuousData()->GetMeanLog(idxpiS), idxpiS + L, 1);
		process->PiecewiseTranslation(GetContinuousData()->GetMeanLog(idxpiNpiS), idxpiNpiS + L, 1);
		process->PiecewiseTranslation(GetContinuousData()->GetMeanLog(idxgentime), idxgentime + L, 1);
		
		// condition the multivariate process
		// on the matrix of quantitative traits.
        // offset of L == 1 if adaptive component, 0 otherwise
        for (int i=0; i<Ncont; i++)	{
            process->SetAndClamp(contdata,L+i,i,contdatatype);
        }

        cerr << "DFE structural parameters\n";

		//create the combination factors 
        // structural parameters of the DFE 
        // these parameters determine the response of piN/piS and dN/dS to changes in Ne:
        // omega_0 = kappa_1 N_e^{-beta}
        // piN/piS = kappa_2 N_e^{-beta}
        //
        // if adaptive component included:
        // omega = omega_0 + omega_a (adaptive version of the model; otherwise, omega_a set to 0)
		
        // shape parameter of DFE
		beta = new Normal(Zero, One);

        // kappa1 and kappa2
        // normally, kappa1 and kappa2 are linked: two different functions of beta and of the mean s of the DFE
        // but if polymorphism data and multiple sequence alignments are not based on same coding sequences
        // then these two variables may not match
        // in that case, kappa1 and kapp2 are considered as 2 independent degrees of freedom of the model
		logkappa1 = new Normal(Zero, One);
		if (!sameseq) {
            normal_logkappa2 = new Normal(Zero, One);
            constrained_logkappa2 = 0;
            logkappa2 = normal_logkappa2;
        }
        else    {
            normal_logkappa2 = 0;
            constrained_logkappa2 = new ConstrainedLogKappa2(beta, logkappa1);
            logkappa2 = constrained_logkappa2;
        }

        // reasonable initial values
		beta->setval(0.2);
		logkappa1->setval(0.9);
		if (!sameseq) {normal_logkappa2->setval(0.4);}
		
		// create the node-trees for omega0, Ne, u, dS, Ne, and tau (all in log)
        // based on the structural parameters and the 3-dim process (piS, piNpiS, gentime)
        // specifically:
        // log u = log piS
		
        // log omega_0 = log piN/piS + log kappa_1 - log kappa_2
        nodeneutralomegatree = new NeutralOmegaLinearCombinationNodeTree(process, logkappa1, logkappa2, idxpiNpiS);

        // using : log Ne -1/beta * (log piN/piS - log kappa_2)
        // log u = log piS - log N_e - log 4.0
        //       = log piS + 1/beta * log piN/piS - 1/beta * log kappa_2 - log 4.0
        nodeutree = new ULinearCombinationNodeTree(process, beta, logkappa2, idxpiS, idxpiNpiS);

        // log Ne = log piS - log u - log 4.0
        // (note that it could also have been computed based on relation given above)
		nodeNetree = new NeLinearCombinationNodeTree(process, nodeutree, idxpiS); 

        // log dS = log u - log tau
		nodesynratetree = new SynrateLinearCombinationNodeTree(process, nodeutree, absrootage, idxgentime);
		
        // if adaptive component included, then compute omega = omega_0 + omega_a
        // check that: sum should be in natural units
		if (!noadapt) {
            cerr << "check node omega tree\n";
            exit(1);
            nodeomegatree = new OmegaLinearCombinationNodeTree(process, nodeneutralomegatree, 0);
        }

        // account for ancestral polymorphism
        if (shiftages)  {
            shiftedchronogram = new ShiftedChronogram(chronogram, nodeNetree, process, absrootage, idxgentime);
            shiftedbranchtimetree = new BranchTimeTree(shiftedchronogram, One);
        }
        else    {
            shiftedchronogram = 0;
            shiftedbranchtimetree = branchtimetree;
        }
	
		// create the branch-specific mean dS and mean dN/dS
		
		// the times given by the chronogram with the rate 
		synratetree = new MeanExpTree(nodesynratetree, shiftedbranchtimetree, INTEGRAL, false);

		// create the dN/dS on each branch, nased on the second entry of the multivariate process
		if (!noadapt) {
            omegatree = new MeanExpTree(nodeomegatree, shiftedbranchtimetree, MEAN, false);
        }
        else    {
            omegatree = new MeanExpTree(nodeneutralomegatree, shiftedbranchtimetree, MEAN, false);
        }
		
		// create u on each branch, nased on the third entry of the multivariate process

		cerr << "codon matrices\n";
		// create a GTR nucleotide matrix
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);

		// create a series of codon matrices, one for each branch of the chronogram
		// each codon matrix is parameterized by the global nucleotide matrix, combined with the branch-specific dN/dS
		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

        cerr << "phylo process\n";
		// create a phylo process based on this array of branch specific matrices
		// and condition it on the multiple alignment codondata
		// pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
		pathconjtree = new MGCodonBranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
		phyloprocess = new PathConjugatePhyloProcess(pathconjtree);


		cerr << "unfold model\n";
		phyloprocess->Unfold();
		if (sample)	{
            cerr << "sample\n";
			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}

		cerr << "register\n";
		// register model
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(beta);
		RootRegister(logkappa1);
		RootRegister(logkappa2);
		RootRegister(absrootage);
		if (MeanChi)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		Register();

		MakeScheduler();
		if (sample)	{
            cerr << "update\n";
			Update();
		}
        cerr << "ok\n";
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~BranchOmegaMultivariateModel() {}
	
	// accessors
	Tree* GetTree() {return tree;}
	
	SequenceAlignment* GetData()	{
		return codondata;
	}

	MeanExpTree* GetSynRateTree() {return synratetree;}
	MeanExpTree* GetOmegaTree() {return omegatree;}


	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	
	SynrateLinearCombinationNodeTree*  GetSynrateNodeTree() {return nodesynratetree;}

	NeutralOmegaLinearCombinationNodeTree*  GetNeutralOmegaNodeTree() {return nodeneutralomegatree;}
	
	OmegaLinearCombinationNodeTree*  GetOmegaNodeTree() {return nodeomegatree;}
	
	NeLinearCombinationNodeTree* GetNeNodeTree() {return nodeNetree;}
	
	ULinearCombinationNodeTree* GetUNodeTree() {return nodeutree;}
	
	Chronogram* GetChronogram() {return chronogram;}

    LengthTree* GetLengthTree() {return branchtimetree;}

	bool isCalibrated()	{
		return iscalib;
	}

	CalibratedChronogram* GetCalibratedChronogram()	{
		if (! iscalib)	{
			cerr << "error : calibrated chronogram does not exist under uncalibrated model\n";
			exit(1);
		}
		return dynamic_cast<CalibratedChronogram*>(chronogram);
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

	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}
	

	CovMatrix* GetCovMatrix() {return sigma;}

	// probability computation

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += chronogram->GetLogProb();
		if (chronoprior)	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		total += beta->GetLogProb();
		total += logkappa1->GetLogProb();
		if (!sameseq) {total += normal_logkappa2->GetLogProb();}
		return total;
	}

	double GetLogLikelihood()	{
		double ret = pathconjtree->GetLogProb();
		return ret;
	}

	// MCMC schedule
	virtual void MakeScheduler()	{

		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");

		for (int i=0; i<nrep; i++)	{

			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			if (isCalibrated())	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");

				if (chronoprior)	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
			}

			int n = taxonset->GetNtaxa() * 100;
			scheduler.Register(new MultiVariatePropagateMove(process,1,0.1,0.1),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.5,0.5),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.9,0.9),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.01,0.99,0.99),n,"propmove");

			scheduler.Register(new SimpleMove(process,1),150,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),150,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),150,"multinormal");

			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.01,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.001,10),1,"conjugate sigma - process");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"cov matrix scaling factors");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"cov matrix scaling factors");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"cov matrix scaling factors");

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
			
			scheduler.Register(new SimpleMove(beta,0.003),10,"beta");
			scheduler.Register(new SimpleMove(beta,0.0008),10,"beta");
			scheduler.Register(new SimpleMove(beta,0.0003),10,"beta");
			
			scheduler.Register(new SimpleMove(logkappa1,0.15),10,"logkappa1");
			scheduler.Register(new SimpleMove(logkappa1,0.08),10,"logkappa1");
			scheduler.Register(new SimpleMove(logkappa1,0.01),10,"logkappa1");
			
			if (!sameseq) {
				scheduler.Register(new SimpleMove(normal_logkappa2,0.05),10,"logkappa2");
				scheduler.Register(new SimpleMove(normal_logkappa2,0.008),10,"logkappa2");
				scheduler.Register(new SimpleMove(normal_logkappa2,0.001),10,"logkappa2");

				scheduler.Register(new BetaKappaMove(beta,normal_logkappa2,1,11),10,"betakappa2");
				scheduler.Register(new BetaKappaMove(beta,normal_logkappa2,0.1,11),10,"betakappa2");
				scheduler.Register(new BetaKappaMove(beta,normal_logkappa2,0.01,11),10,"betakappa2");
			}
		}
	}

	void drawSample()	{
		if (chronoprior)	{
			Chi->Sample();
			Chi2->Sample();
		}
		chronogram->Sample();
		DiagArray->Sample();
		sigma->Sample();
		// sigma->SetIdentity();
		process->Sample();
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		beta->Sample();
		logkappa1->Sample();
		if (!sameseq) {normal_logkappa2->Sample();}
		relrate->Sample();
		stationary->Sample();
		phyloprocess->Sample();
	}

    void UpdateLengthTree() {
        branchtimetree->specialUpdate();
        if (shiftages)  {
            shiftedchronogram->specialUpdate();
            shiftedbranchtimetree->specialUpdate();
        }
    }

    int GetShiftNerr()  {
        return shiftedbranchtimetree->GetNerr();
    }

	// summary statistics
	double GetTotalLength()	{
		return GetSynRateTree()->GetTotal();
	}

	double GetMeanOmega()	{
		return omegatree->GetMean();
	}

	double GetBeta() const	{
		return beta->val();
	}

	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tdS\tomega";
		if (isCalibrated())	{
			os << "\trootage";
		}
		os << "\tbeta";
		os << "\tlogkappa1";
		os << "\tlogkappa2";
        for (int k=0; k<Ncont+L; k++)   {
            os << "\tmean_" << k;
        }
        for (int k=0; k<Ncont+L; k++)   {
            os << "\troot_" << k;
        }
		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}
		os << "\tstatent";
		os << "\trrent";
		if (chronoprior)	{
			os << "\tchronologprior";
			os << "\tchi\tchi2";
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetTotalLength();
		os << '\t' << GetMeanOmega();
		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		os << '\t' << beta->val();
		os << '\t' << logkappa1->val();
		os << '\t' << logkappa2->val();
		for (int k=0; k<Ncont+L; k++)	{
            os << '\t' << process->GetMean(k);
        }

		for (int k=0; k<Ncont+L; k++)	{
            os << '\t' << process->GetVal(tree->GetRoot(),k);
        }

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}
		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->val().GetEntropy();
		if (chronoprior)	{
			os << '\t' << chronogram->GetLogProb();
			os << '\t' << Chi->val() << '\t' << Chi2->val();
		}

		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
			if (chronoprior)	{
				os << *Chi << '\t' << *Chi2 << '\n';
			}
		}
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		os << *beta << '\n';
		os << *logkappa1 << '\n';
		if (!sameseq) {os << *normal_logkappa2 << '\n';}
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
			if (chronoprior)	{
				is >> *Chi >> *Chi2;
			}
		}
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		is >> *beta;
		is >> *logkappa1;
		if (!sameseq) {is >> *normal_logkappa2;}
		is >> *relrate;
		is >> *stationary;
	}
};

#endif
