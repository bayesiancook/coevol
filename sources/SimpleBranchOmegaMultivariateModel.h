
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include <iostream>

#include <string> 

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "CalibratedChronogram.h"
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
#include "LinearCombinationNodeTree.h"
#include "MeanChronogram.h"
#include "AuxCoevol.h"


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

	Chronogram* chronogram;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* sigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	ConjugateMultiVariateTreeProcess* process;

	
	//parameter of the linear combination
	
	Normal* gamma;
	Normal* beta;	
	
	Var<PosReal>* absrootage;
	
	double* synrateslope;
	double* omegaslope;
	double* u_Neslope;
	double* adaptative_omegaslope;
	double* oppositealphaslope;
	

	SynrateLinearCombinationNodeTree* nodesynratetree;
	OmegaLinearCombinationNodeTree* nodeomegatree;
	U_NeLinearCombinationNodeTree* nodeu_Netree;
	Adaptative_omegaLinearCombinationNodeTree* nodeadaptative_omegatree;
	OppositeAlphaLinearCombinationNodeTree* nodeoppositealphatree;

	MeanExpTreeFromMultiVariate* mutratetree;
	MeanExpTreeFromMultiVariate* Netree;
	MeanExpTree* synratetree;
	MeanExpTree* omegatree;
	//MeanExpTree* Netree;
	//MeanExpTree* adaptative_omegatree;

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

	bool withNe;
	bool clamptree;
	bool meanexp;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int nrep;

	int df;

	bool iscalib;

	public:


	BranchOmegaMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, double priorsigma, int indf, int contdatatype, bool inwithNe, bool inclamptree, bool inmeanexp, int innrep, bool sample=true, GeneticCodeType type=Universal)	{

		withNe = inwithNe;

		// if we want the divergence times to be fixed
		clamptree = inclamptree;

		// arithmetic / geodesic averages (see Poujol et al)
		meanexp = inmeanexp;

		// number of components of the Brownian process corresponding to rates

		// here L = 1: u
		L = 1;

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

		cerr << "new chrono\n";
		if (calibfile != "None")	{
			iscalib = true;
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);
			chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
		}
		else	{
			iscalib = false;
			chronogram = new Chronogram(tree,One);
		}
		if (clamptree)	{
			chronogram->Clamp();
		}
	
		cerr << "ok\n";
		
		if (iscalib) {
			absrootage = GetCalibratedChronogram()->GetScale();
		}
		else {
			absrootage = new Const<PosReal>(rootage);	
		}
		
		
		cerr << "checking rootage" << absrootage->val();	

		// Ncont : number of quantitative traits
		// L : number of substitution parameters coevolving with traits (typically, 2: dS and dN/dS).
		// create an array of positive variables kappa_i, i=1..Ncont + L
		double mindiag = 0.001;
		double maxdiag = 1000;
		DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		// create a diagonal matrix, with the kappa_i along the diagonal
		sigmaZero = new SigmaZero(DiagArray);

		// create covariance matrix
		// from an inverse wishart of parameter sigma0
		// and df degrees of freedom
		df = Ncont + L + indf;
		cerr << "sigma\n";
		sigma = new ConjugateInverseWishart(sigmaZero, df);

		// create a multivariate brownian process (of dimension Ncont + L)
		// along the chronogram, and with covariance matrix sigma
		cerr << "process\n";
		process = new ConjugateMultiVariateTreeProcess(sigma,chronogram);
		process->Reset();
		
		if (! withNe) {
			process->PiecewiseTranslation(-21, 0, 1);
		}
		else {
			process->PiecewiseTranslation(18, 0, 1);
		}
		
		// condition the multivariate process
		// on the matrix of quantitative traits.
		// note the offset here : first trait corresponds to entry L+1 of the process, etc.

		// this is because the first L entries of the process correspond to the substitution variables (mut)

		for (int i=0; i<Ncont; i++)	{
			process->SetAndClamp(contdata,L+i,i,contdatatype);
		}

		// just for numerical stability of the starting point
		//for (int l=0; l<L; l++)	{
		//	process->CutOff(1,l);
		//}

		//create the combination factors 
		
		gamma = new Normal(Zero, One);
		beta = new Normal(Zero, One);

		//create the slopes which define the linear combination
		
		synrateslope = new double[L + Ncont];
		omegaslope = new double[L + Ncont];
		u_Neslope = new double[L + Ncont];
		adaptative_omegaslope = new double[L + Ncont];
		oppositealphaslope = new double[L + Ncont];
		
		CreateSynrateSlope();
		CreateOmegaSlope();
		CreateU_NeSlope();
		CreateAdaptative_omegaSlope();
		CreateOppositeAlphaSlope();
		

		// create the node tree obtained from the linear combinations
		
		nodesynratetree = new SynrateLinearCombinationNodeTree(process, absrootage, synrateslope, withNe);
		nodeomegatree = new OmegaLinearCombinationNodeTree(process, gamma, beta, omegaslope, withNe);
		nodeu_Netree = new U_NeLinearCombinationNodeTree(process, u_Neslope); 
		nodeadaptative_omegatree = new Adaptative_omegaLinearCombinationNodeTree(process, nodeomegatree, adaptative_omegaslope);
		nodeoppositealphatree = new OppositeAlphaLinearCombinationNodeTree(process, nodeomegatree, oppositealphaslope);
		
		// create the branch lengths resulting from combining

		// the times given by the chronogram with the rate 
		synratetree = new MeanExpTree(nodesynratetree, chronogram, INTEGRAL, false);

		// create the dN/dS on each branch, nased on the second entry of the multivariate process
		omegatree = new MeanExpTree(nodeomegatree, chronogram, MEAN, false);

		//Netree = new MeanExpTree(nodeNetree, chronogram, MEAN, false);
		
		//adaptative_omegatree = new MeanExpTree(nodeadaptative_omegatree, chronogram, MEAN, false);
		
		// create u on each branch, nased on the third entry of the multivariate process
		if (!withNe) {mutratetree = new MeanExpTreeFromMultiVariate(process,0,MEAN,false,meanexp);}
		if (withNe) {Netree = new MeanExpTreeFromMultiVariate(process,0,MEAN,false,meanexp);}

		cerr << "matrix\n";

		// create a GTR nucleotide matrix
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);

		// create a series of codon matrices, one for each branch of the chronogram
		// each codon matrix is parameterized by the global nucleotide matrix, combined with the branch-specific dN/dS
		matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

		// create a phylo process based on this array of branch specific matrices
		// and condition it on the multiple alignment codondata
		// pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
		pathconjtree = new MGCodonBranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
		phyloprocess = new PathConjugatePhyloProcess(pathconjtree);


		cerr << "unfold\n";
		phyloprocess->Unfold();
		cerr << "sample\n";
		if (sample)	{
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
		RootRegister(gamma);
		RootRegister(beta);
		RootRegister(absrootage);
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~BranchOmegaMultivariateModel() {
		DeleteSynrateSlope();
		DeleteOmegaSlope();
		DeleteU_NeSlope();
		DeleteAdaptative_omegaSlope();
		DeleteOppositeAlphaSlope();
		}
	
	void CreateSynrateSlope() {
		string cha("generation_time");
		string cha2("piS");
		if (!withNe) {
			synrateslope[0] = 1;
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					synrateslope[i+1] = -1;
				}	
				else {
					synrateslope[i+1] = 0;
				}
			}
		}
		else {
			synrateslope[0] = -1;
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					synrateslope[i+1] = -1;
				}	
				else if (GetContinuousData()->GetCharacterName(i)==cha2) {
					synrateslope[i+1] = 1;
				}	
				else {
					synrateslope[i+1] = 0;
				}
			}
		}
	}			

		 
	void DeleteSynrateSlope() {
		delete synrateslope;
	}
	
	void CreateOmegaSlope() {
		string cha("piS");
		if (!withNe) {
			omegaslope[0] = 1;
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					omegaslope[i+1] = -1;
				}	
				else {
					omegaslope[i+1] = 0;
				}
			}
		}
		else {
			omegaslope[0] = -1;			
			for (int i=0; i<Ncont; i++) {
				omegaslope[i+1] = 0;
			}
		}
	}		
		 
	
	void DeleteOmegaSlope() {
		delete omegaslope;
	}
	
	void CreateU_NeSlope() {
		string cha("piS");
		u_Neslope[0] = -1;
		for (int i=0; i<Ncont; i++) {
			if (GetContinuousData()->GetCharacterName(i)==cha) {
				u_Neslope[i+1] = 1;
			}	
			else {
				u_Neslope[i+1] = 0;
			}
		}
	}				
		
	void DeleteU_NeSlope() {
		delete u_Neslope;
	}
	
	
	void CreateAdaptative_omegaSlope() {
		string cha("piNpiS");
		adaptative_omegaslope[0] = 0;
		for (int i=0; i<Ncont; i++) {
			if (GetContinuousData()->GetCharacterName(i)==cha) {
				adaptative_omegaslope[i+1] = -1;
			}	
			else {
				adaptative_omegaslope[i+1] = 0;
			}
		}
	}	
	
	void DeleteAdaptative_omegaSlope() {
		delete adaptative_omegaslope;
	}	
	
	
		
	void CreateOppositeAlphaSlope() {
		string cha("piNpiS");
		oppositealphaslope[0] = 0;
		for (int i=0; i<Ncont; i++) {
			if (GetContinuousData()->GetCharacterName(i)==cha) {
				oppositealphaslope[i+1] = 1;
			}	
			else {
				oppositealphaslope[i+1] = 0;
			}
		}
	}	
	
	void DeleteOppositeAlphaSlope() {
		delete oppositealphaslope;
	}			
		
		
	// accessors
	Tree* GetTree() {return tree;}

	SequenceAlignment* GetData()	{
		return codondata;
	}

	//MeanExpTree* GetAdaptative_omegaTree() {return adaptative_omegatree;}
	//MeanExpTree* GetNeTree() {return Netree;}
	MeanExpTree* GetSynRateTree() {return synratetree;}
	MeanExpTree* GetOmegaTree() {return omegatree;}
	MeanExpTreeFromMultiVariate* GetMutTree() {return mutratetree;}
	MeanExpTreeFromMultiVariate* GetNeTree() {return Netree;}


	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	
	SynrateLinearCombinationNodeTree*  GetSynrateNodeTree() {return nodesynratetree;}
	
	OmegaLinearCombinationNodeTree*  GetOmegaNodeTree() {return nodeomegatree;}
	
	U_NeLinearCombinationNodeTree* GetNeNodeTree() {return nodeu_Netree;}
	
	Adaptative_omegaLinearCombinationNodeTree* GetAdaptative_omegaNodeTree() {return nodeadaptative_omegatree;}
	
	OppositeAlphaLinearCombinationNodeTree* GetOppositeAlphaNodeTree() {return nodeoppositealphatree;}


	Chronogram* GetChronogram() {return chronogram;}

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
		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		total += gamma->GetLogProb();
		total += beta->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = pathconjtree->GetLogProb();
		return ret;
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

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

			scheduler.Register(new SimpleMove(DiagArray,10),10,"kappa");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"kappa");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"kappa");

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
			
			scheduler.Register(new SimpleMove(gamma,0.1),10,"gamma");
			scheduler.Register(new SimpleMove(gamma,0.01),10,"gamma");
			scheduler.Register(new SimpleMove(gamma,0.001),10,"gamma");
			
			scheduler.Register(new SimpleMove(beta,0.3),10,"beta");
			scheduler.Register(new SimpleMove(beta,0.1),10,"beta");
			scheduler.Register(new SimpleMove(beta,0.03),10,"beta");
			
		}
	}

	void drawSample()	{
		chronogram->Sample();
		DiagArray->Sample();
		sigma->Sample();
		// sigma->SetIdentity();
		process->Sample();
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		gamma->Sample();
		beta->Sample();
		relrate->Sample();
		stationary->Sample();
		phyloprocess->Sample();
	}

	// summary statistics
	double GetTotalLength()	{
		return GetSynRateTree()->GetTotal();
	}

	double GetMeanOmega()	{
		return omegatree->GetMean();
	}

	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tdS\tomega";
		if (isCalibrated())	{
			os << "\trootage";
		}
		os << "\tgamma";
		os << "\tbeta";
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
		os << '\t' << gamma->val();
		os << '\t' << beta->val();
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
		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		os << *gamma << '\n';
		os << *beta << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		is >> *gamma;
		is >> *beta;
		is >> *relrate;
		is >> *stationary;
	}
};

#endif
