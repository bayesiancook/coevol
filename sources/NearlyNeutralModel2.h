
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include <iostream>

#include <string> 
#include <cmath>

// #include <boost/math/special_functions/zeta.hpp>

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "CalibratedChronogram.h"
#include "BDChronogram.h"
#include "BDCalibratedChronogram.h"
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

class GammaBetaMove : public MCUpdate, public Mnode	{

	public:

	GammaBetaMove(Rvar<Real>* ingamma, Rvar<Real>* inbeta, double intuning, double ina)	{
		gamma = ingamma;
		beta = inbeta;
		tuning = intuning;
		a = ina;
		gamma->Register(this);
		beta->Register(this);
	}

	double Move(double tuning_modulator = 1)	{

		double acc = 1.0;
		if ((!gamma->isClamped()) && (! beta->isClamped()))	{
			
			Corrupt(true);
			double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
			gamma->setval(gamma->val() + u);
			beta->setval(beta->val() + a*u);
			double logratio = Update();
			acc = (log(Random::Uniform()) < logratio);
			if (! acc)	{
				Corrupt(false);
				Restore();
			}
		}
		return acc;

	}

	private:
	Rvar<Real>* gamma;
	Rvar<Real>* beta;
	double tuning;
	double a;
};

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
	
	Const<PosReal>* K;

	double meanchi;
	double meanchi2;
	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Rvar<PosReal>* Chi;
	Rvar<PosReal>* Chi2;

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
	Normal* beta2;	
	
	Var<PosReal>* absrootage;
	
	double* neutralomegaslope;
	double* omegaslope;
	double* uslope;
	double* synrateslope;
	double* Neslope;
		

	NeutralOmegaLinearCombinationNodeTree* nodeneutralomegatree;
	OmegaLinearCombinationNodeTree* nodeomegatree;
	ULinearCombinationNodeTree* nodeutree;
	SynrateLinearCombinationNodeTree* nodesynratetree;
	NeLinearCombinationNodeTree* nodeNetree;
	

	MeanExpTree* omegatree;
	MeanExpTree* synratetree;

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

	bool clamptree;
	bool meanexp;
	
	bool sameseq;
	
	bool noadapt;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int nrep;

	int df;

	bool iscalib;
	int chronoprior;

	public:

    /*
	double NeutralityIndexFactor(int nind, int jmax)	{


        double denom = 0;
        for (int k=1; k<2*nind; k++)    {
            denom += 1.0 / k;
        }

		double total = 0;
        for (int j=2; j<=jmax; j++)	{
            double z = boost::math::zeta(j) / j;
            double num = 0;
            for (int k=1; k<2*nind; k++)    {
                num += exp(Random::logGamma(2*nind+1) - Random::logGamma(k+1) - Random::logGamma(2*nind-k+1) + Random::logGamma(k) + Random::logGamma(2*nind-k+j) - Random::logGamma(2*nind+j));
            }
            total += z * num / denom;
        }
		cerr << "NI factor : " << total << '\n';
		cerr << "NI(0.2)   : " << 1 + 0.2*total << '\n';
		return total;
	}
    */


	BranchOmegaMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double priorsigma, int indf, int contdatatype, bool insameseq, bool innoadapt, bool inclamptree, bool inmeanexp, int innrep, bool sample=true, GeneticCodeType type=Universal)	{

		// if we want the divergence times to be fixed
		clamptree = inclamptree;
		chronoprior = inchronoprior;

		// arithmetic / geodesic averages (see Poujol et al)
		meanexp = inmeanexp;
		
		sameseq = insameseq;
		
		noadapt = innoadapt;

		// number of components of the Brownian process corresponding to rates

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
		
		if (sameseq) {	 
			cerr << "same seq: still to be implemented\n";
			exit(1);
			/*
			int i(2);
			double k(0);
			double tempo(0);
			do {
				tempo = (boost::math::zeta (i) * factorielle(1999) * factorielle(i) * factorielle(2001)) / (factorielle(2000 + i) * factorielle(1999) * i);
				k += tempo;
				i++;
			}while (tempo > pow(10, -7));
		
			K = new Const<PosReal>(k);
			*/
		}


		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;
		meanchi = -1;
		meanchi2 = -1;

		cerr << "new chrono\n";
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
	
		cerr << "ok\n";
		
		if (iscalib) {
			absrootage = GetCalibratedChronogram()->GetScale();
		}
		else {
			absrootage = new Const<PosReal>(rootage);	
		}
		
		
		cerr << "checking rootage \t" << absrootage->val() << "\n";	

		// Ncont : number of quantitative traits
		// L : number of substitution parameters coevolving with traits (typically, 2: dS and dN/dS).
		// create an array of positive variables kappa_i, i=1..Ncont + L
		double mindiag = 0.01;
		double maxdiag = 100;
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
		
		int indice1;
		int indice2;
		int indice3;
		string char1("piS");
		string char2("piNpiS");
		string char3("generation_time");
		for (int i=0; i<Ncont; i++) {
			if (GetContinuousData()->GetCharacterName(i)==char1) {
				indice1 = i;
			}	
			else if (GetContinuousData()->GetCharacterName(i)==char2) {
				indice2 = i;
			}	
			else if (GetContinuousData()->GetCharacterName(i)==char3) {
				indice3 = i;
			}	
		}	
		
		process->PiecewiseTranslation(GetContinuousData()->GetMeanLog(indice1), indice1+L, 1);
		process->PiecewiseTranslation(GetContinuousData()->GetMeanLog(indice2), indice2+L, 1);
		process->PiecewiseTranslation(GetContinuousData()->GetMeanLog(indice3), indice3+L, 1);
		
		// condition the multivariate process
		// on the matrix of quantitative traits.
		// note the offset here : first trait corresponds to entry L+1 of the process, etc.

		// this is because the first L entries of the process correspond to the substitution variables (mut)

        for (int i=0; i<Ncont; i++)	{
            process->SetAndClamp(contdata,L+i,i,contdatatype);
        }

		//create the combination factors 
		
		gamma = new Normal(Zero, One);
		beta = new Normal(Zero, One);
		if (!sameseq) {beta2 = new Normal(Zero, One);}
		
		gamma->setval(0.2);
		beta->setval(0.9);
		if (!sameseq) {beta2->setval(0.4);}
		
		//create the slopes which define the linear combination
		
		neutralomegaslope = new double[L + Ncont];
		if (!noadapt) {omegaslope = new double[L + Ncont];}
		uslope = new double[L + Ncont];
		synrateslope = new double[L + Ncont];
		Neslope = new double[L + Ncont];
		
		CreateNeutralOmegaSlope();
		if (!noadapt) {CreateOmegaSlope();}
		CreateUSlope();
		CreateSynrateSlope();
		CreateNeSlope();


		// create the node tree obtained from the linear combinations
		
		if (sameseq) {nodeneutralomegatree = new NeutralOmegaLinearCombinationNodeTree(process, gamma, K, neutralomegaslope,sameseq);}
		if (!sameseq) {nodeneutralomegatree = new NeutralOmegaLinearCombinationNodeTree(process, beta, beta2, neutralomegaslope,sameseq);}
		if (!noadapt) {nodeomegatree = new OmegaLinearCombinationNodeTree(process, nodeneutralomegatree, omegaslope);}
		if (sameseq) {nodeutree = new ULinearCombinationNodeTree(process, gamma, beta, K, uslope,sameseq);}
		if (!sameseq) {nodeutree = new ULinearCombinationNodeTree(process, gamma, beta2, uslope, sameseq);}
		nodesynratetree = new SynrateLinearCombinationNodeTree(process, nodeutree, absrootage, synrateslope);
		nodeNetree = new NeLinearCombinationNodeTree(process, nodeutree, Neslope); 
		
		// create the branch lengths resulting from combining
		
		// the times given by the chronogram with the rate 
		synratetree = new MeanExpTree(nodesynratetree, chronogram, INTEGRAL, false);

		// create the dN/dS on each branch, nased on the second entry of the multivariate process
		if (!noadapt) {omegatree = new MeanExpTree(nodeomegatree, chronogram, MEAN, false);}
		if (noadapt) {omegatree = new MeanExpTree(nodeneutralomegatree, chronogram, MEAN, false);}
		
		// create u on each branch, nased on the third entry of the multivariate process

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
		if (!sameseq) {RootRegister(beta2);}
		RootRegister(absrootage);
		if (MeanChi)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
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
	~BranchOmegaMultivariateModel() {
		DeleteNeutralOmegaSlope();
		if (!noadapt) {DeleteOmegaSlope();}
		DeleteUSlope();
		DeleteSynrateSlope();
		DeleteNeSlope();
		}
	
	
	void CreateUSlope() {
		if (sameseq) {
			if (!noadapt) {
				string cha1("piS");
				string cha2("piNpiS");
				uslope[0] = 0;
				for (int i = 0; i<Ncont; i++) {
					if (GetContinuousData()->GetCharacterName(i)==cha1) {
						uslope[i+1] = 1;
					}	
					else if (GetContinuousData()->GetCharacterName(i)==cha2) {
						uslope[i+1] = 1000;
					}
					else {
						uslope[i+1] = 0;
					}
				}
			}
			else {
				string cha1("piS");
				string cha2("piNpiS");
				for (int i = 0; i<Ncont; i++) {
					if (GetContinuousData()->GetCharacterName(i)==cha1) {
						uslope[i] = 1;
					}	
					else if (GetContinuousData()->GetCharacterName(i)==cha2) {
						uslope[i] = 1000;
					}
					else {
						uslope[i] = 0;
					}
				}
			}
		}	
		if (!sameseq) {
			if (!noadapt) {
				string cha1("piS");
				string cha2("piNpiS");
				uslope[0] = 0;
				for (int i = 0; i<Ncont; i++) {
					if (GetContinuousData()->GetCharacterName(i)==cha1) {
						uslope[i+1] = 1;
					}	
					else if (GetContinuousData()->GetCharacterName(i)==cha2) {
						uslope[i+1] = 1000;
					}
					else {
						uslope[i+1] = 0;
					}
				}
			}
			else {
				string cha1("piS");
				string cha2("piNpiS");
				for (int i = 0; i<Ncont; i++) {
					if (GetContinuousData()->GetCharacterName(i)==cha1) {
						uslope[i] = 1;
					}	
					else if (GetContinuousData()->GetCharacterName(i)==cha2) {
						uslope[i] = 1000;
					}
					else {
						uslope[i] = 0;
					}
				}
			}	
		}	
	}	
	
	void DeleteUSlope() {
		delete uslope;
	}	
			
	
	void CreateSynrateSlope() {
		if (!noadapt) {
			string cha("generation_time");
			synrateslope[0] = 0;
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
			string cha("generation_time");
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					synrateslope[i] = -1;
				}	
				else {
					synrateslope[i] = 0;
				}
			}
		}	
	}	
				

		 
	void DeleteSynrateSlope() {
		delete synrateslope;
	}
	
	
	void CreateNeutralOmegaSlope() {
		if (!noadapt) {
			string cha("piNpiS");
			neutralomegaslope[0] = 0;
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					neutralomegaslope[i+1] = 1;
				}	
				else {
					neutralomegaslope[i+1] = 0;
				}
			}
		}
		else {
			string cha("piNpiS");
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					neutralomegaslope[i] = 1;
				}	
				else {
					neutralomegaslope[i] = 0;
				}
			}
		}
	}	
	 
	
	void DeleteNeutralOmegaSlope() {
		delete neutralomegaslope;
	}
	
	
	void CreateOmegaSlope() {
		omegaslope[0] = 1;
		for (int i=0; i<Ncont; i++) {
			omegaslope[i+1] = 0;
		}
	}
	
	void DeleteOmegaSlope() {
		delete omegaslope;
	}
	
	
	void CreateNeSlope() {
		if (!noadapt) {
			string cha("piS");
			Neslope[0] = 0;
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					Neslope[i+1] = 1;
				}	
				else {
					Neslope[i+1] = 0;
				}
			}
		}
		else {				
			string cha("piS");
			for (int i=0; i<Ncont; i++) {
				if (GetContinuousData()->GetCharacterName(i)==cha) {
					Neslope[i] = 1;
				}	
				else {
					Neslope[i] = 0;
				}
			}
		}
	}	
		
	void DeleteNeSlope() {
		delete Neslope;
	}
	

		
	// accessors
	Tree* GetTree() {return tree;}
	
	Var<Real>* GetGamma() {return gamma;}
	
	Const<PosReal>* GetK() {return K;}

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
	
	
	double* GetNeutralOmegaSlope() {return neutralomegaslope;}
	double* GetUSlope() {return uslope;}
	double* GetSynrateSlope() {return synrateslope;}
	double* GetNeSlope() {return Neslope;}


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
		if (chronoprior)	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		total += gamma->GetLogProb();
		total += beta->GetLogProb();
		if (!sameseq) {total += beta2->GetLogProb();}
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
			
			scheduler.Register(new SimpleMove(gamma,0.003),10,"gamma");
			scheduler.Register(new SimpleMove(gamma,0.0008),10,"gamma");
			scheduler.Register(new SimpleMove(gamma,0.0003),10,"gamma");
			
			scheduler.Register(new SimpleMove(beta,0.15),10,"beta");
			scheduler.Register(new SimpleMove(beta,0.08),10,"beta");
			scheduler.Register(new SimpleMove(beta,0.01),10,"beta");
			
			if (!sameseq) {
				scheduler.Register(new SimpleMove(beta2,0.05),10,"beta2");
				scheduler.Register(new SimpleMove(beta2,0.008),10,"beta2");
				scheduler.Register(new SimpleMove(beta2,0.001),10,"beta2");

				scheduler.Register(new GammaBetaMove(gamma,beta2,1,11),10,"gammabeta2");
				scheduler.Register(new GammaBetaMove(gamma,beta2,0.1,11),10,"gammabeta2");
				scheduler.Register(new GammaBetaMove(gamma,beta2,0.01,11),10,"gammabeta2");
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
		gamma->Sample();
		beta->Sample();
		if (!sameseq) {beta2->Sample();}
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
		if (!sameseq) {os << "\tbeta2";}
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
			os << "\tdelta\tkappa";
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
		os << '\t' << gamma->val();
		os << '\t' << beta->val();
		if (!sameseq) {os << '\t' << beta2->val();}
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
		os << *gamma << '\n';
		os << *beta << '\n';
		if (!sameseq) {os << *beta2 << '\n';}
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
		is >> *gamma;
		is >> *beta;
		if (!sameseq) {is >> *beta2;}
		is >> *relrate;
		is >> *stationary;
	}
};

#endif
