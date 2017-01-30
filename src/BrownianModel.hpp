#ifndef BROWNIANMODEL_H
#define	BROWNIANMODEL_H


#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "GTRModel.h"
#include "BrownianProcess.h"
#include "CommutativeBrownianProcess.h"
#include "EvolutionNucBrownianProcess.h"
#include "EvolutionGCBrownianProcess.h"
#include "BrownianMove.h"
#include "ContinuousData.h"
#include "CodonSequenceAlignment.h"
#include "MGCodonTransitionMatrix.h"
#include "Chronogram.h"
#include "GeneralConjugatePath.h"
#include "ConjugateInverseWishart.h"
#include "Jeffreys.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "GeneralTransitionConjugatePath.h"

#include "EvolutionSegmentedNucProcess.h"
#include "EvolutionSegmentedGCProcess.h"
#include "EvolutionSegmentedCompoProcess.h"
#include "RandomSegmentedBranchSitePath.h"
#include "GeneralSegmentedConjugatePath.h"
#include "EvolutionSegmentedOmegaProcess.h"
#include "EvolutionOmegaBrownianProcess.h"

class BrownianModel : public ProbModel
{
	public:
		BrownianModel(string datafile, string treefile, string contdatafile, bool inpathconjugate, bool inmapSegment, int innSegments, Segmentation insegm, bool incommut, int ingc, bool inomega, int inconjsigma, bool infixtimes, bool inparal, bool sample = true, GeneticCodeType type=Universal) 	{

			pathconjugate = inpathconjugate;
			mapSegment = inmapSegment;
			conjsigma = inconjsigma;
			nSegments = innSegments;
			segm = insegm;
			commut = incommut;
			compo = false;
			gc = false;
			if (ingc == 2)	{
				compo = true;
				cerr << "compo model\n";
			}
			else if (ingc == 1)	{
				gc = true;
			}
			// gc = ingc;
			omega = inomega;
			fixtimes = infixtimes;
			paral = inparal;
  

			if(paral) {
				#ifdef _OPENMP
					#pragma omp parallel for
					for (int i=0; i<omp_get_num_threads(); i++)	{
						cerr << "Parallel computing enabled" << endl;
						cerr << "proc: " << omp_get_thread_num() << '\t' << omp_get_num_threads() << '\n';
					}
				#else
					cerr << "Parallel computing not enabled. Please check OpenMP is in your makefile options" << endl;
					paral = false;
				#endif
			}

			cerr << "nseg : " << nSegments << '\n';
			BrownianBridge::setNTreeSegments(nSegments);
			BrownianBridge::setSegmentation(segm);

			// fetch data from file
			data = new FileSequenceAlignment(datafile);
			Nstate = data->GetNstate();	// # states (20 for amino acids)

			if (omega)	{
				if (Nstate != Nnuc)	{
					cerr << "error : omega model only on nucleotide codon sequence alignments\n";
					exit(1);
				}
				codondata = new CodonSequenceAlignment(data, true, type);
			}
			else {
				codondata = 0;
			}
			Nsite = data->GetNsite();	// # columns
			codonStateSpace = new CodonStateSpace(type);

			L = 1;
			if(gc)	{
				L++;
			}
			else if (compo)	{
				L+= Nstate -1;
			}
			if(omega)	{
				L++;
			}

			taxonset = data->GetTaxonSet();

			// get tree from file (newick format)
			tree = new Tree(treefile);

			// check whether tree and data fits together
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

			// ----------
			// construction of the graph
			// ----------

			Zero = new Const<Real>(0);
			One = new Const<PosReal>(1);

			chronogram = new Chronogram(tree,One);


		  
			// create an array of positive variables kappa_i, i=1..Ncont + L
			double mindiag = 0.001;
			double maxdiag = 1000;
			DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
			DiagArray->setval(1.0);
			// create a diagonal matrix, with the kappa_i along the diagonal
			sigmaZero = new SigmaZero(DiagArray);

			// create covariance matrix
			// from an inverse wishart of parameter sigma0
			cerr << "sigma\n";
			if(conjsigma)	{
				conjugatesigma = new ConjugateInverseWishart(sigmaZero, Ncont + L);
				sigma = conjugatesigma;
			}
			else	{
				sigma = new InverseWishartMatrix(sigmaZero, Ncont + L);
			}

			sigma->SetAtSigmaZero();

			/*
			ifstream sis("cov0");
			sis >> *sigma;
			cout << *sigma << '\n';
			sigma->Clamp();
			*/

			//brownian process
			cerr << "process\n";
			if (conjsigma == 2)	{
				conjugateBrownianProcess = new ConjugateBrownianProcess(chronogram,conjugatesigma);
				brownianProcess = conjugateBrownianProcess;
			}
			else	{
				brownianProcess = new BrownianProcess(chronogram,sigma);
			}
			// condition the multivariate process
			// on the matrix of quantitative traits.
			// note the offset here : first trait corresponds to entry L+1 of the process, etc.
			// this is because the first L entries of the process correspond to the substitution variables (dS, dN/dS)
			for (int i=0; i<Ncont; i++)	{
				brownianProcess->GetInstantProcess()->SetAndClamp(contdata,L+i,i);
			}

			// just for numerical stability of the starting point
			for (int l=0; l<L; l++)	{
					brownianProcess->GetInstantProcess()->CutOff(1,l);
			}

			// substitution matrix

			cerr << "SUBS PROCESS\n";
			relrate = new Dirichlet(Nstate*(Nstate-1)/2);
			// if ((Nstate == Naa) && (rrtype == "lg"))	{
			if (Nstate == Naa)	{
				double total = 0;
				for (int i=0; i<Nstate*(Nstate-1)/2; i++)	{
					(*relrate)[i] = LG_RR[i];
					total += (*relrate)[i];
				}
				for (int i=0; i<Nstate*(Nstate-1)/2; i++)	{
					(*relrate)[i] /= total;
				}
				relrate->Clamp();
			}
			else	{
				relrate->setuniform();
			}

			stationary = new Dirichlet(Nstate);
			stationary->setuniform();

			if(commut) {
				matrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary, true);
				commutativeProcess = new CommutativeBrownianProcess(brownianProcess);
			}
			else if (compo)	{
				cerr << "change of basis\n";
				changeofbasis = new SumConstrainedMapping(Nstate);
				if(mapSegment) {
					cerr << "evolution segmented process\n";
					evolutionSegmentedProcess = new EvolutionSegmentedCompoProcess(brownianProcess, relrate, stationary, changeofbasis, 1, paral);
					cerr << "esp ok\n";
				}
				else {
					// ??
					cerr << "compo: only under mapSegment option\n";
					exit(1);
					// evolutionProcess = new EvolutionGCBrownianProcess(brownianProcess, relrate, rootGCrate, 1);
				}
			}
			else if(gc) {
				rootGCrate = new Beta(One, One);
				if(mapSegment) {
					 evolutionSegmentedProcess = new EvolutionSegmentedGCProcess(brownianProcess, relrate, rootGCrate, paral);
				}
				else {
						evolutionProcess = new EvolutionGCBrownianProcess(brownianProcess, relrate, rootGCrate, 1);
				}
			}
			else {
				matrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary, true);
				if(omega) {
					if(paral) {																																						  
						#pragma omp parallel																																				  
						{
						// int nThread = omp_get_num_threads();	  
						int nThread = 1;
							codonMatrix = new MGOmegaCodonTransitionMatrix_OMP(matrix, codonStateSpace, nThread);
							//codonMatrix = new MGOmegaCodonTransitionMatrix(matrix, codonStateSpace);

						}																																									  
					}
					else {		 
						 codonMatrix = new MGOmegaCodonTransitionMatrix(matrix, codonStateSpace);
					}
					if(mapSegment) {
						evolutionSegmentedProcess = new EvolutionSegmentedOmegaProcess(brownianProcess, codonMatrix, paral);
					}
					else {
						evolutionProcess = new EvolutionOmegaBrownianProcess(brownianProcess, matrix, codonStateSpace, 1);
					}
				}
				else {
					if(mapSegment) {
						evolutionSegmentedProcess = new EvolutionSegmentedNucProcess(brownianProcess, matrix, paral);
					}
					else {
						cerr << "simple model\n";
						cerr << "non commut\n";
						evolutionProcess = new EvolutionNucBrownianProcess(brownianProcess, matrix, -1, -1);
					}
				}
			}

			cerr << "PHYLO process\n";

			// a phylogenetic process
			if(commut) {
				if (pathconjugate)	{
					commutPathconjtree = new OneMatrixPathConjugateTree(commutativeProcess,matrix,data);
					phyloprocess = new PathConjugatePhyloProcess(commutPathconjtree);
				}
				else	{
					phyloprocess = new OneMatrixPhyloProcess(commutativeProcess,matrix,data);
				}
			}
			else if(mapSegment) {
				if(pathconjugate) {
					segmentedPathconjtree = new BranchMatrixSegmentedPathConjugateTree(chronogram, evolutionSegmentedProcess, data);
					phyloprocess = new SegmentedPathConjugatePhyloProcess(segmentedPathconjtree);
				}
				else {
					phyloprocess = new BrownianPhyloProcess(chronogram, data, evolutionSegmentedProcess);
				}
			}
			else {
				pathconjtree = new BranchMatrixTransitionPathConjugateTree(chronogram, evolutionProcess, data);
				phyloprocess = new TransitionPathConjugatePhyloProcess(pathconjtree);
			}

			cerr << "Mean of brownian process :" << GetMeanRho() << " - Var :" << GetVarRho() << endl;
			cerr << "unfold\n";
			phyloprocess->Unfold();

			if (sample)	{
					phyloprocess->Sample();
			}

			cerr << "register\n";
			RootRegister(Zero);
			RootRegister(One);
			RootRegister(relrate);
			RootRegister(stationary);

			Register();


			MakeScheduler();
			if(sample)
				Update();

			TraceHeader(cerr);
			Trace(cerr);

			cerr << "model created\n";
		}
		virtual ~BrownianModel();



		BrownianProcess* GetBrownianProcess()	{
			return brownianProcess;
		}

		Tree* GetTree()	{
			return tree;
		}
		int GetNcont() {
			return Ncont;
		}

		SequenceAlignment* GetData()	{
		if (codondata)	{
			return codondata;
		}
		return data;
	}

		CovMatrix* GetSigma() {
			return sigma;
		}

		double GetLogProb()	{
			return GetLogPrior() + GetLogLikelihood();
		}

		Chronogram* GetChronogram()	{
			return chronogram;
		}

		double GetLogPrior();
		double GetLogLikelihood();

		void drawSample();

		void MakeScheduler();
		//void ReadScheduler(string filename);

	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true, false);
		return 1;
	}
	*/

	double GetMeanRho();
		double GetIntegralRho();
	double GetVarRho();
		double GetMeanOmega();
	double GetVarOmega();
		double GetMeanGC();
	double GetVarGC();
		double GetRootGC();


	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os);

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os);
	void Details(ostream& os);

	void ToStream(ostream& os);
	void FromStream(istream& is);

	void ToStreamShort(ostream& os);
	void FromStreamShort(istream& is);

	void PrintEntries(ostream& os);

	void test() {
		  brownianProcess->GetPureBrownianProcess()->CheckLength();
	}

	void PostPredSimu(string name, bool resamplecov, bool resampleprocess, int nsegment);

	void Load(string simufile);

	protected:
	private:

		// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* data;
		CodonSequenceAlignment* codondata;
		CodonStateSpace *codonStateSpace;
	TaxonSet* taxonset;
		ContinuousData* contdata;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;
		//Number of continuous characters
		int Ncont;
		//Number of evolution characters
		int L;
	// ---------
	// the random variables of the model
	// ---------
	Const<Real>* Zero;
	Const<PosReal>* One;

	// chronogram
	Chronogram* chronogram;

	// autocorrelated process
	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	InverseWishartMatrix* sigma;
	ConjugateInverseWishart* conjugatesigma;

	BrownianProcess* brownianProcess;
	ConjugateBrownianProcess* conjugateBrownianProcess;
	int nSegments;	//Number of segments in each brownian bridge
		Segmentation segm;
		bool commut;
		int gc;
		bool omega;
		CommutativeBrownianProcess *commutativeProcess;
		BranchValPtrTree<RandomTransitionMatrix>* evolutionProcess;
		EvolutionSegmentedProcess* evolutionSegmentedProcess;

	// compo model
	bool compo;
	SumConstrainedMapping* changeofbasis;

	// substitution matrix is relrate * stationary
	Dirichlet* relrate;
	Dirichlet* stationary;
		Rvar<UnitReal> *rootGCrate;
	GTRRandomSubMatrixWithNormRates* matrix;
		MGOmegaCodonTransitionMatrix* codonMatrix;

	// phylo process
	PathConjugateTree* commutPathconjtree;
	TransitionPathConjugateTree* pathconjtree;
	SegmentedPathConjugateTree* segmentedPathconjtree;
	PhyloProcess* phyloprocess;


	// OneMatrixPhyloProcess* phyloprocess;
	bool pathconjugate;
		bool mapSegment;
		int conjsigma;
		bool fixtimes;
		bool paral;


};


#endif	/* BROWNIANMODEL_H */

