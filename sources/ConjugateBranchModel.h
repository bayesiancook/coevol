
#ifndef BRANCHOMEGAMULTI
#define BRANCHOMEGAMULTI

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "MG3OmegaCodonSubMatrix.h"
#include "MGOmega3CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"

#include "GCProcess.h"

#include "GeneralConjugatePath.h"
#include "CodonConjugatePath.h"

#include "Jeffreys.h"
#include "MatrixTree.h"

class RootConstrainedGammaTree : public GammaTree	{

	public:

	RootConstrainedGammaTree(Tree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, bool inwithroot = false)  : GammaTree(intree) 	{
		SetWithRoot(inwithroot);
		alpha = inalpha;
		beta = inbeta;
		rootbranchval = new Gamma(alpha,beta);
		rootval = 0;
		if (WithRoot())	{
			rootval = new Gamma(alpha,beta);
			SetBranchVal(0,rootval);
		}
		for (const Link* from=GetRoot()->Next(); from!=GetRoot(); from=from->Next())	{
			SetBranchVal(from->GetBranch(),rootbranchval);
			RecursiveCreate(from->Out());
		}
	}

	~RootConstrainedGammaTree()	{
		for (const Link* from=GetRoot()->Next(); from!=GetRoot(); from=from->Next())	{
			RecursiveDelete(from->Out());
		}
		delete rootval;
		delete rootbranchval;
	}

	protected:

	Gamma* rootbranchval;
	Gamma* rootval;

};

class RootConstrainedBetaTree : public BetaTree  {

	public:

	RootConstrainedBetaTree(Tree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, bool inwithroot = false)  : BetaTree(intree) {
		SetWithRoot(inwithroot);
		alpha = inalpha;
		beta = inbeta;
		rootbranchval = new Beta(alpha,beta);
		rootval = 0;
		if (WithRoot())	{
			rootval = new Beta(alpha,beta);
			SetBranchVal(0,rootval);
		}
		for (const Link* from=GetRoot()->Next(); from!=GetRoot(); from=from->Next())	{
			SetBranchVal(from->GetBranch(),rootbranchval);
			RecursiveCreate(from->Out());
		}
	}

	~RootConstrainedBetaTree()	{
		for (const Link* from=GetRoot()->Next(); from!=GetRoot(); from=from->Next())	{
			RecursiveDelete(from->Out());
		}
		delete rootval;
		delete rootbranchval;
	}

	protected:

	Beta* rootbranchval;
	Beta* rootval;

};


class UGam : public Rvar<PosReal>	{

	public:

	UGam(Var<PosReal>* intime, Var<PosReal>* inalpha, Var<PosReal>* inbeta)	{
		// mean = alpha / beta
		// var = alpha / beta^2

		// mean' = time * alpha / beta
		// var' = time^2 * alpha / beta^2
		// thus beta' = beta / time
		time = intime;
		alpha = inalpha;
		beta = inbeta;
		Register(time);
		Register(alpha);
		Register(beta);
		Sample();
	}

	protected:

	void drawSample()	{
		double v = Random::Gamma(alpha->val(),beta->val() / time->val());
		if (v < 1e-20)	{
			v = 1e-20;
		}
		setval(v);
	}

	double logProb()	{
		return -Random::logGamma(alpha->val()) + (alpha->val())*log(beta->val()/time->val()) - *this * (beta->val()/time->val()) - (1 - (alpha->val()))*log(*this);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	Var<PosReal>* time;

};

class UGamTree : public BranchProcess<PosReal> {

	public:

	// UGamTree(Tree* intree)	: BranchProcess<PosReal>(intree) {}

	UGamTree(LengthTree* intree, Var<PosReal>* inalpha, Var<PosReal>* inbeta, bool inwithroot = false)  : BranchProcess<PosReal>(intree->GetTree()) {
		SetWithRoot(inwithroot);
		alpha = inalpha;
		beta = inbeta;
		lengthtree = intree;
		RecursiveCreate(GetRoot());
	}

	~UGamTree()	{
		RecursiveDelete(GetRoot());
	}

	double GetTotal()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total;

	}

	double GetMean()	{
		int n = 0;
		double total = GetMean(GetRoot(),n);
		return total / n;

	}

	double GetVar()	{
		int n = 0;
		double mean = GetMean(GetRoot(),n);
		n = 0;
		double meansquare = GetMeanSquare(GetRoot(),n);
		mean /= n;
		meansquare /= n;
		return meansquare - mean * mean;
	}

	protected:

	double GetMean(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += this->GetBranchVal(link->GetBranch())->val();
			n++;
			total += GetMean(link->Out(),n);
		}
		return total;
	}

	double GetMeanSquare(const Link* from, int& n)	{
		double total = 0;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = this->GetBranchVal(link->GetBranch())->val();
			total += tmp * tmp;
			n++;
			total += GetMeanSquare(link->Out(),n);
		}
		return total;
	}

	Rvar<PosReal>* CreateBranchVal(const Link* link)	{
		return new UGam(lengthtree->GetBranchVal(link->GetBranch()),alpha,beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
	LengthTree* lengthtree;

};

class ConjugateOmegaTreeMove : public MCUpdate, public Mnode {

	GammaTree* tree;
	MGCodonBranchMatrixPathConjugateTree* pathconjtree;
	Var<PosReal>* alpha;
	Var<PosReal>* beta;

	public:

	const Link* GetRoot() {return tree->GetRoot();}

	ConjugateOmegaTreeMove(GammaTree* intree, MGCodonBranchMatrixPathConjugateTree* inpathconjtree)	{
		tree = intree;
		pathconjtree = inpathconjtree;
		alpha = tree->GetAlpha();
		beta = tree->GetBeta();
		tree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(false);
		pathconjtree->ComputeTotSuffStat();
		LocalResampleOmega(GetRoot()->Next());
		for (const Link* from=GetRoot()->Next(); from!=GetRoot(); from=from->Next())	{
			RecursiveResampleOmega(from->Out());
		}
		Update();
		return 1;
	}

	void RecursiveResampleOmega(const Link* from)	{
		/*
		if ((! from->isRoot())	{
			LocalResampleOmega(from);
		}
		*/

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			LocalResampleOmega(link);
			RecursiveResampleOmega(link->Out());
		}
	}

	void LocalResampleOmega(const Link* from)	{
		double currentomega = tree->GetBranchVal(from->GetBranch())->val();
		double tmp = Random::Gamma(alpha->val() + pathconjtree->GetBranchMGCodonPathConjugate(from->GetBranch())->GetTotNonSynCount(), beta->val() + pathconjtree->GetBranchMGCodonPathConjugate(from->GetBranch())->GetTotNonSynBeta() / currentomega);
		tree->GetBranchVal(from->GetBranch())->setval(tmp);

	}

};

class BranchModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* nucdata;
	CodonSequenceAlignment* codondata;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;
	Const<PosReal>* PriorMu;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	int chronoprior;
	double meanchi;
	double meanchi2;
	bool iscalib;
	// 0 : uniform;
	// 1 : bd;

	// chronogram
	Chronogram* chronogram;

	GammaTree* syngammatree;
	LengthTree* lengthtree;

	// Jeffreys* sigma;
	Jeffreys* sigmaalpha;
	Jeffreys* sigmabeta;
	LengthTree* synratetree;
	// LogNormalTreeProcess* logntree;
	UGamTree* gamsynratetree;

	Jeffreys* gcalpha;
	Jeffreys* gcbeta;
	BetaTree* gctree;

	Jeffreys* omegaalpha;
	Jeffreys* omegabeta;
	GammaTree* omegatree;

	Jeffreys* omegatsalpha;
	Jeffreys* omegatsbeta;
	GammaTree* omegatstree;
	Jeffreys* omegatv0alpha;
	Jeffreys* omegatv0beta;
	GammaTree* omegatv0tree;
	Jeffreys* omegatvgcalpha;
	Jeffreys* omegatvgcbeta;
	GammaTree* omegatvgctree;

	GCStatTree* stattree;
	NucMatrixTree* nucmatrixtree;

	BetaTree* gctree1;
	BetaTree* gctree2;
	BetaTree* gctree3;
	GCStatTree* stattree1;
	GCStatTree* stattree2;
	GCStatTree* stattree3;
	NucMatrixTree* nucmatrixtree1;
	NucMatrixTree* nucmatrixtree2;
	NucMatrixTree* nucmatrixtree3;

	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	// for both
	BranchValPtrTree<RandomSubMatrix>* matrixtree;
	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;
	// BranchMatrixPhyloProcess* phyloprocess;

	double mappingfreq;

	bool clamptree;

	int gc;
	int codonmodel;

	// total number of substitution parameters modelled as non homogeneous
	// int L;

	bool conjpath;
	bool normalise;
	int nrep;
	int gcindex;

	int clampsuffstat;
	string suffstatfile;

	bool fullconj;

	public:

	bool Unconstrained()	{
		return (syngammatree != 0);
	}

	SequenceAlignment* GetData()	{
		if (codonmodel)	{
			return codondata;
		}
		return nucdata;
	}

	bool isCalibrated()	{
		return iscalib;
	}

	BranchModel(string datafile, string treefile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, int ingc, int incodonmodel, int inconjpath, bool infullconj, double inmappingfreq, bool inclamptree, bool innormalise, int innrep, string insuffstatfile, bool sample, GeneticCodeType type)	{


		// analyse options

		mappingfreq = inmappingfreq;
		fullconj = infullconj;

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clamptree = inclamptree;

		codonmodel = incodonmodel;
		gc = ingc;
		if (gc == 3)	{
			if (! codonmodel)	{
				cerr << "error : can apply gc3 formalism only under codon model\n";
				cerr << "to activate codon model, use -dsom option\n";
				exit(1);
			}
		}

		if (inconjpath == -1)	{
			conjpath = true;
		}
		else	{
			conjpath = inconjpath;
		}
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 30 : 1;
		}
		normalise = innormalise;


		// get data from file
		nucdata = new FileSequenceAlignment(datafile);

		if (codonmodel)	{
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		cerr << "tree and data ok\n";
		cerr << '\n';

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		syngammatree = 0;
		iscalib = false;
		if (calibfile == "Unconstrained")	{
			PriorMu = new Const<PosReal>(10);
			syngammatree = new GammaTree(tree,One,PriorMu);
			lengthtree = syngammatree;
			synratetree = syngammatree;
		}
		else	{
			PriorMu = new Const<PosReal>(1);
			if (calibfile != "None")	{

				if (calibfile != "None")	{
					iscalib = true;
					double a = rootage * rootage / rootstdev / rootstdev;
					double b = rootage / rootstdev / rootstdev;
					CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

					if (chronoprior == 0)	{
						chronogram = new CalibratedChronogram(tree,PriorMu,a,b,calibset);
					}
					else {
						cerr << "BD\n";
						MeanChi = new Const<PosReal>(meanchi);
						MeanChi2 = new Const<PosReal>(meanchi2);
						Chi = new Exponential(MeanChi,Exponential::MEAN);
						Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
						chronogram = new BDCalibratedChronogram(tree,PriorMu,Chi,Chi2,a,b,calibset,chronoprior);
					}
				}
				else	{
					cerr << "new chrono\n";
					chronogram = new Chronogram(tree,PriorMu);
					cerr << "ok\n";
				}
			}
			else	{
				cerr << "new chrono\n";
				chronogram = new Chronogram(tree,One);
				cerr << "ok\n";
			}
			/*
			sigma = new Jeffreys(0.001,1000,Zero);
			sigma->setval(1);
			*/
			sigmaalpha = new Jeffreys(0.001,1000,Zero);
			sigmaalpha->setval(1);
			sigmabeta = new Jeffreys(0.001,1000,Zero);
			sigmabeta->setval(1);
			gamsynratetree = new UGamTree(chronogram,sigmaalpha,sigmabeta);
			synratetree = gamsynratetree;
			// logntree = new LogNormalTreeProcess(chronogram,sigma,INTEGRAL);
			// synratetree = logntree;
			/*
			cerr << "??? in branch process synratetree \n";
			exit(1);
			*/
			lengthtree = chronogram;
		}

		if (clamptree)	{
			chronogram->Clamp();
		}


		CreateSubstitutionProcess();

		if (phyloprocess)	{
			phyloprocess->Unfold();
		}
		if (sample)	{
			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}

		cerr << "root register\n";
		// register model
		RootRegister(Zero);
		RootRegister(One);
		RootRegister(PriorMu);
		if (! Unconstrained())	{
			if (chronoprior == 1)	{
				RootRegister(MeanChi);
				RootRegister(MeanChi2);
			}
		}
		RootRegister(relrate);
		if (!gc)	{
			RootRegister(stationary);
		}
		Register();

		cerr << "scheduler\n";
		MakeScheduler();
		if (sample)	{
			cerr << "update\n";
			Update();
			cerr << "ok\n";
		}
	}

	~BranchModel() {}

	void CreateSubstitutionProcess()	{

		phyloprocess = 0;

		nucmatrix = 0;
		stationary = 0;

		omegaalpha = 0;
		omegabeta = 0;
		omegatree = 0;

		gcalpha = 0;
		gcbeta = 0;

		gctree = 0;
		stattree = 0;
		nucmatrixtree = 0;

		gctree1 = 0;
		gctree2 = 0;
		gctree3 = 0;
		stattree1 = 0;
		stattree2 = 0;
		stattree3 = 0;
		nucmatrixtree1 = 0;
		nucmatrixtree2 = 0;
		nucmatrixtree3 = 0;



		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);

		if (codonmodel == 3)	{

			cerr << "omega trees\n";
			omegatsalpha = new Jeffreys(0.001,1000,Zero);
			omegatsbeta = new Jeffreys(0.001,1000,Zero);
			omegatsalpha->setval(10);
			omegatsbeta->setval(100);
			omegatstree = new RootConstrainedGammaTree(tree,omegatsalpha,omegatsbeta,false);

			omegatv0alpha = new Jeffreys(0.001,1000,Zero);
			omegatv0beta = new Jeffreys(0.001,1000,Zero);
			omegatv0alpha->setval(10);
			omegatv0beta->setval(100);
			omegatv0tree = new RootConstrainedGammaTree(tree,omegatv0alpha,omegatv0beta,false);

			omegatvgcalpha = new Jeffreys(0.001,1000,Zero);
			omegatvgcbeta = new Jeffreys(0.001,1000,Zero);
			omegatvgcalpha->setval(10);
			omegatvgcbeta->setval(100);
			omegatvgctree = new RootConstrainedGammaTree(tree,omegatvgcalpha,omegatvgcbeta,false);

			cerr << "ok\n";

			if (gc)	{
				cerr << "error : gc not yet compatible with omega3\n";
				exit(1);
			}
			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
			cerr << "matrix tree\n";
			matrixtree = new Omega3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatstree, omegatv0tree, omegatvgctree, One);
			cerr << "ok\n";
			// make substitution mappings
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				if (clampsuffstat)	{
					cerr << "read suffstat\n";
					pathconjtree->ReadFromFile(suffstatfile);
					cerr << "ok\n";
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else	{
				pathconjtree = 0;
				phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, codondata);
			}
		}
		else if (codonmodel)	{

			omegaalpha = new Jeffreys(0.001,1000,Zero);
			omegabeta = new Jeffreys(0.001,1000,Zero);
			omegaalpha->setval(10);
			omegabeta->setval(100);
			omegatree = new RootConstrainedGammaTree(tree,omegaalpha,omegabeta,false);

			if (gc == 3)	{
				gcindex = 1;
				gcalpha = new Jeffreys(0.001,1000,Zero);
				gcbeta = new Jeffreys(0.001,1000,Zero);
				gcalpha->setval(10);
				gcbeta->setval(10);
				gctree1 = new RootConstrainedBetaTree(tree,gcalpha,gcbeta,true);
				gctree2 = new RootConstrainedBetaTree(tree,gcalpha,gcbeta,true);
				gctree3 = new RootConstrainedBetaTree(tree,gcalpha,gcbeta,true);
				stattree1 = new GCStatTree(gctree1,gctree1->GetBranchVal(0));
				stattree2 = new GCStatTree(gctree2,gctree2->GetBranchVal(0));
				stattree3 = new GCStatTree(gctree3,gctree3->GetBranchVal(0));
				nucmatrixtree1 = new GTRGCNucMatrixTree(relrate,stattree1,normalise);
				nucmatrixtree2 = new GTRGCNucMatrixTree(relrate,stattree2,normalise);
				nucmatrixtree3 = new GTRGCNucMatrixTree(relrate,stattree3,normalise);
				matrixtree = new GC3MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree1, nucmatrixtree2, nucmatrixtree3, omegatree, One);
			}
			else if (gc == 1)	{
				gcindex = 1;
				gcalpha = new Jeffreys(0.001,1000,Zero);
				gcbeta = new Jeffreys(0.001,1000,Zero);
				gcalpha->setval(10);
				gcbeta->setval(10);
				gctree = new RootConstrainedBetaTree(tree,gcalpha,gcbeta,true);
				stattree = new GCStatTree(gctree,gctree->GetBranchVal(0));
				nucmatrixtree = new GTRGCNucMatrixTree(relrate,stattree,normalise);
				matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);
			}
			else	{
				stationary = new Dirichlet(Nnuc);
				nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
				matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);
			}
			// make substitution mappings
			if (conjpath)	{
				if (gc)	{
					pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				}
				else	{
					pathconjtree = new MGCodonBranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
					// pathconjtree = new BranchMatrixPathConjugateTree(synratetree, matrixtree, codondata);
				}
				if (clampsuffstat)	{
					cerr << "read suffstat\n";
					pathconjtree->ReadFromFile(suffstatfile);
					cerr << "ok\n";
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else	{
				pathconjtree = 0;
				phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, codondata);
			}
		}
		else	{

			if (gc)	{
				gcindex = 1;
				gcalpha = new Jeffreys(0.001,1000,Zero);
				gcbeta = new Jeffreys(0.001,1000,Zero);
				gcalpha->setval(10);
				gcbeta->setval(10);
				gctree = new RootConstrainedBetaTree(tree,gcalpha,gcbeta,true);
				stattree = new GCStatTree(gctree,gctree->GetBranchVal(0));
				nucmatrixtree = new GTRGCNucMatrixTree(relrate,stattree,normalise);
			}
			else	{
				stationary = new Dirichlet(Nnuc);
				nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
			}
			// make substitution mappings
			if (conjpath)	{
				if (gc)	{
					pathconjtree = new BranchMatrixPathConjugateTree(synratetree, nucmatrixtree, nucdata);
				}
				else	{
					pathconjtree = new OneMatrixPathConjugateTree(synratetree,nucmatrix,nucdata);
				}
				if (clampsuffstat)	{
					cerr << "create suffstat\n";
					pathconjtree->ReadFromFile(suffstatfile);
					cerr << "ok\n";
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
				}
			}
			else	{
				pathconjtree = 0;
				if (gc)	{
					phyloprocess = new BranchMatrixPhyloProcess(synratetree, nucmatrixtree, nucdata);
				}
				else	{
					phyloprocess = new OneMatrixPhyloProcess(synratetree, nucmatrix, nucdata);
				}
			}
		}

		if (mappingfreq == -1)	{
			if (codonmodel)	{
				mappingfreq = 0.2;
			}
			else	{
				mappingfreq = 1;
			}
		}
	}

	Tree* GetTree() {return tree;}

	BranchVarTree<PosReal>* GetSynRateTree() {
		return synratetree;
	}

	GammaTree* GetSynGammaTree() {
		return syngammatree;
	}

	BetaTree* GetGCTree() {return gctree;}
	BetaTree* GetGCTree1() {return gctree1;}
	BetaTree* GetGCTree2() {return gctree2;}
	BetaTree* GetGCTree3() {return gctree3;}
	GammaTree* GetOmegaTree() {
		if (! Codon())	{
			cerr << "error : get omegatree not defined\n";
			exit(1);
		}
		return omegatree;
	}

	GammaTree* GetOmegaTsTree() {
		if (! Codon3())	{
			cerr << "error : get omegatree not defined\n";
			exit(1);
		}
		return omegatstree;
	}

	GammaTree* GetOmegaTv0Tree() {
		if (! Codon3())	{
			cerr << "error : get omegatree not defined\n";
			exit(1);
		}
		return omegatv0tree;
	}

	GammaTree* GetOmegaTvGCTree() {
		if (! Codon3())	{
			cerr << "error : get omegatree not defined\n";
			exit(1);
		}
		return omegatvgctree;
	}

	bool Codon() {return (codonmodel==1);}
	bool Codon3() {return (codonmodel == 3);}
	bool isGCActivated() {return (gc == 1);}
	bool isGC3Activated() {return (gc == 3);}

	Chronogram* GetChronogram() {return chronogram;}
	LengthTree* GetLengthTree() {return lengthtree;}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		if (Unconstrained())	{
			total += syngammatree->GetLogProb();
		}
		else	{
			if (chronoprior)	{
				total += Chi->GetLogProb();
				total += Chi2->GetLogProb();
			}
			total += chronogram->GetLogProb();
			/*
			total += sigma->GetLogProb();
			total += GetSynRateTree()->GetLogProb();
			*/
			total += sigmaalpha->GetLogProb();
			total += sigmabeta->GetLogProb();
			total += gamsynratetree->GetLogProb();
		}

		total += relrate->GetLogProb();

		if (codonmodel == 3)	{
			total += omegatsalpha->GetLogProb();
			total += omegatsbeta->GetLogProb();
			total += omegatstree->GetLogProb();
			total += omegatv0alpha->GetLogProb();
			total += omegatv0beta->GetLogProb();
			total += omegatv0tree->GetLogProb();
			total += omegatvgcalpha->GetLogProb();
			total += omegatvgcbeta->GetLogProb();
			total += omegatvgctree->GetLogProb();
		}
		else if (codonmodel)	{
			total += omegaalpha->GetLogProb();
			total += omegabeta->GetLogProb();
			total += omegatree->GetLogProb();
		}

		if (gc == 3)	{
			total += gcalpha->GetLogProb();
			total += gcbeta->GetLogProb();
			total += gctree1->GetLogProb();
			total += gctree2->GetLogProb();
			total += gctree3->GetLogProb();
		}
		else if (gc)	{
			total += gcalpha->GetLogProb();
			total += gcbeta->GetLogProb();
			total += gctree->GetLogProb();
		}
		else	{
			total += stationary->GetLogProb();
		}
		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		if (clampsuffstat)	{
			ret = pathconjtree->GetLogProb();
		}
		else	{
			ret = phyloprocess->GetLogProb();
		}
		return ret;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			if (! clampsuffstat)	{
				scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree,mappingfreq),1,"mapping + sufficient stat");
				// scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
			}
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		for (int i=0; i<nrep; i++)	{
			if (Unconstrained())	{
				scheduler.Register(new SimpleMove(syngammatree,1),10,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.1),10,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.01),10,"syngamtree");
			}
			else if (! clamptree)	{
				if (chronoprior)	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

				if (isCalibrated())	{
					scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
					scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
					scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
				}
			}

			/*
			scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");

			scheduler.Register(new SimpleMove(logntree,10),10,"lognormal");
			scheduler.Register(new SimpleMove(logntree,1),10,"lognormal");
			scheduler.Register(new SimpleMove(logntree,0.1),10,"lognormal");
			scheduler.Register(new SimpleMove(logntree,0.01),10,"lognormal");
			*/

			scheduler.Register(new SimpleMove(sigmaalpha,10),100,"sigmaalpha");
			scheduler.Register(new SimpleMove(sigmaalpha,1),100,"sigmaalpha");
			scheduler.Register(new SimpleMove(sigmaalpha,0.1),100,"sigmaalpha");

			scheduler.Register(new SimpleMove(sigmabeta,10),100,"sigmabeta");
			scheduler.Register(new SimpleMove(sigmabeta,1),100,"sigmabeta");
			scheduler.Register(new SimpleMove(sigmabeta,0.1),100,"sigmabeta");

			scheduler.Register(new SimpleMove(gamsynratetree,10),10,"ugam");
			scheduler.Register(new SimpleMove(gamsynratetree,1),10,"ugam");
			scheduler.Register(new SimpleMove(gamsynratetree,0.1),10,"ugam");
			scheduler.Register(new SimpleMove(gamsynratetree,0.01),10,"ugam");

			if (gc == 3)	{
				scheduler.Register(new SimpleMove(gcalpha,10),100,"gcalpha");
				scheduler.Register(new SimpleMove(gcalpha,1),100,"gcalpha");
				scheduler.Register(new SimpleMove(gcalpha,0.1),100,"gcalpha");

				scheduler.Register(new SimpleMove(gcbeta,10),100,"gcbeta");
				scheduler.Register(new SimpleMove(gcbeta,1),100,"gcbeta");
				scheduler.Register(new SimpleMove(gcbeta,0.1),100,"gcbeta");

				scheduler.Register(new SimpleMove(gctree1,10),10,"gctree");
				scheduler.Register(new SimpleMove(gctree1,1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree1,0.1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree1,0.01),10,"gctree");

				scheduler.Register(new SimpleMove(gctree2,10),10,"gctree");
				scheduler.Register(new SimpleMove(gctree2,1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree2,0.1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree2,0.01),10,"gctree");

				scheduler.Register(new SimpleMove(gctree3,10),10,"gctree");
				scheduler.Register(new SimpleMove(gctree3,1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree3,0.1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree3,0.01),10,"gctree");
			}
			else if (gc)	{
				scheduler.Register(new SimpleMove(gcalpha,10),100,"gcalpha");
				scheduler.Register(new SimpleMove(gcalpha,1),100,"gcalpha");
				scheduler.Register(new SimpleMove(gcalpha,0.1),100,"gcalpha");

				scheduler.Register(new SimpleMove(gcbeta,10),100,"gcbeta");
				scheduler.Register(new SimpleMove(gcbeta,1),100,"gcbeta");
				scheduler.Register(new SimpleMove(gcbeta,0.1),100,"gcbeta");

				scheduler.Register(new SimpleMove(gctree,10),10,"gctree");
				scheduler.Register(new SimpleMove(gctree,1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree,0.1),10,"gctree");
				scheduler.Register(new SimpleMove(gctree,0.01),10,"gctree");
			}
			else	{
				scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
				scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
				scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
				scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
			}

			if (codonmodel == 3)	{
				scheduler.Register(new SimpleMove(omegatsalpha,10),100,"omalpha");
				scheduler.Register(new SimpleMove(omegatsalpha,1),100,"omalpha");
				scheduler.Register(new SimpleMove(omegatsalpha,0.1),100,"omalpha");

				scheduler.Register(new SimpleMove(omegatsbeta,10),100,"ombeta");
				scheduler.Register(new SimpleMove(omegatsbeta,1),100,"ombeta");
				scheduler.Register(new SimpleMove(omegatsbeta,0.1),100,"ombeta");

				scheduler.Register(new SimpleMove(omegatstree,10),10,"omtree");
				scheduler.Register(new SimpleMove(omegatstree,1),10,"omtree");
				scheduler.Register(new SimpleMove(omegatstree,0.1),10,"omtree");
				scheduler.Register(new SimpleMove(omegatstree,0.01),10,"omtree");

				scheduler.Register(new SimpleMove(omegatv0alpha,10),100,"omalpha");
				scheduler.Register(new SimpleMove(omegatv0alpha,1),100,"omalpha");
				scheduler.Register(new SimpleMove(omegatv0alpha,0.1),100,"omalpha");

				scheduler.Register(new SimpleMove(omegatv0beta,10),100,"ombeta");
				scheduler.Register(new SimpleMove(omegatv0beta,1),100,"ombeta");
				scheduler.Register(new SimpleMove(omegatv0beta,0.1),100,"ombeta");

				scheduler.Register(new SimpleMove(omegatv0tree,10),10,"omtree");
				scheduler.Register(new SimpleMove(omegatv0tree,1),10,"omtree");
				scheduler.Register(new SimpleMove(omegatv0tree,0.1),10,"omtree");
				scheduler.Register(new SimpleMove(omegatv0tree,0.01),10,"omtree");

				scheduler.Register(new SimpleMove(omegatvgcalpha,10),100,"omalpha");
				scheduler.Register(new SimpleMove(omegatvgcalpha,1),100,"omalpha");
				scheduler.Register(new SimpleMove(omegatvgcalpha,0.1),100,"omalpha");

				scheduler.Register(new SimpleMove(omegatvgcbeta,10),100,"ombeta");
				scheduler.Register(new SimpleMove(omegatvgcbeta,1),100,"ombeta");
				scheduler.Register(new SimpleMove(omegatvgcbeta,0.1),100,"ombeta");

				scheduler.Register(new SimpleMove(omegatvgctree,10),10,"omtree");
				scheduler.Register(new SimpleMove(omegatvgctree,1),10,"omtree");
				scheduler.Register(new SimpleMove(omegatvgctree,0.1),10,"omtree");
				scheduler.Register(new SimpleMove(omegatvgctree,0.01),10,"omtree");
			}
			else if (codonmodel)	{
				scheduler.Register(new SimpleMove(omegaalpha,10),100,"omalpha");
				scheduler.Register(new SimpleMove(omegaalpha,1),100,"omalpha");
				scheduler.Register(new SimpleMove(omegaalpha,0.1),100,"omalpha");

				scheduler.Register(new SimpleMove(omegabeta,10),100,"ombeta");
				scheduler.Register(new SimpleMove(omegabeta,1),100,"ombeta");
				scheduler.Register(new SimpleMove(omegabeta,0.1),100,"ombeta");

				if (fullconj)	{
					MGCodonBranchMatrixPathConjugateTree* tmp = dynamic_cast<MGCodonBranchMatrixPathConjugateTree*> (pathconjtree);
					if (! tmp)	{
						cerr << "error: conjugate move called on non conjugate pathtree\n";
						exit(1);
					}
					scheduler.Register(new ConjugateOmegaTreeMove(omegatree,tmp),1,"conjomtree");
				}
				else	{
					scheduler.Register(new SimpleMove(omegatree,10),10,"omtree");
					scheduler.Register(new SimpleMove(omegatree,1),10,"omtree");
					scheduler.Register(new SimpleMove(omegatree,0.1),10,"omtree");
					scheduler.Register(new SimpleMove(omegatree,0.01),10,"omtree");
				}
			}

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
		}
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";

		if (Unconstrained())	{
			syngammatree->Sample();
		}
		else	{
			if (chronoprior)	{
				Chi->Sample();
				Chi2->Sample();
			}
			chronogram->Sample();
			/*
			sigma->Sample();
			logntree->Sample();
			*/
			sigmaalpha->Sample();
			sigmabeta->Sample();
			gamsynratetree->Sample();
		}

		relrate->Sample();

		if (codonmodel == 3)	{
			omegatsalpha->Sample();
			omegatsbeta->Sample();
			omegatstree->Sample();
			omegatv0alpha->Sample();
			omegatv0beta->Sample();
			omegatv0tree->Sample();
			omegatvgcalpha->Sample();
			omegatvgcbeta->Sample();
			omegatvgctree->Sample();
		}
		else if (codonmodel)	{
			omegaalpha->Sample();
			omegabeta->Sample();
			omegatree->Sample();
		}

		if (gc == 3)	{
			gcalpha->Sample();
			gcbeta->Sample();
			gctree1->Sample();
			gctree2->Sample();
			gctree3->Sample();
		}
		else if (gc)	{
			gcalpha->Sample();
			gcbeta->Sample();
			gctree->Sample();
		}
		else	{
			stationary->Sample();
		}

		if (! clampsuffstat)	{
			phyloprocess->Sample();
		}

		cerr << "ok\n";
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	double GetMeanSynRate()	{
		if (Unconstrained())	{
			return syngammatree->GetMean();
		}
		return gamsynratetree->GetMean();
		// return logntree->GetMeanRate();
	}

	double GetTotalTime()	{
		if (Unconstrained())	{
			return syngammatree->GetTotal();
		}
		return chronogram->GetTotalTime();
	}

	double GetMeanOmega()	{
		if (codonmodel != 1)	{
			cerr << "error: get mean omega\n";
			exit(1);
		}
		return omegatree->GetMean();
	}

	double GetMeanOmegaTs()	{
		if (codonmodel != 3)	{
			cerr << "error: get mean omega\n";
			exit(1);
		}
		return omegatstree->GetMean();
	}

	double GetMeanOmegaTv0()	{
		if (codonmodel != 3)	{
			cerr << "error: get mean omega\n";
			exit(1);
		}
		return omegatv0tree->GetMean();
	}

	double GetMeanOmegaTvGC()	{
		if (codonmodel != 3)	{
			cerr << "error: get mean omega\n";
			exit(1);
		}
		return omegatvgctree->GetMean();
	}

	double GetMeanGC1()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree1->GetMean();
	}

	double GetMeanGC2()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree2->GetMean();
	}

	double GetMeanGC3()	{
		if (gc != 3)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree3->GetMean();
	}

	double GetMeanGC()	{
		if (! gc)	{
			return (*stationary)[1] + (*stationary)[2];
		}
		if (gc != 1)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree->GetMean();
	}

	double GetVarGC()	{
		if (gc != 1)	{
			cerr << "error : get mean gc \n";
			exit(1);
		}
		return gctree->GetVar();
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		if (codonmodel == 3)	{
			os << "\tsynrate";
			os << "\tomegats\talphats\tbetats";
			os << "\tomegatv0\talphatv0\tbetatv0";
			os << "\tomegatvgc\talphatvgc\tbetatvgc";
		}
		else if (codonmodel)	{
			os << "\tsynrate\tomega\talpha\tbeta";
		}
		else	{
			os << "\trate";
		}
		if (! Unconstrained())	{
			os << "\talpha\tbeta";
		}
		if (gc == 3)	{
			os << "\tmeangc1\tmeangc2\tmeangc3\talpha\tbeta";
		}
		else if (gc)	{
			os << "\tmeangc\talpha\tbeta";
		}
		else	{
			os << "\tgc\tstatent";
		}
		os << "\trrent";

		if (isCalibrated())	{
			os << "\trootage";
		}
		if (chronoprior == 1)	{
			os << "\tp1\tp2";
			os << "\tnumerror";
		}
		os << '\n';
	}

	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		if (codonmodel == 3)	{
			os << '\t' << GetMeanOmegaTs() << '\t' << *omegatsalpha << '\t' << *omegatsbeta;
			os << '\t' << GetMeanOmegaTv0() << '\t' << *omegatv0alpha << '\t' << *omegatv0beta;
			os << '\t' << GetMeanOmegaTvGC() << '\t' << *omegatvgcalpha << '\t' << *omegatvgcbeta;
		}
		else if (codonmodel)	{
			os << '\t' << GetMeanOmega() << '\t' << *omegaalpha << '\t' << *omegabeta;
		}


		if (! Unconstrained())	{
			os << '\t' << sigmaalpha->val() << '\t' << sigmabeta->val();
		}

		if (gc == 3)	{
			os << '\t' << GetMeanGC1();
			os << '\t' << GetMeanGC2();
			os << '\t' << GetMeanGC3();
			os << '\t' << *gcalpha << '\t' << *gcbeta;
		}
		else if (gc)	{
			os << '\t' << GetMeanGC();
			os << '\t' << *gcalpha << '\t' << *gcbeta;
		}
		else	{
			os << '\t' << GetMeanGC();
			os << '\t' << stationary->val().GetEntropy();
		}

		os << '\t' << relrate->val().GetEntropy();

		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		if (chronoprior == 1)	{
			os << '\t' << *Chi << '\t' << *Chi2;
			os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}

		os << '\n';
		os.flush();
	}

	void PrintSuffStat(ostream& os)	{
		pathconjtree->PrintSuffStat(os);
	}

	void PrintSynNonSyn(string name)	{
		if (Codon())	{
			if (! pathconjtree)	{
				cerr << "error in branchmodel syn nonsyn: only under conjugate path\n";
				exit(1);
			}
			pathconjtree->ActivateSufficientStatistic();
			ofstream os((name + ".nans").c_str(), ios_base::app);
			pathconjtree->PrintSynNonSyn(os);
			os.flush();
			/*
			ofstream osgc((name + ".nansgc").c_str(), ios_base::app);
			pathconjtree->PrintSynNonSyn(osgc,1,2);
			osgc.flush();
			ofstream osat((name + ".nansat").c_str(), ios_base::app);
			pathconjtree->PrintSynNonSyn(osat,0,3);
			osat.flush();
			*/
		}
	}

	void ToStream(ostream& os)	{
		if (Unconstrained())	{
			os << *syngammatree << '\n';
		}
		else	{
			os << *chronogram << '\n';
			if (isCalibrated())	{
				os << *GetCalibratedChronogram()->GetScale() << '\n';
			}
			if (chronoprior)	{
				os << *Chi << '\t' << *Chi2 << '\n';
			}
			os << *sigmaalpha << '\n';
			os << *sigmabeta << '\n';
			os << *gamsynratetree << '\n';
			/*
			os << *sigma << '\n';
			os << *logntree << '\n';
			*/
		}


		os << *relrate << '\n';

		if (codonmodel == 3)	{
			os << *omegatsalpha << '\n';
			os << *omegatsbeta << '\n';
			os << *omegatstree << '\n';
			os << *omegatv0alpha << '\n';
			os << *omegatv0beta << '\n';
			os << *omegatv0tree << '\n';
			os << *omegatvgcalpha << '\n';
			os << *omegatvgcbeta << '\n';
			os << *omegatvgctree << '\n';
		}
		else if (codonmodel)	{
			os << *omegaalpha << '\n';
			os << *omegabeta << '\n';
			os << *omegatree << '\n';
		}

		if (gc == 3)	{
			os << *gcalpha << '\n';
			os << *gcbeta << '\n';
			os << *gctree1 << '\n';
			os << *gctree2 << '\n';
			os << *gctree3 << '\n';
		}
		else if (gc)	{
			os << *gcalpha << '\n';
			os << *gcbeta << '\n';
			os << *gctree << '\n';
		}
		else	{
			os << *stationary << '\n';
		}
	}

	void FromStream(istream& is)	{
		if (Unconstrained())	{
			is >> *syngammatree;
		}
		else	{
			is >> *chronogram;
			if (isCalibrated())	{
				is >> *GetCalibratedChronogram()->GetScale();
			}
			if (chronoprior)	{
				is >> *Chi >> *Chi2;
			}
			/*
			is >> *sigma;
			is >> *logntree;
			*/
			is >> *sigmaalpha;
			is >> *sigmabeta;
			is >> *gamsynratetree;
		}

		is >> *relrate;

		if (codonmodel == 3)	{
			is >> *omegatsalpha;
			is >> *omegatsbeta;
			is >> *omegatstree;
			is >> *omegatv0alpha;
			is >> *omegatv0beta;
			is >> *omegatv0tree;
			is >> *omegatvgcalpha;
			is >> *omegatvgcbeta;
			is >> *omegatvgctree;
		}
		else if (codonmodel)	{
			is >> *omegaalpha;
			is >> *omegabeta;
			is >> *omegatree;
		}

		if (gc == 3)	{
			is >> *gcalpha;
			is >> *gcbeta;
			is >> *gctree1;
			is >> *gctree2;
			is >> *gctree3;
		}
		else if (gc)	{
			is >> *gcalpha;
			is >> *gcbeta;
			is >> *gctree;
		}
		else	{
			is >> *stationary;
		}
	}
};

#endif
