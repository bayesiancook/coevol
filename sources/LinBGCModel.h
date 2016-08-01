
#ifndef LINBGC
#define LINBGC

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSequenceAlignment.h"

#include "Chronogram.h"

#include "BranchProcess.h"
#include "NodeProcess.h"
#include "MatrixTree.h"

#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"
#include "MeanExpTree.h"
#include "Normal.h"

// #include "GCProcess.h"

#include "GeneralConjugatePath.h"

#include "Jeffreys.h"
#include "Partition.h"

#include "NonRevCodonSubMatrix.h"

#include "ConjugateInverseWishart.h"

class ConjugateMultivariateNormalIIDArray : public IIDArray<RealVector>	{

	public:
	ConjugateMultivariateNormalIIDArray(int insize, MultiNormalSemiConjugate* insigma) : IIDArray<RealVector>(insize)	{
		sigma = insigma;
		Create();
	}

	void Reset()	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->SetAtZero();
		}
	}

	double GetMean(int k)	{
		double mean = 0;
		for (int i=0; i<GetSize(); i++)	{
			mean += (*GetVal(i))[k];
		}
		mean /= GetSize();
		return mean;
	}

	double GetVar(int k)	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<GetSize(); i++)	{
			double tmp = (*GetVal(i))[k];
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= GetSize();
		var /= GetSize();
		var -= mean * mean;
		return var;
	}

	double GetMeanVar()	{
		double mean = 0;
		for (int k=0; k<GetVal(0)->GetDim(); k++)	{
			mean += GetVar(k);
		}
		return mean / GetVal(0)->GetDim();
	}

	protected:

	Rvar<RealVector>* CreateVal(int site)	{
		return new ConjugateMultivariateNormal(sigma);
	}

	MultiNormalSemiConjugate* sigma;
};

class ConjugateInverseWishartIIDArray : public IIDArray<CovMatrix>	{

	public:
	ConjugateInverseWishartIIDArray(int insize, SigmaZero* incenter, int indf) : IIDArray<CovMatrix>(insize)	{
		center = incenter;
		df = indf;
		Create();
	}

	ConjugateInverseWishart* GetConjugateInverseWishart(int site)	{
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(GetVal(site));
		if (!tmp)	{
			cerr << "error in conj inv wishart iid array\n";
			exit(1);
		}
		return tmp;
	}

	protected:

	Rvar<CovMatrix>* CreateVal(int site)	{
		return new ConjugateInverseWishart(center,df);
	}

	SigmaZero* center;
	int df;
};


class ResampleSigmaMove : public MCUpdate	{

	public:

	ResampleSigmaMove(ConjugateInverseWishart* insigma)	{
		sigma = insigma;
	}

	double Move(double tuning_modulator){
		sigma->Integrate();
		sigma->Resample();
		return 1;
	}

	private:

	ConjugateInverseWishart* sigma;
};

class ResampleSigmaIIDArrayMove : public MCUpdate	{

	public:

	ResampleSigmaIIDArrayMove(ConjugateInverseWishartIIDArray* inarray)	{
		array = inarray;
	}

	double Move(double tuning_modulator){
		for (int i=0; i<array->GetSize(); i++)	{
			array->GetConjugateInverseWishart(i)->Integrate();
			array->GetConjugateInverseWishart(i)->Resample();
		}
		return 1;
	}

	private:

	ConjugateInverseWishartIIDArray* array;

};

class GeneBranchRate : public Dvar<PosRealVector>	{

	public:

	GeneBranchRate(Var<PosRealVector>* insynrate, Var<PosRealVector>* inomega, Var<RealVector>* inx, Var<RealVector>* iny, Var<RealVector>* inz)	{

		setval(PosRealVector(2*insynrate->GetDim()));
		bkvalue = PosRealVector(2*insynrate->GetDim());
		synrate = insynrate;
		omega = inomega;
		x = inx;
		y = iny;
		z = inz;
		Register(synrate);
		Register(omega);
		Register(x);
		Register(y);
		Register(z);
		specialUpdate();
	}

	protected:

	void specialUpdate()	{
		for (int i=0; i<GetDim()/2; i++)	{
			(*this)[i] = (*synrate)[i] * exp((*x)[i] + (*y)[i] + (*z)[i]);
			(*this)[GetDim()/2 + i] = (*omega)[i] * exp((*x)[GetDim()/2 + i] + (*y)[GetDim()/2 + i] + (*z)[GetDim()/2 + i]);
		}
	}

	Var<PosRealVector>* synrate;
	Var<PosRealVector>* omega;
	Var<RealVector>* x;
	Var<RealVector>* y;
	Var<RealVector>* z;

};

class NonRevCodonMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:

	NonRevCodonMatrixTree(CodonStateSpace* incodonstatespace, LengthTree* intree, BranchPartition* inpartition, GeneBranchRate** inrate, Var<UnitReal>* inrootgc, bool innormalise = false, int indiscn = 10)	{
		SetWithRoot(true);
		tree = intree;
		partition = inpartition;
		rate = inrate;
		rootgc = inrootgc;
		normalise = innormalise;
		discn = indiscn;
		codonstatespace = incodonstatespace;
		RecursiveCreate(GetRoot());
	}


	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomRootCodonSubMatrix(codonstatespace,rootgc,false);
		}
		return new RandomNonRevCodonSubMatrix(codonstatespace,tree->GetBranchVal(link->GetBranch()),rate[partition->GetAlloc(link->GetBranch())-1],false,discn);
	}

	Tree* GetTree() {return tree->GetTree();}

	private:

	LengthTree* tree;
	BranchPartition* partition;
	GeneBranchRate** rate;
	Var<UnitReal>* rootgc;
	bool normalise;
	int discn;
	CodonStateSpace* codonstatespace;
};


/*
class PartitionMatrixTree : public BranchValPtrTree<RandomSubMatrix>	{

	public:

	PartitionMatrixTree(Tree* intree, BranchPartition* inpartition, RandomSubMatrix** inmatrixarray, RandomSubMatrix* inrootmatrix)	{
		SetWithRoot(true);
		tree = intree;
		partition = inpartition;
		matrixarray = inmatrixarray;
		rootmatrix = inrootmatrix;
		RecursiveCreate(GetRoot());
	}

	~PartitionMatrixTree()	{
		// nothing because there was no matrix allocation upon creating the object
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return rootmatrix;
		}
		return matrixarray[partition->GetAlloc(link->GetBranch())-1];
	}

	Tree* GetTree() {return tree;}

	private:

	Tree* tree;
	BranchPartition* partition;
	RandomSubMatrix** matrixarray;
	RandomSubMatrix* rootmatrix;
};
*/

class LinBGCModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	FileSequenceAlignment** nucdata;
	CodonSequenceAlignment** codondata;
	TaxonSet* taxonset;

	int Ngene;
	int* Nsite;

	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	// int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	// chronogram
	Chronogram* chronogram;
	Jeffreys* synsigma;
	LogNormalTreeProcess* synratetree;

	BranchPartition* partition;
	int Nmat;

	int L; // (6 independent nucleotide pairs, because of strand symmetry)
	JeffreysIIDArray* BranchDiagArray;
	JeffreysIIDArray* GeneDiagArray;
	JeffreysIIDArray* BranchGeneDiagArray;
	SigmaZero* BranchSigmaZero;
	SigmaZero* GeneSigmaZero;
	SigmaZero* BranchGeneSigmaZero;
	ConjugateInverseWishart* branchsigma;
	ConjugateInverseWishart* genesigma;
	ConjugateInverseWishartIIDArray* branchgenesigma;

	Jeffreys* synrate0;
	Jeffreys* omega0;
	IIDExp* synrate;
	IIDExp* omega;
	ConjugateMultivariateNormalIIDArray* brancheffect;
	ConjugateMultivariateNormalIIDArray* geneeffect;
	ConjugateMultivariateNormalIIDArray** branchgeneeffect;

	// root effects
	Jeffreys* rootalpha;
	Jeffreys* rootbeta;
	BetaIIDArray* rootgc;
	// RandomSubMatrix** rootmatrix;

	GeneticCodeType codetype;
	CodonStateSpace* codonstatespace;
	// RandomSubMatrix*** genebranchcodonmatrix;

	GeneBranchRate*** genebranchrate;

	NonRevCodonMatrixTree** matrixtree;
	// PartitionMatrixTree** matrixtree;

	// phylo process
	BranchMatrixPathConjugateTree** pathconjtree;
	PhyloProcess** phyloprocess;

	bool priorsampling;
	int clampsuffstat;
	string suffstatfile;

	bool conjpath;

	bool normalise;


	int withbranch;
	int withgene;
	int withbranchgene;

	int nrep;
	int discn;

	int df;

	public:

	SequenceAlignment* GetData(int gene)	{
		return codondata[gene];
	}

	LinBGCModel(string datafile, string treefile, double priorsigma, int indf, int inconjpath, bool innormalise, int innrep, string insuffstatfile, string mix, int inwithbranch, int inwithgene, int inwithbranchgene, int indiscn, GeneticCodeType type = Universal, bool sample=true)	{

		discn = indiscn;
		df = indf;

		withbranch = inwithbranch;
		withgene = inwithgene;
		withbranchgene = inwithbranchgene;

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		codetype = type;
		codonstatespace = new CodonStateSpace(codetype);

		// get data from file

		ifstream is(datafile.c_str());
		is >> Ngene;
		Nsite = new int[Ngene];
		nucdata = new FileSequenceAlignment*[Ngene];
		codondata = new CodonSequenceAlignment*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			string filename;
			is >> filename;
			cerr << filename << '\n';
			nucdata[gene] = new FileSequenceAlignment(filename);
			codondata[gene] = new CodonSequenceAlignment(nucdata[gene], true, type);
			Nsite[gene] = GetData(gene)->GetNsite();
			cerr << filename << '\t' << Nsite[gene] << '\n';
		}

		priorsampling = false;

		if (inconjpath == -1)	{
			conjpath = true;
		}
		else if (inconjpath == 2)	{
			conjpath = false;
			priorsampling = true;
		}
		else	{
			conjpath = inconjpath;
		}
		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 10 : 1;
		}

		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		taxonset = nucdata[0]->GetTaxonSet();

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

		chronogram = new Chronogram(tree,One);
		chronogram->Clamp();

		double min = 0.001;
		double max = 1000;

		synsigma = new Jeffreys(min,max,Zero);
		synratetree = new LogNormalTreeProcess(chronogram,synsigma,INTEGRAL);
		synratetree->Reset();

		if (mix == "branch")	{
			partition = new BranchPartition(tree,true);
		}
		else if (mix != "None")	{
			partition = new BranchPartition(tree,mix);
		}
		else	{
			partition = new BranchPartition(tree);
		}

		Nmat = partition->GetNComponent() - 1;
		cerr << "partition  " << Nmat << '\n';

		L = 6;

		double mindiag = 0.001;
		double maxdiag = 1000;
		BranchDiagArray = new JeffreysIIDArray(2*L,mindiag,maxdiag,Zero);
		GeneDiagArray = new JeffreysIIDArray(2*L,mindiag,maxdiag,Zero);
		BranchGeneDiagArray = new JeffreysIIDArray(2*L,mindiag,maxdiag,Zero);

		cerr << priorsigma << '\n';
		if (priorsigma == -1)	{
			BranchDiagArray->setval(1.0);
			GeneDiagArray->setval(1.0);
			BranchGeneDiagArray->setval(1.0);
		}
		else	{
			BranchDiagArray->ClampAt(priorsigma);
			GeneDiagArray->ClampAt(priorsigma);
			BranchGeneDiagArray->ClampAt(priorsigma);
		}

		BranchSigmaZero = new SigmaZero(BranchDiagArray);
		GeneSigmaZero = new SigmaZero(GeneDiagArray);
		BranchGeneSigmaZero = new SigmaZero(BranchGeneDiagArray);

		cerr << "branch sigma \n";
		branchsigma = new ConjugateInverseWishart(BranchSigmaZero,2*L+df);
		cerr << "gene sigma\n";
		genesigma = new ConjugateInverseWishart(GeneSigmaZero,2*L+df);
		cerr << "sigma iid\n";
		branchgenesigma = new ConjugateInverseWishartIIDArray(Nmat,BranchGeneSigmaZero,2*L+df);

		cerr << "basal rates\n";
		// L  basal rates and L basal dN/dS
		// exponential variables
		synrate0 = new Jeffreys(min,max,Zero);
		omega0 = new Jeffreys(min,max,Zero);
		synrate0->setval(0.1);
		omega0->setval(0.1);
		synrate = new IIDExp(L,synrate0);
		// synrate->SetAtOne();
		for (int i=0; i<L; i++)	{
			(*synrate)[i] = 0.1;
		}
		omega = new IIDExp(L,omega0);
		// omega->SetAtOne();
		for (int i=0; i<L; i++)	{
			(*omega)[i] = 0.2;
		}

		cerr << "normal effects\n";
		// multivariate normal effects
		brancheffect = new ConjugateMultivariateNormalIIDArray(Nmat,branchsigma);
		brancheffect->Reset();
		geneeffect = new ConjugateMultivariateNormalIIDArray(Ngene,genesigma);
		geneeffect->Reset();
		branchgeneeffect = new ConjugateMultivariateNormalIIDArray*[Nmat];
		for (int mat=0; mat<Nmat; mat++)	{
			branchgeneeffect[mat] = new ConjugateMultivariateNormalIIDArray(Ngene,branchgenesigma->GetConjugateInverseWishart(mat));
			branchgeneeffect[mat]->Reset();
		}

		// gene branch rates
		// gene branch matrices
		// genebranchcodonmatrix =  new RandomSubMatrix**[Ngene];
		genebranchrate = new GeneBranchRate**[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			// genebranchcodonmatrix[gene] = new RandomSubMatrix*[Nmat];
			genebranchrate[gene] = new GeneBranchRate*[Nmat];
			for (int mat=0; mat<Nmat; mat++)	{
				genebranchrate[gene][mat] = new GeneBranchRate(synrate,omega,brancheffect->GetVal(mat),geneeffect->GetVal(gene),branchgeneeffect[mat]->GetVal(gene));
				// genebranchcodonmatrix[gene][mat] = new RandomNonRevCodonSubMatrix(codonstatespace,genebranchrate[gene][mat],false,discn);
			}
		}

		// root eq frequencies
		rootalpha = new Jeffreys(min,max,Zero);
		rootbeta = new Jeffreys(min,max,Zero);
		rootgc = new BetaIIDArray(Ngene,rootalpha,rootbeta);
		/*
		rootmatrix = new RandomSubMatrix*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			rootmatrix[gene] = new RandomRootCodonSubMatrix(codonstatespace,rootgc->GetVal(gene),false);
		}
		*/

		matrixtree = new NonRevCodonMatrixTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			rootgc->GetVal(gene)->setval(0.4);
			matrixtree[gene] = new NonRevCodonMatrixTree(codonstatespace,synratetree,partition,genebranchrate[gene],rootgc->GetVal(gene),normalise,discn);
		}

		/*
		matrixtree = new PartitionMatrixTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			matrixtree[gene] = new PartitionMatrixTree(tree,partition,genebranchcodonmatrix[gene],rootmatrix[gene]);
		}
		*/

		cerr << "phylo processes\n";
		cerr << discn << '\n';

		if (conjpath)	{
			pathconjtree = new BranchMatrixPathConjugateTree*[Ngene];
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				pathconjtree[gene] = new BranchMatrixPathConjugateTree(synratetree, matrixtree[gene], GetData(gene));
				phyloprocess[gene] = new PathConjugatePhyloProcess(pathconjtree[gene]);
			}
		}
		else	{
			pathconjtree = 0;
			if (priorsampling)	{
				phyloprocess = 0;
			}
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene] = new BranchMatrixPhyloProcess(synratetree, matrixtree[gene], GetData(gene));
			}
		}

		if (phyloprocess)	{
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene]->Unfold();
			}
		}
		if (sample)	{
			if (phyloprocess)	{
				for (int gene=0; gene<Ngene; gene++)	{
					phyloprocess[gene]->Sample();
				}
			}
		}

		// register model
		RootRegister(Zero);
		RootRegister(One);
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}

	}

	Tree* GetTree() {return tree;}

	Chronogram* GetChronogram() {
		return chronogram;
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		total += synsigma->GetLogProb();
		total += synratetree->GetLogProb();

		total += synrate0->GetLogProb();
		total += synrate->GetLogProb();
		total += omega0->GetLogProb();
		total += omega->GetLogProb();

		total += BranchDiagArray->GetLogProb();
		total += GeneDiagArray->GetLogProb();
		total += BranchGeneDiagArray->GetLogProb();
		total += branchsigma->GetLogProb();
		total += genesigma->GetLogProb();
		total += branchgenesigma->GetLogProb();

		total += brancheffect->GetLogProb();
		total += geneeffect->GetLogProb();
		for (int mat=0; mat<Nmat; mat++)	{
			total += branchgeneeffect[mat]->GetLogProb();
		}

		total += rootalpha->GetLogProb();
		total += rootbeta->GetLogProb();
		total += rootgc->GetLogProb();

		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		if (priorsampling)	{
			return 0;
		}
		else if (clampsuffstat)	{
			for (int gene=0; gene<Ngene; gene++)	{
				ret += pathconjtree[gene]->GetLogProb();
			}
		}
		else	{
			for (int gene=0; gene<Ngene; gene++)	{
				ret += phyloprocess[gene]->GetLogProb();
			}
		}
		return ret;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			if (! clampsuffstat)	{
				for (int gene=0; gene<Ngene; gene++)	{
					scheduler.Register(new DSemiConjugateMappingMove(phyloprocess[gene],pathconjtree[gene]),1,"mapping + sufficient stat");
				}
			}
		}
		else	{
			if (phyloprocess)	{
				for (int gene=0; gene<Ngene; gene++)	{
					scheduler.Register(new SimpleMove(phyloprocess[gene],1),1,"mapping");
				}
			}
		}

		for (int i=0; i<nrep; i++)	{

			/*
			scheduler.Register(new SimpleMove(synratetree,10),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,1),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,0.1),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,0.01),10,"synratetree");

			scheduler.Register(new SimpleMove(synsigma,10),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,1),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,0.1),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,0.01),100,"syn sigma");
			*/

			// scheduler.Register(new SimpleMove(synrate,10),100,"synrate");
			scheduler.Register(new SimpleMove(synrate,1),10,"synrate");
			scheduler.Register(new SimpleMove(synrate,0.1),10,"synrate");
			scheduler.Register(new SimpleMove(synrate,0.01),10,"synrate");

			scheduler.Register(new SimpleMove(synrate0,10),100,"synrate0");
			scheduler.Register(new SimpleMove(synrate0,1),100,"synrate0");
			scheduler.Register(new SimpleMove(synrate0,0.1),100,"synrate0");
			scheduler.Register(new SimpleMove(synrate0,0.01),100,"synrate0");

			// scheduler.Register(new SimpleMove(omega,10),100,"omega");
			scheduler.Register(new SimpleMove(omega,1),10,"omega");
			scheduler.Register(new SimpleMove(omega,0.1),10,"omega");
			scheduler.Register(new SimpleMove(omega,0.01),10,"omega");

			scheduler.Register(new SimpleMove(omega0,10),100,"omega0");
			scheduler.Register(new SimpleMove(omega0,1),100,"omega0");
			scheduler.Register(new SimpleMove(omega0,0.1),100,"omega0");
			scheduler.Register(new SimpleMove(omega0,0.01),100,"omega0");

			scheduler.Register(new SimpleMove(rootgc,1),10,"root bgc");
			scheduler.Register(new SimpleMove(rootgc,0.1),10,"root bgc");

			// HYPER ROOT
			scheduler.Register(new SimpleMove(rootalpha,1),100,"root alpha");
			scheduler.Register(new SimpleMove(rootalpha,0.1),100,"root alpha");

			scheduler.Register(new SimpleMove(rootbeta,1),100,"root beta");
			scheduler.Register(new SimpleMove(rootbeta,0.1),100,"root beta");

			// EFFECTS
			if (withbranch)	{
				scheduler.Register(new SimpleMove(brancheffect,1),10,"branch effect");
				scheduler.Register(new SimpleMove(brancheffect,0.1),10,"branch effect");
			}

			if (withgene)	{
				scheduler.Register(new SimpleMove(geneeffect,1),10,"gene effect");
				scheduler.Register(new SimpleMove(geneeffect,0.1),10,"gene effect");
			}

			if (withbranchgene)	{
				for (int mat=0; mat<Nmat; mat++)	{
					scheduler.Register(new SimpleMove(branchgeneeffect[mat],1),10,"branch gene effect");
					scheduler.Register(new SimpleMove(branchgeneeffect[mat],0.1),10,"branch gene effect");
				}
			}

			// CONJ + HYPER
			if (withbranch)	{
				scheduler.Register(new ResampleSigmaMove(branchsigma), 1, "branch sigma");
			}
			if (withgene)	{
				scheduler.Register(new ResampleSigmaMove(genesigma), 1, "gene sigma");
			}
			if (withbranchgene)	{
				scheduler.Register(new ResampleSigmaIIDArrayMove(branchgenesigma), 1, "branch gene sigma");
			}

			if (withbranch)	{
				scheduler.Register(new SimpleMove(BranchDiagArray,10),10,"branch sigma0");
				scheduler.Register(new SimpleMove(BranchDiagArray,1),10,"branch sigma0");
				scheduler.Register(new SimpleMove(BranchDiagArray,0.1),10,"branch sigma0");
			}

			if (withgene)	{
				scheduler.Register(new SimpleMove(GeneDiagArray,10),10,"gene sigma0");
				scheduler.Register(new SimpleMove(GeneDiagArray,1),10,"gene sigma0");
				scheduler.Register(new SimpleMove(GeneDiagArray,0.1),10,"gene sigma0");
			}

			if (withbranchgene)	{
				scheduler.Register(new SimpleMove(BranchGeneDiagArray,10),10,"branch gene sigma0");
				scheduler.Register(new SimpleMove(BranchGeneDiagArray,1),10,"branch gene sigma0");
				scheduler.Register(new SimpleMove(BranchGeneDiagArray,0.1),10,"branch gene sigma0");
			}

			// NON CONJ
			/*
			scheduler.Register(new SimpleMove(branchsigma,10),100,"branch sigma");
			scheduler.Register(new SimpleMove(branchsigma,1),100,"branch sigma");
			scheduler.Register(new SimpleMove(branchsigma,0.1),100,"branch sigma");
			scheduler.Register(new SimpleMove(branchsigma,0.01),100,"branch sigma");

			scheduler.Register(new SimpleMove(genesigma,10),100,"gene sigma");
			scheduler.Register(new SimpleMove(genesigma,1),100,"gene sigma");
			scheduler.Register(new SimpleMove(genesigma,0.1),100,"gene sigma");
			scheduler.Register(new SimpleMove(genesigma,0.01),100,"gene sigma");

			for (int mat=0; mat<Nmat; mat++)	{
				scheduler.Register(new SimpleMove(branchgenesigma,10),100,"branch gene sigma");
				scheduler.Register(new SimpleMove(branchgenesigma,1),100,"branch gene sigma");
				scheduler.Register(new SimpleMove(branchgenesigma,0.1),100,"branch gene sigma");
				scheduler.Register(new SimpleMove(branchgenesigma,0.01),100,"branch gene sigma");
			}
			*/

		}
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,false);
		// scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";

		/*
		synsigma->Sample();
		synratetree->Sample();
		*/

		synrate0->Sample();
		synrate->Sample();
		omega0->Sample();
		omega->Sample();

		BranchDiagArray->Sample();
		GeneDiagArray->Sample();
		BranchGeneDiagArray->Sample();
		branchsigma->Sample();
		genesigma->Sample();
		branchgenesigma->Sample();

		brancheffect->Sample();
		geneeffect->Sample();
		for (int mat=0; mat<Nmat; mat++)	{
			branchgeneeffect[mat]->Sample();
		}

		rootalpha->Sample();
		rootbeta->Sample();
		rootgc->Sample();

		cerr << "ok\n";
	}

	int GetNgene()	{
		return Ngene;
	}

	LogNormalTreeProcess* GetSynRateTree()	{
		return synratetree;
	}

	/*
	void UpdateGeneProcess(int gene)	{
		if (geneprocesstype == 4)	{
			argeneprocess[gene]->specialUpdate();
			mixedgeneprocess[gene]->specialUpdate();
		}
		else if (! geneprocesstype || (geneprocesstype == 3))	{
			argeneprocess[gene]->specialUpdate();
		}
	}

	void UpdateGeneProcesses()	{
		for (int gene=0; gene<Ngene; gene++)	{
			UpdateGeneProcess(gene);
		}
	}

	void UpdateBranchProductProcesses()	{
		for (int gene=0; gene<Ngene; gene++)	{
			branchproductprocess[gene]->specialUpdate();
		}
	}

	void UpdateBGC()	{
		GetBGCProcess()->specialUpdate();
		UpdateGeneProcesses();
		UpdateBranchProductProcesses();
	}
	*/

	double GetMeanSynRate()	{
		return GetSynRateTree()->GetMeanRate();
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		/*
		os << "\tsynrate";
		os << "\tsynsigma";
		os << "\tsynrate0";
		os << "\tomega0";
		*/
		os << "\tat2ta\tcg2gc\tat2cg\tcg2at\tat2gc\tcg2ta";
		os << "\tomat2ta\tomcg2gc\tomat2cg\tomcg2at\tomat2gc\tomcg2ta";
		os << "\tmeangc\tvargc";
		os << "\tbranch";
		os << "\tgene";
		for (int mat=0; mat<Nmat; mat++)	{
			os << "\tbranchgene_" << mat;
		}
		/*
		for (int i=0; i<2*L; i++)	{
			for (int j=i+1; j<2*L; j++)	{
				os << "\tsigma_" << i << "_" << j;
			}
		}
		*/

		/*
		os << '\t' << "rootalpha";
		os << '\t' << "rootbeta";
		*/
		os << '\n';
	}

	void Trace(ostream& os)	{

		// os.precision(10);
		os << GetLogPrior() << '\t' << GetLogLikelihood();
		/*
		os << '\t' << GetMeanSynRate();
		os << '\t' << synsigma->val();
		os << '\t' << synrate0->val();
		os << '\t' << omega0->val();
		*/
		for (int i=0; i<L; i++)	{
			os << '\t' << (*synrate)[i];
		}
		for (int i=0; i<L; i++)	{
			os << '\t' << (*omega)[i];
		}
		os << '\t' << rootgc->GetMean();
		os << '\t' << rootgc->GetVar();
		os << '\t' << brancheffect->GetMeanVar();
		os << '\t' << geneeffect->GetMeanVar();
		for (int mat=0; mat<Nmat; mat++)	{
			os << '\t' << branchgeneeffect[mat]->GetMeanVar();
		}
		/*
		CovMatrix* mat = branchgenesigma->GetVal(0);
		for (int i=0; i<2*L; i++)	{
			for (int j=i+1; j<2*L; j++)	{
				os << '\t' << (*mat)[i][j];
			}
		}
		*/

		/*
		os << '\t' << rootalpha->val();
		os << '\t' << rootbeta->val();
		*/
		os << '\n';
		os.flush();
	}

	void TraceSigmaHeader(ostream& os)	{
		int first = 1;
		for (int k=0; k<Nmat; k++)	{
			for (int i=0; i<2*L; i++)	{
				for (int j=i+1; j<2*L; j++)	{
					if (! first)	{
						os << '\t';
					}
					else	{
						first = 0;
					}
					os << k << "_" << i << "_" << j;
				}
			}
		}
	}

	void TraceSigma(ostream& os)	{

		int first = 1;
		for (int k=0; k<Nmat; k++)	{
			CovMatrix* mat = branchgenesigma->GetVal(k);
			for (int i=0; i<2*L; i++)	{
				for (int j=i+1; j<2*L; j++)	{
					if (! first)	{
						os << '\t';
					}
					else	{
						first = 0;
					}
					os << (*mat)[i][j] / sqrt((*mat)[i][i] * (*mat)[j][j]);
				}
			}
		}
		os << '\n';
	}

	void ToStream(ostream& os)	{
		os << *synsigma << '\n';
		os << *synratetree << '\n';

		os << *synrate0 << '\n';
		os << *synrate << '\n';
		os << *omega0 << '\n';
		os << *omega << '\n';

		os << *BranchDiagArray << '\n';
		os << *GeneDiagArray << '\n';
		os << *BranchGeneDiagArray << '\n';
		os << *branchsigma << '\n';
		os << *genesigma << '\n';
		os << *branchgenesigma << '\n';

		os << *brancheffect << '\n';
		os << *geneeffect << '\n';
		for (int mat=0; mat<Nmat; mat++)	{
			os << *branchgeneeffect[mat] << '\n';
		}

		os << *rootalpha << '\n';
		os << *rootbeta << '\n';
		os << *rootgc << '\n';
	}



	void FromStream(istream& is)	{
		is >> *synsigma;
		is >> *synratetree;

		is >> *synrate0;
		is >> *synrate;
		is >> *omega0;
		is >> *omega;

		is >> *BranchDiagArray;
		is >> *GeneDiagArray;
		is >> *BranchGeneDiagArray;
		is >> *branchsigma;
		is >> *genesigma;
		is >> *branchgenesigma;

		is >> *brancheffect;
		is >> *geneeffect;
		for (int mat=0; mat<Nmat; mat++)	{
			is >> *branchgeneeffect[mat];
		}

		is >> *rootalpha;
		is >> *rootbeta;
		is >> *rootgc;
	}
};

#endif
