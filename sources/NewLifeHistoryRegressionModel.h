
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
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
#include "CoalCalibratedChronogram.h"
// #include "InterpolatedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "MultiVarNormal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GCProcess.h"

#include "GeneralConjugatePath.h"

#include "Jeffreys.h"

#include "Partition.h"

#include "SplitLengthTree.h"
#include "SplitChronogram.h"
#include "SplitMultiVariateMove.h"
#include "MultiVariatePropagateMove.h"

#include "PartitionMultiVariateTreeProcess.h"

#include "WhiteNoise.h"


class IIDNormalArray : public IIDArray<RealVector>	{

	public:

	IIDNormalArray(int insize, int indim, Var<Real>* inmean, Var<PosReal>* invar, Var<PosReal>* invar0 = 0) : IIDArray<RealVector>(insize)	{
		mean = inmean;
		var = invar;
		dim = indim;
		if (invar0)	{
			var0 = invar0;
		}
		else	{
			var0 = var;
		}
		Create();
	}

	int GetDim()	{
		return dim;
	}

	IIDNormal* GetIIDNormal(int i)	{
		return dynamic_cast<IIDNormal*>(GetVal(i));
	}

	void SetAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetIIDNormal(i)->SetAtZero();
		}
	}

	void ClampAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetIIDNormal(i)->ClampAtZero();
		}
	}

	protected:

	Rvar<RealVector>* CreateVal(int mat)	{
		if (!mat)	{
			return new IIDNormal(dim,mean,var);
		}
		return new IIDNormal(dim,mean,var0);
	}

	private:

	int dim;
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<PosReal>* var0;
};

class DLinRegContInstantValue : public Dvar<RealVector>	{

	public:

	DLinRegContInstantValue(Var<PosReal>* indate, Var<RealVector>* inb, Var<RealVector>* indrift1, Var<PosReal>* inphi1, Var<RealVector>* indrift2 = 0, Var<PosReal>* inphi2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0) {
		date = indate;
		b = inb;
		drift1 = indrift1;
		phi1 = inphi1;
		drift2 = indrift2;
		phi2 = inphi2;
		agescale = inagescale;
		kt = inkt;

		Register(date);
		Register(b);
		Register(drift1);
		Register(phi1);
		Register(drift2);
		Register(phi2);
		Register(agescale);
	}


	double GetTrend(double t, int index)	{
		if (t >= kt)	{
			return (*drift1)[index] * exp(-phi1->val() * t);
		}
		return (*drift1)[index] * exp(-phi1->val() * t) + (*drift2)[index] * (1 - exp(-phi2->val() * (kt - t)));

	}

	void specialUpdate()	{
		for (int i=0; i<GetDim(); i++)	{
			double t = date->val() * agescale->val();
			for (int i=0; i<GetDim() ; i++) {
				double f = GetTrend(t,i);
				(*this)[i] = (*b)[i] + f;
			}
		}
	}

	protected:

	Var<PosReal>* date;
	Var<RealVector>* b;
	Var<RealVector>* drift1;
	Var<PosReal>* phi1;
	Var<RealVector>* drift2;
	Var<PosReal>* phi2;
	Var<PosReal>* agescale;
	double kt;

};

class DLinRegCont: public NodeValPtrTree<Dvar<RealVector> > {

	public:

	DLinRegCont(NodeBranchVarTree<PosReal,PosReal>* inchrono, NodeVarTree<RealVector>* intree, Var<RealVector>* indrift1 = 0, Var<PosReal>* inphi1 = 0, Var<RealVector>* indrift2 = 0, Var<PosReal>* inphi2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0)	{
		chrono = inchrono;
		tree = intree;
		drift1 = indrift1;
		phi1 = inphi1;
		drift2 = indrift2;
		phi2 = inphi2;
		agescale = inagescale;
		kt = inkt;
		RecursiveCreate(GetRoot());
		specialUpdate();
	}

	~DLinRegCont()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree(){
		return tree->GetTree();
	}

	DLinRegContInstantValue* GetInstantValue(const Node* node)	{
		DLinRegContInstantValue* m = dynamic_cast<DLinRegContInstantValue*> (GetNodeVal(node));
		return m;
	}

	double GetExpVal(const Link* link, int index)	{
		return exp((*GetNodeVal(link->GetNode()))[index]);
	}

	double GetVal(const Link* link, int index)	{
		return (*GetNodeVal(link->GetNode()))[index];
	}

	double GetMean(int index)	{
		int n = 0;
		double tmp = RecursiveGetMean(GetRoot(),index,n);
		return tmp / n;
	}

	double GetMeanExp(int index)	{
		int n = 0;
		double tmp = RecursiveGetMeanExp(GetRoot(),index,n);
		return tmp / n;
	}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	protected:

	void RecursiveSpecialUpdate(const Link* from)	{
		GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	double RecursiveGetMean(const Link* from, int index, int& tot)	{
		double total = GetVal(from,index);
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMean(link->Out(), index, tot);
		}
		return total;
	}

	double RecursiveGetMeanExp(const Link* from, int index, int& tot)	{
		double total = GetVal(from,index);
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMeanExp(link->Out(), index, tot);
		}
		return total;
	}

	Dvar<RealVector>* CreateNodeVal(const Link* from)	{
		return new DLinRegContInstantValue(chrono->GetNodeVal(from->GetNode()), tree->GetNodeVal(from->GetNode()), drift1, phi1, drift2, phi2, agescale, kt);
	}

	NodeBranchVarTree<PosReal,PosReal>* chrono;
	NodeVarTree<RealVector>* tree;
	Var<RealVector>* drift1;
	Var<PosReal>* phi1;
	Var<RealVector>* drift2;
	Var<PosReal>* phi2;
	Var<PosReal>* agescale;
	double kt;
};

class DLinRegSubBranchMean: public Dvar<PosReal>	{

	public:

	DLinRegSubBranchMean(bool inmean, Var<PosReal>* intime, Var<PosReal>* indate, Var<RealVector>* inup1, Var<RealVector>* indown1, Var<RealVector>* inup2, Var<RealVector>* indown2, int inindex, Var<RealVector>* inreg, Var<RealVector>* indrift1, Var<PosReal>* inphi1, Var<RealVector>* indrift2 = 0, Var<PosReal>* inphi2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0, Var<RealVector>* inpulse = 0, Var<PosReal>* inpulsephi = 0) {
		mean = inmean;
		time = intime;
		date = indate;
		up1 = inup1;
		down1 = indown1;
		up2 = inup2;
		down2 = indown2;
		index = inindex;
		reg = inreg;
		drift1 = indrift1;
		phi1 = inphi1;
		drift2 = indrift2;
		phi2 = inphi2;
		agescale = inagescale;
		kt = inkt;
		pulse = inpulse;
		pulsephi = inpulsephi;

		Register(time);
		Register(date);
		Register(up1);
		Register(down1);
		Register(up2);
		Register(down2);
		Register(reg);
		Register(drift1);
		Register(phi1);
		Register(drift2);
		Register(phi2);
		Register(agescale);
		Register(pulse);
		Register(pulsephi);
	}

	double GetTrend(double t, int i)	{
		if (t >= kt)	{
			return (*drift1)[i] * exp(-phi1->val() * t);
		}
		return (*drift1)[i] * exp(-phi1->val() * t) + (*drift2)[i] * (1 - exp(-phi2->val() * (kt - t))) + (*pulse)[i] * exp(-pulsephi->val() * (kt-t));

	}

	double GetRate(double tdown, double tup, double t)	{
		double total = ((t-tdown) * (*up2)[index] + (tup-t) * (*down2)[index]) / (tup - tdown);
		for (int i=0; i<reg->GetDim(); i++)	{
			double tmp1 = ((t-tdown) * (*up1)[i] + (tup-t) * (*down1)[i]) / (tup - tdown);
			double tmp2 = GetTrend(t,i);
			total += (*reg)[i] * (tmp1 + tmp2);
		}
		return exp(total);
	}

	/*
	double getmean(double tdown, double tup, double N)	{
		double total = 0;
		for (int i=0; i<=N; i++)	{
			double t = tup + (tdown-tup) * ((double) i) / N;
			double tmp = GetRate(tdown,tup,t);
			if ((i==0) || (i==N))	{
				total += 0.5 * tmp;
			}
			else	{
				total += tmp;
			}
		}
		return total / N;
	}
	*/

	double getmean(double tdown, double tup, double N)	{
		double total = 0;
		double du = (*up2)[index];
		double dd = (*down2)[index];
		for (int i=0; i<=N; i++)	{
			double t = tdown + (tup-tdown) * ((double) i) / N;
			double u = (t-tdown) / (tup-tdown);
			double v = (tup-t) / (tup-tdown);
			if (tup == tdown)	{
				cerr << "equal times\n";
				exit(1);
			}

			double temp1 = u*du + v*dd;

			double temp2 = 0;
			double expo1 = exp(-phi1->val() * t);
			double expo2 = exp(-phi2->val() * (kt-t));
			double expo3 = exp(-pulsephi->val() * (kt-t));
			for (int j=0; j<reg->GetDim(); j++)	{
				double tmp1 = u * (*up1)[j] + v * (*down1)[j];
				double tmp2 = (*drift1)[j] * expo1;
				if (t<=kt)	{
					tmp2 += (*drift2)[j] * (1-expo2) + (*pulse)[j] * expo3;
				}
				temp2 += (*reg)[j] * (tmp1 + tmp2);
			}

			double temp = exp(temp1 + temp2);
			if ((i==0) || (i==N))	{
				total += 0.5 * temp;
			}
			else	{
				total += temp;
			}
		}
		return total / N;
	}

	void specialUpdate()	{
		if (time)	{
			double tdown = date->val() * agescale->val();
			double tup = tdown  + time->val() * agescale->val();
			if (mean)	{
				setval(getmean(tdown,tup,10));
			}
			else	{
				setval(time->val() * getmean(tdown,tup,10));
			}
		}
		else	{
			setval(1.0);
		}
	}

	protected:

	int index;
	bool mean;
	Var<PosReal>* time;
	Var<PosReal>* date;
	Var<RealVector>* up1;
	Var<RealVector>* down1;
	Var<RealVector>* up2;
	Var<RealVector>* down2;
	Var<RealVector>* reg;
	Var<RealVector>* drift1;
	Var<PosReal>* phi1;
	Var<RealVector>* drift2;
	Var<PosReal>* phi2;
	Var<PosReal>* agescale;
	double kt;
	Var<RealVector>* pulse;
	Var<PosReal>* pulsephi;

};

class DLinRegSub: public BranchValPtrTree<Dvar<PosReal> > {

	public:

	DLinRegSub(bool inwithroot, NodeBranchVarTree<PosReal,PosReal>* inchrono, NodeVarTree<RealVector>* inb1, NodeVarTree<RealVector>* inb2, int inindex, Var<RealVector>* inreg, Var<RealVector>* indrift1 = 0, Var<PosReal>* inphi1 = 0, Var<RealVector>* indrift2 = 0, Var<PosReal>* inphi2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0, Var<RealVector>* inpulse = 0, Var<PosReal>* inpulsephi = 0)	{
		SetWithRoot(inwithroot);
		chrono = inchrono;
		b1 = inb1;
		b2 = inb2;
		index = inindex;
		reg = inreg;
		drift1 = indrift1;
		phi1 = inphi1;
		drift2 = indrift2;
		phi2 = inphi2;
		agescale = inagescale;
		kt = inkt;
		pulse = inpulse;
		pulsephi = inpulsephi;
		RecursiveCreate(GetRoot());
		specialUpdate();
	}

	~DLinRegSub()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree(){
		return b1->GetTree();
	}

	DLinRegSubBranchMean* GetInstantValue(const Branch* branch)	{
		DLinRegSubBranchMean* m = dynamic_cast<DLinRegSubBranchMean*> (GetBranchVal(branch));
		return m;
	}

	double GetVal(const Link* link)	{
		return GetBranchVal(link->GetBranch())->val();
	}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	double GetTotal()	{
		double tmp = RecursiveGetTotal(GetRoot());
		return tmp;
	}

	void CutOff(double cutoff)	{
		RecursiveCutOff(GetRoot(),cutoff);
	}

	protected:

	double RecursiveGetTotal(const Link* from)	{
		double total = 0;
		if (WithRoot() || (! from->isRoot()))	{
			total = GetVal(from);
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetTotal(link->Out());
		}
		return total;
	}

	void RecursiveSpecialUpdate(const Link* from)	{
		if (WithRoot() || (! from->isRoot()))	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	void RecursiveCutOff(const Link* from, double cutoff)	{
		if (WithRoot() || (! from->isRoot()))	{
			if (GetBranchVal(from->GetBranch())->val() > cutoff)	{
				GetBranchVal(from->GetBranch())->setval(cutoff);
			}
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCutOff(link->Out(),cutoff);
		}
	}

	Dvar<PosReal>* CreateBranchVal(const Link* from)	{
		return new DLinRegSubBranchMean(WithRoot(),chrono->GetBranchVal(from->GetBranch()), chrono->GetNodeVal(from->GetNode()), b1->GetNodeVal(from->Out()->GetNode()), b1->GetNodeVal(from->GetNode()), b2->GetNodeVal(from->Out()->GetNode()), b2->GetNodeVal(from->GetNode()), index, reg, drift1, phi1, drift2, phi2, agescale, kt, pulse, pulsephi);
	}

	int index;
	NodeBranchVarTree<PosReal,PosReal>* chrono;
	NodeVarTree<RealVector>* b1;
	NodeVarTree<RealVector>* b2;
	Var<RealVector>* reg;
	Var<RealVector>* drift1;
	Var<PosReal>* phi1;
	Var<RealVector>* drift2;
	Var<PosReal>* phi2;
	Var<PosReal>* agescale;
	double kt;
	Var<RealVector>* pulse;
	Var<PosReal>* pulsephi;
};

template <class R, class D> class DSemiConjugateMove : public MCUpdate	{

	public:

	DSemiConjugateMove(R* inrandom, D* indsemi, double intuning, int inn) : random(inrandom), dsemi(indsemi), tuning(intuning), n(inn) {}

	double Move(double tuning_modulator = 1)	{
		dsemi->ActivateSufficientStatistic();
		double total = 0;
		for (int i=0; i<n; i++)	{
			total += random->Move(tuning* tuning_modulator);
		}
		total /= n;
		dsemi->InactivateSufficientStatistic();
		return total;
	}

	protected:

	R* random;
	D* dsemi;
	double tuning;
	int n;
};

class DSemiConjugateMappingMove : public MCUpdate	{

	public:

	DSemiConjugateMappingMove(PhyloProcess* inprocess, PathConjugateTree* inpathconjtree) : process(inprocess), pathconjtree(inpathconjtree) {}

	double Move(double tuning_modulator=1)	{
		pathconjtree->InactivateSufficientStatistic();
		process->Move(1);
		pathconjtree->ActivateSufficientStatistic();
		return 1;
	}

	protected:

	PhyloProcess* process;
	PathConjugateTree* pathconjtree;
};


class MatrixTree : public BranchValPtrTree<RandomSubMatrix>	{


	public:

	MatrixTree(CodonStateSpace* instatespace, RandomSubMatrix* innucmatrix, LengthTree* inomegatree, Var<PosReal>* inrootomega) {
		SetWithRoot(true);
		omegatree = inomegatree;
		nucmatrix = innucmatrix;
		statespace = instatespace;
		rootomega = inrootomega;
		RecursiveCreate(GetRoot());
	}

	~MatrixTree()	{
		RecursiveDelete(GetRoot());
	}

	protected:

	RandomSubMatrix* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,rootomega);
		}
		return new RandomMGOmegaCodonSubMatrix(statespace,nucmatrix,omegatree->GetBranchVal(link->GetBranch()));
	}

	Tree* GetTree() {return omegatree->GetTree();}

	private:

	CodonStateSpace* statespace;
	LengthTree* omegatree;
	RandomSubMatrix* nucmatrix;
	Var<PosReal>* rootomega;

};


class BranchMatrixPhyloProcess : public PhyloProcess	{


	protected:

	public:

	BranchMatrixPhyloProcess(LengthTree* intree, BranchValPtrTree<RandomSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PhyloProcess(intree,indata)	{
		matrixtree = inmatrixtree;
	}

	virtual RandomBranchSitePath* 	CreateRandomBranchSitePath(const Link* link, int site)	{
		return  new RandomBranchSitePath(this,tree->GetBranchLength(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	protected:
	BranchValPtrTree<RandomSubMatrix>* matrixtree;
};


class LifeHistoryRegressionModel : public ProbModel {

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
	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	int Ninterpol;

	int chronoprior;
	double meanchi;
	double meanchi2;
	// 0 : uniform;
	// 1 : bd;
	// 2 : bd with cauchy proper lower bounds
	// 3 : rbd with cauchy proper lower bounds
	// 4 : coal

	double T0;
	Const<PosReal>* N0;
	Const<PosReal>* N1;
	Const<PosReal>* N2;
	Const<PosReal>* coald;
	Const<PosReal>* coalr;

	// chronogram
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;

	NodeBranchVarTree<PosReal,PosReal>* finalchrono;

	JeffreysIIDArray* ContDiagArray;
	JeffreysIIDArray* SubDiagArray;
	SigmaZero* ContSigmaZero;
	SigmaZero* SubSigmaZero;

	Rvar<CovMatrix>* contsigma;
	Rvar<CovMatrix>* subsigma;

	IIDNormalArray* reg;
	Gamma* phi;

	MultiVarNormal* drift;
	Gamma* driftphi;

	MultiVarNormal* drift2;
	Gamma* driftphi2;
	double kt;

	MultiVarNormal* pulse;
	Gamma* pulsephi;

	Const<RealVector>* controotmean;
	Const<PosRealVector>* controotvar;

	Const<RealVector>* subrootmean;
	Const<PosRealVector>* subrootvar;

	MultiVariateTreeProcess* contprocess;
	MultiVariateTreeProcess* subprocess;

	DLinRegCont* loglifehistory;
	DLinRegSub* synratetree;
	DLinRegSub* nonsynratetree;
	BranchValPtrTree<Dvar<PosReal> >* omegatree;

	// nucleotide mutation matrix is relrate * stationary
	Dirichlet* relrate;

	// for homogeneous model
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	// for both
	BranchValPtrTree<RandomSubMatrix>* matrixtree;
	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool clamproot;
	bool clamptree;

	int omegaratiotree;
	// 0 : rate model
	// 1 : dS omega
	// 2 : dS dN

	// total number of substitution parameters modeled as non homogeneous
	int L;

	bool conjpath;
	bool priorsampling;

	bool normalise;

	int nrep;

	int df;

	string bounds;
	bool withdrift;
	bool withexpdrift;
	bool withdoubleexpdrift;

	int withpulse;

	int clampsuffstat;
	string suffstatfile;

	bool autoregressive;

	public:

	SequenceAlignment* GetData()	{
		if (omegaratiotree)	{
			return codondata;
		}
		return nucdata;
	}

	SplitTree* GetSplitTree()	{
		SplitTree* tmp = dynamic_cast<SplitTree*>(splittree);
		if (! tmp)	{
			cerr << "error in GetSplitTree : null pointer\n";
			cerr << splittree << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	Var<PosReal>* GetRootScale()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale();
		}
		return mu;
	}

	bool Split()	{
		return Ninterpol != 1;
	}

	LifeHistoryRegressionModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double inN0, double inN1, double inN2, double inT0, double priorsigma, int indf, bool inclampdiag, bool inautoregressive, int inconjpath, int contdatatype, int inomegaratiotree, bool inclamproot, bool inclamptree, bool innormalise, int innrep, int inNinterpol, int inwithdrift, int inwithpulse, double inkt, string insuffstatfile, string rootfile, bool sample=true, GeneticCodeType type=Universal)	{

		Ninterpol = inNinterpol;
		kt = inkt;

		autoregressive = inautoregressive;
		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		if (inwithdrift == 3)	{
			withdrift = true;
			withdoubleexpdrift = true;
			withexpdrift = false;
		}
		else if (inwithdrift == 2)	{
			withdrift = true;
			withexpdrift = true;
			withdoubleexpdrift = false;
		}
		else if (inwithdrift == 1)	{
			withdrift = true;
			withexpdrift = false;
			withdoubleexpdrift = false;
		}
		else {
			withdrift = false;
			withdoubleexpdrift = false;
			withexpdrift = false;
		}

		withpulse = inwithpulse;

		if (withpulse || withdoubleexpdrift)	{
			if (kt == 0)	{
				cerr << "error : should specify a time boundary (KT) with double exp or pulse\n";
				exit(1);
			}
		}
		df = indf;

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		N0 = new Const<PosReal>(PosReal(inN0));
		N1 = new Const<PosReal>(PosReal(inN1));
		N2 = new Const<PosReal>(PosReal(inN2));
		T0 = inT0;
		coald = new Const<PosReal>(PosReal(meanchi));
		coalr = new Const<PosReal>(PosReal(meanchi2));

		clampdiag = inclampdiag;
		clamproot = inclamproot;
		clamptree = inclamptree;
		omegaratiotree = inomegaratiotree;

		if (omegaratiotree)	{
			L = 2;
		}
		else	{
			L = 1;
		}

		// get data from file

		nucdata = new FileSequenceAlignment(datafile);


		if (omegaratiotree)	{
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

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
			nrep = conjpath ? 30 : 1;
		}
		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		taxonset = nucdata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		cerr << "before : " << tree->GetFullSize(tree->GetRoot()) << '\n';
		cerr << Ninterpol << '\n';
		if (Split())	{
			/*
			if (withexpdrift || withdoubleexpdrift)	{
				cerr << "error : split and exp drift not yet compatible\n";
				exit(1);
			}
			cerr << "subdivide\n";
			// tree->Subdivide(tree->GetRoot(),Ninterpol);
			*/
			splittree = new SplitTree(tree,Ninterpol);
		}
		else	{
			splittree = tree;
		}
		cerr << "after  : " << splittree->GetFullSize(splittree->GetRoot()) << '\n';

		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			cerr << "error : life history model requires continuous data\n";
			exit(1);
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

		RootAlpha = 0;
		RootBeta = 0;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		PriorMu = new Const<PosReal>(1);
		mu = new Gamma(One,PriorMu);
		mu->ClampAt(1);

		if (calibfile != "None")	{
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,mu,RootAlpha,RootBeta,calibset);
			}
			else if (chronoprior == 4)	{
				chronogram = new CoalCalibratedChronogram(tree,mu,coald,coalr,N0,N1,N2,T0,calibset);
			}
			else {
				cerr << "BD\n";
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,RootAlpha,RootBeta,calibset,chronoprior);
			}
		}
		else	{
			chronogram = new Chronogram(tree,mu);
		}

		if (clamptree)	{
			chronogram->Clamp();
		}

		if (Split())	{
			finalchrono = new SplitChronogram(chronogram,GetSplitTree());
		}
		else	{
			finalchrono = chronogram;
		}

		double mindiag = 0.001;
		double maxdiag = 1000;

		ContDiagArray = new JeffreysIIDArray(Ncont,mindiag,maxdiag,Zero);
		SubDiagArray = new JeffreysIIDArray(L,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			ContDiagArray->setval(1.0);
			SubDiagArray->setval(1.0);
		}
		else	{
			ContDiagArray->ClampAt(priorsigma);
			SubDiagArray->ClampAt(priorsigma);
		}

		ContSigmaZero = new SigmaZero(ContDiagArray);
		SubSigmaZero = new SigmaZero(SubDiagArray);

		cerr << "sigma\n";
		if (clampdiag)	{
			contsigma = new DiagonalCovMatrix(ContSigmaZero, ContSigmaZero->GetDim()+df);
		}
		else	{
			contsigma = new InverseWishartMatrix(ContSigmaZero, ContSigmaZero->GetDim()+df);
		}
		subsigma = new DiagonalCovMatrix(SubSigmaZero, SubSigmaZero->GetDim()+df);

		controotmean = new Const<RealVector>(RealVector(Ncont));
		controotvar = new Const<PosRealVector>(PosRealVector(Ncont));
		subrootmean = new Const<RealVector>(RealVector(L));
		subrootvar = new Const<PosRealVector>(PosRealVector(L));
		if (rootfile != "None")	{
			ifstream is(rootfile.c_str());
			cerr << "root hyper for cont\n";
			for (int i=0; i<Ncont; i++)	{
				double mean, var;
				is >> mean >> var;
				cerr << i << '\t' << mean << '\t' << var << '\n';
				(*controotmean)[i] = mean;
				(*controotvar)[i] = var;
			}
			cerr << "root hyper for sub\n";
			for (int i=0; i<L; i++)	{
				double mean, var;
				is >> mean >> var;
				cerr << i << '\t' << mean << '\t' << var << '\n';
				(*subrootmean)[i] = mean;
				(*subrootvar)[i] = var;
			}
		}
		else	{
			cerr << "default root\n";
			for (int i=0; i<Ncont; i++)	{
				(*controotmean)[i] = 0;
				(*controotvar)[i] = 2;
			}
			for (int i=0; i<L; i++)	{
				(*subrootmean)[i] = 0;
				(*subrootvar)[i] = 2;
			}
		}

		cerr << "drift\n";

		drift = new MultiVarNormal(Zero, ContDiagArray);
		drift->SetAtZero();
		driftphi = new Gamma(One,One);

		drift2 = new MultiVarNormal(Zero, ContDiagArray);
		drift2->SetAtZero();
		driftphi2 = new Gamma(One,One);

		pulse = new MultiVarNormal(Zero, SubDiagArray);
		pulsephi = new Gamma(One,One);
		pulse->SetAtZero();
		if (! withpulse)	{
			pulse->ClampAtZero();
			pulsephi->ClampAt(1);
		}

		if (withdoubleexpdrift)	{
			// nothing
		}
		else if (withexpdrift)	{
			drift2->ClampAtZero();
			driftphi2->ClampAt(1);
		}
		else if (withdrift)	{
			driftphi->ClampAt(1);
			drift2->ClampAtZero();
			driftphi2->ClampAt(1);
		}
		else	{
			drift->ClampAtZero();
			driftphi->ClampAt(1);
			drift2->ClampAtZero();
			driftphi2->ClampAt(1);
		}

		cerr << "regarray\n";
		if (autoregressive)	{
			phi = new Gamma(One,One);
			reg = new IIDNormalArray(L,Ncont+1,Zero,One,One);
		}
		else	{
			phi = 0;
			reg = new IIDNormalArray(L,Ncont,Zero,One,One);
		}
		reg->SetAtZero();
		if (clampdiag)	{
			reg->ClampAtZero();
		}

		contprocess = new MultiVariateTreeProcess(contsigma, finalchrono, 0, 0, controotmean, controotvar);
		cerr << "set and clamp\n";
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				contprocess->SetAndClamp(contdata,i,i,contdatatype);
			}
		}

		subprocess = new MultiVariateTreeProcess(subsigma, chronogram, 0, 0, subrootmean, subrootvar);
		for (int l=0; l<L; l++)	{
			subprocess->CutOff(1,l);
		}

		// loglifetraits = new DLinRegCont();

		CreateSubstitutionProcess();

		if (phyloprocess)	{
			phyloprocess->Unfold();
			if (sample)	{
				phyloprocess->Sample();
			}
		}

		// register model
		RootRegister(PriorMu);
		RootRegister(Zero);
		RootRegister(One);
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		if (chronoprior == 4)	{
			RootRegister(N0);
			RootRegister(N1);
			RootRegister(N2);
			RootRegister(coald);
			RootRegister(coalr);
		}
		else if (chronoprior)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		RootRegister(relrate);
		RootRegister(stationary);
		RootRegister(controotmean);
		RootRegister(controotvar);
		RootRegister(subrootmean);
		RootRegister(subrootvar);
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

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~LifeHistoryRegressionModel() {}

	void CreateSubstitutionProcess()	{

		synratetree = new DLinRegSub(false,finalchrono,contprocess,subprocess,0,reg->GetVal(0),drift,driftphi,drift2,driftphi2,GetRootScale(),kt,pulse,pulsephi);
		// synratetree->CutOff(1);
		nonsynratetree = 0;
		if (omegaratiotree == 0)	{
			omegatree = 0;
		}
		else if (omegaratiotree == 1)	{
			DLinRegSub* tmp = new DLinRegSub(true,finalchrono,contprocess,subprocess,1,reg->GetVal(1),drift,driftphi,drift2,driftphi2,GetRootScale(),kt,pulse,pulsephi);
			// tmp->CutOff(1);
			omegatree = tmp;
		}
		else	{
			nonsynratetree = new DLinRegSub(false,finalchrono,contprocess,subprocess,1,reg->GetVal(1),drift,driftphi,drift2,driftphi2,GetRootScale(),kt,pulse,pulsephi);
			// nonsynratetree->CutOff(1);
			omegatree = new RatioTree(nonsynratetree,synratetree);
		}

		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);

		if (! omegaratiotree)	{
			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);

			// make substitution mappings
			if (conjpath)	{
				pathconjtree = new OneMatrixPathConjugateTree(synratetree,nucmatrix,GetData());

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
				if (priorsampling)	{
					phyloprocess = 0;
				}
				phyloprocess = new OneMatrixPhyloProcess(synratetree, nucmatrix, GetData());
			}
		}
		else	{

			stationary = new Dirichlet(Nnuc);
			nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,normalise);
			matrixtree = new MatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrix, omegatree, One);

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
				if (priorsampling)	{
					phyloprocess = 0;
				}
				else	{
					phyloprocess = new BranchMatrixPhyloProcess(synratetree, matrixtree, codondata);
				}
			}
		}
	}

	Tree* GetTree() {return tree;}
	Tree* GetFineGrainedTree() {return splittree;}

	void UpdateOmegaTree()	{
		if (omegaratiotree == 2)	{
			nonsynratetree->specialUpdate();
			((RatioTree*) omegatree)->specialUpdate();
		}
		else if (omegaratiotree == 1)	{
			((DLinRegSub*) omegatree)->specialUpdate();
		}
		else	{
			cerr << "error : update omega tree called under simple rate model\n";
			exit(1);
		}
	}

	DLinRegSub* GetSynRateTree() {
		DLinRegSub* tmp = dynamic_cast<DLinRegSub*>(synratetree);
		if (! tmp)	{
			cerr << "error in get synratetree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	DLinRegSub* GetOmegaTree() {
		if (! omegaratiotree)	{
			cerr << "error : get omega tree: model is not codon\n";
			exit(1);
		}
		DLinRegSub* tmp = dynamic_cast<DLinRegSub*>(omegatree);
		if (! tmp)	{
			cerr << "error in get omegatree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	MultiVariateTreeProcess* GetContProcess() {return contprocess;}
	MultiVariateTreeProcess* GetSubProcess() {return subprocess;}
	Chronogram* GetChronogram() {
		return chronogram;
	}
	LengthTree* GetLengthTree() {return finalchrono;}
	SplitLengthTree* GetSplitLengthTree()	{
		SplitLengthTree* tmp = dynamic_cast<SplitLengthTree*>(finalchrono);
		if (! tmp)	{
			cerr << "error in GetSplitLengthTree : null pointer\n";
			cerr << finalchrono << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}

	CovMatrix* GetContMatrix() {return contsigma;}
	CovMatrix* GetSubMatrix() {return subsigma;}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double GetPhi()	{
		if (phi)	{
			return phi->val();
		}
		return 0;
	}

	double Update(bool check = false)	{
		double ret = ProbModel::Update();
		/*
		if (phyloprocess)	{
			phyloprocess->Sample();
			ret = ProbModel::Update();
		}
		*/
		return ret;
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();

		total += ContDiagArray->GetLogProb();
		total += SubDiagArray->GetLogProb();
		total += contsigma->GetLogProb();
		total += subsigma->GetLogProb();
		total += reg->GetLogProb();
		if (autoregressive)	{
			total += phi->GetLogProb();
		}
		total += drift->GetLogProb();
		total += driftphi->GetLogProb();
		total += drift2->GetLogProb();
		total += driftphi2->GetLogProb();
		total += pulse->GetLogProb();
		total += pulsephi->GetLogProb();

		total += contprocess->GetLogProb();
		total += subprocess->GetLogProb();

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		if (priorsampling)	{
			return 0;
		}
		else if (clampsuffstat)	{
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
				scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
			}
		}
		else	{
			if (phyloprocess)	{
				scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
			}
		}

		vector <SplitMultiVariateNodeMove*> nodesplitarray;
		vector <SplitMultiVariateBranchMove*> branchsplitarray;
		if (Split())	{
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,10));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,0.1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(contprocess,0.01));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,10));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,0.1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,contprocess,0.01));
		}

		for (int i=0; i<nrep; i++)	{
			if (! clamptree)	{
				if ((chronoprior >= 1) && (chronoprior <= 3))	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
			}

			if (RootAlpha && (! clamptree))	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
			}

			if (Split())	{

				scheduler.Register(nodesplitarray[0],10,"split node process");
				scheduler.Register(nodesplitarray[1],10,"split node process");
				scheduler.Register(nodesplitarray[2],10,"split node process");
				scheduler.Register(nodesplitarray[3],10,"split node process");
				scheduler.Register(branchsplitarray[0],10,"split node process");
				scheduler.Register(branchsplitarray[1],10,"split node process");
				scheduler.Register(branchsplitarray[2],10,"split node process");
				scheduler.Register(branchsplitarray[3],10,"split node process");

			}

			scheduler.Register(new SimpleMove(drift,10),10,"drift");
			scheduler.Register(new SimpleMove(drift,1),10,"drift");
			scheduler.Register(new SimpleMove(drift,0.1),10,"drift");
			scheduler.Register(new SimpleMove(drift,0.01),10,"drift");

			scheduler.Register(new SimpleMove(driftphi,10),10,"drift phi");
			scheduler.Register(new SimpleMove(driftphi,1),10,"drift phi");
			scheduler.Register(new SimpleMove(driftphi,0.1),10,"drift phi");
			scheduler.Register(new SimpleMove(driftphi,0.01),10,"drift phi");

			scheduler.Register(new SimpleMove(drift2,10),10,"drift2");
			scheduler.Register(new SimpleMove(drift2,1),10,"drift2");
			scheduler.Register(new SimpleMove(drift2,0.1),10,"drift2");
			scheduler.Register(new SimpleMove(drift2,0.01),10,"drift2");

			scheduler.Register(new SimpleMove(driftphi2,10),10,"drift phi2");
			scheduler.Register(new SimpleMove(driftphi2,1),10,"drift phi2");
			scheduler.Register(new SimpleMove(driftphi2,0.1),10,"drift phi2");
			scheduler.Register(new SimpleMove(driftphi2,0.01),10,"drift phi2");

			scheduler.Register(new SimpleMove(pulse,10),10,"pulse");
			scheduler.Register(new SimpleMove(pulse,1),10,"pulse");
			scheduler.Register(new SimpleMove(pulse,0.1),10,"pulse");
			scheduler.Register(new SimpleMove(pulse,0.01),10,"pulse");

			scheduler.Register(new SimpleMove(pulsephi,10),10,"pulse phi");
			scheduler.Register(new SimpleMove(pulsephi,1),10,"pulse phi");
			scheduler.Register(new SimpleMove(pulsephi,0.1),10,"pulse phi");
			scheduler.Register(new SimpleMove(pulsephi,0.01),10,"pulse phi");

			scheduler.Register(new SimpleMove(contprocess,10),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(contprocess,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(subprocess,10),10,"sub multinormal");
			scheduler.Register(new SimpleMove(subprocess,1),10,"sub multinormal");
			scheduler.Register(new SimpleMove(subprocess,0.1),10,"sub multinormal");
			scheduler.Register(new SimpleMove(subprocess,0.01),10,"submultinormal");

			scheduler.Register(new SimpleMove(contsigma,10),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,0.1),100,"cont sigma");
			scheduler.Register(new SimpleMove(contsigma,0.01),100,"cont sigma");

			scheduler.Register(new SimpleMove(subsigma,10),100,"sub sigma");
			scheduler.Register(new SimpleMove(subsigma,1),100,"sub sigma");
			scheduler.Register(new SimpleMove(subsigma,0.1),100,"sub sigma");
			scheduler.Register(new SimpleMove(subsigma,0.01),100,"sub sigma");

			scheduler.Register(new SimpleMove(ContDiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(ContDiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(ContDiagArray,0.1),10,"theta");

			scheduler.Register(new SimpleMove(SubDiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(SubDiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(SubDiagArray,0.1),10,"theta");

			scheduler.Register(new SimpleMove(reg,10),10,"reg");
			scheduler.Register(new SimpleMove(reg,1),10,"reg");
			scheduler.Register(new SimpleMove(reg,0.1),10,"reg");
			scheduler.Register(new SimpleMove(reg,0.01),10,"reg");

			if (autoregressive)	{
				scheduler.Register(new SimpleMove(phi,10),100,"phi");
				scheduler.Register(new SimpleMove(phi,1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.01),100,"phi");
			}

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
		}
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "sample\n";

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			Chi->Sample();
			Chi2->Sample();
		}
		chronogram->Sample();

		ContDiagArray->Sample();
		SubDiagArray->Sample();

		contsigma->Sample();
		subsigma->Sample();
		reg->Sample();
		phi->Sample();
		// sigma->SetIdentity();
		drift->Sample();
		driftphi->Sample();
		drift2->Sample();
		driftphi2->Sample();
		pulse->Sample();
		pulsephi->Sample();

		contprocess->Sample();
		subprocess->Sample();

		GetSynRateTree()->specialUpdate();
		if (omegaratiotree)	{
			UpdateOmegaTree();
		}

		relrate->Sample();
		stationary->Sample();

		if (phyloprocess)	{
			phyloprocess->Sample();
		}

		cerr << "ok\n";
	}

	double GetRootAge()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	/*
	double GetMaxdS()	{
		if (! omegaratiotree)	{
			cerr << "error : get max dS called under simple rate model\n";
			exit(1);
		}
		return GetSynRateTree()->GetMax();
	}

	double GetMaxdN()	{
		if (! omegaratiotree)	{
			cerr << "error : get max dN called under simple rate model\n";
			exit(1);
		}
		if (! nonsynratetree)	{
			cerr << "error : cannot call getmaxdN\n";
			exit(1);
		}
		return nonsynratetree->GetMax();
	}

	double GetMaxOmega()	{
		if (! omegaratiotree)	{
			cerr << "error : get max omega called under simple rate model\n";
			exit(1);
		}
		return ((MeanExpTreeFromMultiVariate*) omegatree)->GetMax();
	}
	*/

	double GetMeanSynRate()	{
		return GetSynRateTree()->GetTotal();
	}

	double GetTotalTime()	{
		return chronogram->GetTotalTime();
	}

	double GetMeanOmega()	{
		if (! omegaratiotree)	{
			cerr << "error : get mean omega called under simple rate model\n";
			exit(1);
		}
		if (omegaratiotree)	{
			return ((RatioTree*) omegatree)->GetMean();
		}
		else	{
			return ((MeanExpTreeFromMultiVariate*) omegatree)->GetMean();
		}
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		if (omegaratiotree)	{
			os << "\tsynrate\tomega";
		}
		else	{
			os << "\trate";
		}

		// os << "\trootleft\trootright";

		for (int k=0; k<Ncont; k++)	{
			for (int l=k+1; l<Ncont; l++)	{
				os << '\t' << "cont_" << k << '_' << l;
			}
		}
		for (int k=0; k<L; k++)	{
			for (int l=0; l<reg->GetDim(); l++)	{
				os << '\t' << "reg_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "cont_" << k << '_' << k;
		}
		for (int k=0; k<L; k++)	{
			os << '\t' << "sub_" << k << '_' << k;
		}
		if (autoregressive)	{
			os << '\t' << "phi";
		}
		if (withdrift)	{
			os << "\tdim";
			for (int k=0; k<Ncont; k++)	{
				os << '\t' << "drift_" << k;
			}
			if (withdoubleexpdrift)	{
				os << '\t' << "phi";
				os << '\t' << "dim";
				for (int k=0; k<Ncont; k++)	{
					os << '\t' << "drift2_" << k;
				}
				os << '\t' << "phi2";
			}
			else if (withexpdrift)	{
				os << '\t' << "phi";
			}
		}
		if (withpulse)	{
			os << "\tdim";
			for (int k=0; k<Ncont; k++)	{
				os << '\t' << "pulse_" << k;
			}
			os << '\t' << "pulsephi";
		}
		if (RootAlpha)	{
			os << "\trootage";
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tp1\tp2";
		}

		os << "\tdim";
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "root_" << k;
		}

		os << "\tdim";
		for (int k=0; k<L; k++)	{
			os << '\t' << "root_" << k;
		}

		os << "\tstatent";
		os << "\trrent";

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tnumerror";
		}

		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		if (omegaratiotree)	{
			os << '\t' << GetMeanOmega();
		}

		for (int k=0; k<Ncont; k++)	{
			for (int l=k+1; l<Ncont; l++)	{
				os << '\t' << (*contsigma)[k][l];
			}
		}
		for (int k=0; k<L; k++)	{
			for (int l=0; l<reg->GetDim(); l++)	{
				os << '\t' << (*(*reg)[k])[l];
			}
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << (*contsigma)[k][k];
		}
		for (int k=0; k<L; k++)	{
			os << '\t' << (*subsigma)[k][k];
		}
		if (autoregressive)	{
			os << '\t' << phi->val();
		}
		if (withdrift)	{
			os << '\t' << *drift;
			if (withdoubleexpdrift)	{
				os << '\t' << *driftphi;

				os << '\t' << *drift2;
				os << '\t' << *driftphi2;
			}
			else if (withexpdrift)	{
				os << '\t' << *driftphi;
			}
		}
		if (withpulse)	{
			os << '\t' << *pulse;
			os << '\t' << *pulsephi;
		}
		if (RootAlpha)	{
			os << '\t' << GetRootAge();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}
		os << '\t' << GetContProcess()->GetMultiNormal(GetFineGrainedTree()->GetRoot())->val();
		os << '\t' << GetSubProcess()->GetNodeVal(GetFineGrainedTree()->GetRoot()->GetNode())->val();

		os << '\t' << stationary->val().GetEntropy();

		os << '\t' << relrate->val().GetEntropy();

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << BDCalibratedChronogram::NumErrorCount;
		}

		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *mu << '\n';
		os << *chronogram << '\n';
		if (RootAlpha)	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		os << *ContDiagArray << '\n';
		os << *SubDiagArray << '\n';
		os << *contsigma << '\n';
		os << *subsigma << '\n';
		os << *reg<< '\n';
		if (autoregressive)	{
			os << *phi << '\n';
		}
		os << *drift << '\n';
		if (withdoubleexpdrift)	{
			os << *driftphi << '\n';

			os << *driftphi2 << '\n';
			os << *drift2 << '\n';
		}
		else if (withexpdrift)	{
			os << *driftphi << '\n';
		}
		if (withpulse)	{
			os << *pulse << '\n';
			os << *pulsephi << '\n';
		}
		os << '\n';
		os << *contprocess << '\n';
		os << *subprocess << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is)	{
		is >> *mu;
		is >> *chronogram;
		if (RootAlpha)	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			is >> *Chi >> *Chi2;
		}

		is >> *ContDiagArray;
		is >> *SubDiagArray;
		is >> *contsigma;
		is >> *subsigma;
		is >> *reg;
		if (autoregressive)	{
			is >>  *phi;
		}
		is >> *drift;
		if (withdoubleexpdrift)	{
			is >> *driftphi;

			is >> *driftphi2;
			is >> *drift2;
		}
		else if (withexpdrift)	{
			is >> *driftphi;
		}
		if (withpulse)	{
			is >> *pulse;
			is >> *pulsephi;
		}
		is >> *contprocess;
		is >> *subprocess;
		is >> *relrate;
		is >> *stationary;
	}

	void GetNucMatrix(ifstream& is)	{
		is >> *relrate;
		is >> *stationary;
		cerr << *relrate << '\n';
		cerr << *stationary << '\n';
	}
};

#endif
