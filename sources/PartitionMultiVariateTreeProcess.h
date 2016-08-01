
#ifndef PARTMV_H
#define PARTMV_H

#include "ConjugateMultiVariateTreeProcess.h"
#include "Partition.h"
#include "ValArray.h"


class PartitionMultiVariateTreeProcess : public virtual MultiVariateTreeProcess	{

	public:

	PartitionMultiVariateTreeProcess() {}

	PartitionMultiVariateTreeProcess(VarArray<CovMatrix>* insigmaarray, LengthTree* intree, BranchPartition* inparttree = 0, LengthTree* inscaletree = 0, VarArray<RealVector>* indriftarray = 0, VarArray<PosReal>* indriftphiarray = 0, NodeBranchVarTree<PosReal,PosReal>* inchrono = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0, VarArray<RealVector>* indriftarray2 = 0, VarArray<PosReal>* indriftphiarray2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0, GenericTimeLine* intimeline = 0, Var<Real>* inalpha = 0) {

		sigma = 0;
		drift = 0;
		rootmean = inrootmean;
		rootvar = inrootvar;
		sigmaarray = insigmaarray;
		driftarray = indriftarray;
		driftphiarray = indriftphiarray;
		driftarray2 = indriftarray2;
		driftphiarray2 = indriftphiarray2;
		agescale = inagescale;
		kt = inkt;
		chrono = inchrono;
		tree = intree;
		parttree = inparttree;
		scaletree = inscaletree;
		timeline = intimeline;
		alpha = inalpha;
		RecursiveCreate(GetRoot());
	}

	~PartitionMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}

	int GetDim()	{
		return sigmaarray->GetVal(0)->GetDim();
	}

	int GetNmatrix()	{
		return sigmaarray->GetSize();
	}

	int GetNdrift()	{
		return driftarray->GetSize();
	}

	protected:

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			int alloc = parttree->GetAlloc(link->GetBranch());
			Var<CovMatrix>* sigma = 0;
			if (sigmaarray->GetSize() == 1)	{
				sigma = sigmaarray->GetVal(0);
			}
			else	{
				sigma = sigmaarray->GetVal(alloc);
			}
			Var<RealVector>* drift = 0;
			if (driftarray)	{
				drift = driftarray->GetVal(alloc);
			}
			Var<PosReal>* driftphi = 0;
			Var<RealVector>* drift2 = 0;
			Var<PosReal>* driftphi2 = 0;
			Var<PosReal>* date = 0;
			if (driftphiarray)	{
				driftphi = driftphiarray->GetVal(alloc);
				date = chrono->GetNodeVal(link->GetNode());
				if (! date)	{
					cerr << "error in paertition mv proc: mull date\n";
					cerr << link << '\t' << link->Next() << '\t' << link->Next()->Next() << '\t' << link->Next()->Next()->Next() << '\n';
					exit(1);
				}
				if (driftphiarray2)	{
					drift2 = driftarray2->GetVal(alloc);
					driftphi2 = driftphiarray2->GetVal(alloc);
				}
			}
			return new MultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), scaletree ? scaletree->GetBranchVal(link->GetBranch()): 0, drift, driftphi,date, drift2, driftphi2, agescale, kt, timeline, alpha);
		}
		else{
			Var<CovMatrix>* sigma = sigmaarray->GetVal(0);
 			return new MultiNormal(0,sigma,rootmean,rootvar);
		}
	}

	VarArray<CovMatrix>* sigmaarray;
	VarArray<RealVector>* driftarray;
	VarArray<PosReal>* driftphiarray;
	VarArray<RealVector>* driftarray2;
	VarArray<PosReal>* driftphiarray2;
	// Var<PosReal>* agescale;
	double kt;
	// NodeBranchVarTree<PosReal,PosReal>* chrono;
	BranchPartition* parttree;
	GenericTimeLine* timeline;
	Var<Real>* alpha;
};


class ConjugatePartitionMultiVariateTreeProcess : public PartitionMultiVariateTreeProcess, ConjugateMultiVariateTreeProcess	{

	public:

	ConjugatePartitionMultiVariateTreeProcess(VarArray<CovMatrix>* insigmaarray, LengthTree* intree, BranchPartition* inparttree = 0, LengthTree* inscaletree = 0, VarArray<RealVector>* indriftarray = 0, VarArray<PosReal>* indriftphiarray = 0, NodeBranchVarTree<PosReal,PosReal>* inchrono = 0, Var<RealVector>* inrootmean =  0, Var<PosRealVector>* inrootvar = 0, VarArray<RealVector>* indriftarray2 = 0, VarArray<PosReal>* indriftphiarray2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0, GenericTimeLine* intimeline = 0, Var<Real>* inalpha = 0) {
		sigma = 0;
		drift = 0;
		rootmean = inrootmean;
		rootvar = inrootvar;
		sigmaarray = insigmaarray;
		driftarray = indriftarray;
		driftphiarray = indriftphiarray;
		driftarray2 = indriftarray2;
		driftphiarray2 = indriftphiarray2;
		agescale = inagescale;
		kt = inkt;
		chrono = inchrono;
		tree = intree;
		parttree = inparttree;
		scaletree = inscaletree;
		timeline = intimeline;
		alpha = inalpha;

		RecursiveCreate(GetRoot());
	}

	~ConjugatePartitionMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}

	ConjugateInverseWishart* GetConjugateInverseWishart(int mat)	{
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(sigmaarray->GetVal(mat));
		if (! tmp)	{
			cerr << "error in conjugate partition multivariate tree process: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			int alloc = parttree->GetAlloc(link->GetBranch());
			ConjugateInverseWishart* sigma = 0;
			if (sigmaarray->GetSize() == 1)	{
				sigma = GetConjugateInverseWishart(0);
			}
			else	{
				sigma = GetConjugateInverseWishart(alloc);
			}
			Var<RealVector>* drift = 0;
			if (driftarray)	{
				drift = driftarray->GetVal(alloc);
			}
			Var<PosReal>* driftphi = 0;
			Var<RealVector>* drift2 = 0;
			Var<PosReal>* driftphi2 = 0;
			Var<PosReal>* date = 0;
			if (driftphiarray)	{
				driftphi = driftphiarray->GetVal(0);
				date = chrono->GetNodeVal(link->GetNode());
				if (! date)	{
					cerr << "error in paertition mv proc: mull date\n";
					exit(1);
				}
				if (driftphiarray2)	{
					drift2 = driftarray2->GetVal(alloc);
					driftphi2 = driftphiarray2->GetVal(alloc);
				}
			}
			// return new ConjugateMultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), scaletree ? scaletree->GetBranchVal(link->GetBranch()): 0, drift, driftphi, date);
			return new ConjugateMultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), scaletree ? scaletree->GetBranchVal(link->GetBranch()): 0, drift, driftphi,date, drift2, driftphi2, agescale, kt, timeline,alpha);
		}
		else{
			ConjugateInverseWishart* sigma = GetConjugateInverseWishart(0);
 			// return new ConjugateMultiNormal(sigma);
 			return new ConjugateMultiNormal(0,sigma,rootmean,rootvar);
		}
	}
};


class ConjugatePartitionMultiVariateMove : public MCUpdate {

	ConjugatePartitionMultiVariateTreeProcess* tree;
	int nmove;
	double tuning;

	public:

	ConjugatePartitionMultiVariateMove(ConjugatePartitionMultiVariateTreeProcess* intree, double intuning, int innmove){
		tree = intree;
		tuning = intuning;
		nmove =innmove;
	}

	ConjugateInverseWishart* GetConjugateInverseWishart(int mat)	{
		return tree->GetConjugateInverseWishart(mat);
	}

	int GetNmatrix()	{
		return tree->GetNmatrix();
	}

	double Move(double tuning_modulator){
		for (int mat=0; mat<GetNmatrix(); mat++)	{
			GetConjugateInverseWishart(mat)->Integrate();
		}
		double r = 0;
		for(int i=0; i<nmove; i++){
			r += tree->Move(tuning * tuning_modulator);
		}
		for (int mat=0; mat<GetNmatrix(); mat++)	{
			GetConjugateInverseWishart(mat)->Resample();
		}
		return nmove ? r/nmove : nmove;
	}
};

#endif

