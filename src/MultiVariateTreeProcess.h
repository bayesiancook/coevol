#ifndef MULTIVARIATETREEPROCESS
#define MULTIVARIATETREEPROCESS

#include "InverseWishartMatrix.h"
#include "ValTree.h"
#include "ContinuousData.h"
#include "AncestralData.h"
#include "MatrixAlgebra.h"



class MultiVariateTreeProcess : public MCMC , public NodeValPtrTree<Rvar<RealVector> > {

	protected:

	Var<CovMatrix>* sigma;
	LengthTree* tree;
	LengthTree* scaletree;
	Var<RealVector>* drift;
	Var<RealVector>* rootmean;
	Var<PosRealVector>* rootvar;

	Var<PosReal>* agescale;
	NodeBranchVarTree<PosReal,PosReal>* chrono;
	GlobalScalingFunction* scalefunction;

	public:


	MultiVariateTreeProcess()	{
		scaletree = 0;
		drift = 0;
		scalefunction = 0;
	}


	MultiVariateTreeProcess(Var<CovMatrix>* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0){
		sigma = insigma;
		tree = intree;
		scaletree = inscaletree;
		drift = indrift;
		rootmean = inrootmean;
		rootvar = inrootvar;
		scalefunction = 0;
		RecursiveCreate(GetRoot());
	}

	MultiVariateTreeProcess(Var<CovMatrix>* insigma, NodeBranchVarTree<PosReal,PosReal>* inchrono, Var<PosReal>* inagescale, GlobalScalingFunction* inscalefunction, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0){
		sigma = insigma;
		tree = inchrono;
		chrono = inchrono;
		scaletree = 0;
		drift = 0;
		rootmean = inrootmean;
		rootvar = inrootvar;
		scalefunction = inscalefunction;
		agescale = inagescale;
		RecursiveCreate(GetRoot());
	}

	~MultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	Tree* GetTree(){
		return tree->GetTree();
	}

	LengthTree* GetLengthTree(){
		return tree;
	}

	double GetBranchLength(const Branch* branch)	{
		if (! branch)	{
			return 0;
		}
		return tree->GetBranchVal(branch)->val();
	}

	void SetBounds(const BoundSet& inboundset, int offset)	{
		for (map<const Node*,BoundForMultiNormal>::const_iterator i = inboundset.GetNodeMap().begin(); i!=inboundset.GetNodeMap().end(); i++)	{
			GetMultiNormal(i->first)->SetBound(i->second,offset);
		}
		BoundCutOff();
	}

	MultiNormal* GetMultiNormal(const Link* link)	{
		MultiNormal* m = dynamic_cast<MultiNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	MultiNormal* GetMultiNormal(const Node* node)	{
		MultiNormal* m = dynamic_cast<MultiNormal*> (GetNodeVal(node));
		return m;
	}

	double GetExpVal(const Link* link, int index)	{
		return exp((*GetNodeVal(link->GetNode()))[index]);
	}

	double GetVal(const Link* link, int index)	{
		return (*GetNodeVal(link->GetNode()))[index];
	}

	virtual int GetDim(){
		return sigma->GetDim();
	}


	void drawSample()	{
		RecursivedrawSample(this->GetRoot());
	}

	int GetNnode()	{
		return GetTree()->GetNnode();
	}

	int GetNinternalNode()	{
		return GetTree()->GetNinternalNode();
	}

	void RegisterNode(DAGnode* node)	{
		RecursiveRegister(node,GetRoot());
	}

	double* GetNodeVals(double* ptr)	{
		return GetNodeVals(GetRoot(),ptr);
	}

	double* GetNodeVals(const Link* from, double* ptr)	{
		for (int i=0; i<GetDim(); i++)	{
			*ptr = GetNodeVal(from->GetNode())->val()[i];
			ptr++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ptr = GetNodeVals(link->Out(), ptr);
		}
		return ptr;
	}

	double* SetNodeVals(double* ptr)	{
		return SetNodeVals(GetRoot(),ptr);
	}

	double* SetNodeVals(const Link* from, double* ptr)	{
		for (int i=0; i<GetDim(); i++)	{
			GetNodeVal(from->GetNode())->val()[i] = *ptr;
			ptr++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ptr = SetNodeVals(link->Out(), ptr);
		}
		return ptr;
	}

	double GetMin(int index)	{
		double current = (*GetNodeVal(GetRoot()->GetNode()))[index];
		double min = RecursiveGetMin(GetRoot(),index,current);
		return min;
	}

	double RecursiveGetMin(const Link* from, int index, double current)	{
		double tmp = (*GetNodeVal(from->GetNode()))[index];
		if (current > tmp)	{
			current = tmp;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMin(link->Out(),index,current);
			if (current > tmp)	{
				current = tmp;
			}
		}
		return current;
	}

	double GetMax(int index)	{
		double current = (*GetNodeVal(GetRoot()->GetNode()))[index];
		double max = RecursiveGetMax(GetRoot(),index,current);
		return max;
	}

	double RecursiveGetMax(const Link* from, int index, double current)	{
		double tmp = (*GetNodeVal(from->GetNode()))[index];
		if (current < tmp)	{
			current = tmp;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveGetMax(link->Out(),index,current);
			if (current < tmp)	{
				current = tmp;
			}
		}
		return current;
	}

	void GetEmpiricalCovariance(double** cov)	{

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				cov[i][j] =0;
			}
		}
		int n = 0;
		RecursiveGetEmpiricalCovariance(GetRoot(),cov,n);
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				cov[i][j] /= n;
			}
		}
	}

	void RecursiveGetEmpiricalCovariance(const Link* from, double** cov, int& n)	{
		if (! from->isRoot())	{
			double temp[GetDim()];
			for (int i=0; i<GetDim(); i++)	{
				temp[i] = ((*GetNodeVal(from->GetNode()))[i] - (*GetNodeVal(from->Out()->GetNode()))[i]) / sqrt(GetBranchLength(from->GetBranch()));
				for (int i=0; i<GetDim(); i++)	{
					for (int j=0; j<GetDim(); j++)	{
						cov[i][j] += temp[i] * temp[j];
					}
				}
			}
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetEmpiricalCovariance(link->Out(),cov,n);
		}
	}

	void SampleFixRoot()	{
		RecursiveSampleFixRoot(this->GetRoot());
	}

	void ClampLeaves(int index)	{
		RecursiveClampLeaves(GetRoot(),index);
	}

	void RecursiveClampLeaves(const Link* from, int index)	{
		if (from->isLeaf())	{
			MultiNormal* m = GetMultiNormal(from);
			m->Clamp(index);
		}
		else	{
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				RecursiveClampLeaves(link->Out(),index);
			}
		}
	}

	void ClampTransversal(int index, double d)	{
		RecursiveClampTransversal(GetRoot(),index,d);
	}

	void RecursiveClampTransversal(const Link* from, int index, double d)	{
		MultiNormal* m = GetMultiNormal(from);
		(*m)[index] = d;
		m->Clamp(index);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveClampTransversal(link->Out(),index,d);
		}
	}

	void BoundCutOff()	{
		RecursiveBoundCutOff(GetRoot());
	}

	void RecursiveBoundCutOff(const Link* from)	{
		GetMultiNormal(from)->CheckBounds();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveBoundCutOff(link->Out());
		}
	}

	void LowerCutOff(double cutoff, int index)	{
		RecursiveLowerCutOff(GetRoot(),cutoff,index);
	}

	void RecursiveLowerCutOff(const Link* from, double cutoff, int index)	{
		if ((! GetMultiNormal(from)->ClampVector[index]) && ((*GetNodeVal(from->GetNode()))[index] < cutoff))	{
			(*GetNodeVal(from->GetNode()))[index] = cutoff;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveLowerCutOff(link->Out(),cutoff,index);
		}
	}

	void CutOff(double cutoff, int index)	{
		RecursiveCutOff(GetRoot(),cutoff,index);
	}

	void RecursiveCutOff(const Link* from, double cutoff, int index)	{
		if ((! GetMultiNormal(from)->ClampVector[index]) && (fabs((*GetNodeVal(from->GetNode()))[index]) > cutoff))	{
			(*GetNodeVal(from->GetNode()))[index] = cutoff;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCutOff(link->Out(),cutoff,index);
		}
	}

	void GetLeafData(ContinuousData* data, int offset)	{
		RecursiveGetLeafData(GetRoot(),data,offset);
	}

	void RecursiveGetLeafData(const Link* from, ContinuousData* data, int offset)	{
		if (from->isLeaf())	{
			for (int i=0; i<data->GetNsite(); i++)	{
				if (data->Data[data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName())][i] != -1)	{
					data->Data[data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName())][i] = exp(GetVal(from,i+offset));
				}
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetLeafData(link->Out(),data,offset);
		}
	}

	void GetAncData(AncestralData* data, int offset)	{
		RecursiveGetAncData(GetRoot(),data,offset);
	}

	void RecursiveGetAncData(const Link* from, AncestralData* data, int offset)	{
		if (data->Exists(from))	{
			for (int i=0; i<data->GetNsite(); i++)	{
				int pos = i+offset;
				double tmp = GetMultiNormal(from)->val()[pos];
				if (data->GetState(from,i) != -1)	{
					data->SetState(from, i, exp(tmp));
				}
			}
		}
		else	{
			cerr << "did not find : " << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\n';
			exit(1);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetAncData(link->Out(),data,offset);
		}
	}

	double PiecewiseTranslationMove(double tuning, int index, int k){
		int n = 0;
		double tot = RecursivePiecewiseTranslationMove(this->GetRoot(),tuning,index,k,n);
		return tot / n;
	}

	double PiecewiseTranslation(double u, int index, int k){
		return RecursivePiecewiseTranslation(this->GetRoot(),u,index,k);
	}

	double RecursivePiecewiseTranslation(const Link* from, double u, int index, int k)	{
		double total = GetMultiNormal(from)->PiecewiseTranslation(u, index, k);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursivePiecewiseTranslation(link->Out(),u,index,k);
		}
		return total;
	}

	void RootNeighborPiecewiseTranslation(double u, int index, int k)	{
		const Link* from = GetRoot();
		GetMultiNormal(from)->PiecewiseTranslation(u, index, k);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			GetMultiNormal(link->Out())->PiecewiseTranslation(u, index, k);
		}
	}

	double Move(double tuning){
		int n = 0;
		double tot = RecursiveMove(this->GetRoot(),tuning,n);
		return tot / n;
	}

	double SegmentMove(double tuning, int imin, int size){
		int n = 0;
		double tot = RecursiveSegmentMove(this->GetRoot(),tuning,imin,size,n);
		return tot / n;
	}

	double GetLogProb()	{
		return RecursiveGetLogProb(this->GetRoot());
	}


	double GetMean(int index)	{
		int n = 0;
		double tmp = RecursiveGetMean(GetRoot(),index,n);
		return tmp / n;
	}

	double GetLeafMean(int index)	{
		int n = 0;
		double tmp = RecursiveGetLeafMean(GetRoot(),index,n);
		return tmp / n;
	}

	double GetMeanExp(int index)	{
		int n = 0;
		double tmp = RecursiveGetMeanExp(GetRoot(),index,n);
		return tmp / n;
	}

	// index : column in the continuous data
	// pos : entry to be clamped in the multivariate normal process
	void SetAndClamp(ContinuousData* data, int pos, int index = 0, int type = 0){
		int k = 0;
		int n = 0;
		RecursiveSetAndClamp(GetRoot(), data, pos, index, type, k, n);
		cerr << data->GetCharacterName(index) << " : " << n-k << " out of " << n << " missing\n";
	}

	void SetAndClamp(AncestralData* data, int offset = 0, int type = 0){
		RecursiveSetAndClamp(GetRoot(), data, offset,  type);
	}

	void SetLeafStates(Normal*** leafstates, ContinuousData* data, Var<PosReal>* variance, int pos, int index=0, int type = 0){
		RecursiveSetLeafStates(GetRoot(), leafstates, data, variance, pos, index, type);
	}

	void ClampRoot()	{
		GetNodeVal(GetRoot()->GetNode())->SetAtZero();
		GetNodeVal(GetRoot()->GetNode())->Clamp();
	}

	void ClampAll()	{
		RecursiveClampAll(GetRoot());
	}

	void RecursiveClampAll(const Link* from)	{
		GetNodeVal(from->GetNode())->Clamp();
		for (int i=0; i<GetDim(); i++)	{
			GetMultiNormal(from->GetNode())->Clamp(i);
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveClampAll(link->Out());
		}
	}

	string Construct(){
		return RecursiveConstruct(GetRoot());
	}

	void Reset()	{
		RecursiveReset(GetRoot());
	}

	void Add(MultiVariateTreeProcess* inprocess)	{
		RecursiveAdd(inprocess,GetRoot());
	}

	void ScalarMultiplication(double inscal)	{
		RecursiveScalarMultiplication(inscal,GetRoot());
	}

	void RecursiveRegister(DAGnode* node, const Link* from)	{
		GetMultiNormal(from)->Register(node);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(node,link->Out());
		}
	}

	void RecursiveDeregister(DAGnode* node, const Link* from)	{
		GetMultiNormal(from)->DeregisterFrom(node);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeregister(node,link->Out());
		}
	}

	void CheckIdentity()	{
		CheckIdentity(GetRoot());
	}

	void CheckIdentity(const Link* from)	{
		if (!GetMultiNormal(from)->CheckIdentity())	{
			if (from->Out()->isLeaf())	{
				cerr << "link out is leaf\n";
			}
			if (from->isRoot())	{
				cerr << "link is root\n";
			}
			exit(1);
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			CheckIdentity(link->Out());

		}
	}

	void Unclamp(int index)	{
		RecursiveUnclamp(GetRoot(),index);
	}

	void ML(int index)	{
		MLBackward(GetRoot(),index);
		MLForward(GetRoot(),0,index);
	}

	protected :

	void RecursivedrawSample(const Link* from)	{
		GetMultiNormal(from)->Sample();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivedrawSample(link->Out());
		}
	}

	void RecursiveSampleFixRoot(const Link* from)	{
		if (! from->isRoot())	{
			GetMultiNormal(from)->drawSampleUnclamped();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSampleFixRoot(link->Out());
		}
	}

	void RecursiveReset(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
		GetMultiNormal(from)->SetAtZero();
	}

	void RecursiveAdd(MultiVariateTreeProcess* inprocess, const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(inprocess,link->Out());
		}
		GetMultiNormal(from)->Add(*(inprocess->GetMultiNormal(from)));
	}

	void RecursiveScalarMultiplication(double inscal, const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveScalarMultiplication(inscal,link->Out());
		}
		GetMultiNormal(from)->ScalarMultiplication(inscal);
	}

	double RecursiveMove(const Link* from, double tuning, int& count)	{
		double total = GetMultiNormal(from)->Move(tuning);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveMove(link->Out(),tuning,count);
		}
		return total;
	}

	double RecursiveSegmentMove(const Link* from, double tuning, int imin, int size, int& count)	{
		double total = GetMultiNormal(from)->SegmentMove(tuning,imin,size);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveSegmentMove(link->Out(),tuning,imin,size,count);
		}
		return total;
	}

	double RecursivePiecewiseTranslationMove(const Link* from, double tuning, int index, int k, int& count)	{
		double total = GetMultiNormal(from)->PiecewiseTranslationMove(tuning, index, k);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursivePiecewiseTranslationMove(link->Out(),tuning,index,k,count);
		}
		return total;
	}

	double RecursiveGetLogProb(const Link* from)	{
		double total = GetMultiNormal(from)->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLogProb(link->Out());
		}
		return total;
	}

	double RecursiveGetMean(const Link* from, int index, int& tot)	{
		double total = (*GetMultiNormal(from))[index];
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMean(link->Out(), index, tot);
		}
		return total;
	}

	double RecursiveGetLeafMean(const Link* from, int index, int& tot)	{
		double total = 0;
		if (from->isLeaf())	{
			total += (*GetMultiNormal(from))[index];
			tot++;
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLeafMean(link->Out(), index, tot);
		}
		return total;
	}

	double RecursiveGetMeanExp(const Link* from, int index, int& tot)	{
		double total = exp((*GetMultiNormal(from))[index]);
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMeanExp(link->Out(), index, tot);
		}
		return total;
	}


	void RecursiveUnclamp(const Link* from, int index){
		if (from->isLeaf())	{
			GetMultiNormal(from)->UnClamp(index);
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveUnclamp(link->Out(), index);
		}
	}

	void RecursiveSetLeafStates(const Link* from, Normal*** leafstates, ContinuousData* data, Var<PosReal>* variance, int pos, int index, int type)	{

		if (from->isLeaf())	{
			int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			if (tax != -1)	{
				double tmp = data->GetState(tax, index);
				if (tmp != -1)	{
					MultiNormal* m = GetMultiNormal(from);
					leafstates[tax][index] = new Normal(m,variance,pos);
					double tmp = data->GetState(tax, index);
					if (type == 1)	{
						if ((tmp <= 0) || (tmp >= 1))	{
							cerr << "error in recursive clamp at: cont data not between 0 and 1\n";
							exit(1);
						}
						leafstates[tax][index]->ClampAt(log(tmp / (1 - tmp)));
						m->val()[pos] = log(tmp/(1-tmp)) + 0.01 * Random::sNormal();
					}
					else if (type == 2)	{
						leafstates[tax][index]->ClampAt(tmp);
						m->val()[pos] = tmp + 0.01 * Random::sNormal();
					}
					else	{
						if (tmp <= 0)	{
							cerr << "error in recursive clamp at: negative cont data\n";
							exit(1);
						}
						leafstates[tax][index]->ClampAt(log(tmp));
						m->val()[pos] = log(tmp) + 0.01 * Random::sNormal();
					}
				}
				else	{
					leafstates[tax][index] = 0;
				}
			}
			else	{
				cerr << "set leaf states: " << from->GetNode()->GetName() << " not found\n";
			}

		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetLeafStates(link->Out(), leafstates, data, variance, pos, index, type);
		}
	}

	void RecursiveSetAndClamp(const Link* from, ContinuousData* data, int pos, int index, int type, int& k, int & n){
		if(from->isLeaf()){
			n++;
			int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			// int tax = data->GetTaxonSet()->GetTaxonIndexWithIncompleteName(from->GetNode()->GetName());
			// int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
			// cerr << "set and clamp : " << from->GetNode()->GetName() << '\t' << tax << '\t' << pos << '\t' << index << '\t';
			if (tax != -1)	{
				double tmp = data->GetState(tax, index);
				// cerr << tmp << '\t';
				if (tmp != -1)	{
					k++;
					MultiNormal* m = GetMultiNormal(from);
					if (type == 2)	{
						m->val()[pos] = tmp;
					}
					else if (type == 1)	{ // logit
						if ((tmp <= 0) || (tmp >= 1))	{
							cerr << "error in recursive clamp at: cont data not between 0 and 1\n";
							exit(1);
						}
						m->val()[pos] = log(tmp / (1 - tmp));
					}
					else	{
						if (tmp <= 0)	{
							cerr << "error in recursive clamp at: negative cont data\n";
							exit(1);
						}
						m->val()[pos] = log(tmp);
					}
					m->Clamp(pos);
					// cerr << (*m)[pos];
				}
			}
			else	{
				cerr << "set and clamp : " << from->GetNode()->GetName() << " not found\n";
			}
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetAndClamp(link->Out(), data, pos, index, type,k,n);
		}
	}

	void RecursiveSetAndClamp(const Link* from, AncestralData* data, int offset, int type){
		if (data->Exists(from))	{
			for (int i=0; i<data->GetNsite(); i++)	{
				int pos = i+offset;
				double tmp = data->GetState(from, i);
				MultiNormal* m = GetMultiNormal(from);
				if (type == 2)	{
					m->val()[pos] = tmp;
				}
				else if (type == 1)	{ // logit
					if (tmp <= 0)	{
						cerr << "error in recursive clamp at: negative cont data\n";
						exit(1);
					}
					m->val()[pos] = log(tmp / (1 - tmp));
				}
				else	{
					if (tmp <= 0)	{
						cerr << "error in recursive clamp at: negative cont data\n";
						exit(1);
					}
					m->val()[pos] = log(tmp);
				}
				m->Clamp(pos);
			}
		}
		else	{
			cerr << "did not find : " << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\n';
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetAndClamp(link->Out(), data, offset, type);
		}
	}

	/*
	void BackwardSetEmpirical(const Link* from)	{
		MultiNormal* m = GetMultiNormal(from);
		if (from->isLeaf())	{
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{

		}
	}
	*/

	void MLBackward(const Link* from, int index)	{

		double nu = 1.0 / GetBranchLength(from->GetBranch());

		if (from->isLeaf())	{
			if (GetMultiNormal(from)->ClampVector[index])	{
				MLmu[from] = (*GetMultiNormal(from))[index];
				MLk[from] = 0;
				MLm[from] = nu;
			}
			else	{
				MLmu[from] = 0;
				MLk[from] = 0;
				MLm[from] = 0;
			}
		}
		else	{
			double k0 = 0;
			double mu0 = 0;
			for(const Link* link=from->Next(); link!=from; link=link->Next())	{
				MLBackward(link->Out(), index);
				k0 += MLm[link->Out()];
				mu0 += MLm[link->Out()] * MLmu[link->Out()];
			}
			mu0 /= k0;
			MLmu[from] = mu0;
			MLk[from] = k0;

			if (! from->isRoot())	{
				double m = nu * k0 / (nu + k0);
				MLm[from] = m;
			}
			else	{
				MLm[from] = 0;
			}
		}
		// cerr <<  MLmu[from] << '\t' << MLk[from] << '\t' << MLm[from] << '\t' << from->isLeaf() << '\t' << from->isRoot() << '\t' << GetMultiNormal(from)->ClampVector[index] << '\t' << (*GetMultiNormal(from))[index] << '\n';
	}

	void MLForward(const Link* from, const Link* parent, int index)	{
		if (from->isRoot())	{
			if (! GetMultiNormal(from)->ClampVector[index])	{
				(*GetMultiNormal(from))[index] = MLmu[from];
				GetMultiNormal(from)->Clamp(index);
			}
		}
		else	{
			double x0 = (*GetMultiNormal(parent))[index];
			double k = MLk[from];
			double mu = MLmu[from];
			double nu = 1.0 / GetBranchLength(from->GetBranch());
			double alpha = (k*mu + nu*x0) / (k + nu);
			if (! GetMultiNormal(from)->ClampVector[index])	{
				(*GetMultiNormal(from))[index] = alpha;
				GetMultiNormal(from)->Clamp(index);
			}
		}
		// cerr << (*GetMultiNormal(from))[index] << '\n';
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			MLForward(link->Out(), from, index);
		}
	}

	map<const Link*, double> MLmu;
	map<const Link*, double> MLk;
	map<const Link*, double> MLm;


	// useless ??
	string RecursiveConstruct(Link* from){
		ostringstream retour;
		if(!from->isLeaf()){
			retour << "(";
			for(Link* link=from->Next(); link!=from; link=link->Next()){
				retour << RecursiveConstruct(link->Out());
				if(link->Next() != from){
					retour << ",";
				}
			}
		}
		retour << ")" << from->GetNode()->GetName();
		for(int i=0; i<GetDim(); i++){
			retour << "_" << exp(GetMultiNormal(from)->val()[i]);
		}
		if(from->isRoot()){
			retour << ";";
		}
		return retour.str();
	}




	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			if (scalefunction)	{
				return new MultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), chrono->GetNodeVal(link->GetNode()), agescale, scalefunction);
			}
			return new MultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), scaletree ? scaletree->GetBranchVal(link->GetBranch()): 0, drift);
		}
		else{
 			// return new MultiNormal(sigma);
 			return new MultiNormal(0,sigma,rootmean,rootvar);
		}
	}

};

class MultiVariateSegmentMove : public MCUpdate {

	MultiVariateTreeProcess* tree;
	double tuning;
	int imin;
	int size;

	public:

	MultiVariateSegmentMove(MultiVariateTreeProcess* intree, double intuning, int inimin, int insize){
		tree = intree;
		tuning = intuning;
		imin = inimin;
		size = insize;
	}

	double Move(double tuning_modulator){
		return tree->SegmentMove(tuning,imin,size);
	}
};

class MultiVariateWholeTreePiecewiseTranslationMove : public MCUpdate, public Mnode {

	MultiVariateTreeProcess* tree;
	double tuning;
	int index;
	int k;

	public:

	MultiVariateWholeTreePiecewiseTranslationMove(MultiVariateTreeProcess* intree, double intuning, int inindex, int ink){
		tree = intree;
		tuning = intuning;
		index = inindex;
		k = ink;
		tree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		tree->PiecewiseTranslation(u,index,k);
		double logratio = Update();
		// cerr << u << '\t' << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

#endif

