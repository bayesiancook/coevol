#ifndef MIXMULTIVARIATETREEPROCESS
#define MIXMULTIVARIATETREEPROCESS

#include "InverseWishartMatrix.h"
#include "ValTree.h"
#include "ContinuousData.h"
#include "AncestralData.h"



class MixMultiVariateTreeProcess : public MCMC , public NodeValPtrTree<Rvar<RealVector> > {

	protected:

	VarArray<CovMatrix>* sigma;
	LengthTree* tree;
	BranchPartition* partition;

	public:


	MixMultiVariateTreeProcess()	{}


	MixMultiVariateTreeProcess(VarArray<CovMatrix>* insigma, LengthTree* intree, BranchPartition* inpartition){
		sigma = insigma;
		tree = intree;
		scaletree = inscaletree;
		RecursiveCreate(GetRoot());
	}

	~MixMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	Tree* GetTree(){
		return tree->GetTree();
	}


	LengthTree* GetLengthTree(){
		return tree;
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

	int GetDim(){
		return sigma->GetDim();
	}


	void drawSample()	{
		RecursivedrawSample(this->GetRoot());
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

	double Move(double tuning){
		int n = 0;
		double tot = RecursiveMove(this->GetRoot(),tuning,n);
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
		cerr << "character : " << index << '\t' << k << '\t' << n << '\n';
	}

	void SetAndClamp(AncestralData* data, int offset = 0, int type = 0){
		RecursiveSetAndClamp(GetRoot(), data, offset,  type);
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

	void Add(MixMultiVariateTreeProcess* inprocess)	{
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


	protected :

	void RecursivedrawSample(const Link* from)	{
		GetMultiNormal(from)->Sample();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivedrawSample(link->Out());
		}
	}

	void RecursiveReset(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
		GetMultiNormal(from)->SetAtZero();
	}

	void RecursiveAdd(MixMultiVariateTreeProcess* inprocess, const Link* from)	{
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
			return new MultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), scaletree ? scaletree->GetBranchVal(link->GetBranch()): 0);
		}
		else{
 			return new MultiNormal(sigma);
		}
	}

};


class MixMultiVariateWholeTreePiecewiseTranslationMove : public MCUpdate, public Mnode {

	MixMultiVariateTreeProcess* tree;
	double tuning;
	int index;
	int k;

	public:

	MixMultiVariateWholeTreePiecewiseTranslationMove(MixMultiVariateTreeProcess* intree, double intuning, int inindex, int ink){
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

