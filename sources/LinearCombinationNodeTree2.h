#ifndef NODELINEARCOMBINATIONTREE2_H
#define NODELINEARCOMBINATIONTREE2_H

#include "ValTree.h"


#include <cmath>
#include <sstream>

class NeutralOmegaLinearCombination : public Dvar<Real> {
	
	public :
	
	NeutralOmegaLinearCombination(Var<RealVector>* inx, Var<Real>* inlogkappa1, Var<Real>* inlogkappa2, int inidxpiNpiS) {
		x = inx;
        logkappa1 = inlogkappa1;
        logkappa2 = inlogkappa2;
        idxpiNpiS = inidxpiNpiS;
		Register(x);
		Register(logkappa1);
		Register(logkappa2);
		specialUpdate();
	}
			
    // since:
    // omega_0 = kappa_1 N_e^{-beta}
    // piN/piS = kappa_2 N_e^{-beta}
    //
    // -> log omega_0 = log piN/piS + log kappa_1 - log kappa_2
	void specialUpdate() {
        setval((*x)[idxpiNpiS] + logkappa1->val() - logkappa2->val());
	}		
	
	private :
	
	Var<RealVector>* x;
	Var<Real>* logkappa1;
	Var<Real>* logkappa2;
    int idxpiNpiS;
};		


class OmegaLinearCombination : public Dvar<Real> {
	
	public :
	
    // process x at index idx, gives log omega_a
    // log omega_0 given as a separate real variable
	OmegaLinearCombination(Var<RealVector>* inx, Var<Real>* inlogom0, int inidx)    {
		x = inx;
		logom0 = inlogom0;
        idx = inidx;
		Register(x);
		Register(logom0);
		specialUpdate();
	}
			
    // omega = omega_0 + omega_a
	void specialUpdate() {
        setval(exp(logom0->val()) + exp((*x)[idx]));
	}	
	
	private :
	
	Var<RealVector>* x;
	Var<Real>* logom0;
    int idx;
};	

class ULinearCombination : public Dvar<Real> {
	
	public :
	
	ULinearCombination(Var<RealVector>* inx, Var<Real>* inbeta, Var<Real>* inlogkappa2, int inidxpiS, int inidxpiNpiS) {
		x = inx;
        beta = inbeta;
        logkappa2 = inlogkappa2;
        idxpiS = inidxpiS;
        idxpiNpiS = inidxpiNpiS;
		Register(x);
		Register(beta);
        Register(logkappa2);
		specialUpdate();
	}
			
    // using : log Ne -1/beta * (log piN/piS - log kappa_2)
    // log u = log piS - log N_e - log 4.0
    //       = log piS + 1/beta * log piN/piS - 1/beta * log kappa_2 - log 4.0
	void specialUpdate() {
        setval((*x)[idxpiS] + 1/beta->val()*(*x)[idxpiNpiS] - logkappa2->val()/beta->val() - log(4.0));
	}
	
	private :
	
	Var<RealVector>* x;
	Var<Real>* beta;
	Var<Real>* logkappa2;
    int idxpiS;
    int idxpiNpiS;
};


class NeLinearCombination : public Dvar<Real> {
	
	public :
	
	NeLinearCombination(Var<RealVector>* inx, Var<Real>* inlogu, int inidxpiS) {
		x = inx;
		logu = inlogu;
        idxpiS = inidxpiS;
		Register(x);
		Register(logu);
		specialUpdate();
	}
			
		
	void specialUpdate() {
        setval((*x)[idxpiS] - logu->val() - log(4.0));
	}
	
	private :			
	
	Var<RealVector>* x;
	Var<Real>* logu;
	int idxpiS;
		
};		

class SynrateLinearCombination : public Dvar<Real> {
	
	public :
	
	SynrateLinearCombination(Var<RealVector>* inx, Var<Real>* inlogu, Var<PosReal>* inrootage, int inidxgentime)    {
		x = inx;
        logu = inlogu;
		rootage = inrootage;
        idxgentime = inidxgentime;
		Register(x);
		Register(logu);
		Register(rootage);
		specialUpdate();
	}
		
    // dS is per tree depth, which is itself given in myr
    // log dS = log u - log tau + log(rootage*365.10^6)
    // generation time is given in days
    // 365.10^6: # days per myr
	void specialUpdate() {
        setval(logu->val() - (*x)[idxgentime] + log(rootage->val()) + log(365.0) + 6*log(10.0));
	}
	
	private :			
	
	Var<RealVector>* x;
	Var<Real>* logu;
	Var<PosReal>* rootage;
    int idxgentime;
};		


class NeutralOmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	NeutralOmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* inlogkappa1, Var<Real>* inlogkappa2, int inidxpiNpiS)  {
		process	= inprocess;
        logkappa1 = inlogkappa1;
        logkappa2 = inlogkappa2;
        idxpiNpiS = inidxpiNpiS;
		RecursiveCreate(GetRoot());
	}
	
	~NeutralOmegaLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :
	
	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<Real>* CreateNodeVal(const Link* link){
        return new NeutralOmegaLinearCombination(process->GetNodeVal(link->GetNode()), logkappa1, logkappa2, idxpiNpiS);
	}
		
	
	NodeVarTree<RealVector>* process;
	Var<Real>* logkappa1;
    Var<Real>* logkappa2;
    int idxpiNpiS;
};	

class OmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	OmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeneutralomegatree, int inidx)  {
		process	= inprocess;
		nodeneutralomegatree = innodeneutralomegatree;
        idx = inidx;
		RecursiveCreate(GetRoot());
	}
	
	~OmegaLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}

	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
		
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		return new OmegaLinearCombination(process->GetNodeVal(link->GetNode()), nodeneutralomegatree->GetNodeVal(link->GetNode()), idx);
	}
	
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeneutralomegatree;
    int idx;
};				

class ULinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	ULinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* inbeta, Var<Real>* inlogkappa2, int inidxpiS, int inidxpiNpiS)    {
		process	= inprocess;
        beta = inbeta;
        logkappa2 = inlogkappa2;
        idxpiS = inidxpiS;
        idxpiNpiS = inidxpiNpiS;
		RecursiveCreate(GetRoot());
	}
		
	~ULinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	Dvar<Real>* CreateNodeVal(const Link* link){
        return new ULinearCombination(process->GetNodeVal(link->GetNode()), beta, logkappa2, idxpiS, idxpiNpiS);
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<Real>* beta;
	Var<Real>* logkappa2;
    int idxpiS;
    int idxpiNpiS;
};						

class NeLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	NeLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeutree, int inidxpiS)   {
		process	= inprocess;
		nodeutree = innodeutree;
        idxpiS = inidxpiS;
		RecursiveCreate(GetRoot());
	}
	
	~NeLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		return new NeLinearCombination(process->GetNodeVal(link->GetNode()), nodeutree->GetNodeVal(link->GetNode()), idxpiS);
	}
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeutree;
    int idxpiS;
};

class SynrateLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	SynrateLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeutree, Var<PosReal>* inrootage, int inidxgentime) {
		process	= inprocess;
		nodeutree = innodeutree;
		rootage = inrootage;
        idxgentime = inidxgentime;
		RecursiveCreate(GetRoot());
	}
	
	~SynrateLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		return new SynrateLinearCombination(process->GetNodeVal(link->GetNode()), nodeutree->GetNodeVal(link->GetNode()), rootage, idxgentime);
	}
	
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeutree;
	Var<PosReal>* rootage;
    int idxgentime;
};

#endif	


