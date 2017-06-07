#ifndef NODELINEARCOMBINATIONTREE2_H
#define NODELINEARCOMBINATIONTREE2_H

#include "ValTree.h"


#include <cmath>
#include <sstream>


class NeutralOmegaLinearCombination : public Dvar<Real> {
	
	public :
	
	
	NeutralOmegaLinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Const<PosReal>* inK, double* inneutralomegaslope, bool insameseq) {
		x = inx;
		gamma = ingamma;
		K = inK;
		neutralomegaslope = inneutralomegaslope;
		sameseq = insameseq;
		Register(x);
		Register(gamma);
		specialUpdate();
	}
	
	NeutralOmegaLinearCombination(Var<RealVector>* inx, Var<Real>* inbeta, Var<Real>* inbeta2, double* inneutralomegaslope, bool insameseq) {
		x = inx;
		beta = inbeta;
		beta2 = inbeta2;
		neutralomegaslope = inneutralomegaslope;
		sameseq = insameseq;
		Register(x);
		Register(beta);
		Register(beta2);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		if (sameseq) {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				a+= (*x)[i] * neutralomegaslope[i] * ( 1 / ( 1 + *gamma * *K ) );
			}
			setval(a);	
		}
		if (!sameseq) {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				a+= (*x)[i] * neutralomegaslope[i];
				//cerr << (*x)[i];
				//cerr << '\t';
				//cerr << a;
				//cerr << '\n';
			}
			a+= *beta - *beta2;
			//cerr << a;
			//cerr << '\n';
			setval(a);		
		}
	}		
	
	private :
	
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Const<PosReal>* K;
	Var<Real>* beta;
	Var<Real>* beta2;
	double* neutralomegaslope;
	bool sameseq;
	
};		



class OmegaLinearCombination : public Dvar<Real> {
	
	public :
	
	OmegaLinearCombination(Var<RealVector>* inx, Var<Real>* ino, double* inomegaslope) {
		x = inx;
		o = ino;
		omegaslope = inomegaslope;
		Register(x);
		Register(o);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= exp((*x)[i]) * omegaslope[i];
			//cerr << (*x)[i];
			//cerr << '\n';
		}
		a+= exp(*o);
		//cerr << a;
		//cerr << '\n';
		setval(log(a));	
	}	
	
	private :
	
	
	Var<RealVector>* x;
	Var<Real>* o;
	double* omegaslope;
	
};	




class ULinearCombination : public Dvar<Real> {
	
	public :
	
	ULinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta, Const<PosReal>* inK, double* inuslope, bool insameseq) {
		x = inx;
		gamma = ingamma;
		beta = inbeta,
		K = inK;
		uslope = inuslope;
		sameseq = insameseq;
		Register(beta);
		Register(gamma);
		Register(x);
		specialUpdate();
	}
	
	ULinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta2, double* inuslope, bool insameseq) {
		x = inx;
		gamma = ingamma;
		beta2 = inbeta2,
		uslope = inuslope;
		sameseq = insameseq;
		Register(beta2);
		Register(gamma);
		Register(x);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		if (sameseq) {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				if (uslope[i] == 1000) {
					a+= ( 1 / *gamma ) * (*x)[i];
				}
				else {
					a+= (*x)[i] * uslope[i];
				}
			}	
			a += ( 1 / *gamma ) * ( *beta - log(*K) - *gamma * log(4) );	
			setval(a);	
		}
		if (!sameseq) {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				if (uslope[i] == 1000) {
					a+= (*x)[i] / *gamma;
				}
				else {
					a+= (*x)[i] * uslope[i];
				}
				//cerr << (*x)[i];
				//cerr << "\t";
				//cerr << a;
				//cerr << "\n";
			}	
			a -= ( *beta2 / *gamma ) + log(4);	
			setval(a);	
			//cerr << a;
			//cerr << "\n";
		}	
	}
	
	private :
	
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Var<Real>* beta;
	Var<Real>* beta2;
	Const<PosReal>* K;
	double* uslope;
	bool sameseq;
	
};


class SynrateLinearCombination : public Dvar<Real> {
	
	public :
	
	SynrateLinearCombination(Var<RealVector>* inx, Var<Real>* ino, Var<PosReal>* inrootage, double* insynrateslope) {
		x = inx;
		o = ino;
		rootage = inrootage;
		synrateslope = insynrateslope;
		Register(x);
		Register(o);
		Register(rootage);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= (*x)[i] * synrateslope[i];
			//cerr << (*x)[i];
			//cerr << '\t';
			//cerr << a;
			//cerr << '\n'; 
		}
		a+= *o;
		//cerr << a;
		//cerr << "\n";
		a += log(rootage->val() * 365 * pow(10, 6));
		//cerr << rootage->val();
		//cerr << "\n";
		//cerr << a;
		//cerr << '\n';
		setval(a);	
	}
		
					
	
	private :			
	
	Var<RealVector>* x;
	Var<Real>* o;
	double* synrateslope;
	Var<PosReal>* rootage;
		
};		
		
		
		
	
class NeLinearCombination : public Dvar<Real> {
	
	public :
	
	NeLinearCombination(Var<RealVector>* inx, Var<Real>* ino, double* inNeslope) {
		x = inx;
		o = ino;
		Neslope = inNeslope;
		Register(x);
		Register(o);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= (*x)[i] * Neslope[i];
		}
		a -= *o;
		a -= log(4);
		setval(a);	
	}
		
					
	
	private :			
	
	Var<RealVector>* x;
	Var<Real>* o;
	double* Neslope;
		
};		



class NeutralOmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	NeutralOmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Const<PosReal>* inK, double* inneutralomegaslope, bool insameseq) {
		process	= inprocess;
		gamma = ingamma;
		K = inK;
		neutralomegaslope = inneutralomegaslope;
		sameseq = insameseq;
		RecursiveCreate(GetRoot());
	}
	
	NeutralOmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* inbeta, Var<Real>* inbeta2, double* inneutralomegaslope, bool insameseq) {
		process	= inprocess;
		beta = inbeta;
		beta2 = inbeta2;
		neutralomegaslope = inneutralomegaslope;
		sameseq = insameseq;
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
		if ((! from->isRoot()))	{
			GetNodeVal(from->GetNode())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<Real>* CreateNodeVal(const Link* link){
		if (sameseq) {
			return new NeutralOmegaLinearCombination(process->GetNodeVal(link->GetNode()), gamma, K, neutralomegaslope, sameseq);
		}
		if (!sameseq) {
			return new NeutralOmegaLinearCombination(process->GetNodeVal(link->GetNode()), beta, beta2, neutralomegaslope, sameseq);
		}
	}
		
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Const<PosReal>* K;
	Var<Real>* beta;
	Var<Real>* beta2;
	double* neutralomegaslope;
	bool sameseq;
	
};	





class OmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	OmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeneutralomegatree, double* inomegaslope) {
		process	= inprocess;
		nodeneutralomegatree = innodeneutralomegatree;
		omegaslope = inomegaslope;
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
		if ((! from->isRoot()))	{
			GetNodeVal(from->GetNode())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		return new OmegaLinearCombination(process->GetNodeVal(link->GetNode()), nodeneutralomegatree->GetNodeVal(link->GetNode()), omegaslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeneutralomegatree;
	double* omegaslope;
	
};				
			

		
class ULinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	ULinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Var<Real>* inbeta, Const<PosReal>* inK, double* inuslope, bool insameseq) {
		process	= inprocess;
		gamma = ingamma;
		beta = inbeta;
		K = inK;
		uslope = inuslope;
		sameseq = insameseq;
		RecursiveCreate(GetRoot());
	}
		
	ULinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Var<Real>* inbeta2, double* inuslope, bool insameseq) {
		process	= inprocess;
		gamma = ingamma;
		beta2 = inbeta2;
		uslope = inuslope;
		sameseq = insameseq;
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
		if ((! from->isRoot()))	{
			GetNodeVal(from->GetNode())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		if (sameseq) {
			return new ULinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta, K, uslope, sameseq);
		}
		if (!sameseq) {
			return new ULinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta2, uslope, sameseq);
		}	
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Var<Real>* beta;
	Var<Real>* beta2;
	Const<PosReal>* K;
	double* uslope;
	bool sameseq;
	
};						
	
	
	
class SynrateLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	SynrateLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeutree, Var<PosReal>* inrootage, double* insynrateslope) {
		process	= inprocess;
		nodeutree = innodeutree;
		rootage = inrootage;
		synrateslope = insynrateslope;
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
		if ((! from->isRoot()))	{
			GetNodeVal(from->GetNode())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		return new SynrateLinearCombination(process->GetNodeVal(link->GetNode()), nodeutree->GetNodeVal(link->GetNode()), rootage, synrateslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeutree;
	Var<PosReal>* rootage;
	double* synrateslope;
	
};
		
		
	
class NeLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	NeLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeutree, double* inNeslope) {
		process	= inprocess;
		nodeutree = innodeutree;
		Neslope = inNeslope;
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
		if ((! from->isRoot()))	{
			GetNodeVal(from->GetNode())->specialUpdate();
		}
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}
	
	
	Dvar<Real>* CreateNodeVal(const Link* link){
		return new NeLinearCombination(process->GetNodeVal(link->GetNode()), nodeutree->GetNodeVal(link->GetNode()), Neslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeutree;
	double* Neslope;
	
};
		
			
			
			
#endif	
	
	
	
	
	
