#ifndef NODELINEARCOMBINATIONTREE_H
#define NODELINEARCOMBINATIONTREE_H

#include "ValTree.h"


#include <cmath>
#include <sstream>

class SynrateLinearCombination : public Dvar<Real> {
	
	public :
	
	SynrateLinearCombination(Var<RealVector>* inx, Var<PosReal>* inrootage, double* insynrateslope) {
		x = inx;
		rootage = inrootage;
		synrateslope = insynrateslope;
		Register(x);
		Register(rootage);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= (*x)[i] * synrateslope[i];
		}
		a += -1 * log(rootage->val() * 365 * pow(10, 6));
		setval(a);	
	}	
	
	private :			
	
	Var<RealVector>* x;
	double* synrateslope;
	Var<PosReal>* rootage;
		
};		
		
	

class OmegaLinearCombination : public Dvar<Real> {
	
	public :
	
	OmegaLinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta, double* inomegaslope) {
		x = inx;
		gamma = ingamma;
		beta = inbeta;
		omegaslope = inomegaslope;
		Register(x);
		Register(gamma);
		Register(beta);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= (*x)[i] * omegaslope[i] * *gamma;
		}
		a += *beta + *gamma * log(4);	
		setval(a);	
	}	
	
	private :
	
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Var<Real>* beta;
	double* omegaslope;
	
};		

class NeLinearCombination : public Dvar<Real> {
	
	public :
	
	NeLinearCombination(Var<RealVector>* inx, double* inNeslope) {
		x = inx;
		Neslope = inNeslope;
		Register(x);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= (*x)[i] * Neslope[i];
		}
		a += log(4);	
		setval(a);	
	}	
	
	private :
	
	
	Var<RealVector>* x;
	double* Neslope;
	
};


class Adaptative_omegaLinearCombination : public Dvar<Real> {
	
	public :
	
	Adaptative_omegaLinearCombination(Var<RealVector>* inx, Var<Real>* ino, double* inadaptative_omegaslope) {
		x = inx;
		o = ino;
		adaptative_omegaslope = inadaptative_omegaslope;
		Register(x);
		Register(o);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= exp((*x)[i]) * adaptative_omegaslope[i];
		}
		a+= exp(*o);
		setval(log(a));	
	}	
	
	private :
	
	
	Var<RealVector>* x;
	Var<Real>* o;
	double* adaptative_omegaslope;
	
};			
			
	
class SynrateLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	SynrateLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<PosReal>* inrootage, double* insynrateslope) {
		process	= inprocess;
		rootage = inrootage;
		synrateslope = insynrateslope;
		RecursiveCreate(GetRoot());
	}
	
	~SynrateLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	Dvar<Real>* CreateNodeVal(const Link* link){
		return new SynrateLinearCombination(process->GetNodeVal(link->GetNode()), rootage, synrateslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<PosReal>* rootage;
	double* synrateslope;
	
};		
		
	
class OmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	OmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Var<Real>* inbeta, double* inomegaslope) {
		process	= inprocess;
		gamma = ingamma;
		beta = inbeta;
		omegaslope = inomegaslope;
		RecursiveCreate(GetRoot());
	}
	
	~OmegaLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	Dvar<Real>* CreateNodeVal(const Link* link){
		return new OmegaLinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta, omegaslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Var<Real>* beta;
	double* omegaslope;
	
};	
	
		
class NeLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	NeLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, double* inNeslope) {
		process	= inprocess;
		Neslope = inNeslope;
		RecursiveCreate(GetRoot());
	}
	
	~NeLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	Dvar<Real>* CreateNodeVal(const Link* link){
		return new NeLinearCombination(process->GetNodeVal(link->GetNode()), Neslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	double* Neslope;
	
};			


class Adaptative_omegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	Adaptative_omegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, NodeVarTree<Real>* innodeomegatree, double* inadaptative_omegaslope) {
		process	= inprocess;
		nodeomegatree = innodeomegatree;
		adaptative_omegaslope = inadaptative_omegaslope;
		RecursiveCreate(GetRoot());
	}
	
	~Adaptative_omegaLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	Dvar<Real>* CreateNodeVal(const Link* link){
		return new Adaptative_omegaLinearCombination(process->GetNodeVal(link->GetNode()), nodeomegatree->GetNodeVal(link->GetNode()), adaptative_omegaslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	NodeVarTree<Real>* nodeomegatree;
	double* adaptative_omegaslope;
	
};				
			
#endif	
	
	
	
	
	
