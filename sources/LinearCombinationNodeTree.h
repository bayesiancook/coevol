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
		a += -1 * log(rootage->val() * 365);
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
		double* p;
		*p = SynrateLinearCombination(process->GetNodeVal(link->GetNode()), rootage, synrateslope);
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
		double* p;
		*p = OmegaLinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta, omegaslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Var<Real>* beta;
	double* omegaslope;
	
};		
			
#endif	
	
	
	
	
	
