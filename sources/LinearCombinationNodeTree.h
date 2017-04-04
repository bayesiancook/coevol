#ifndef NODELINEARCOMBINATIONTREE_H
#define NODELINEARCOMBINATIONTREE_H

#include "ValTree.h"


#include <cmath>
#include <sstream>

class SynrateLinearCombination : public Dvar<Real> {
	
	public :
	
	SynrateLinearCombination(Var<RealVector>* inx, Var<PosReal>* inrootage, double* insynrateslope, bool inwithNe) {
		x = inx;
		rootage = inrootage;
		synrateslope = insynrateslope;
		withNe = inwithNe;
		Register(x);
		Register(rootage);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		if (!withNe) {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				a+= (*x)[i] * synrateslope[i];
			}
			a += log(rootage->val() * 365 * pow(10, 6));
			setval(a);	
		}
		else {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				a+= (*x)[i] * synrateslope[i];
			}
			a += log(rootage->val() * 365 * pow(10, 6) / 4);
			setval(a);	
		}
	}	
			
			
	
	private :			
	
	Var<RealVector>* x;
	double* synrateslope;
	Var<PosReal>* rootage;
	bool withNe;
		
};		
		
	

class OmegaLinearCombination : public Dvar<Real> {
	
	public :
	
	OmegaLinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta, double* inomegaslope, bool inwithNe) {
		x = inx;
		gamma = ingamma;
		beta = inbeta;
		omegaslope = inomegaslope;
		withNe = inwithNe;
		Register(x);
		Register(gamma);
		Register(beta);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		if (!withNe) {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				a+= (*x)[i] * omegaslope[i] * *gamma;
			}
			a += *beta + *gamma * log(4);	
			setval(a);	
		}
		else {
			double a(0);
			for (int i=0; i<x->GetDim(); i++) {
				a+= (*x)[i] * omegaslope[i] * *gamma;
			}
			a += *beta;	
			setval(a);	
		}
	}		
	
	private :
	
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Var<Real>* beta;
	double* omegaslope;
	bool withNe;
	
};		

class U_NeLinearCombination : public Dvar<Real> {
	
	public :
	
	U_NeLinearCombination(Var<RealVector>* inx, double* inu_Neslope) {
		x = inx;
		u_Neslope = inu_Neslope;
		Register(x);
		specialUpdate();
	}
			
		
	void specialUpdate() {
		double a(0);
		for (int i=0; i<x->GetDim(); i++) {
			a+= (*x)[i] * u_Neslope[i];
		}
		a -= log(4);	
		setval(a);	
	}	
	
	private :
	
	
	Var<RealVector>* x;
	double* u_Neslope;
	
};


			
	
class SynrateLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	SynrateLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<PosReal>* inrootage, double* insynrateslope, bool inwithNe) {
		process	= inprocess;
		rootage = inrootage;
		synrateslope = insynrateslope;
		withNe = inwithNe;
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
		return new SynrateLinearCombination(process->GetNodeVal(link->GetNode()), rootage, synrateslope, withNe);
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<PosReal>* rootage;
	double* synrateslope;
	bool withNe;
	
};
		
		
	
class OmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	OmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Var<Real>* inbeta, double* inomegaslope, bool inwithNe) {
		process	= inprocess;
		gamma = ingamma;
		beta = inbeta;
		omegaslope = inomegaslope;
		withNe = inwithNe;
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
		return new OmegaLinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta, omegaslope, withNe);
	}
	
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Var<Real>* beta;
	double* omegaslope;
	bool withNe;
	
};	
	
		
class U_NeLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	U_NeLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, double* inu_Neslope) {
		process	= inprocess;
		u_Neslope = inu_Neslope;
		RecursiveCreate(GetRoot());
	}
	
	~U_NeLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :

	Dvar<Real>* CreateNodeVal(const Link* link){
		return new U_NeLinearCombination(process->GetNodeVal(link->GetNode()), u_Neslope);
	}
	
	
	NodeVarTree<RealVector>* process;
	double* u_Neslope;
	
};			


#endif	
	
	
	
	
	
