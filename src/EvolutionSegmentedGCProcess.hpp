/*
 * File:   EvolutionSegmentedGCProcess.h
 * Author: bhorvill
 *
 * Created on May 29, 2013, 4:17 PM
 */

#ifndef EVOLUTIONSEGMENTEDGCPROCESS_H
#define	EVOLUTIONSEGMENTEDGCPROCESS_H

#include "EvolutionSegmentedProcess.h"
#include "Chrono.h"
#include "BrownianMove.h"

class GCMatrix : virtual public EvolutionSegmentedMatrix {
	public :


		GCMatrix(RandomBrownianPath* bridge, Var<RealVector>* up, Var<RealVector>* down, Profile *inrelrate, int k, int ingc)  :
			TransitionMatrix(4), EvolutionSegmentedMatrix(bridge, up, down, 4, k) {
			relrate = inrelrate;
			rStationary = 0;
			gc = ingc;
		  
			ComputeArray();
		}

			//Implements TransitionMatrix::ComputeArray()
			virtual void ComputeR(double lengthLeft, double l) {

			  
				double length = brownianBridge->getLength();

				//Compute the values at the sub-branch extremities, without accounting the brownian bridge (using Thales theorem)					 
				double r1 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);
				double gc1 = up[gc] + lengthLeft/length * (down[gc]-up[gc]);
				lengthLeft+=brownianBridge->getSegmentLength(k);
				double r2 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);
				double gc2 = up[gc] + lengthLeft/length * (down[gc]-up[gc]);

				//Compute the values of the integrated rates (r and omega), and the interval length
				double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, POS_R)) + exp(r2 + brownianBridge->getValue(k+1, POS_R)) );

				double gamma1 = 1.0/(1.0+exp(-(gc1 + brownianBridge->getValue(k, gc))));
				double gamma2 = 1.0/(1.0+exp(-(gc2 + brownianBridge->getValue(k+1, gc))));
				double gamma = (gamma1+gamma2)/2;

				for(int i=0; i<GetNstate(); i++) {
					double tot = 0;
					for(int j=0; j<GetNstate(); j++) {
						if(i!=j) {
	
							R[i][j] = rm*l*RelativeRate(i,j) * ((j == 1 || j == 2) ? gamma/2 : (1-gamma)/2);
							tot+=R[i][j];
						}
					}
					R[i][i]=-tot;
				}


			}

			double RelativeRate(int i, int j) {return (*relrate)[rrindex(i,j,GetNstate())];}

			static int rrindex(int i, int j, int nstate)	{
					return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
			}

			//Compute the stationary from the beginning
			virtual void ComputeStationary() {

				for(int i=0; i<GetNstate(); i++) {
					stationary[i] = 0.25;
				}		  
		}

	private :
	  
		Profile *relrate;			//The 4*4 matrix which gives the nucleotidic transition rate dependantly of the relative rate and the stationary.
		Profile* rStationary;
		int gc;							 //The position of gamma in the multi-dimensional values
};


class RootGCMatrix : virtual public TransitionMatrix, virtual public RandomTransitionMatrix {
	public :
	  RootGCMatrix(Var<Profile>* instationary, int Nstate) : TransitionMatrix(Nstate, 0, instationary->GetArray()) {
		  Register(instationary);

	  }

	  virtual void ComputeArray() {
	  }
	  virtual void ComputeStationary() {
	  }
	  virtual void SetParameters() {
	  }

	private :
};


class EvolutionSegmentedGCBranch : public EvolutionSegmentedBranch {
  
	public :

	EvolutionSegmentedGCBranch(RandomBrownianPath* inbridge, Var<RealVector>* inup, Var<RealVector>* indown, Var<Profile>* inrelrate, int nState, bool paral) :
		EvolutionSegmentedBranch(inbridge, inup, indown, nState, paral) {
		relrate = inrelrate;
		Register(relrate);
		CreateAllP();
	}

	EvolutionSegmentedGCBranch(Var<RealVector>* inup, Var<Profile>* instationary) :
		EvolutionSegmentedBranch(inup, instationary->GetDim()) {
		rStationary = instationary;
		Register(rStationary);
		globalTransMatrix = CreateRootMatrix();
	}

	virtual EvolutionSegmentedMatrix* CreateP(int i) {
		return new GCMatrix(bridge, up, down, relrate, i, 1);
	}

	virtual RandomTransitionMatrix* CreateRootMatrix() { 
		return new RootGCMatrix(rStationary, rStationary->GetDim());
	}

	Var<Profile> *relrate;
	Var<Profile> *rStationary;

};

class RootGCStationary : public Dvar<Profile> {
	public :
		RootGCStationary(Var<UnitReal> *ingamma) {
			dim = 4;
			profile = new double[dim];
			gamma = ingamma;
			Register(gamma);
			specialUpdate();

		}

		void specialUpdate() {
			for(int i=0; i<dim; i++)
				profile[i] = ((i == 1 || i == 2) ? gamma->val()/2 : (1-gamma->val())/2);
		}

		Var<UnitReal> *gamma;
};


class EvolutionSegmentedGCProcess : public EvolutionSegmentedProcess {

public :

	EvolutionSegmentedGCProcess(BrownianProcess *inprocess, Var<Profile>* inrelrate, Var<UnitReal>* rootGCrate, bool paral) :
		EvolutionSegmentedProcess(inprocess, paral) {
		relrate = inrelrate;
		rStationary = new RootGCStationary(rootGCrate);
		RecursiveCreate(GetRoot());
	}

	virtual EvolutionSegmentedBranch* CreateBranchVal(const Link* from) {
	   // grep the instant value associated to this node
		MultiNormal* vup = process->GetInstantProcess()->GetMultiNormal(from);

		if(from->isRoot())
			return new EvolutionSegmentedGCBranch(vup, rStationary);
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = process->GetInstantProcess()->GetMultiNormal(from->Out());

		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = process->GetPureBrownianProcess()->GetRandomBrownianPath(from);

		return new EvolutionSegmentedGCBranch(brownianPath, vup, vdown, relrate, rStationary->GetDim(), paral);
	}



	Var<Profile> *relrate;
	Var<Profile> *rStationary;


};


#endif	/* EVOLUTIONSEGMENTEDGCPROCESS_H */

