/*
 * File:   EvolutionSegmentedNUCProcess.h
 * Author: bhorvill
 *
 * Created on May 29, 2013, 4:17 PM
 */

#ifndef EVOLUTIONSEGMENTEDNUCPROCESS_H
#define	EVOLUTIONSEGMENTEDNUCPROCESS_H

#include "EvolutionSegmentedProcess.h"
#include "GTRSubMatrix.h"



class NucMatrix : virtual public EvolutionSegmentedMatrix {
	public :


		NucMatrix(RandomBrownianPath* bridge, Var<RealVector>* up, Var<RealVector>* down, GTRRandomSubMatrixWithNormRates* innucMatrix, int k)  :
			TransitionMatrix(4), EvolutionSegmentedMatrix(bridge, up, down, 4, k) {
			nucMatrix = innucMatrix;
		  
			ComputeArray();
		}


		//Implements TransitionMatrix::ComputeArray()
		virtual void ComputeR(double lengthLeft, double l) {

			double length = brownianBridge->getLength();

			//Compute the values at the sub-branch extremities, without accounting the brownian bridge (using Thales theorem)					 
			double r1 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);
			lengthLeft+=brownianBridge->getSegmentLength(k);
			double r2 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);

			//Compute the values of the integrated rate, and the interval length
			double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, POS_R)) + exp(r2 + brownianBridge->getValue(k+1, POS_R)) );

			for(int i=0; i<GetNstate(); i++) {
				double tot = 0;
				for(int j=0; j<GetNstate(); j++) {
					if(i!=j) {
						R[i][j] = rm*l*nucMatrix->GetRow(i)[j];
						tot+=R[i][j];
					}
				}
				R[i][i]=-tot;
			}


		}


		//Compute the stationary from the beginning
		virtual void ComputeStationary() {
			for(int i=0; i<GetNstate(); i++)
				stationary[i] = nucMatrix->GetStationary()[i];
		}

	private :
	  
		GTRRandomSubMatrixWithNormRates* nucMatrix;
};


class RootNucMatrix : virtual public TransitionMatrix, virtual public RandomTransitionMatrix {
	public :
	  RootNucMatrix(GTRRandomSubMatrixWithNormRates* innucMatrix, int Nstate) : TransitionMatrix(Nstate, 0, new double[Nstate]) {
		  nucMatrix = innucMatrix;
		  Register(nucMatrix);
	  }

	  virtual void ComputeArray() {
	  }
	  virtual void ComputeStationary() {
		  for(int i=0; i<GetNstate(); i++) {
			   stationary[i] = nucMatrix->GetStationary()[i];
		   }
	  }
	  virtual void SetParameters() {
	  }

	private :
	  GTRRandomSubMatrixWithNormRates* nucMatrix;


};


class EvolutionSegmentedNucBranch : public EvolutionSegmentedBranch {
  
	public :

	EvolutionSegmentedNucBranch(RandomBrownianPath* inbridge, Var<RealVector>* inup, Var<RealVector>* indown, GTRRandomSubMatrixWithNormRates* innucMatrix, bool paral) :
		EvolutionSegmentedBranch(inbridge, inup, indown, innucMatrix->GetNstate(), paral) {
		nucMatrix = innucMatrix;
		Register(nucMatrix);
		CreateAllP();
	}

	EvolutionSegmentedNucBranch(Var<RealVector>* inup, GTRRandomSubMatrixWithNormRates* innucMatrix) :
		EvolutionSegmentedBranch(inup, innucMatrix->GetNstate()) {
		nucMatrix = innucMatrix;
		Register(nucMatrix);
		globalTransMatrix = CreateRootMatrix();
	}

	virtual EvolutionSegmentedMatrix* CreateP(int i) {
		return new NucMatrix(bridge, up, down, nucMatrix, i);
	}

	virtual RandomTransitionMatrix* CreateRootMatrix() { 
		return new RootNucMatrix(nucMatrix, nucMatrix->GetNstate());
	}

	GTRRandomSubMatrixWithNormRates* nucMatrix;

};


class EvolutionSegmentedNucProcess : public EvolutionSegmentedProcess {

public :

	EvolutionSegmentedNucProcess(BrownianProcess *inprocess, GTRRandomSubMatrixWithNormRates* innucMatrix, bool paral) :
		EvolutionSegmentedProcess(inprocess, paral) {
		nucMatrix = innucMatrix;
		RecursiveCreate(GetRoot());
	}

	virtual EvolutionSegmentedBranch* CreateBranchVal(const Link* from) {
	   // grep the instant value associated to this node
		MultiNormal* vup = process->GetInstantProcess()->GetMultiNormal(from);

		if(from->isRoot())
			return new EvolutionSegmentedNucBranch(vup, nucMatrix);
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = process->GetInstantProcess()->GetMultiNormal(from->Out());

		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = process->GetPureBrownianProcess()->GetRandomBrownianPath(from);

		return new EvolutionSegmentedNucBranch(brownianPath, vup, vdown, nucMatrix, paral);
	}


private:
	GTRRandomSubMatrixWithNormRates* nucMatrix;


};


#endif	/* EVOLUTIONSEGMENTEDNUCPROCESS_H */

