
#ifndef EVOLUTIONSEGMENTEDOMEGAPROCESS_H
#define	EVOLUTIONSEGMENTEDOMEGAPROCESS_H

#include "EvolutionSegmentedProcess.h"
#include "MGCodonTransitionMatrix.h"



#ifdef _OPENMP
#include <omp.h>
#endif


class OmegaMatrix : public EvolutionSegmentedMatrix {
	public :

		OmegaMatrix(RandomBrownianPath* bridge, Var<RealVector>* up, Var<RealVector>* down, MGOmegaCodonTransitionMatrix* incodonMatrix, int k, int inposomega)  :
			TransitionMatrix(incodonMatrix->GetNstate()), EvolutionSegmentedMatrix(bridge, up, down, incodonMatrix->GetNstate(), k) {

			codonMatrix = incodonMatrix;
			posomega = inposomega;


			tmp = new double*[GetNstate()];
			for(int i=0; i<GetNstate(); i++) {
				tmp[i] = new double[GetNstate()];
				for(int j=0; j<GetNstate(); j++) {
					tmp[i][j] = 0;
				}
			}
			order = 2;
			ComputeArray();
		}

		//Compute the stationary from the beginning
		virtual void ComputeStationary() {

			for(int i=0; i<GetNstate(); i++) {
				stationary[i] = codonMatrix->GetStationary()[i];
			}		  
		}

		virtual void ComputeR(double lengthLeft, double l) {

			double length = brownianBridge->getLength();
		  
			//Compute the values at the sub-branch extremities, without accounting the brownian bridge (using Thales theorem)					 
			double r1 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);
			double omega1 = up[posomega] + lengthLeft/length * (down[posomega]-up[posomega]);
			lengthLeft+=l;
			double r2 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);
			double omega2 = up[posomega] + lengthLeft/length * (down[posomega]-up[posomega]);

			//Compute the values of the integrated rates (r and omega), and the interval length
			double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, POS_R)) + exp(r2 + brownianBridge->getValue(k+1, POS_R)) );
			double omegam = 0.5 * ( exp(omega1 + brownianBridge->getValue(k, posomega)) + exp(omega2 + brownianBridge->getValue(k+1, posomega)) );


			if(GetOrder() > 1) {
				for(int i = 0; i< GetNstate(); i++) {
					for(int j=0; j< GetNstate(); j++) {
						R[i][j] = 0;
					}
				}
			}
			codonMatrix->ComputeValues(R, rm*l, omegam);

					  
		}

		virtual void PowerOf2(double** y, int z) {
			if(z == 1)
				return;
			codonMatrix->square(y, tmp);
			#ifdef _OPENMP
			//#pragma omp parallel for
			#endif
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++)
					y[i][j] = tmp[i][j];
			}
			EvolutionSegmentedMatrix::PowerOf2(y, z/2);
		}

	private :

		MGOmegaCodonTransitionMatrix *codonMatrix;

		double** tmp;

		int posomega;						  //The position of omega in the multi-dimensional values

};

class RootOmegaMatrix : virtual public TransitionMatrix, virtual public RandomTransitionMatrix {
	public :
	  RootOmegaMatrix(MGCodonTransitionMatrix* inm) : TransitionMatrix(inm->GetNstate(), 0, new double[inm->GetNstate()]) {
		  m = inm;
		  Register(m);
	  }

	   virtual void ComputeArray() {
	   }
	   virtual void ComputeStationary() {
		   for(int i=0; i<GetNstate(); i++) {
			   stationary[i] = m->GetStationary()[i];
		   }
		}

	   virtual void SetParameters() {
	   }

	private :
		MGCodonTransitionMatrix* m;
};

class EvolutionSegmentedOmegaBranch : public EvolutionSegmentedBranch {
  
	public :

	EvolutionSegmentedOmegaBranch(RandomBrownianPath* inbridge, Var<RealVector>* inup, Var<RealVector>* indown, MGOmegaCodonTransitionMatrix* incodonMatrix, bool paral) :
		EvolutionSegmentedBranch(inbridge, inup, indown, incodonMatrix->GetNstate(), paral) {
		codonMatrix = incodonMatrix;
		Register(codonMatrix);
		CreateAllP();
	}

	EvolutionSegmentedOmegaBranch(Var<RealVector>* inup, MGOmegaCodonTransitionMatrix* incodonMatrix) :
		EvolutionSegmentedBranch(inup, incodonMatrix->GetNstate()) {
		codonMatrix = incodonMatrix;
		Register(codonMatrix);
		globalTransMatrix = CreateRootMatrix();
	}

	virtual EvolutionSegmentedMatrix* CreateP(int i) {
		return new OmegaMatrix(bridge, up, down, codonMatrix, i, 1);
	}

	virtual RandomTransitionMatrix* CreateRootMatrix() {
		return new RootOmegaMatrix(codonMatrix);
	}


private :

	MGOmegaCodonTransitionMatrix* codonMatrix;



};


class EvolutionSegmentedOmegaProcess : public EvolutionSegmentedProcess {

public :

	EvolutionSegmentedOmegaProcess(BrownianProcess *inprocess, MGOmegaCodonTransitionMatrix* incodonMatrix, bool paral) :
		EvolutionSegmentedProcess(inprocess, paral) {
		codonMatrix = incodonMatrix;
		RecursiveCreate(GetRoot());
	}

	virtual EvolutionSegmentedBranch* CreateBranchVal(const Link* from) {
	   // grep the instant value associated to this node
		MultiNormal* vup = process->GetInstantProcess()->GetMultiNormal(from);

		if(from->isRoot())
			return new EvolutionSegmentedOmegaBranch(vup, codonMatrix);
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = process->GetInstantProcess()->GetMultiNormal(from->Out());

		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = process->GetPureBrownianProcess()->GetRandomBrownianPath(from);

		return new EvolutionSegmentedOmegaBranch(brownianPath, vup, vdown, codonMatrix, paral);
	}


	MGOmegaCodonTransitionMatrix* codonMatrix;


};




#endif	/* EVOLUTIONSEGMENTEDOMEGAPROCESS_H */

