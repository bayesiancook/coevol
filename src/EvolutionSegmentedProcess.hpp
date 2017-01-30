#ifndef EVOLUTIONSEGMENTEDPROCESS_H
#define	EVOLUTIONSEGMENTEDPROCESS_H

#include "DAGnode.h"
#include "BrownianBridge.h"
#include "SubMatrix.h"
#include "RandomBrownianPath.h"
#include "BrownianProcess.h"
#include "ValTree.h"


#ifdef _OPENMP
#include <omp.h>
#endif


class EvolutionSegmentedMatrix : virtual public TransitionMatrix {


	public :

		EvolutionSegmentedMatrix(RandomBrownianPath* inbridge, Var<RealVector>* inup, Var<RealVector>* indown, int dim, int ink)  :
			TransitionMatrix(dim) {

			brownianBridge = inbridge;
			up = inup->val().GetArray();
			down = indown->val().GetArray();

			k = ink;
			order = 1;
			nbOrders = new int[11];
			for(int i=0; i<11; i++)
				nbOrders[i] = 0;
		  
			tmp = new double*[GetNstate()];
			for(int i=0; i<GetNstate(); i++) {
				tmp[i] = new double[GetNstate()];
			}
		}

		 //Implements TransitionMatrix::ComputeArray()
		virtual void ComputeArray() {

			 //ChronoMove::chronoMatrix1.Start();


			double lengthLeft = 0;

			for(int i=0; i<k; i++) {
				lengthLeft+=brownianBridge->getSegmentLength(i);
			}

			/* ChronoMove::chronoMatrix1.Stop();
	 
			ChronoMove::chronoMatrix2.Start();*/
			ComputeR(lengthLeft, brownianBridge->getSegmentLength(k));
			/*ChronoMove::chronoMatrix2.Stop();
			ChronoMove::chronoMatrix3.Start();*/
			approachExponential();
			//ChronoMove::chronoMatrix3.Stop();


		}

		void UpdateMatrix() {
			ComputeArrayAndStat();
		}

		virtual void ComputeR(double lengthLeft, double l) = 0;
		virtual void SetParameters() {}

		//Some accessors
	   void setUp(const double* inup) {
		   up = inup;
	   }

	   void setDown(const double* indown) {
		   down = indown;
	   }
	  
	   int GetOrder() {
		   return order;
	   }
	  
	   void ResetNbOrders() {
		   for(int i=0; i<11; i++)
			   nbOrders[i] = 0;
	   }
	  
	   int GetNbOrder(int o) {
		   return nbOrders[o];
	   }
	 
	  
	   void approachExponential() {

				order = 1;
				int z = 0;

				double maxdiag = 0; //Max of diagonal coefficients
				for(int i=0; i<GetNstate(); i++) {
					if(maxdiag< -R[i][i]) {
						maxdiag = -R[i][i];
					}
				}

		 
				while(z<10 && 10*maxdiag > order) {
					order*=2;
					z++;
				}
				nbOrders[z]++;

			   /* if(order>1) {
					ChronoMove::order = true;
					cerr << order << " " << brownianBridge->getLength() << " " << up[0] << " " << down[0] << " " << up[1] << " " << down[1] << endl;
				}   */

				if(order==1) {
					 for(int i=0; i<GetNstate(); i++)
						 R[i][i]+=1;
					 return;

				}
		  

				for(int i=0; i<GetNstate(); i++)
					for(int j=0; j<GetNstate(); j++)
						R[i][j] = R[i][j] / order + (i==j ? 1 : 0);
				//y is now an approximation of exp(R/z)

		  

				PowerOf2(R, order);

			}


		virtual void PowerOf2(double** y, int z) {
			if(z == 1)
				return;

			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					tmp[i][j] = 0;
					for(int k=0; k<GetNstate(); k++) {
						tmp[i][j]+= y[i][k]*y[k][j];
					}
				}
			}
			#ifdef _OPENMP
			//#pragma omp parallel for
			#endif
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++)
					y[i][j] = tmp[i][j];
			}
			EvolutionSegmentedMatrix::PowerOf2(y, z/2);
		}




	protected :
		BrownianBridge *brownianBridge;		//The pure multi-variate brownian bridge of the branch
		const double *up, *down;			//The values at the extremities of the branch

		const static int POS_R=0;		//The position of r in the multi-dimensional values

		int k;							  //Segment position

		double** tmp;					   //A temporary matrix

		int order;			  //The exponential approximation order. It is a power of 2
		int* nbOrders;		  //Stores the number of occurences of order 1, 2, .., 1024 (size 11)
};


class SegmentedBranchGlobalMatrix : virtual public TransitionMatrix, virtual public RandomTransitionMatrix {
	 public :
		 SegmentedBranchGlobalMatrix(Dnode *parent, double** inR, int Nstate) : TransitionMatrix(Nstate, inR, 0) {
			Register(parent);
		}
		virtual void ComputeArray() {
		}
		virtual void ComputeStationary() {
		}

		virtual void SetParameters() {
		}  
};

class EvolutionSegmentedBranch : public Dnode {

	protected:

		RandomBrownianPath* bridge;
		Var<RealVector>* up;
		Var<RealVector>* down;
		EvolutionSegmentedMatrix** P;
		RandomTransitionMatrix* globalTransMatrix;
		double*** cumul;
		bool cumulUpdated;
		int nSegment;
		int nState;
		bool paral;


	public :

		EvolutionSegmentedBranch(RandomBrownianPath* inbridge, Var<RealVector>* inup, Var<RealVector>* indown, int innState, bool inparal) {
					
			SetName("Evolution Segmented Branch");

			bridge = inbridge;
			up = inup;
			nState = innState;
			down = indown;
			paral = inparal;
			Register(bridge);
			Register(up);
			Register(down);
			Register(bridge->getUp());
			Register(bridge->getDown());


			nSegment = bridge->getNSegments();
			P = new EvolutionSegmentedMatrix*[nSegment];

			cumul = new double**[nSegment+1];
			for(int k=0; k<nSegment+1; k++) {
				cumul[k] = new double*[nState];
				for(int i=0; i<nState; i++)
					cumul[k][i] = new double[nState];
			}
			globalTransMatrix = new  SegmentedBranchGlobalMatrix(this, cumul[nSegment], nState);


			cumulUpdated = false;

	 
		}

		EvolutionSegmentedBranch(Var<RealVector>* inup, int innState) {
			SetName("Evolution Root Branch");
			nState = innState;

			bridge = 0;
			up = inup;
			down = 0;
			paral = false;
			Register(up);


			nSegment = 0;
			P = 0;
			cumul = 0;
			cumulUpdated = false;

			//globalTransMatrix = new GlobalTransMatrix(this, rStationary, nState);
		}

		~EvolutionSegmentedBranch() {
			for(int i=0; i<nSegment; i++) {
				delete P[i];
			}
			delete[] P;
			for(int i=0; i<nSegment; i++) {
				for(int j=0; j<nState; j++) {
					delete[] cumul[i][j];
				}
				delete[] cumul[i];
			}
			delete[] cumul;
			delete globalTransMatrix;
		}

		virtual void CreateAllP() {
			for(int i=0; i<nState; i++)
				for(int j=0; j<nState; j++)
					cumul[0][i][j] = (i==j ? 1 : 0);

			for(int k=0; k<nSegment; k++) {
				P[k] = CreateP(k);
				for(int i=0; i<nState; i++)
				   for(int j=0; j<nState; j++) {
					   cumul[k+1][i][j] = 0;
					   for(int l=0; l<nState; l++)
								cumul[k+1][i][j] += cumul[k][i][l]*P[k]->GetR()[l][j];
				   }
			}
			cumulUpdated = true;

		}
		virtual EvolutionSegmentedMatrix* CreateP(int i) = 0;
		virtual RandomTransitionMatrix* CreateRootMatrix() = 0;

		void specialUpdate() {

			if(paral) {
				#ifdef _OPENMP
				// #pragma omp parallel for
				#endif
				for(int i=0; i<GetNsegment(); i++)
					P[i]->UpdateMatrix();
			}
			else {
				for(int i=0; i<GetNsegment(); i++)
					P[i]->UpdateMatrix();
			}
			cumulUpdated = false;
		}

		virtual void updateCumul() {

			for(int k=0; k<nSegment; k++) {
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(int i=0; i<GetNstate(); i++) {
					for(int j=0; j<GetNstate(); j++) {
						cumul[k+1][i][j] = 0;
						for(int h =0; h<GetNstate(); h++) {
							cumul[k+1][i][j] += cumul[k][i][h] * P[k]->GetR()[h][j];
						}
					}
				}
			}

			cumulUpdated = true;
		}

		RandomBrownianPath* GetBridge() { return bridge;}
		TransitionMatrix* GetP(int i) { return P[i]; }
		EvolutionSegmentedMatrix** GetP() {return P;}
		double** GetCumul(int i) {
			if(!cumulUpdated)
				updateCumul();
			return cumul[i];
		}
		RandomTransitionMatrix* GetGlobalTransMatrix() {

			if(!cumulUpdated)
				updateCumul();
			return globalTransMatrix;
		}

		const double* GetStationary() {
			return GetGlobalTransMatrix()->GetStationary();
		}

		bool isRoot() {return !P;}
		int GetNstate() {return nState;}
		int GetNsegment() {return nSegment;}

		void ComputeNbOrders(int* tab) {
			if(!isRoot()) {
				for(int i=0; i<GetNsegment(); i++) {
					for(int j=0; j<11; j++) {
						tab[j] += P[i]->GetNbOrder(j);
					}
					P[i]->ResetNbOrders();
				}
			}
		}

		//Var<Profile> *GetRStationary() { return rStationary;}
 


};



class EvolutionSegmentedProcess : public BranchValPtrTree<EvolutionSegmentedBranch> {

protected:
	BrownianProcess *process;
	int *nbOrders;
	bool paral;
public :
	EvolutionSegmentedProcess(BrownianProcess *inprocess, bool inparal) {
		SetWithRoot(true);
	process = inprocess;
		nbOrders = new int[11];
		paral = inparal;
		for(int i=0; i<11; i++)
			nbOrders[i] = 0;
	}

	 // various accessors
	virtual Tree* GetTree() {
		return process->GetTree();
	}

	Link* GetRoot() { return GetTree()->GetRoot(); }

	void ComputeNbOrders() {

		for(int i=0; i<11; i++) {
			nbOrders[i] = 0;
		}
		ComputeNbOrders(GetRoot());
	}

	void ComputeNbOrders(const Link* from) {

		for(const Link* link = from->Next(); link!=from; link = link->Next()) {
			GetBranchVal(link->GetBranch())->ComputeNbOrders(nbOrders);
			ComputeNbOrders(link->Out());

		}
	}

	int *GetNbOrders() {

		ComputeNbOrders();
		return nbOrders;
	}

};


#endif	/* EVOLUTIONSEGMENTEDPROCESS_H */

