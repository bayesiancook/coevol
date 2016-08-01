#ifndef EVOLUTIONOMEGABROWNIANPROCESS_H
#define	EVOLUTIONOMEGABROWNIANPROCESS_H


#include "RandomSubMatrix.h"
#include "ValTree.h"
#include "BrownianProcess.h"
#include "InverseWishartMatrix.h"


class SegmentOmegaMatrix : public SubMatrix {
	public :
		SegmentOmegaMatrix(SubMatrix* innucMatrix, double inomega, CodonStateSpace *instateSpace) : SubMatrix(instateSpace->GetNstate()) {
			nucMatrix = innucMatrix;
			omega = inomega;
			stateSpace = instateSpace;
		}

		void ComputeArray(int i) {



			double tot = 0;
			for(int j=0; j<GetNstate(); j++) {
				if(i!=j) {
					double s = 0;
					int diff = stateSpace->GetDifferingPosition(i, j);
					if(diff<3)
						s = nucMatrix->GetRow(stateSpace->GetCodonPosition(diff, i))[stateSpace->GetCodonPosition(diff, j)];
					if(stateSpace->Synonymous(i, j))
						s *= omega;

					Q[i][j] = s;
					tot+=s;
				}
			}
			Q[i][i]=-tot;
			if(Q[i][i] == 0) {
				cerr << "Q[i][i] = 0" << endl;
				exit(0);
			}


		}

		void ComputeStationary() {
			for(int i=0; i<GetNstate(); i++)
				mStationary[i] = nucMatrix->GetStationary()[stateSpace->GetCodonPosition(0, i)]*nucMatrix->GetStationary()[stateSpace->GetCodonPosition(1, i)]*nucMatrix->GetStationary()[stateSpace->GetCodonPosition(2, i)];
		}


	private :
		SubMatrix* nucMatrix;
		double omega;
		CodonStateSpace *stateSpace;
};


class EvolutionOmegaMatrix : public TransitionMatrix {

	public :

		//Constructor of a transition matrix associated to a non-root link
		EvolutionOmegaMatrix(RandomBrownianPath* inbridge, const double* inup, const double* indown, SubMatrix* innucMatrix, CodonStateSpace* instateSpace, int inposomega)  :
			TransitionMatrix(instateSpace->GetNstate()) {
			nucMatrix = innucMatrix;
			stateSpace = instateSpace;

			brownianBridge = inbridge;
			dim = brownianBridge->getDim();
			up = inup;
			down = indown;

			posomega = inposomega;

			tmp = new double*[GetNstate()];
				for(int i=0; i<GetNstate(); i++) {
					tmp[i] = new double[GetNstate()];
				}
			root = false;

			ComputeArray();
			ComputeStationary();
		}

		//Special constructor for the root
		EvolutionOmegaMatrix(const double* inup, int indim, SubMatrix* innucMatrix, CodonStateSpace* instateSpace, int inposomega) :
		TransitionMatrix(instateSpace->GetNstate()) {

			nucMatrix = innucMatrix;
			stateSpace = instateSpace;

			up = inup;
			dim = indim;

			down = 0;
			brownianBridge = 0;

			posomega = inposomega;

			tmp = new double*[GetNstate()];
				for(int i=0; i<GetNstate(); i++) {
					tmp[i] = new double[GetNstate()];
			}
			root = true;
			ComputeArray();
			ComputeStationary();

		}

		~EvolutionOmegaMatrix() {
			for(int i=0; i<GetNstate(); i++)
				delete[] tmp[i];
			delete[] tmp;
		}


		//Implements TransitionMatrix::ComputeArray()
		virtual void ComputeArray() {

			bool exactExp = false;

			//Set the array to identity
				for(int i=0; i<GetNstate(); i++)
					for(int j=0; j<GetNstate(); j++)
						R[i][j]= (i==j ? 1 : 0);

			if(!root) {

				int nSegments = brownianBridge->getNSegments();
				double length = brownianBridge->getLength();
				double lengthLeft = 0;

				double** Rinter = new double*[GetNstate()];
					for(int i=0; i<GetNstate(); i++)
						Rinter[i] = new double[GetNstate()];

				for(int k=0; k<nSegments; k++) {

					//Compute the values at the sub-branch extremities, without accounting the brownian bridge (using Thales theorem)					 
					double r1 = up[r] + lengthLeft/length * (down[r]-up[r]);
					double omega1 = up[posomega] + lengthLeft/length * (down[posomega]-up[posomega]);
					lengthLeft+=brownianBridge->getSegmentLength(k);
					double r2 = up[r] + lengthLeft/length * (down[r]-up[r]);
					double omega2 = up[posomega] + lengthLeft/length * (down[posomega]-up[posomega]);

					//Compute the values of the integrated rates (r and omega), and the interval length
					double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, r)) + exp(r2 + brownianBridge->getValue(k+1, r)) );
					double omegam = 0.5 * ( exp(omega1 + brownianBridge->getValue(k, posomega)) + exp(omega2 + brownianBridge->getValue(k+1, posomega)) );


					SubMatrix* segsub = new SegmentOmegaMatrix(nucMatrix, omegam, stateSpace);
					if(exactExp) {
						try {
						   segsub->ComputeExponential(rm*brownianBridge->getSegmentLength(k), Rinter);
						}
						catch (int e) {
							cerr << r1 << " " << r2 << " " << brownianBridge->getValue(k, r) << " " << brownianBridge->getValue(k+1, r) << endl;
							cerr << up[r] << " " << down[r] << " " << length << " " << brownianBridge->getSegmentLength(k) << endl;
							cerr << brownianBridge->getUnitVariance();
							exit(0);
						}
					}
					else {  
						segsub->ApproachExponential(rm*brownianBridge->getSegmentLength(k), Rinter);
					}
					delete segsub;
					multiply(Rinter);

				}

				for(int i=0; i<GetNstate(); i++)
					delete[] Rinter[i];
				delete[] Rinter;

			}
		}

		//Compute the stationary from the beginning
		virtual void ComputeStationary() {
			for(int i=0; i<GetNstate(); i++) {
				stationary[i] = nucMatrix->GetStationary()[stateSpace->GetCodonPosition(0, i)]*nucMatrix->GetStationary()[stateSpace->GetCodonPosition(1, i)]*nucMatrix->GetStationary()[stateSpace->GetCodonPosition(2, i)];
			}
		}

		void multiply(double** Rinter) {
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					tmp[i][j] = 0;
					for(int k=0; k<GetNstate(); k++) {
						tmp[i][j]+= Rinter[i][k]*R[k][j];
					}
				}
			}

			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++)
					R[i][j] = tmp[i][j];
			}

		}

		//Some accessors
		void setUp(const double* inup) {
			up = inup;
		}

		void setDown(const double* indown) {
			down = indown;
		}

		//Statistics methods
		double meanR() {
			return 0;
		}



	protected :
		BrownianBridge *brownianBridge;		//The pure multi-variate brownian bridge of the branch
		const double *up, *down;			//The values at the extremities of the branch
		int dim;				//The variables dimension

		SubMatrix *nucMatrix;		//The 4*4 matrix which gives the nucleotidic transition rate relatively to the relative rate and the stationary.
		CodonStateSpace *stateSpace;

		const static int r=0;		//The position of r in the multi-dimensional values
		int posomega;						  //The position of omega in the multi-dimensional values

		bool root;
		double** tmp;   //An auxilliary matrix for mathematical operations (matricial multiplications, powers...)
};






class RandomEvolutionOmegaMatrix : public EvolutionOmegaMatrix, public RandomTransitionMatrix {

	public :

	RandomEvolutionOmegaMatrix(RandomBrownianPath *inbrownianPath,MultiNormal* inup, MultiNormal* indown, GTRRandomSubMatrixWithNormRates* innucMatrix, CodonStateSpace *instateSpace, int posomega) :
		EvolutionOmegaMatrix(inbrownianPath, inup->GetArray(), indown->GetArray(), innucMatrix, instateSpace, posomega) {
			brownianPath = inbrownianPath;
			nucMatrix = innucMatrix;
			stateSpace = instateSpace;
			up = inup;
			down = indown;
			//stationary = 0;

			SetName("Non commutative transition matrix for nucleotidic sequences");
			Register(up);
			Register(down);
			Register(brownianPath);
			Register(brownianPath->getUp());
			Register(brownianPath->getDown());
			Register(nucMatrix);

		}

		RandomEvolutionOmegaMatrix(MultiNormal* inup, GTRRandomSubMatrixWithNormRates* innucMatrix, CodonStateSpace *instateSpace, int posomega) :
		EvolutionOmegaMatrix(inup->GetArray(), inup->GetDim(), innucMatrix, instateSpace, posomega) {
			up = inup;
			nucMatrix = innucMatrix;
			stateSpace = instateSpace;

			SetName("Non commutative transition matrix for nucleotidic sequences");

			Register(up);
			Register(nucMatrix);


		}

	void SetParameters() {

		}

	void ComputeArray() {EvolutionOmegaMatrix::ComputeArray();}
	void ComputeStationary() {EvolutionOmegaMatrix::ComputeStationary();}

	private :

	RandomBrownianPath *brownianPath;
	MultiNormal *up, *down;
	Var<Profile> *relrate;
		Var<Profile> *rStationary;

		GTRRandomSubMatrixWithNormRates *nucMatrix;
		CodonStateSpace *stateSpace;


};



class EvolutionOmegaBrownianProcess : public BranchValPtrTree<RandomTransitionMatrix> {

public:
	EvolutionOmegaBrownianProcess(BrownianProcess *inprocess, GTRRandomSubMatrixWithNormRates* innucMatrix, CodonStateSpace *instateSpace, int inposomega) : brownianProcess(inprocess), posomega(inposomega) {
		SetWithRoot(true);
		nucMatrix = innucMatrix;
		stateSpace = instateSpace;
		RecursiveCreate(GetRoot());
	}

	 virtual ~EvolutionOmegaBrownianProcess() {
		RecursiveDelete(GetRoot());
	}


	 virtual RandomTransitionMatrix* CreateBranchVal(const Link* link) {
		// grep the instant value associated to this node
		MultiNormal* vup = brownianProcess->GetInstantProcess()->GetMultiNormal(link);
		if(link->isRoot())
			return new RandomEvolutionOmegaMatrix(vup, nucMatrix, stateSpace, posomega);
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = brownianProcess->GetInstantProcess()->GetMultiNormal(link->Out());
		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = brownianProcess->GetPureBrownianProcess()->GetRandomBrownianPath(link);

		// make the new transition matrix, and return the pointer
		return new RandomEvolutionOmegaMatrix(brownianPath, vup, vdown, nucMatrix, stateSpace, posomega);
	}

	void specialUpdate() {
		RecursiveSpecialUpdate(GetRoot());
	}

	void RecursiveSpecialUpdate(const Link* from) {
	  
		for(Link* link=from->Next(); link!=from; link=link->Next()) {
			RecursiveSpecialUpdate(link->Out());
			GetBranchVal(link->GetBranch())->specialUpdate();
		}
	}

	// various accessors
	virtual Tree* GetTree() {
		return brownianProcess->GetTree();
	}


private:
	BrownianProcess* brownianProcess;
	GTRRandomSubMatrixWithNormRates *nucMatrix;
	CodonStateSpace *stateSpace;
	int posomega;
};


#endif	/* EVOLUTIONOMEGABROWNIANPROCESS_H */

