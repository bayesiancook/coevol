
#ifndef BROWNIANGCSTATIONARY_H
#define	BROWNIANGCSTATIONARY_H

#include "RandomSubMatrix.h"
#include "ValTree.h"
#include "BrownianProcess.h"
#include "InverseWishartMatrix.h"

class RootStationary : public Dvar<Profile> {
	public :
		RootStationary(Var<UnitReal> *ingamma) {
			this->dim = 4;
			this->profile = new double[4];
			gamma = ingamma;
			Register(gamma);
			specialUpdate();

		}

		void specialUpdate() {
			for(int i=0; i<4; i++)
				this->GetArray()[i] = ((i == 1 || i == 2) ? gamma->val()/2 : (1-gamma->val())/2);
		}

		Var<UnitReal> *gamma;
};


class SegmentGCMatrix : public SubMatrix {
	public :
		SegmentGCMatrix(Profile *inrelrate, double ingamma) : SubMatrix(4) {
			relrate = inrelrate;
			gamma = ingamma;
		}

		void ComputeArray(int i) {
			double tot = 0;
			for(int j=0; j<GetNstate(); j++) {
				if(i!=j) {
					Q[i][j] = RelativeRate(i,j) * ((j == 1 || j == 2) ? gamma/2 : (1-gamma)/2);
					tot+=Q[i][j];
				}
			}
			Q[i][i]=-tot;

		}

		void ComputeStationary() {
			for(int i=0; i<GetNstate(); i++)
				mStationary[i] = ((i == 1 || i == 2) ? gamma/2 : (1-gamma)/2);
		}

		double RelativeRate(int i, int j) {return (*relrate)[rrindex(i,j,GetNstate())];}

	static int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}


	private :
		Profile *relrate;
		double gamma;
};

class EvolutionGCMatrix : public TransitionMatrix {

	public :

			//Constructor of a transition matrix associated to a non-root link
			EvolutionGCMatrix(BrownianBridge *inbrownianBridge, const double* inup, const double* indown, Profile *inrelrate, int ingc) :
					TransitionMatrix(4) {
					relrate = inrelrate;
					rStationary = 0;

					brownianBridge = inbrownianBridge;
					dim = brownianBridge->getDim();
					up = inup;
					down = indown;

					gc = ingc;
					tmp = new double*[GetNstate()];
					for(int i=0; i<GetNstate(); i++) {
						tmp[i] = new double[GetNstate()];
					}
					root = false;

					ComputeArray();
					ComputeStationary();
			}


			//Special constructor for the root
			EvolutionGCMatrix(const double* inup, int indim, Profile *inrelrate, Profile *instationary, int ingc) :
					TransitionMatrix(4) {

					relrate = inrelrate;
					rStationary = instationary;
					up = inup;
					dim = indim;

					down = 0;
					brownianBridge = 0;

					gc = ingc;

					tmp = new double*[GetNstate()];
					for(int i=0; i<GetNstate(); i++) {
						tmp[i] = new double[GetNstate()];
					}
					root = true;
					ComputeArray();
					ComputeStationary();



			}


			~EvolutionGCMatrix() {
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
						double gc1 = up[gc] + lengthLeft/length * (down[gc]-up[gc]);
						lengthLeft+=brownianBridge->getSegmentLength(k);
						double r2 = up[r] + lengthLeft/length * (down[r]-up[r]);
						double gc2 = up[gc] + lengthLeft/length * (down[gc]-up[gc]);

						//Compute the values of the integrated rates (r and omega), and the interval length
						double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, r)) + exp(r2 + brownianBridge->getValue(k+1, r)) );

						double gamma1 = 1.0/(1.0+exp(-(gc1 + brownianBridge->getValue(k, gc))));
						double gamma2 = 1.0/(1.0+exp(-(gc2 + brownianBridge->getValue(k+1, gc))));
						double gamma = (gamma1+gamma2)/2;

						SubMatrix* segsub = new SegmentGCMatrix(relrate, gamma);
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

				if(root) {

					for(int i=0; i<GetNstate(); i++) {
						stationary[i] = rStationary->GetArray()[i];
					}
				}
				else {
					for(int i=0; i<GetNstate(); i++) {
						stationary[i] = 0.25;
					}		  
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

			Profile *relrate;			//The 4*4 matrix which gives the nucleotidic transition rate dependantly of the relative rate and the stationary.
			Profile *rStationary;			   //root stationary

			const static int r=0;		//The position of r in the multi-dimensional values
			int gc;		//The position of omega in the multi-dimensional values

			bool root;
			double** tmp;   //An auxilliary matrix for mathematical operations (matricial multiplications, powers...)
};






class RandomEvolutionGCMatrix : public EvolutionGCMatrix, public RandomTransitionMatrix {

	public :

	RandomEvolutionGCMatrix(RandomBrownianPath *inbrownianPath,MultiNormal* inup, MultiNormal* indown,Var<Profile> *inrelrate, int  gc) :
		EvolutionGCMatrix(inbrownianPath, inup->GetArray(), indown->GetArray(), inrelrate, gc) {
			brownianPath = inbrownianPath;
			up = inup;
			down = indown;
			relrate = inrelrate;
			//stationary = 0;

			SetName("Non commutative transition matrix for nucleotidic sequences");
			Register(up);
			Register(down);
			Register(brownianPath);
			Register(brownianPath->getUp());
			Register(brownianPath->getDown());
			Register(relrate);

		}

		RandomEvolutionGCMatrix(MultiNormal* inup, Var<Profile> *inrelrate, Var<Profile> *instationary, int gc) :
		EvolutionGCMatrix(inup->GetArray(), inup->GetDim(), inrelrate, instationary, gc) {
			up = inup;
			relrate = inrelrate;
			rStationary = instationary;

			SetName("Non commutative transition matrix for nucleotidic sequences");

			Register(up);
			Register(relrate);
			Register(rStationary);


		}

	void SetParameters() {

		}

	void ComputeArray() {EvolutionGCMatrix::ComputeArray();}
	void ComputeStationary() {EvolutionGCMatrix::ComputeStationary();}

	private :

	RandomBrownianPath *brownianPath;
	MultiNormal *up, *down;
	Var<Profile> *relrate;
		Var<Profile> *rStationary;


};



class EvolutionGCBrownianProcess : public BranchValPtrTree<RandomTransitionMatrix> {

public:
	EvolutionGCBrownianProcess(BrownianProcess *inprocess, Var<Profile>* inrelrate, Var<UnitReal>* rootGCrate, int ingc) : brownianProcess(inprocess), relrate(inrelrate), gc(ingc) {
		SetWithRoot(true);
		rStationary = new RootStationary(rootGCrate);
		RecursiveCreate(GetRoot());
	}

	 virtual ~EvolutionGCBrownianProcess() {
		RecursiveDelete(GetRoot());
	}


	 virtual RandomTransitionMatrix* CreateBranchVal(const Link* link) {
		// grep the instant value associated to this node
		MultiNormal* vup = brownianProcess->GetInstantProcess()->GetMultiNormal(link);
		if(link->isRoot())
			return new RandomEvolutionGCMatrix(vup, relrate, rStationary, gc);
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = brownianProcess->GetInstantProcess()->GetMultiNormal(link->Out());
		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = brownianProcess->GetPureBrownianProcess()->GetRandomBrownianPath(link);

		// make the new transition matrix, and return the pointer
		return new RandomEvolutionGCMatrix(brownianPath, vup, vdown, relrate, gc);
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
	Var<Profile> *relrate;
	Var<Profile> *rStationary;
	int gc;
};


#endif	/* BROWNIANGCSTATIONARY_H */

