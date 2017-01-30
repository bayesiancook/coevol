#ifndef EVOLUTIONNUCBROWNIANPROCESS_H
#define EVOLUTIONNUCBROWNIANPROCESS_H

#include "RandomSubMatrix.h"
#include "ValTree.h"
#include "BrownianProcess.h"
#include "InverseWishartMatrix.h"


class EvolutionNucMatrix : public TransitionMatrix {

	public :

			//Constructor of a transition matrix associated to a non-root link
			EvolutionNucMatrix(BrownianBridge *inbrownianBridge, const double* inup, const double* indown, SubMatrix *innucMatrix, int ingc=-1) :
					TransitionMatrix(innucMatrix->GetNstate()) {
					nucMatrix = innucMatrix;

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
			}


			//Special constructor for the root
			EvolutionNucMatrix(const double* inup, int indim, SubMatrix *innucMatrix, int ingc = -1) :
					TransitionMatrix(innucMatrix->GetNstate()) {

					nucMatrix = innucMatrix;
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

			}


			~EvolutionNucMatrix() {
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
						lengthLeft+=brownianBridge->getSegmentLength(k);
						double r2 = up[r] + lengthLeft/length * (down[r]-up[r]);

						//Compute the values of the integrated rates (r and omega), and the interval length
						double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, r)) + exp(r2 + brownianBridge->getValue(k+1, r)) );
						if(exactExp) {
							try {
								nucMatrix->ComputeExponential(rm*brownianBridge->getSegmentLength(k), Rinter);
							}
							catch (int e) {
								cerr << r1 << " " << r2 << " " << brownianBridge->getValue(k, r) << " " << brownianBridge->getValue(k+1, r) << endl;
								cerr << up[r] << " " << down[r] << " " << length << " " << brownianBridge->getSegmentLength(k) << endl;
								cerr << brownianBridge->getUnitVariance();
								exit(0);
							}
						}
						else {
							//approachExp(rm*brownianBridge->getSegmentLength(k), Rinter);
							nucMatrix->ApproachExponential(rm*brownianBridge->getSegmentLength(k), Rinter);
						}
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
						stationary[i] = nucMatrix->Stationary(i);
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

		   /* void pow(double** y, int z) {
				if(z == 1)
					return;
				for(int i=0; i<GetNstate(); i++) {
					for(int j=0; j<GetNstate(); j++) {
						tmp[i][j] = 0;
						for(int k=0; k<GetNstate(); k++) {
							tmp[i][j]+= y[i][k]*y[k][j];
						}
					}
				}
				for(int i=0; i<GetNstate(); i++) {
					for(int j=0; j<GetNstate(); j++)
						y[i][j] = tmp[i][j];
				}
				pow(y, z/2);
			}*/

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

		   /* void approachExp(double l, double** y) {
				int z = 1;  //Reducing coefficient

				double maxdiag = 0; //Max of diagonal coefficients
				for(int i=0; i<GetNstate(); i++)
					if(maxdiag< (*nucMatrix)(i,i)*l)
						maxdiag = (*nucMatrix)(i,i)*l;

				while(z<1024 && maxdiag > z)
					z*=2;


				for(int i=0; i<GetNstate(); i++)
					for(int j=0; j<GetNstate(); j++)
						y[i][j] = l* (*nucMatrix)(i,j) / z + (i==j ? 1 : 0);
				//y is now an approximation of exp(NucMatrix/z)

				this->pow(y, z);

				for(int i=0; i<GetNstate(); i++) {
					double tot = 0;
					for(int j=0; j<GetNstate(); j++)
						tot+=y[i][j];
					if(fabs(tot-1)>10e-8) {
						cerr << "Error in approachExp : row does not sum to 1 " << endl;
						cerr << "Sum : " << tot;
						exit(1);
					}
				}
			}*/


	protected :
			BrownianBridge *brownianBridge;		//The pure multi-variate brownian bridge of the branch
			const double *up, *down;			//The values at the extremities of the branch
			int dim;				//The variables dimension

			SubMatrix *nucMatrix;			//The 4*4 matrix which gives the nucleotidic transition rate dependantly of the relative rate and the stationary.

			const static int r=0;		//The position of r in the multi-dimensional values
			int gc;		//The position of omega in the multi-dimensional values

			bool root;
			double** tmp;   //An auxilliary matrix for mathematical operations (matricial multiplications, powers...)
};






class RandomEvolutionNucMatrix : public EvolutionNucMatrix, public RandomTransitionMatrix {

	public :

	RandomEvolutionNucMatrix(RandomBrownianPath *inbrownianPath,MultiNormal* inup, MultiNormal* indown,RandomSubMatrix *innucMatrix) :
		EvolutionNucMatrix(inbrownianPath, inup->GetArray(), indown->GetArray(), innucMatrix) {
			brownianPath = inbrownianPath;
			up = inup;
			down = indown;
			nucMatrix = innucMatrix;

			SetName("Non commutative transition matrix for nucleotidic sequences");

			Register(up);
			Register(down);
			Register(brownianPath);
			Register(brownianPath->getUp());
			Register(brownianPath->getDown());
			Register(nucMatrix);

		}

		RandomEvolutionNucMatrix(MultiNormal* inup, RandomSubMatrix *innucMatrix) :
		EvolutionNucMatrix(inup->GetArray(), inup->GetDim(), innucMatrix) {
			up = inup;
			nucMatrix = innucMatrix;

			SetName("Non commutative transition matrix for nucleotidic sequences");

			Register(up);
			Register(nucMatrix);

		}

	void SetParameters() {

		}

	void ComputeArray() {EvolutionNucMatrix::ComputeArray();}
	void ComputeStationary() {EvolutionNucMatrix::ComputeStationary();}

	private :

	RandomBrownianPath *brownianPath;
	MultiNormal *up, *down;
	RandomSubMatrix *nucMatrix;


};



class EvolutionNucBrownianProcess : public BranchValPtrTree<RandomTransitionMatrix> {

public:
	EvolutionNucBrownianProcess(BrownianProcess *inprocess, RandomSubMatrix* innucMatrix, int ingc, int intt) : brownianProcess(inprocess), nucMatrix(innucMatrix), gc(ingc) {
		SetWithRoot(true);
		RecursiveCreate(GetRoot());
	}

	 virtual ~EvolutionNucBrownianProcess() {
		RecursiveDelete(GetRoot());
	}


	 virtual RandomTransitionMatrix* CreateBranchVal(const Link* link) {

		// grep the instant value associated to this node
		MultiNormal* vup = brownianProcess->GetInstantProcess()->GetMultiNormal(link);
		if(link->isRoot())
			return new RandomEvolutionNucMatrix(vup, nucMatrix);
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = brownianProcess->GetInstantProcess()->GetMultiNormal(link->Out());
		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = brownianProcess->GetPureBrownianProcess()->GetRandomBrownianPath(link);

		// make the new transition matrix, and return the pointer
		return new RandomEvolutionNucMatrix(brownianPath, vup, vdown, nucMatrix);
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
	RandomSubMatrix* nucMatrix;
	int gc;
};



#endif

