#ifndef MGCODONTRANSITIONMATRIX_H
#define	MGCODONTRANSITIONMATRIX_H

#include "SubMatrix.h"
#include "CodonStateSpace.h"
#include "GTRSubMatrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

class MGCodonTransitionMatrix : public TransitionMatrix, public RandomTransitionMatrix {

	public :
		MGCodonTransitionMatrix(GTRRandomSubMatrixWithNormRates *innucMatrix, CodonStateSpace *instateSpace) : TransitionMatrix(instateSpace->GetNstate()) {
			nucMatrix = innucMatrix;
			stateSpace = instateSpace;

			Register(nucMatrix);

			syn = new bool*[GetNstate()];
			for(int i=0; i<GetNstate(); i++) {
				syn[i] = new bool[GetNstate()];
				for(int j=0; j<GetNstate(); j++) {
					syn[i][j] = stateSpace->Synonymous(i,j);
					R[i][j] = 0;
				}
			}
		}

		virtual void ComputeArray() {
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					if(i!=j) {
						int diff = stateSpace->GetDifferingPosition(i, j);
						if(diff < 3)
							R[i][j] = nucMatrix->GetRow(stateSpace->GetCodonPosition(diff, i))[stateSpace->GetCodonPosition(diff, j)];
					}
				}
			}
		}

	virtual void ComputeStationary() {
			double tot = 0;
			for(int i=0; i<GetNstate(); i++) {
				stationary[i] = nucMatrix->GetStationary()[stateSpace->GetCodonPosition(0, i)]*nucMatrix->GetStationary()[stateSpace->GetCodonPosition(1, i)]*nucMatrix->GetStationary()[stateSpace->GetCodonPosition(2, i)];
				tot+=stationary[i];
			}
			for(int i=0; i<GetNstate(); i++) {
				stationary[i] /= tot;
			}
		}

		bool isSyn(int i, int j) {
			return syn[i][j];
		}

		virtual void SetParameters() {

		}

		CodonStateSpace* GetStateSpace() {
			return stateSpace;
		}

		SubMatrix* GetNucMatrix() {
			return nucMatrix;
		}


	protected :
		GTRRandomSubMatrixWithNormRates* nucMatrix;
		CodonStateSpace *stateSpace;

		bool** syn;


};

class MGOmegaCodonTransitionMatrix : public MGCodonTransitionMatrix {

	public :
		MGOmegaCodonTransitionMatrix(GTRRandomSubMatrixWithNormRates *innucMatrix, CodonStateSpace *instateSpace) : MGCodonTransitionMatrix(innucMatrix, instateSpace) {
			rows = new int[nVal];
			columns = new int[nVal];

			int k=0;
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					if(i==j || stateSpace->GetDifferingPosition(i,j)<3) {
						rows[k] = i;
						columns[k] = j;
						k++;
					}
				}
			}


			ki = new int[nFactor];
			kj = new int[nFactor];
			k = 0;

			for(int i = 0; i<nVal; i++) {
				for(int j=0; j<nVal; j++) {
					if(columns[i] == rows[j]) {
						ki[k] = i;
						kj[k] = j;
						k++;
					}
				}
			}


		  

		}

		virtual void ComputeValues(double** tab, double r, double omega) {

			int currentrow = 0;
			double tot = 0;

			for(int k = 0; k < nVal; k++) {
				int i = rows[k], j = columns[k];
				if(i!=currentrow) {
					tab[currentrow][currentrow] = -tot;
					tot=0;
					currentrow = i;
				}
				if(rows[k]!=columns[k]) {
					if(isSyn(i, j)) {
						tab[i][j] = r*R[i][j];
					}
					else
						tab[i][j] = r*omega*R[i][j];
					tot += tab[i][j];
				}

				tab[currentrow][currentrow] = -tot;
			}
					 
		}

	  
		void square(double** input, double** output) {
 
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					output[i][j] = 0;
				}
			}
			for(int k = 0; k <nFactor; k++) {
				output[rows[ki[k]]][columns[kj[k]]] += input[rows[ki[k]]][columns[ki[k]]]*input[rows[kj[k]]][columns[kj[k]]];
			}
		}

	protected :

		int *rows, *columns;
		int *ki, *kj;

		const static int nVal = 587;
		const static int nFactor = 5673;

};



class MGOmegaCodonTransitionMatrix_OMP : public MGOmegaCodonTransitionMatrix {

public :
	MGOmegaCodonTransitionMatrix_OMP(GTRRandomSubMatrixWithNormRates *innucMatrix, CodonStateSpace *instateSpace, int innThreads) : MGOmegaCodonTransitionMatrix(innucMatrix, instateSpace) {
		nThreads = innThreads;

		pSyn = new double**[nThreads];
		pR = new double**[nThreads];
		pnVal = new double[nThreads];
		pRows = new int*[nThreads];
		pColumns = new int*[nThreads];

		for(int k=0; k<nThreads; k++) {
			pSyn[k] = new double*[GetNstate()];
			pR[k] = new double*[GetNstate()];
			pnVal[k] = nVal;
			pRows[k] = new int[nVal];
			pColumns[k] = new int[nVal];

			for(int i=0; i<GetNstate(); i++) {
				pSyn[k][i] = new double[GetNstate()];
				pR[k][i] = new double[GetNstate()];
				for(int j=0; j<GetNstate(); j++) {
					pSyn[k][i][j] = syn[i][j];
					pR[k][i][j] = R[i][j];
				}
			}

			for(int i=0; i<nVal; i++) {
				pRows[k][i] = rows[i];
				pColumns[k][i] = columns[i];
			}
		}

	}

   virtual void ComputeValues(double** tab, double r, double omega) {

	   int thr = 1;
	   // int thr = omp_get_thread_num();
	  
	 
		int currentrow = 0;
		double tot = 0;

		for(int k = 0; k < pnVal[thr]; k++) {
			int i = pRows[thr][k], j = pColumns[thr][k];
			if(i!=currentrow) {
				tab[currentrow][currentrow] = -tot;
				tot=0;
				currentrow = i;
			}
			if(pRows[thr][k]!=pColumns[thr][k]) {
				if(pSyn[thr][i][j]) {
					tab[i][j] = r*pR[thr][i][j];
				}
				else
					tab[i][j] = r*omega*pR[thr][i][j];
				tot += tab[i][j];
			}
			tab[currentrow][currentrow] = -tot;
		}

	}
  
   void updatepR() {
	   for(int k=0; k<nThreads; k++) {		  
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					pR[k][i][j] = R[i][j];
				}
			}
	   }

   }
  
	virtual void ComputeArray() {
		MGCodonTransitionMatrix::ComputeArray();
		for(int k=0; k<nThreads; k++) {
			for(int i=0; i<GetNstate(); i++) {
				for(int j=0; j<GetNstate(); j++) {
					pR[k][i][j] = R[i][j];
				}
			}
		}
	}



protected :
	int nThreads;
	double*** pSyn;
	double*** pR;
	double* pnVal;
	int** pRows, **pColumns;


};

#endif	/* MGCODONTRANSITIONMATRIX_H */

