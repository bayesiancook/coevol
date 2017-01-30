
#ifndef SPLITAAMAT
#define SPLITAAMAT

#include "SimilarityMatrix.h"

class SplitAAMatrix {

	public:

	SplitAAMatrix(CodonStateSpace* instatespace)	{
		statespace = instatespace;
		NonCTNearest = new int*[Naa];
		// NonCTNearest = new int*[statespace->GetNstate()];
		// CTtype = new int*[statespace->GetNstate()];
		for (int i=0; i<Naa; i++)	{
		// for (int i=0; i<statespace->GetNstate(); i++)	{
			NonCTNearest[i] = new int[Naa];
			// NonCTNearest[i] = new int[statespace->GetNstate()];
			// CTtype[i] = new int[statespace->GetNstate()];
		}
		MakeArrays();
	}

	bool IsNonCTNearest(int a, int b)	{
		return NonCTNearest[a][b];
	}

	void Print(ostream& os, SimilarityMatrix* mat)	{
		for (int i=0; i<Naa; i++)	{
			for (int j=0; j<Naa; j++)	{
				if (i != j)	{
					if (NonCTNearest[i][j] != -1)	{
						os << statespace->GetProteinStateSpace()->GetState(i) << '\t' << statespace->GetProteinStateSpace()->GetState(j) << '\t' << mat->isRadical(i,j) << '\t' << NonCTNearest[i][j] << '\n';
					}
				}
			}
		}
	}

	private:

	void MakeArrays()	{
		/*
		for (int i=0; i<statespace->GetNstate(); i++)	{
			for (int j=0; j<statespace->GetNstate(); j++)	{
		*/
		for (int i=0; i<Naa; i++)	{
			for (int j=0; j<Naa; j++)	{
				if (i != j)	{
					NonCTNearest[i][j] = statespace->IsNonCTNearest(i,j);
				}
				else	{
					NonCTNearest[i][j] = -1;
				}
			}
		}
	}

	int** NonCTNearest;
	CodonStateSpace* statespace;
};

#endif

