#include "AminoAcidOmegaSubMatrix.h"


///	AAOmegaSubMatrix

void AAOmegaSubMatrix :: ComputeArray(int n){
	double total = 0;
	for (int i=0; i< Naa; i++){
		if(n!=i){
			Q[n][i] = RelativeRate(n,i) * mStationary[i];
			if(mSimilarityMatrix->isRadical(n,i))
				Q[n][i] *= GetOmega();
			total += Q[n][i];
		}
	}
	Q[n][n] = -total;
}

void RandomAAOmegaSubMatrix :: SetParameters(){
	for (int i=0; i<relrate->GetDim(); i++)	{
		rescaledrelrate[i] = (*relrate)[i] * relrate->GetDim();
	}
	SetRelativeRate(rescaledrelrate);
	CopyStationary(stat->val().GetArray());
	setOmega(randomOmega->val());
}

void SplitAAOmegaSubMatrix :: ComputeArray(int n){
	double total = 0;
	for (int i=0; i< Naa; i++){
		if(n!=i){
			Q[n][i] = RelativeRate(n,i) * mStationary[i];
			if (mSplitMat->IsNonCTNearest(n,i) == 1)	{
				if(mSimilarityMatrix->isRadical(n,i))	{
					Q[n][i] *= GetOmegaTv();
				}
			}
			else	{
				if (splittype && (mSplitMat->IsNonCTNearest(n,i) == -1))	{
					Q[n][i] = 0;
				}
				else	{
					Q[n][i] *= GetTsTv();
					if(mSimilarityMatrix->isRadical(n,i))	{
						Q[n][i] *= GetOmegaTs();
					}
				}
			}
			total += Q[n][i];
		}
	}
	Q[n][n] = -total;
}


void RandomSplitAAOmegaSubMatrix :: SetParameters(){
	for (int i=0; i<relrate->GetDim(); i++)	{
		rescaledrelrate[i] = (*relrate)[i] * relrate->GetDim();
	}
	SetRelativeRate(rescaledrelrate);
	CopyStationary(stat->val().GetArray());
	setTsTv(randomTsTv->val());
	if (randomOmegaTs)	{
		setOmegaTs(randomOmegaTs->val());
	}
	else	{
		setOmegaTs(1);
	}
	setOmegaTv(randomOmegaTv->val());
}

///      AADoubleOmegaSubMatrix

void AADoubleOmegaSubMatrix::ComputeArray(int n)	{

	double total = 0;
	for (int i=0; i<GetNstate(); i++)	{
		if (i!= n){
			int pos = GetDifferingPosition(n,i);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,n);
				int b = GetCodonPosition(pos,i);
				if (a == b)	{
					cerr << GetCodonStateSpace()->GetState(n) << '\t' << GetCodonStateSpace()->GetState(i) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[n][i] = (*NucMatrix)(a,b);
				if (! Synonymous(n,i))	{
					Q[n][i] *= GetOmegadNdS();
				}
				if(mSimilarityMatrix->isRadical(statespace->Translation(n),statespace->Translation(i)))	{
					Q[n][i] *= GetOmegaKrKc();
				}
			}
			else	{
				Q[n][i] = 0;
			}
			total += Q[n][i];
		}
	}
	Q[n][n] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

void AADoubleOmegaSubMatrix::SetOmega(double inomegadNdS, double inomegaKrKc){
		omegadNdS = inomegadNdS;
		omegaKrKc = inomegaKrKc;
}

void RandomAADoubleOmegaSubMatrix::SetParameters(){
		SetNucMatrix(matrix);
		SetOmega(RandomOmegadNdS->val(), RandomOmegaKrKc->val());
}

///	AAregSubMatrix

void AAregSubMatrix :: ComputeArray(int n){
	double total = 0;
	for (int i=0; i<Naa; i++){
		if(n!=i){
			Q[n][i] = exp(GetSlope(n,i) * log(GetOmega()) + GetIntercept(n,i)) * mStationary[i];
			total += Q[n][i];
		}
	}
	Q[n][n] = -total;
}


void AAregSubMatrix :: CopyStationary(const double* instat) {
	for (int i=0; i< Naa ; i++){
		mStationary[i] = instat[i];
	}
}

void AAregSubMatrix :: CopySlopeIntercept(double* inslope, double* inintercept, int dimension){
	for (int i=0; i< dimension; i++){
		mSlope [i] = inslope[i];
		mIntercept [i] = inintercept[i];
	}
}

void AAregSubMatrix :: setOmega(double inomega) {
	omega = inomega;
}


void RandomAAregSubMatrix :: SetParameters(){
	setOmega(randomOmega->val());
	CopyStationary(stat->val().GetArray());
	CopySlopeIntercept(slope->val().GetArray(), intercept->val().GetArray(), slope->GetDim());
}


