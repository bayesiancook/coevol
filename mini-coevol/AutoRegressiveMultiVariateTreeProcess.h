#ifndef ARMULTIVARIATETREEPROCESS
#define ARMULTIVARIATETREEPROCESS

#include "MultiVariateTreeProcess.h"
#include "Chronogram.h"
#include "ContinuousData.h"


// Vector of reals draw from a Wishart distribution.
class AutoRegressiveMultiNormal : public MultiNormal	{

	public:

	AutoRegressiveMultiNormal(Var<CovMatrix>* insigma, Var<RealVector>* inmean, Var<PosReal>* inphi, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0) : MultiNormal(insigma,inup,intime)	{
		phi = inphi;
		mean = inmean;
		Register(phi);
		Register(mean);
		Sample();
	}


	~AutoRegressiveMultiNormal(){
	}

	void drawSample(){
		if(isRoot()){
			double mu = 1.0 / sqrt(2 * phi->val());
			double* t = new double[GetDim()];
			sigma->drawVal(t);
			for ( int i=0; i <GetDim(); i++){
				if(!ClampVector[i]){
					val()[i] = (*mean)[i] + t[i] * mu;
				}
			}
		}
		else{
			double* t = new double[GetDim()];
			sigma->drawVal(t);
			double expo = exp(-phi->val() * time->val());
			double mu = sqrt ((1 - exp(-2 * phi->val() * time->val())) / 2 / phi->val());
			for ( int i=0; i <GetDim(); i++){
				if(!ClampVector[i]){
					val()[i] = t[i] * mu + up->val()[i] * expo + (*mean)[i] * (1 - expo);
				}
			}
			delete[] t;
		}
	}

	double logProb() {
		if(isRoot()){
			//calcul de transposé de X * sigma * X
			double tXSX = 0;
			double* dval = new double[GetDim()];
			double mu = 1.0 / sqrt(2 * phi->val());
			for (int i=0; i<GetDim() ; i++) {
				dval[i] =  (val()[i] - (*mean)[i]) / mu;
			}
			for (int i=0; i<GetDim() ; i++) {
				tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
				for (int j=0; j<i ; j++) {
					tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
				}
			}
			double d = -0.5 * (sigma->GetLogDeterminant() + tXSX) -  GetDim() * log(mu);
			return d;
		}
		else{
			//calcul de transposé de X * sigma * X
			double tXSX = 0;
			double* dval = new double[GetDim()];
			double expo = exp(-phi->val() * time->val());
			double mu = sqrt ((1 - exp(-2 * phi->val() * time->val())) / 2 / phi->val());
			for (int i=0; i<GetDim() ; i++) {
				dval[i] =  (val()[i] - up->val()[i] * expo - (*mean)[i] * (1 - expo))/mu;
			}
			for (int i=0; i<GetDim() ; i++) {
				tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
				for (int j=0; j<i ; j++) {
					tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
				}
			}
			double d = -0.5 * (sigma->GetLogDeterminant() + tXSX) -  GetDim() * log(mu);
			delete[] dval;
			return d;
		}
	}

	protected:

	Var<PosReal>* phi;
	Var<RealVector>* mean;

};

class AutoRegressiveMultiVariateTreeProcess : public virtual MultiVariateTreeProcess	{

	protected:

	Var<PosReal>* phi;
	Var<RealVector>* mean;

	public:


	AutoRegressiveMultiVariateTreeProcess()	{}


	AutoRegressiveMultiVariateTreeProcess(Var<CovMatrix>* insigma, Var<RealVector>* inmean, Var<PosReal>* inphi, LengthTree* intree) : MultiVariateTreeProcess() {
		sigma = insigma;
		phi = inphi;
		mean = inmean;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~AutoRegressiveMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	AutoRegressiveMultiNormal* GetAutoRegressiveMultiNormal(const Link* link)	{
		AutoRegressiveMultiNormal* m = dynamic_cast<AutoRegressiveMultiNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			return new AutoRegressiveMultiNormal(sigma, mean, phi, GetAutoRegressiveMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()));
		}
		else{
 			return new AutoRegressiveMultiNormal(sigma, mean, phi);
		}
	}

};

/*
class SigmaPhiMove : public MCUpdate, public Mnode {

	Rvar<CovMatrix>* sigma;
	Rvar<PosReal>* phi;
	double tuning;

	public:

	SigmaPhiMove(Rvar<CovMatrix>* insigma, Rvar<PosReal>* inphi, double intuning)	{
		sigma = insigma;
		phi = inphi;
		tuning = intuning;
		Register(phi);
		Register(sigma);
	}

	double Move(double tuning_modulator){
		return tree->SegmentMove(tuning,imin,size);
	}
};
*/

#endif
