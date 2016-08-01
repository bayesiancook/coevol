#ifndef BWARMULTIVARIATETREEPROCESS
#define BWARMULTIVARIATETREEPROCESS

#include "MultiVariateTreeProcess.h"
#include "Chronogram.h"
#include "ContinuousData.h"


// Vector of reals draw from a Wishart distribution.
class BrownianAutoRegressiveMultiNormal : public MultiNormal	{

	public:

	BrownianAutoRegressiveMultiNormal(Var<CovMatrix>* insigma, Var<RealVector>* inmean, Var<PosReal>* inphi, Var<Real>* incovgamma, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0) : MultiNormal(insigma,inup,intime)	{

		phi = inphi;
		covgamma = incovgamma;
		mean = inmean;
		Register(phi);
		Register(covgamma);
		Register(mean);
		if (GetDim() != 2)	{
			cerr << "error in Brownian Autoregressive\n";
			exit(1);
		}
		v = new CovMatrix(GetDim());
		Sample();
	}


	~BrownianAutoRegressiveMultiNormal(){
		delete v;
	}

	void drawSample(){
		if(isRoot()){
			for ( int i=0; i <GetDim(); i++){
				if(!ClampVector[i]){
					val()[i] = 0;
				}
			}
		}
		else{
			double alpha = (*mean)[0];
			double beta = (*mean)[1];
			(*v)[0][0] = (*sigma)[0][0] * time->val();
			(*v)[0][1] = (*v)[1][0] = alpha * (*sigma)[0][0] * time->val();
			(*v)[1][1] = alpha * alpha * (*sigma)[0][0] * time->val() + (*sigma)[1][1] / 2 / phi->val() * (1 - exp(-2 * phi->val() * time->val()));
			v->CorruptDiag();
			v->Diagonalise();

			double* t = new double[GetDim()];
			v->drawVal(t);
			if (! ClampVector[0])	{
				(*this)[0] = (*up)[0] + t[0];
			}
			if (! ClampVector[1])	{
				double expo = exp(-phi->val() * time->val());
				(*this)[1] = beta + alpha * (*up)[0] * (1 - expo) + ((*up)[1] - beta) * expo + t[1];
			}
			delete[] t;
		}
	}

	double logProb() {
		if(isRoot()){
			// return 0;
			double gamma = covgamma->val();
			double alpha = (*mean)[0];
			double beta = (*mean)[1];
			double v = ((gamma - alpha) * (gamma - alpha) * (*sigma)[0][0] + (*sigma)[1][1]) / 2 / phi->val();
			// double v = (alpha * alpha * (*sigma)[0][0] + (*sigma)[1][1]) / 2 / phi->val();
			double tmp = (*this)[1] - beta - alpha * (*this)[0];
			return -0.5 * (log(v) + tmp * tmp / v);
		}
		else{
			double gamma = covgamma->val();
			double alpha = (*mean)[0];
			double beta = (*mean)[1];
			(*v)[0][0] = (*sigma)[0][0] * time->val();
			(*v)[0][1] = (*v)[1][0] = (*sigma)[0][0] * (alpha * time->val() - (alpha - gamma) * (1 - exp(-phi->val() * time->val())) / phi->val());
			(*v)[1][1] = (*sigma)[0][0] * (alpha * alpha *time->val() - 2 * alpha * (alpha - gamma) * (1 - exp(-phi->val() * time->val())) / phi->val() + (alpha - gamma) * (alpha - gamma) * (1 - exp(-2 * phi->val() * time->val())) / 2 / phi->val()) + (*sigma)[1][1] / 2 / phi->val() * (1 - exp(-2 * phi->val() * time->val()));
			v->CorruptDiag();
			v->Diagonalise();

			// double expo = 1;
			double expo = exp(-phi->val() * time->val());
			double tmp0 = (*this)[0] - (*up)[0];
			// double tmp1 = (*this)[1] - (*up)[1];
			double tmp1 = ((*this)[1] - beta) - (alpha * (*up)[0] * (1 - expo) + ((*up)[1] - beta) * expo);
			double** inv = v->GetInvMatrix();
			/*
			for (int i=0; i<2; i++)	{
				for (int j=0; j<2; j++)	{
					double tmp = 0;
					for (int k=0; k<2; k++)	{
						tmp += (*v)[i][k] * inv[k][j];
					}
					if (i==j)	{
						tmp -=1;
					}
					if (fabs(tmp) > 1e-7)	{
						cerr << "error: matrix inverse\n";
						cerr << tmp << '\n';
						cerr << '\n';
						for (int ii=0; ii<2; ii++)	{
							for (int jj=0; jj<2; jj++)	{
								cerr << (*v)[ii][jj] << '\t';
							}
							cerr << '\n';
						}
						cerr << '\n';
						for (int ii=0; ii<2; ii++)	{
							for (int jj=0; jj<2; jj++)	{
								cerr << inv[ii][jj] << '\t';
							}
							cerr << '\n';
						}
						cerr << '\n';
						// exit(1);
					}
				}
			}
			*/
			double tXSX = inv[0][0] * tmp0 * tmp0 + 2 * inv[0][1] * tmp0 * tmp1 + inv[1][1] * tmp1 * tmp1;
			double d = -0.5 * (v->GetLogDeterminant() + tXSX);
			return d;
		}
	}

	/*
	double logProb() {
		if(isRoot()){
			// return 0;
			double alpha = (*mean)[0];
			double beta = (*mean)[1];
			double v = (*sigma)[1][1] / 2 / phi->val();
			double tmp = (*this)[1] - beta - alpha * (*this)[0];
			return -0.5 * (log(v) + tmp * tmp / v);
		}
		else{
			double alpha = (*mean)[0];
			double beta = (*mean)[1];
			(*v)[0][0] = (*sigma)[0][0] * time->val();
			(*v)[0][1] = (*v)[1][0] = alpha * (*sigma)[0][0] * time->val();
			// (*v)[1][1] = alpha * alpha * (*sigma)[0][0] * time->val() + (*sigma)[1][1] * time->val();
			(*v)[1][1] = alpha * alpha * (*sigma)[0][0] * time->val() + (*sigma)[1][1] / 2 / phi->val() * (1 - exp(-2 * phi->val() * time->val()));
			v->CorruptDiag();
			v->Diagonalise();

			// double expo = 1;
			double expo = exp(-phi->val() * time->val());
			double tmp0 = (*this)[0] - (*up)[0];
			// double tmp1 = (*this)[1] - (*up)[1];
			double tmp1 = ((*this)[1] - beta) - (alpha * (*up)[0] * (1 - expo) + ((*up)[1] - beta) * expo);
			double** inv = v->GetInvMatrix();
			double tXSX = inv[0][0] * tmp0 * tmp0 + 2 * inv[0][1] * tmp0 * tmp1 + inv[1][1] * tmp1 * tmp1;
			double d = -0.5 * (v->GetLogDeterminant() + tXSX);
			return d;
		}
	}
	*/

	protected:

	Var<PosReal>* phi;
	Var<Real>* covgamma;
	Var<RealVector>* mean;
	CovMatrix* v;

};

class BrownianAutoRegressiveMultiVariateTreeProcess : public virtual MultiVariateTreeProcess	{

	protected:

	Var<PosReal>* phi;
	Var<Real>* covgamma;
	Var<RealVector>* mean;

	public:


	BrownianAutoRegressiveMultiVariateTreeProcess()	{}


	BrownianAutoRegressiveMultiVariateTreeProcess(Var<CovMatrix>* insigma, Var<RealVector>* inmean, Var<PosReal>* inphi, Var<Real>* incovgamma, LengthTree* intree) : MultiVariateTreeProcess() {
		sigma = insigma;
		phi = inphi;
		covgamma = incovgamma;
		mean = inmean;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~BrownianAutoRegressiveMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	BrownianAutoRegressiveMultiNormal* GetBrownianAutoRegressiveMultiNormal(const Link* link)	{
		BrownianAutoRegressiveMultiNormal* m = dynamic_cast<BrownianAutoRegressiveMultiNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			return new BrownianAutoRegressiveMultiNormal(sigma, mean, phi, covgamma, GetBrownianAutoRegressiveMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()));
		}
		else{
 			return new BrownianAutoRegressiveMultiNormal(sigma, mean, phi, covgamma);
		}
	}

};

#endif
