#ifndef TIMELINEMULTIVARIATETREEPROCESS_H
#define TIMELINEMULTIVARIATETREEPROCESS_H

#include "ConjugateMultiVariateTreeProcess.h"
#include "Chronogram.h"
#include "ContinuousData.h"
#include "TimeLine.h"

// Vector of reals draw from a Wishart distribution.
class TimeLineMultiNormal : public MultiNormal	{

	public:

	TimeLineMultiNormal(Var<CovMatrix>* insigma, Var<RealVector>* inglob1 = 0, Var<RealVector>* inglob2 = 0, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0) : MultiNormal(insigma,inup,intime)	{
		glob1 = inglob1;
		glob2 = inglob2;
		if (glob1)	{
			Register(glob1);
		}
		if (glob2)	{
			Register(glob2);
		}
		Sample();
	}


	~TimeLineMultiNormal(){
	}

	void drawSample(){
		if(isRoot()){
			for ( int i=0; i <GetDim(); i++){
				val()[i] = 0;
				// val()[i] = 2*rootmax * (Random::Uniform() - 0.5);
			}
		}
		else{
			double* t = new double[GetDim()];
			sigma->drawVal(t);
			for ( int i=0; i <GetDim(); i++){
				if(!ClampVector[i]){
					double tt = time->val();
					if (scale)	{
						tt *= scale->val();
					}
					val()[i] = t[i] * sqrt(tt) + (*glob2)[i] - (*glob1)[i] + up->val()[i];
				}
			}
			delete[] t;
		}
	}

	double logProb() {
		if(isRoot()){
			double total = 0;
			for (int i=0; i<GetDim(); i++)	{
				if (fabs(val()[i]) < rootmax)	{
					total -= log(2 * rootmax);
				}
				else	{
					total -= 100;
				}
			}
			return total;
		}
		else{
			//calcul de transposé de X * sigma * X
			double tXSX = 0;
			double* dval = new double[GetDim()];
			double tt = scale ? (time->val() * scale->val()) : ((double) time->val());
			double roott = sqrt(tt);
			for (int i=0; i<GetDim() ; i++) {
				dval[i] =  (val()[i] - up->val()[i] - (*glob2)[i] + (*glob1)[i])/roott;
			}
			for (int i=0; i<GetDim() ; i++) {
				tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
				for (int j=0; j<i ; j++) {
					tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
				}
			}
			double d = -0.5 * (sigma->GetLogDeterminant() + tXSX + GetDim() * 2 * log(roott));
			delete[] dval;
			return d;
		}
	}

	protected:

	Var<RealVector>* glob1;
	Var<RealVector>* glob2;

};

class TimeLineMultiVariateTreeProcess : public virtual MultiVariateTreeProcess	{

	protected:

	TimeLine* timeline;

	public:


	TimeLineMultiVariateTreeProcess()	{}


	TimeLineMultiVariateTreeProcess(Var<CovMatrix>* insigma, TimeLine* intimeline, Chronogram* intree) : MultiVariateTreeProcess() {
		sigma = insigma;
		timeline = intimeline;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~TimeLineMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	TimeLineMultiNormal* GetTimeLineMultiNormal(const Link* link)	{
		TimeLineMultiNormal* m = dynamic_cast<TimeLineMultiNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			return new TimeLineMultiNormal(sigma, timeline->GetValFromLink(link->Out()), timeline->GetValFromLink(link), GetTimeLineMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()));
		}
		else{
 			return new TimeLineMultiNormal(sigma);
		}
	}
};

class ConjugateTimeLineMultiNormal : public ConjugateSampling<RealVector>, public TimeLineMultiNormal	{

	public:

	ConjugateTimeLineMultiNormal(MultiNormalSemiConjugate* insigma, Var<RealVector>* inglob1 = 0, Var<RealVector>* inglob2 = 0, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0) :
		Rvar<RealVector>(),
		ConjugateSampling<RealVector>(),
		TimeLineMultiNormal(insigma, inglob1, inglob2, inup, intime)
	{
		conjugate_up.insert(insigma);
		contrast = 0;
		if (! intime)	{
			if (inup)	{
				cerr << "error in con time line mn\n";
				exit(1);
			}
			isroot = true;
		}
		else	{
			isroot = false;
		}
	}

	~ConjugateTimeLineMultiNormal() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		if (up)	{
			MultiNormalSemiConjugate* prior = dynamic_cast<MultiNormalSemiConjugate*>(parent);
			if (! prior)	{
				cout << "cast error in ConjugateMultiNormal::AddSuffStat\n";
				exit(1);
			}
			prior->AddToShape();
			ComputeContrast();
			prior->AddToScale(contrast);
		}
	}

	void ComputeContrast()	{
		if (! contrast)	{
			contrast = new double[GetDim()];
		}
		double tt = scale ? (time->val() * scale->val()) : ((double) time->val());
		double roott = sqrt(tt);
		for (int i=0; i<GetDim(); i++)	{
			contrast[i] = ((*this)[i] - (*glob2)[i] + (*glob1)[i] - (*up)[i]) / roott;
		}
	}

	protected:
	bool isroot;

	private:
	double* contrast;

};

class ConjugateTimeLineMultiVariateTreeProcess : public TimeLineMultiVariateTreeProcess, ConjugateMultiVariateTreeProcess	{

	public:

	ConjugateTimeLineMultiVariateTreeProcess(ConjugateInverseWishart* insigma, TimeLine* intimeline, LengthTree* intree){
		sigma = insigma;
		timeline = intimeline;
		TimeLineMultiVariateTreeProcess::sigma = insigma;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~ConjugateTimeLineMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			if (! tree->GetBranchVal(link->GetBranch()))	{
				cerr << "error in conjugate time line mv tree process\n";
				exit(1);
			}
			return new ConjugateTimeLineMultiNormal(sigma, timeline->GetValFromLink(link->Out()), timeline->GetValFromLink(link), GetTimeLineMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()));
		}
		else{
 			return new ConjugateTimeLineMultiNormal(sigma);
		}
	}

};

#endif
