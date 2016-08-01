
#include "MultiVariateTreeProcess.h"

class PreciseMultiNormal : public MultiNormal	{

	public:


	PreciseMultiNormal(Var<CovMatrix>* insigma, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0, int inresolution = 1000) : MultiNormal(insigma,inup, intime), resolution(inresolution) {

		integral = new double[GetDim()];

	}

	~PreciseMultiVariateNormal()	{
		delete[] integral;
	}

	void drawSample()	{

		for (int i=0; i<GetDim(); i++)	{
			integral[i] = 0;
		}

		if(isRoot()){
			for ( int i=0; i <GetDim(); i++){
				val()[i] = 0;
			}
		}
		else{
			for ( int i=0; i <GetDim(); i++){
				val()[i] = up->val()[i];
			}
			double* t = new double[GetDim()];
			double sdt = sqrt(time->val() / resolution);
			for (int k=0; k<resolution; k++)	{
				sigma->drawVal(t);
				for ( int i=0; i <GetDim(); i++){
					integral[i] += (val()[i] + 0.5 * t[i] * sdt) * sdt;
					val()[i] += t[i] * sdt;
				}
			}
			delete[] t;
		}
	}

	double GetIntegratedValue(int index)	{
		return integral(index);
	}

	private:

	int resolution;
	double* integral;

};


class PreciseMultiVariateTreeProcess : public virtual MultiVariateTreeProcess	{

	public:

	PreciseMultiVariateTreeProcess() {}

	PreciseMultiVariateTreeProcess(ConjugateInverseWishart* insigma, LengthTree* intree, int inresolution){
		sigma = insigma;
		tree = intree;
		resolution = inresolution;
		RecursiveCreate(GetRoot());
	}

	~PreciseMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	PreciseMultiNormal* GetPreciseMultiNormal(const Link* link)	{
		PreciseMultiNormal* m = dynamic_cast<PreciseMultiNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			return new PreciseMultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()),resolution);
		}
		else{
 			return new PreciseMultiNormal(sigma,0,0,resolution);
		}
	}

	int resolution;
};

class PreciseMeanExpFromMultiVariate : public MeanExpFromMultiVariate	{

	public:

	PreciseMeanExpFromMultiVariate(MultiNormal* inup, MultiNormal* indown, Var<PosReal>* intime, int inindex, BranchValType inbval) : MeanExpFromMultiVariate(inup,indownn,intime,inindex,inbval,false) {}

	void specialUpdate(){
		double tmp = indown->GetIntegratedValue(index);
		if (bval == MEAN)	{
			tmp /= intime->val();
		}
		setval(tmp);
	}

};

class PreciseMeanExpTreeFromMultiVariate : public MeanExpTreeFromMultiVariate	{

	public:

	PreciseMeanExpTreeFromMultiVariate(MultiVariateTreeProcess* inprocess, int inindex, BranchValType inbval, bool withroot)	{
		process = inprocess;
		index = inindex;
		bval = inbval;
		meanexp = false;
		SetWithRoot(withroot);
		if ((bval == INTEGRAL) && (WithRoot()))	{
			cerr << "error in mean exp tree : integral is always without root\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~PreciseMeanExpTreeFromMultiVariate()	{
		// RecursiveDelete(GetRoot());
	}

	protected:

	Dvar<PosReal>* CreateBranchVal(const Link* link)	{
		if (link->isRoot())	{
			if (bval == INTEGRAL)	{
				cerr << "error in MeanExpFromMult::CreateBranchVal\n";
				exit(1);
			}
			return new PreciseMeanExpFromMultiVariate(process->GetPreciseMultiNormal(link), process->GetPreciseMultiNormal(link), 0, index, bval);
		}
		return new PreciseMeanExpFromMultiVariate(process->GetPreciseMultiNormal(link), process->GetPreciseMultiNormal(link->Out()), process->GetLengthTree()->GetBranchVal(link->GetBranch()), index, bval);
	}
};

