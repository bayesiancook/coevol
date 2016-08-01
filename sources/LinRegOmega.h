
#ifndef LINREG_H
#define LINREG_H

#include "ValArray.h"
#include "ValTree.h"

template<class T>  class BidimArray : public MCMC {

	public:

	BidimArray(int inN, int inP)	{
		N = inN;
		P = inP;
		array = 0;
	}

	int GetN()	{
		return N;
	}

	int GetP()	{
		return P;
	}

	virtual ~BidimArray()	{
		Delete();
	}

	void RegisterArray(DAGnode* innode)	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				GetCell(i,j)->Register(innode);
			}
		}
	}

	void RegisterArray(DAGnode* innode, int k)	{
		for (int i=0; i<N; i++)	{
			GetCell(i,k)->Register(innode);
		}
	}

	void SetAtRandom(double mean, double delta, int k)	{
		for (int i=0; i<GetN(); i++)	{
			GetCell(i,k)->setval(mean + delta * (Random::Uniform() - 0.5));
		}
	}

	double* GetVals(double* ptr)	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				(*ptr++) = GetCell(i,j)->val();
			}
		}
		return ptr;
	}

	double* SetVals(double* ptr)	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				GetCell(i,j)->setval(*ptr++);
			}
		}
		return ptr;
	}

	void ClampAt(double** inval)	{
		for (int i=0; i<GetN(); i++)	{
			for (int j=0; j<GetP(); j++)	{
				GetCell(i,j)->ClampAt(inval[i][j]);
			}
		}
	}

	void Create()	{
		array = new Rvar<T>**[N];
		for (int i=0; i<N; i++)	{
			array[i] = new Rvar<T>*[P];
			for (int j=0; j<P; j++)	{
				array[i][j] = CreateCell(i,j);
			}
		}
	}

	void Delete()	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				delete array[i][j];
			}
			delete[] array[i];
		}
		delete[] array;
		array = 0;
	}

	Rvar<T>* GetCell(int i, int j)	{
		return array[i][j];
	}

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				total += array[i][j]->GetLogProb();
			}
		}
		return total;
	}

	double Move(double tuning)	{
		double total = 0;
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				total += array[i][j]->Move(tuning);
			}
		}
		return total / N/P;
	}

	double ArrayMove(int* alloc, int myid, double tuning)	{
		double acc = 0;
		double tot = 0;
		for (int i=0; i<N; i++)	{
			if (alloc[i] == myid)	{
				for (int j=0; j<P; j++)	{
					acc += array[i][j]->Move(tuning);
					tot ++;
				}
			}
		}
		return acc / tot;
	}

	void drawSample()	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				array[i][j]->Sample();
			}
		}
	}

	void ToStream(ostream& os)	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				os << *array[i][j] << '\t';
			}
			os << '\n';
		}
	}

	void FromStream(istream& is)	{
		for (int i=0; i<N; i++)	{
			for (int j=0; j<P; j++)	{
				is >> *array[i][j];
			}
		}
	}

	double GetMean(int j)	{
		double mean = 0;
		for (int i=0; i<N; i++)	{
			double tmp = array[i][j]->val();
			mean += tmp;
		}
		mean /= N;
		return mean;
	}

	double GetVar(int j)	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<N; i++)	{
			double tmp = array[i][j]->val();
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= N;
		var /= N;
		var -= mean*mean;
		return var;
	}

	void RegisterChild(DAGnode* innode, int i)	{
		for (int j=0; j<P; j++)	{
			innode->Register(GetCell(i,j));
		}
	}

	void RegisterChildHorizontal(DAGnode* innode, int j)	{
		for (int i=0; i<N; i++)	{
			innode->Register(GetCell(i,j));
		}
	}

	void SetAt(double d, int j)	{
		for (int i=0; i<N; i++)	{
			GetCell(i,j)->setval(d);
		}
	}

	void Clamp(int j)	{
		for (int i=0; i<N; i++)	{
			GetCell(i,j)->Clamp();
		}
	}

	protected:

	virtual Rvar<T>* CreateCell(int i, int j) = 0;

	int N;
	int P;
	Rvar<T>*** array;

};

// first simple version of bidimarray
// only one parameter : lambda
// or a PosRealVector ?

/*
class PrecNormal : public virtual Rvar<Real>	{

	public:

	PrecNormal(Var<Real>* inmean, Var<PosReal>* inprec)	{

		mean = inmean;
		prec = inprec;
		Register(mean);
		Register(prec);
		Sample();
	}

	double logProb()	{
		return -0.5 * log(2 * Pi / prec->val()) -0.5 * (this->val() - mean->val()) * (this->val() - mean->val()) * prec->val();
	}

	protected:

	void drawSample()	{
		setval((Random::sNormal() + mean->val()) / sqrt(prec->val()));
	}

	private:

	Var<Real>* mean;
	Var<PosReal>* prec;

};

class BidimIIDNormal : public BidimArray<Real>	{

	public:

	BidimIIDNormal(int inN, int inP, Var<Real>** inmean, Var<PosReal>** inprec) : BidimArray<Real>(inN,inP)	{
		mean = inmean;
		prec = inprec;
		Create();
	}

	~BidimIIDNormal()	{
		Delete();
	}

	protected:

	Rvar<Real>* CreateCell(int i, int j)	{
		return new PrecNormal(mean[j],prec[j]);
	}

	Var<Real>** mean;
	Var<PosReal>** prec;

};
*/

class BidimIIDNormal : public BidimArray<Real>	{

	public:

	BidimIIDNormal(int inN, int inP, Var<Real>** inmean, Var<PosReal>** invar) : BidimArray<Real>(inN,inP)	{
		mean = inmean;
		var = invar;
		Create();
	}

	~BidimIIDNormal()	{
		Delete();
	}

	/*
	void SetAtRandom(double mean, double delta, int k)	{
		for (int i=0; i<GetN(); i++)	{
			GetCell(i,k)->setval(mean + delta * (Random::Uniform() - 0.5));
		}
	}
	*/

	protected:

	Rvar<Real>* CreateCell(int i, int j)	{
		return new Normal(mean[j],var[j]);
	}

	Var<Real>** mean;
	Var<PosReal>** var;

};

class BidimIIDUniform: public BidimArray<Real>	{

	public:

	BidimIIDUniform(int inN, int inP, double inmin, double inmax, DAGnode* inroot) : BidimArray<Real>(inN,inP)	{
		root = inroot;
		min = inmin;
		max = inmax;
		Create();
	}

	~BidimIIDUniform()	{
		Delete();
	}

	void SetAtRandom(double mean, double delta, int k)	{
		for (int i=0; i<GetN(); i++)	{
			GetCell(i,k)->setval(mean + delta * (Random::Uniform() - 0.5));
		}
	}

	protected:

	Rvar<Real>* CreateCell(int i, int j)	{
		return new Uniform(min,max,root);
	}

	double min;
	double max;
	DAGnode* root;
};

class Bernouilli : public Rvar<Int>	{

	public:

	Bernouilli(Var<UnitReal>* intheta)	{
		theta = intheta;
		Register(theta);
		Sample();
	}

	protected:

	void drawSample()	{
		if (Random::Uniform() < theta->val())	{
			setval(1);
		}
		else	{
			setval(0);
		}
	}

	double logProb()	{
		return (val() == 1) ? log(theta->val()) : log(1-theta->val());
	}

	double ProposeMove(double tuning)	{
		if (val() == 1)	{
			setval(0);
		}
		else	{
			setval(1);
		}
		return 0;
	}

	Var<UnitReal>* theta;
};

class BidimIIDBernouilli : public BidimArray<Int>	{

	public:

	BidimIIDBernouilli(int inN, int inP, Beta** intheta) : BidimArray<Int>(inN,inP)	{
		theta = intheta;
		Create();
	}

	~BidimIIDBernouilli()	{
		Delete();
	}

	protected:

	Rvar<Int>* CreateCell(int i, int j)	{
		return new Bernouilli(theta[j]);
	}

	Beta** theta;

};


class LinRegNormal : public Rvar<Real>	{

	public:

	LinRegNormal(Var<PosReal>* intime, BidimArray<Real>* inregcoef, BidimArray<Int>* intoggle, int ingene, Var<Real>* inup, Var<RealVector>* incovup, Var<RealVector>* incovdown, Var<PosReal>* intau)  {
		regcoef = inregcoef;
		toggle = intoggle;
		gene = ingene;
		up = inup;
		covup = incovup;
		covdown = incovdown;
		time = intime;
		tau = intau;
		if (time)	{
			Register(time);
		}
		Register(tau);
		regcoef->RegisterChild(this,gene);
		if (toggle)	{
			toggle->RegisterChild(this,gene);
		}
		if (up)	{
			Register(up);
		}
		if (covup)	{
			Register(covup);
		}
		if (covdown)	{
			Register(covdown);
		}
		Sample();
	}

	~LinRegNormal()	{
	}

	bool isRoot()	{
		return (up == 0);
	}

	int GetCovariateDim()	{
		if ((! covup) || (! covdown))	{
			cerr << "error in get cov dim\n";
			exit(1);
		}
		return covup->GetDim();
	}

	void ClampAtZero(){
		setval(0);
		Clamp();
	}

	double Translation(double u)	{
		if (! isClamped())	{
			setval(val() + u);
		}
		return 0;
	}

	double logProb()	{
		double total = 0;
		if (! isRoot())	{
			double tmp = 0;
			if (toggle)	{
				for (int j=0; j<GetCovariateDim(); j++)	{
					tmp += toggle->GetCell(gene,j)->val() * regcoef->GetCell(gene,j)->val() * ( (*covdown)[j] - (*covup)[j] );
				}
			}
			else	{
				for (int j=0; j<GetCovariateDim(); j++)	{
					tmp += regcoef->GetCell(gene,j)->val() * ( (*covdown)[j] - (*covup)[j] );
				}
			}
			double tt = time->val() / tau->val();
			double v = val() - up->val() - tmp;
				total -= 0.5 * ( v * v / tt + log(tt));
		}
		return total;
	}

	void drawSample()	{
		if(! isClamped())	{
			if (up)	{
				double tmp = 0;
				if (toggle)	{
					for (int j=0; j<GetCovariateDim(); j++)	{
						tmp += toggle->GetCell(gene,j)->val() * regcoef->GetCell(gene,j)->val() * ( (*covdown)[j] - (*covup)[j] );
					}
				}
				else	{
					for (int j=0; j<GetCovariateDim(); j++)	{
						tmp += regcoef->GetCell(gene,j)->val() * ( (*covdown)[j] - (*covup)[j] );
					}
				}
				double tt = time->val() / tau->val();
				double u = Random::sNormal();
				double v = up->val() + tmp + sqrt(tt) * u;
				setval(v);
			}
			else	{
				setval(0);
			}
		}
	}

	protected:

	BidimArray<Real>* regcoef;
	BidimArray<Int>* toggle;
	int gene;
	Var<PosReal>* tau;
	Var<PosReal>* time;
	Var<Real>* up;
	Var<RealVector>* covup;
	Var<RealVector>* covdown;
};

class LinRegNormalProcess : public MCMC, public NodeValPtrTree<Rvar<Real> > {
// class LinRegNormalProcess : public MCMC, public Additive, public NodeValPtrTree<Rvar<Real> > {

	public:

	LinRegNormalProcess() {
	}

	LinRegNormalProcess(LengthTree* inlengthtree, NodeVarTree<RealVector>* intree, BidimArray<Real>* inregcoef, BidimArray<Int>* intoggle, int ingene, Var<PosReal>* intau)	{
		lengthtree = inlengthtree;
		tree = intree;
		regcoef = inregcoef;
		toggle = intoggle;
		gene = ingene;
		tau = intau;
		if (regcoef->GetP() != tree->GetNodeVal(GetRoot()->GetNode())->GetDim())	{
			cerr << "error in linregnormalprocess: non matching dimension\n";
			cerr << regcoef->GetP() << '\n';
			cerr << tree->GetNodeVal(GetRoot()->GetNode())->GetDim() << '\n';
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~LinRegNormalProcess()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree(){
		return tree->GetTree();
	}

	LengthTree* GetLengthTree(){
		return lengthtree;
	}

	int GetDim()	{
		return regcoef->GetP();
	}

	void drawSample()	{
		RecursivedrawSample(this->GetRoot());
	}

	double* GetNodeVals(double* ptr)	{
		return GetNodeVals(GetRoot(),ptr);
	}

	double* GetNodeVals(const Link* from, double* ptr)	{
		*ptr = GetNodeVal(from->GetNode())->val();
		ptr++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ptr = GetNodeVals(link->Out(), ptr);
		}
		return ptr;
	}

	double* SetNodeVals(double* ptr)	{
		return SetNodeVals(GetRoot(),ptr);
	}

	double* SetNodeVals(const Link* from, double* ptr)	{
		GetNodeVal(from->GetNode())->setval(*ptr);
		ptr++;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ptr = SetNodeVals(link->Out(), ptr);
		}
		return ptr;
	}

	/*
	int ScalarAddition(double d)	{
		Translation(d);
		return 0;
	}
	*/

	double Translation(double u){
		return RecursiveTranslation(this->GetRoot(),u);
	}

	double RecursiveTranslation(const Link* from, double u)	{
		double total = GetLinRegNormal(from)->Translation(u);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveTranslation(link->Out(),u);
		}
		return total;
	}

	LinRegNormal* GetLinRegNormal(const Link* link)	{
		LinRegNormal* m = dynamic_cast<LinRegNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	LinRegNormal* GetLinRegNormal(const Node* node)	{
		LinRegNormal* m = dynamic_cast<LinRegNormal*> (GetNodeVal(node));
		return m;
	}

	double GetExpVal(const Link* link)	{
		return exp(GetNodeVal(link->GetNode())->val());
	}

	double GetVal(const Link* link)	{
		return GetNodeVal(link->GetNode())->val();
	}

	void CutOff(double cutoff)	{
		RecursiveCutOff(GetRoot(),cutoff);
	}

	void RecursiveCutOff(const Link* from, double cutoff)	{
		if (GetNodeVal(from->GetNode())->val() > cutoff)	{
			GetNodeVal(from->GetNode())->setval(cutoff);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCutOff(link->Out(),cutoff);
		}
	}

	void RecursiveRegister(DAGnode* node, const Link* from)	{
		GetLinRegNormal(from)->Register(node);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(node,link->Out());
		}
	}

	void RecursiveDeregister(DAGnode* node, const Link* from)	{
		GetLinRegNormal(from)->DeregisterFrom(node);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDeregister(node,link->Out());
		}
	}

	double Move(double tuning){
		int n = 0;
		double tot = RecursiveMove(this->GetRoot(),tuning,n);
		return tot / n;
	}

	double GetLogProb()	{
		return RecursiveGetLogProb(this->GetRoot());
	}


	double GetMean()	{
		int n = 0;
		double tmp = RecursiveGetMean(GetRoot(),n);
		return tmp / n;
	}

	double GetMeanExp()	{
		int n = 0;
		double tmp = RecursiveGetMeanExp(GetRoot(),n);
		return tmp / n;
	}

	void Reset()	{
		RecursiveReset(GetRoot());
	}

	void Clamp()	{
		RecursiveClamp(GetRoot());
	}

	void ClampRoot(double d)	{
		GetLinRegNormal(GetRoot())->ClampAt(d);
	}

	protected:

	void RecursivedrawSample(const Link* from)	{
		GetLinRegNormal(from)->Sample();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursivedrawSample(link->Out());
		}
	}

	void RecursiveReset(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
		GetLinRegNormal(from)->setval(0);
	}

	void RecursiveClamp(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveClamp(link->Out());
		}
		GetLinRegNormal(from)->Clamp();
	}

	double RecursiveMove(const Link* from, double tuning, int& count)	{
		double total = GetLinRegNormal(from)->Move(tuning);
		count++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveMove(link->Out(),tuning,count);
		}
		return total;
	}

	double RecursiveGetLogProb(const Link* from)	{
		double total = GetLinRegNormal(from)->GetLogProb();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetLogProb(link->Out());
		}
		return total;
	}

	double RecursiveGetMean(const Link* from, int& tot)	{
		double total = GetLinRegNormal(from)->val();
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMean(link->Out(), tot);
		}
		return total;
	}

	double RecursiveGetMeanExp(const Link* from, int& tot)	{
		double total = exp(GetLinRegNormal(from)->val());
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMeanExp(link->Out(), tot);
		}
		return total;
	}

	Rvar<Real>* CreateNodeVal(const Link* from)	{
		if (from->isRoot())	{
			return new LinRegNormal(0,regcoef,toggle,gene,0,tree->GetNodeVal(from->GetNode()),tree->GetNodeVal(from->GetNode()),tau);
		}
		return new LinRegNormal(lengthtree->GetBranchVal(from->GetBranch()),regcoef,toggle,gene,GetNodeVal(from->Out()->GetNode()),tree->GetNodeVal(from->Out()->GetNode()),tree->GetNodeVal(from->GetNode()),tau);
	}

	LengthTree* lengthtree;
	NodeVarTree<RealVector>* tree;
	BidimArray<Real>* regcoef;
	BidimArray<Int>* toggle;
	int gene;
	Var<PosReal>* tau;
};

#endif

