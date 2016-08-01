
#ifndef LINREG_H
#define LINREG_H

class IIDNormalArray : public IIDArray<RealVector>	{

	public:

	IIDNormalArray(int insize, int indim, Var<Real>* inmean, Var<PosReal>* invar, Var<PosReal>* invar0 = 0) : IIDArray<RealVector>(insize)	{
		mean = inmean;
		var = invar;
		dim = indim;
		if (invar0)	{
			var0 = invar0;
		}
		else	{
			var0 = var;
		}
		Create();
	}

	IIDNormal* GetIIDNormal(int i)	{
		return dynamic_cast<IIDNormal*>(GetVal(i));
	}

	int GetDim()	{
		return dim;
	}

	void SetAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetIIDNormal(i)->SetAtZero();
		}
	}

	void ClampAtZero()	{
		for (int i=0; i<GetSize(); i++)	{
			GetIIDNormal(i)->ClampAtZero();
		}
	}

	protected:

	Rvar<RealVector>* CreateVal(int mat)	{
		if (!mat)	{
			return new IIDNormal(dim,mean,var);
		}
		return new IIDNormal(dim,mean,var0);
	}

	private:

	int dim;
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<PosReal>* var0;
};

class LinRegNormal : public Rvar<RealVector>	{

	public:

	LinRegNormal(Var<PosReal>* intime, IIDNormalArray* inregcoef, DiagonalCovMatrix* insigma, Var<RealVector>* inup, Var<RealVector>* incovup, Var<RealVector>* incovdown, Var<PosReal>* inphi = 0)  {
		regcoef = inregcoef;
		sigma = insigma;
		up = inup;
		covup = incovup;
		covdown = incovdown;
		time = intime;
		phi = inphi;
		setval(RealVector(insigma->GetDim()));
		Register(sigma);
		if (time)	{
			Register(time);
		}
		regcoef->RegisterChild(this);
		if (up)	{
			Register(up);
		}
		if (covup)	{
			Register(covup);
		}
		if (covdown)	{
			Register(covdown);
		}
		if (phi)	{
			Register(phi);
			if (!covdown)	{
				cerr << "error : covdown should never be 0\n";
				exit(1);
			}
			if (regcoef->GetVal(0)->GetDim() != covdown->GetDim() + 1)	{
				cerr << "error : non matching dimension\n";
				exit(1);
			}
		}
		ClampVector = new bool[GetDim()];
		for ( int i=0; i <GetDim(); i++){
			UnClamp(i);
		}
		Sample();
	}

	~LinRegNormal()	{
		delete[] ClampVector;
	}

	bool isRoot()	{
		return (up == 0);
	}

	int GetCovariateDim()	{
		if (! covup)	{
			cerr << "error in get cov dim\n";
			exit(1);
		}
		return covup->GetDim();
	}

	/*
	void Clamp(int index){
		ClampVector[index] = true;
	}
	*/

	void ClampAt(double inval, int index){
		val()[index] = inval;
		ClampVector[index] = true;
	}

	void ClampAtZero(){
		for (int i=0; i<GetDim(); i++)	{
			val()[i] = 0;
			ClampVector[i] = true;
		}
	}

	void UnClamp(int index){
		ClampVector[index] = false;
	}

	double PiecewiseTranslation(double u, int index, int n)	{
		for ( int i=0; i <n; i++){
			if (i+index >= GetDim())	{
				cerr << "error in MultiNormal::PiecewiseTranslationMove : " << index << '\t' << n << '\t' << GetDim() << "\n";
				exit(1);
			}
			if ((! isClamped()) && (!ClampVector[index+i]))	{
				(*this)[index + i] += u;
			}
		}
		return 0;
	}

	double logProb()	{
		double total = 0;
		if (phi)	{
			if (isRoot())	{
				for (int i=0; i<GetDim(); i++)	{
					double tmp = 0;
					for (int j=0; j<GetCovariateDim(); j++)	{
						tmp += (*regcoef->GetVal(i))[j] * (*covdown)[j];
					}
					tmp += (*regcoef->GetVal(i))[GetCovariateDim()];
					double tt = (*sigma)[i][i] / 2 / phi->val();
					double v = (*this)[i] - tmp;
					total -= 0.5 * ( v * v / tt + log(tt));
				}
			}
			else	{
				double expo = exp(-phi->val() * time->val());
				for (int i=0; i<GetDim(); i++)	{
					double tmp = 0;
					for (int j=0; j<GetCovariateDim(); j++)	{
						tmp += (*regcoef->GetVal(i))[j] * ( (*covdown)[j] - expo * (*covup)[j] );
					}
					tmp += (*regcoef->GetVal(i))[GetCovariateDim()] * (1-expo);
					double tt = (1 - expo * expo) * (*sigma)[i][i] / 2 / phi->val();
					double v = (*this)[i] - expo * (*up)[i] - tmp;
					total -= 0.5 * ( v * v / tt + log(tt));
				}
			}
		}
		else	{
			if (! isRoot())	{
				for (int i=0; i<GetDim(); i++)	{
					double tmp = 0;
					for (int j=0; j<GetCovariateDim(); j++)	{
						tmp += (*regcoef->GetVal(i))[j] * ( (*covdown)[j] - (*covup)[j] );
					}
					double tt = (*sigma)[i][i] * time->val();
					double v = (*this)[i] - (*up)[i] - tmp;
					total -= 0.5 * ( v * v / tt + log(tt));
				}
			}
		}
		return total;
	}

	void drawSample()	{
		/*
		if (phi)	{
		}
		else	{
		*/
			if (isRoot())	{
				for (int i=0; i<GetDim(); i++)	{
					if(!ClampVector[i]){
						(*this)[i] = 0;
					}
				}
			}
			else	{
				for (int i=0; i<GetDim(); i++)	{
					if(!ClampVector[i]){
						double tmp = 0;
						for (int j=0; j<GetCovariateDim(); j++)	{
							tmp += (*regcoef->GetVal(i))[j] * ( (*covdown)[j] - (*covup)[j] );
						}
						double t = sqrt((*sigma)[i][i] * time->val());
						double u = Random::sNormal();
						double v = (*up)[i] + tmp + t * u;
						(*this)[i] = v;
					}
				}
			}
		// }
	}

	protected:

	bool* ClampVector;
	IIDNormalArray* regcoef;
	DiagonalCovMatrix* sigma;
	Var<PosReal>* time;
	Var<RealVector>* up;
	Var<RealVector>* covup;
	Var<RealVector>* covdown;
	Var<PosReal>* phi;

};

class LinRegNormalProcess : public MCMC, public NodeValPtrTree<Rvar<RealVector> > {

	public:

	LinRegNormalProcess() {
	}

	LinRegNormalProcess(LengthTree* inlengthtree, NodeVarTree<RealVector>* intree, IIDNormalArray* inregcoef, DiagonalCovMatrix* insigma, Var<PosReal>* inphi)	{
		lengthtree = inlengthtree;
		tree = intree;
		regcoef = inregcoef;
		sigma = insigma;
		phi = inphi;
		if (regcoef->GetVal(0)->GetDim() != tree->GetNodeVal(GetRoot()->GetNode())->GetDim())	{
			cerr << "error in linregnormalprocess: non matching dimension\n";
			cerr << regcoef->GetVal(0)->GetDim() << '\n';
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

	int GetCovariateDim()	{
		return regcoef->GetVal(0)->GetDim();
	}

	int GetDim()	{
		return regcoef->GetSize();
	}

	void drawSample()	{
		RecursivedrawSample(this->GetRoot());
	}

	double PiecewiseTranslation(double u, int index, int k){
		return RecursivePiecewiseTranslation(this->GetRoot(),u,index,k);
	}

	double RecursivePiecewiseTranslation(const Link* from, double u, int index, int k)	{
		double total = GetLinRegNormal(from)->PiecewiseTranslation(u, index, k);
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursivePiecewiseTranslation(link->Out(),u,index,k);
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

	double GetExpVal(const Link* link, int index)	{
		return exp((*GetNodeVal(link->GetNode()))[index]);
	}

	double GetVal(const Link* link, int index)	{
		return (*GetNodeVal(link->GetNode()))[index];
	}

	void CutOff(double cutoff, int index)	{
		RecursiveCutOff(GetRoot(),cutoff,index);
	}

	void RecursiveCutOff(const Link* from, double cutoff, int index)	{
		(*GetNodeVal(from->GetNode()))[index] = cutoff;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCutOff(link->Out(),cutoff,index);
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


	double GetMean(int index)	{
		int n = 0;
		double tmp = RecursiveGetMean(GetRoot(),index,n);
		return tmp / n;
	}

	double GetMeanExp(int index)	{
		int n = 0;
		double tmp = RecursiveGetMeanExp(GetRoot(),index,n);
		return tmp / n;
	}

	void Reset()	{
		RecursiveReset(GetRoot());
	}

	void Clamp()	{
		RecursiveClamp(GetRoot());
	}

	void Add(MultiVariateTreeProcess* inprocess)	{
		RecursiveAdd(inprocess,GetRoot());
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
		GetLinRegNormal(from)->SetAtZero();
	}

	void RecursiveClamp(const Link* from)	{
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveClamp(link->Out());
		}
		GetLinRegNormal(from)->Clamp();
	}

	void RecursiveAdd(MultiVariateTreeProcess* inprocess, const Link* from)	{
		cerr << "error in lin reg recursive add\n";
		exit(1);
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(inprocess,link->Out());
		}
		// GetLinRegNormal(from)->Add(*(inprocess->GetLinRegNormal(from)));
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

	double RecursiveGetMean(const Link* from, int index, int& tot)	{
		double total = (*GetLinRegNormal(from))[index];
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMean(link->Out(), index, tot);
		}
		return total;
	}

	double RecursiveGetMeanExp(const Link* from, int index, int& tot)	{
		double total = exp((*GetLinRegNormal(from))[index]);
		tot++;
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			total += RecursiveGetMeanExp(link->Out(), index, tot);
		}
		return total;
	}

	Rvar<RealVector>* CreateNodeVal(const Link* from)	{
		if (from->isRoot())	{
			return new LinRegNormal(0,regcoef,sigma,0,tree->GetNodeVal(from->GetNode()),tree->GetNodeVal(from->GetNode()),phi);
		}
		return new LinRegNormal(lengthtree->GetBranchVal(from->GetBranch()),regcoef,sigma,GetNodeVal(from->Out()->GetNode()),tree->GetNodeVal(from->Out()->GetNode()),tree->GetNodeVal(from->GetNode()),phi);
	}

	DiagonalCovMatrix* sigma;
	LengthTree* lengthtree;
	NodeVarTree<RealVector>* tree;
	IIDNormalArray* regcoef;
	Var<PosReal>* phi;
};

// partitions

class BidimIIDNormalArray 	: public MCMC {

	public:

	BidimIIDNormalArray(int inncomp, int insize, int indim, Var<Real>* inmean, Var<PosReal>* invar, Var<PosReal>* invar0 = 0)	{
		ncomp = inncomp;
		size = insize;
		dim = indim;
		mean = inmean;
		var = invar;
		if (invar0)	{
			var0 = invar0;
		}
		else	{
			var0 = var;
		}
		Create();
	}

	int GetDim()	{
		return dim;
	}

	int GetSize()	{
		return size;
	}

	int GetNcomp()	{
		return ncomp;
	}

	IIDNormalArray* GetIIDNormalArray(int i)	{
		return val[i];
	}

	IIDNormal* GetIIDNormal(int i, int j)	{
		return val[i]->GetIIDNormal(j);
	}

	double Move(double tuning)	{
		double tot = 0;
		for (int i=0; i<ncomp; i++)	{
			tot += val[i]->Move(tuning);
		}
		return tot / ncomp;
	}

	void drawSample()	{
		for (int i=0; i<ncomp; i++)	{
			val[i]->drawSample();
		}
	}

	void SetAtZero()	{
		for (int i=0; i<ncomp; i++)	{
			GetIIDNormalArray(i)->SetAtZero();
		}
	}

	void ClampAtZero()	{
		for (int i=0; i<ncomp; i++)	{
			val[i]->ClampAtZero();
		}
	}

	double GetLogProb()	{
		double tot = 0;
		for (int i=0; i<ncomp; i++)	{
			tot += val[i]->GetLogProb();
		}
		return tot;
	}

	void ToStream(ostream& os)	{
		for (int i=0; i<ncomp; i++)	{
			os << *(val[i]) << '\t';
		}
		os << '\n';
	}

	void FromStream(istream& is)	{
		for (int i=0; i<ncomp; i++)	{
			is >> *(val[i]);
		}
	}

	friend ostream& operator<<(ostream& os, BidimIIDNormalArray& a)	{
		a.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, BidimIIDNormalArray& a)  {
		a.FromStream(is);
		return is;
	}

	protected:

	void Create()	{
		val = new IIDNormalArray*[ncomp];
		val[0] = new IIDNormalArray(size,dim,mean,var);
		for (int i=1; i<ncomp; i++)	{
			val[i] = new IIDNormalArray(size,dim,mean,var0);
		}
	}

	IIDNormalArray** val;
	int ncomp;
	int size;
	int dim;
	Var<Real>* mean;
	Var<PosReal>* var;
	Var<PosReal>* var0;
};

class PartitionLinRegNormalProcess : public LinRegNormalProcess	{

	public:

	PartitionLinRegNormalProcess() {}

	PartitionLinRegNormalProcess(LengthTree* inlengthtree, NodeVarTree<RealVector>* intree, BidimIIDNormalArray* inregcoefarray, VarArray<CovMatrix>* insigmaarray, BranchPartition* inparttree, Var<PosReal>* inphi)	{
		sigmaarray = insigmaarray;
		regcoefarray = inregcoefarray;
		sigma = dynamic_cast<DiagonalCovMatrix*>(sigmaarray->GetVal(0));
		regcoef = regcoefarray->GetIIDNormalArray(0);
		tree = intree;
		lengthtree = inlengthtree;
		parttree = inparttree;
		phi = inphi;
		RecursiveCreate(GetRoot());
	}

	~PartitionLinRegNormalProcess(){
		RecursiveDelete(GetRoot());
	}

	int GetDim()	{
		return sigmaarray->GetVal(0)->GetDim();
	}

	protected:

	Rvar<RealVector>* CreateNodeVal(const Link* from){
		if(!from->isRoot()){
			int alloc = parttree->GetAlloc(from->GetBranch());
			Var<CovMatrix>* sigma = 0;
			if (sigmaarray->GetSize() == 1)	{
				sigma = sigmaarray->GetVal(0);
			}
			else	{
				sigma = sigmaarray->GetVal(alloc);
			}
			DiagonalCovMatrix* diagsigma = dynamic_cast<DiagonalCovMatrix*>(sigma);
			if (! diagsigma)	{
				cerr << "error in particionlinreg : null diag matrix\n";
				cerr << "sigma\n";
				cerr << alloc << '\t' << sigma << '\n';
				exit(1);
			}
			IIDNormalArray* regcoef = 0;
			if (regcoefarray->GetNcomp() == 1)	{
				regcoef = regcoefarray->GetIIDNormalArray(0);
			}
			else	{
				regcoef = regcoefarray->GetIIDNormalArray(alloc);
			}
			return new LinRegNormal(lengthtree->GetBranchVal(from->GetBranch()),regcoef,diagsigma,GetNodeVal(from->Out()->GetNode()),tree->GetNodeVal(from->Out()->GetNode()),tree->GetNodeVal(from->GetNode()),phi);
		}
		else	{
			Var<CovMatrix>* sigma = sigmaarray->GetVal(0);
			DiagonalCovMatrix* diagsigma = dynamic_cast<DiagonalCovMatrix*>(sigma);
			if (! diagsigma)	{
				cerr << "error in particionlinreg : null diag matrix\n";
				cerr << "sigma\n";
				cerr << "root"  << '\t' << sigma << '\n';
				exit(1);
			}
			IIDNormalArray* regcoef = regcoefarray->GetIIDNormalArray(0);
			return new LinRegNormal(0,regcoef,diagsigma,0,tree->GetNodeVal(from->GetNode()),tree->GetNodeVal(from->GetNode()),phi);
		}
	}

	BidimIIDNormalArray* regcoefarray;
	VarArray<CovMatrix>* sigmaarray;
	BranchPartition* parttree;
	Var<PosReal>* phi;
};

#endif

