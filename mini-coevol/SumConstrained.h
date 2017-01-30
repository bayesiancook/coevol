
#ifndef SUMCONS_H
#define SUMCONS_H

#include "RandomTypes.h"
#include "ValTree.h"

class SumConstrainedMapping	{

	public:

	SumConstrainedMapping(int indim)	{

		dim = indim;
		base = new double*[dim];
		for (int i=0; i<dim; i++)	{
			base[i] = new double[dim];
		}
		CreateBasis();
	}

	~SumConstrainedMapping()	{
		for (int i=0; i<dim; i++)	{
			delete[] base[i];
		}
		delete[] base;
	}

	int GetDim() {return dim;}

	void CreateBasis()	{

		int bkseed = Random::GetSeed();
		Random::InitRandom(101);

		for (int i=0; i<dim; i++)	{
			base[0][i] = 1.0 / sqrt(dim);
		}
		for (int i=1; i<dim; i++)	{
			for (int k=0; k<dim; k++)	{
				base[i][k] = Random::Uniform() - 0.5;
			}
			for (int j=0; j<i; j++)	{
				double scal = 0;
				for (int k=0; k<dim; k++)	{
					scal += base[i][k] * base[j][k];
				}
				for (int k=0; k<dim; k++)	{
					base[i][k] -= scal * base[j][k];
				}
				double norm = 0;
				for (int k=0; k<dim; k++)	{
					norm += base[i][k] * base[i][k];
				}
				norm = sqrt(norm);
				for (int k=0; k<dim; k++)	{
					base[i][k] /= norm;
				}
			}
		}

		// check ON:
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				double tmp = 0;
				for (int k=0; k<dim; k++)	{
					tmp += base[i][k] * base[j][k];
				}
				if (i==j)	{
					if (fabs(tmp-1) > 1e-8)	{
						cerr << "error: not norm\n";
						exit(1);
					}
				}
				else	{
					if (fabs(tmp) > 1e-8)	{
						cerr << "error: not ortho\n";
						exit(1);
					}
				}
			}
		}

		Random::InitRandom(bkseed);
	}

	void ToStream(ostream& os)	{
		os << dim << '\n';
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				os << base[i][j] << '\t';
			}
			os << '\n';
		}
	}

	void FromStream(istream& is)	{
		int tempdim;
		is >> tempdim;
		if (tempdim != dim)	{
			cerr << "error : non matching dimension in sum constrained mapping\n";
			exit(1);
		}
		for (int i=0; i<dim; i++)	{
			for (int j=0; j<dim; j++)	{
				is >> base[i][j];
			}
		}
	}

	int dim;
	double** base;

};

class SumConstrainedRealVector : public Dvar<RealVector>	{

	public:

	SumConstrainedRealVector(Var<RealVector>* insource, SumConstrainedMapping* inmapping)	{
		setval(RealVector(insource->GetDim() + 1));
		bkvalue = RealVector(insource->GetDim() + 1);
		source = insource;
		mapping = inmapping;
		temp = new double[GetDim()];
		Register(source);
		specialUpdate();
	}

	~SumConstrainedRealVector()	{
		delete[] temp;
	}

	protected:

	void specialUpdate()	{
		double** b = mapping->base;
		temp[0] = 0;
		for (int i=0; i<source->GetDim(); i++)	{
			temp[i+1] = (*source)[i];
		}
		// double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			double tmp = 0;
			for (int j=0; j<GetDim(); j++)	{
				tmp += b[j][i] * temp[j];
			}
			(*this)[i] = tmp;
			// total += tmp;
		}
		/*
		if (fabs(total) > 1e-8)	{
			cerr << "error in sum constrained special update\n";
			cerr << total << '\n';
			exit(1);
		}
		*/
	}

	Var<RealVector>* source;
	SumConstrainedMapping* mapping;
	double* temp;
};

class SumConstrainedProfile : public Dvar<Profile>	{

	public:

	SumConstrainedProfile(Var<RealVector>* inup, Var<RealVector>* indown, int inoffset, SumConstrainedMapping* inmapping)	{
		setval(Profile(inmapping->GetDim()));
		bkvalue = Profile(inmapping->GetDim());
		up = inup;
		down = indown;
		offset = inoffset;
		mapping = inmapping;
		tempup = new double[GetDim()];
		tempdown = new double[GetDim()];
		Register(up);
		Register(down);
		specialUpdate();
	}

	~SumConstrainedProfile()	{
		delete[] tempup;
		delete[] tempdown;
	}

	protected:

	void specialUpdate()	{
		double** b = mapping->base;
		tempup[0] = 0;
		tempdown[0] = 0;
		for (int i=0; i<GetDim()-1; i++)	{
			if (i+offset > up->GetDim())	{
				cerr << "error in sum constrained : overflow\n";
				exit(1);
			}
			tempup[i+1] = (*up)[i+offset];
			tempdown[i+1] = (*down)[i+offset];
		}
		double total = 0;
		double totup = 0;
		double totdown = 0;
		for (int i=0; i<GetDim(); i++)	{
			double tmpup = 0;
			double tmpdown = 0;
			for (int j=0; j<GetDim(); j++)	{
				tmpup += b[j][i] * tempup[j];
				tmpdown += b[j][i] * tempdown[j];
			}
			totup += tmpup;
			totdown += tmpdown;
			(*this)[i] = exp(tmpup) + exp(tmpdown);
			total += (*this)[i];
		}
		if (fabs(totup) > 1e-6)	{
			cerr << "error : non matching sum\n";
			cerr << "up : " << totup << '\n';
		}
		if (fabs(totdown) > 1e-6)	{
			cerr << "error : non matching sum\n";
			cerr << "down : " << totdown << '\n';
		}
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] /= total;
		}
	}

	Var<RealVector>* up;
	Var<RealVector>* down;
	int offset;
	SumConstrainedMapping* mapping;
	double* tempup;
	double* tempdown;
};

class SumConstrainedStatTree : public BranchValPtrTree<Dvar<Profile> >	{

	public:

	SumConstrainedStatTree(NodeVarTree<RealVector>* inprocess, int inoffset, SumConstrainedMapping* inmapping) {
		process = inprocess;
		offset = inoffset;
		mapping = inmapping;
		SetWithRoot(false);
		RecursiveCreate(GetRoot());
		cerr << "in sum cons tree : " << offset << '\t' << mapping->GetDim() << '\n';
	}

	~SumConstrainedStatTree()	{
		RecursiveDelete(GetRoot());
	}

	Tree* GetTree() {return process->GetTree();}

	void specialUpdate()	{
		RecursiveSpecialUpdate(GetRoot());
	}

	double GetMeanEntropy()	{
		int n = 0;
		double total = GetTotalEntropy(GetRoot(),n);
		return total / n;
	}

	protected:

	void RecursiveSpecialUpdate(const Link* from)	{
		if (! from->isRoot())	{
			GetBranchVal(from->GetBranch())->specialUpdate();
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSpecialUpdate(link->Out());
		}
	}

	double GetTotalEntropy(const Link* from, int& n)	{
		double total = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalEntropy(link->Out(),n);
		}
		if (! from->isRoot())	{
			total += GetBranchVal(from->GetBranch())->GetEntropy();
			n++;
		}
		return total;
	}

	Dvar<Profile>* CreateBranchVal(const Link* link)	{
		return new SumConstrainedProfile(process->GetNodeVal(link->GetNode()), process->GetNodeVal(link->Out()->GetNode()), offset, mapping);
	}

	NodeVarTree<RealVector>* process;
	int offset;
	SumConstrainedMapping* mapping;

};

#endif
