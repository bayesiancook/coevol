#include "ValTree.h"
#include "Var.h"
#include "MultiVariateTreeProcess.h"

class MeanICTree	{

	public:

	MeanICTree(NodeVarTree<RealVector>* inprocess, LengthTree* inlengthtree)	{
		process = inprocess;
		lengthtree = inlengthtree;
		dim = inprocess->GetNodeVal(GetRoot()->GetNode())->GetDim();
		RecursiveCreate(GetRoot());
		gmean = new double[dim];
		gpp = new double[dim];
	}

	MeanICTree(MultiVariateTreeProcess* inprocess)	{
		process = inprocess;
		lengthtree = inprocess->GetLengthTree();
		dim = inprocess->GetNodeVal(GetRoot()->GetNode())->GetDim();
		RecursiveCreate(GetRoot());
		gmean = new double[dim];
		gpp = new double[dim];
	}

	LengthTree* GetLengthTree()	{
		return lengthtree;
	}

	~MeanICTree()	{
		RecursiveDelete(GetRoot());
		delete[] gmean;
		delete[] gpp;
	}

	Tree* GetTree() {return process->GetTree();}

	Link* GetRoot() {return process->GetRoot();}

	int GetDim() {return dim;}

	void Reset()	{
		RecursiveReset(GetRoot());
		size = 0;
		for (int k=0; k<dim; k++)	{
			gpp[k] = 0;
		}
	}

	void Add()	{
		for (int k=0; k<dim; k++)	{
			gmean[k] = 0;
		}
		RecursiveAdd(GetRoot());
		size++;
		for (int k=0; k<dim; k++)	{
			if (gmean[k] > 0)	{
				gpp[k]++;
			}
		}
	}

	void Tabulate(ostream& os, bool leafonly = false)	{
		RecursiveTabulate(os,GetRoot(), leafonly);
	}

	void Normalise()	{
		RecursiveNormalise(GetRoot());
		for (int k=0; k<dim; k++)	{
			gpp[k] /= size;
		}
	}

	double GetGmean(int k)	{
		return gmean[k];
	}


	void DAgostinosKandZ()	{

		double* m1 = new double[GetDim()];
		double* m2 = new double[GetDim()];
		double* m3 = new double[GetDim()];
		double* m4 = new double[GetDim()];

		double n = GetContrastMean(m1);
		GetContrastCentralMoment(m2,m1,2);
		GetContrastCentralMoment(m3,m1,3);
		GetContrastCentralMoment(m4,m1,4);

		cout << '\n';
		cout << "D'Agostino's test\n";
		cout << "i\tomnibus\tskewness\tkurtosis\n";
		cout << '\n';
		for (int i=0; i<GetDim(); i++)	{
			double g1 = m3[i] / sqrt(m2[i] * m2[i] * m2[i]);
			double g2 = m4[i] / m2[i] / m2[i] - 3;

			double mu2 = 6.0 * (n-2) / (n+1) / (n+3);
			double gamma2 = 36.0 * (n-7) * (n*n + 2*n -5) / (n-2) / (n+5) / (n+7) / (n+9);
			double W2 = sqrt(2*gamma2 +4)- 1;
			double W = sqrt(W2);
			double delta = 1.0 / sqrt(log(W));
			double alpha = sqrt(2 / (W2 - 1));
			double t = g1 / alpha / sqrt(mu2);
			double Z1 = delta * log(t + sqrt(t*t+1));

			double mmu1 = -6.0 / (n+1);
			double mmu2 = 24.0 * n * (n-2) * (n-3) / (n+1) / (n+1) / (n+3) / (n+5);
			double ggamma1 = 6.0 * (n*n -5*n + 2) / (n+7) / (n+9) * sqrt(6.0 * (n+3) * (n+5) / n / (n-2) / (n-3));
			double A = 6.0 + 8.0 / ggamma1* (2.0 / ggamma1 + sqrt(1 + 4.0 / ggamma1 / ggamma1));
			double Z2 = sqrt(9.0 * A / 2) * (1 - 2.0 / 9 / A - exp( 1.0 / 3 * log( (1 - 2.0 / A) / (1.0 + (g2 - mmu1) / sqrt(mmu2) * sqrt(2.0 / (A - 4))))));

			double ZZ1 = Z1 * Z1;
			double ZZ2 = Z2 * Z2;
			double KK = ZZ1 + ZZ2;
			cout << i << '\t' << KK << '\t' << ZZ1 << '\t' << ZZ2 << '\n';
		}

		delete[] m1;
		delete[] m2;
		delete[] m3;
		delete[] m4;

	}

	int GetContrastMean(double* m)	{
		for (int i=0; i<GetDim(); i++)	{
			m[i] = 0;
		}
		int n = 0;
		RecursiveGetContrastMean(GetRoot(),m,n);
		for (int i=0; i<GetDim(); i++)	{
			m[i] /= n;
		}
		return n;
	}

	int GetContrastCentralMoment(double* mk, double* m1, int k)	{
		for (int i=0; i<GetDim(); i++)	{
			mk[i] = 0;
		}
		int n = 0;
		RecursiveGetContrastCentralMoment(GetRoot(),mk,m1,k,n);
		for (int i=0; i<GetDim(); i++)	{
			mk[i] /= n;
		}
		return n;
	}

	void RecursiveGetContrastMean(const Link* from, double* m, int& n)	{
		if (! from->isRoot())	{
			double* c = mean[from->GetBranch()];
			for (int i=0; i<GetDim(); i++)	{
				m[i] += c[i];
			}
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetContrastMean(link->Out(),m,n);
		}
	}

	void RecursiveGetContrastCentralMoment(const Link* from, double* mk, const double* m1, int k, int& n)	{
		if (! from->isRoot())	{
			double* c = mean[from->GetBranch()];
			for (int i=0; i<GetDim(); i++)	{
				mk[i] += pow((c[i] - m1[i]),k);
			}
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetContrastCentralMoment(link->Out(),mk,m1,k,n);
		}
	}

	private:

	void RecursiveCreate(const Link* from)	{
		if (! from->isRoot())	{
			mean[from->GetBranch()] = new double[dim];
			var[from->GetBranch()] = new double[dim];
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveCreate(link->Out());
		}
	}

	void RecursiveDelete(const Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveDelete(link->Out());
		}
		if (! from->isRoot())	{
			delete[] var[from->GetBranch()];
			delete[] mean[from->GetBranch()];
		}
	}

	void RecursiveAdd(Link* from)	{
		if (!from->isRoot())	{
			for (int i=0; i<dim; i++)	{
				double tmp = ((process->GetNodeVal(from->GetNode())->val())[i] - (process->GetNodeVal(from->Out()->GetNode())->val())[i]) / sqrt(GetLengthTree()->GetBranchVal(from->GetBranch())->val());
				double tmp2 = ((process->GetNodeVal(from->GetNode())->val())[i] - (process->GetNodeVal(from->Out()->GetNode())->val())[i]);
				mean[from->GetBranch()][i] += tmp;
				var[from->GetBranch()][i] += tmp * tmp;
				/*
				if (tmp > 0)	{
					gmean[i]++;
				}
				*/
				gmean[i] += tmp2;
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(link->Out());
		}
	}

	void RecursiveReset(Link* from)	{
		if (! from->isRoot())	{
			for (int i=0; i<dim; i++)	{
				mean[from->GetBranch()][i] = 0;
				var[from->GetBranch()][i] = 0;
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
	}


	void RecursiveNormalise(Link* from)	{
		if (! from->isRoot())	{
			for (int i=0; i<dim; i++)	{
				mean[from->GetBranch()][i] /= size;
				var[from->GetBranch()][i] /= size;
				var[from->GetBranch()][i] -= mean[from->GetBranch()][i] * mean[from->GetBranch()][i];
				if (var[from->GetBranch()][i] < 0)	{
					var[from->GetBranch()][i] = 0;
				}
			}
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
		}
	}

	void RecursiveTabulate(ostream& os, Link* from, bool leafonly)	{
		if ((! from->isRoot()) && ((! leafonly) || (from->isLeaf())))	{
		os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
			for (int i=0; i<dim; i++)	{
				os << mean[from->GetBranch()][i] << '\t' << sqrt(var[from->GetBranch()][i]) << '\t';
			}
			os << '\n';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulate(os,link->Out(),leafonly);
		}
	}

	double* gmean;
	double* gpp;
	map<const Branch*,double*> mean;
	map<const Branch*,double*> var;
	NodeVarTree<RealVector>* process;
	LengthTree* lengthtree;
	// MultiVariateTreeProcess* process;
	int dim;
	int size;

};

