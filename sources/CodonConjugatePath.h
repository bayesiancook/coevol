
#ifndef MGCODONCONJPATH_H
#define MGCODONCONJPATH_H

#include "GeneralConjugatePath.h"
#include "CodonSubMatrix.h"

class MGCodonPathConjugate : public virtual PathConjugate	{

	public:

	MGCodonPathConjugate(Var<PosReal>* inlength, Var<PosReal>* inrate, RandomSubMatrix* inmatrix, Var<Profile>* instationary) : PathConjugate(inlength,inrate,inmatrix,instationary)	{
		nucflag = false;
		codonmatrix = dynamic_cast<MGCodonSubMatrix*>(matrix);
		if (! codonmatrix)	{
			cerr << "error in mg codon path: null matrix\n";
			exit(1);
		}
		CreateNucSuffStat();
	}

	virtual ~MGCodonPathConjugate()	{
		DeleteNucSuffStat();
	}

	void CreateNucSuffStat()	{
		nucrootcount = new int[Nnuc];
		synpaircount = new int*[Nnuc];
		nonsynpaircount = new int*[Nnuc];
		synpairbeta = new double*[Nnuc];
		nonsynpairbeta = new double*[Nnuc];
		for (int i=0; i<Nnuc; i++)	{
			synpaircount[i] = new int[Nnuc];
			nonsynpaircount[i] = new int[Nnuc];
			synpairbeta[i] = new double[Nnuc];
			nonsynpairbeta[i] = new double[Nnuc];
		}
	}

	void DeleteNucSuffStat()	{
		delete[] nucrootcount;
		for (int i=0; i<Nnuc; i++)	{
			delete[] synpaircount[i];
			delete[] nonsynpaircount[i];
			delete[] synpairbeta[i];
			delete[] nonsynpairbeta[i];
		}
		delete[] synpaircount;
		delete[] nonsynpaircount;
		delete[] synpairbeta;
		delete[] nonsynpairbeta;
	}

	void ComputeTotSuffStat()	{

		ComputeNucSuffStat();
		totsyncount = 0;
		totnonsyncount = 0;
		totsynbeta = 0;
		totnonsynbeta = 0;

		double** synq = codonmatrix->GetSynNucArray();
		double** nonsynq = codonmatrix->GetNonSynNucArray();
		// double omega = codonmatrix->GetOmega();

		for (int i=0; i<Nnuc; i++)	{
			for (int j=0; j<Nnuc; j++)	{
				if (i != j)	{
					totsyncount += synpaircount[i][j];
					totnonsyncount += nonsynpaircount[i][j];
					totsynbeta += synpairbeta[i][j] * synq[i][j];
					totnonsynbeta += nonsynpairbeta[i][j] * nonsynq[i][j];
				}
			}
		}
		totsynbeta *= GetRate() * GetTime();
		totnonsynbeta *= GetRate() * GetTime();
	}

	int GetTotSynCount()	{
		return totsyncount;
	}
	int GetTotNonSynCount()	{
		return totnonsyncount;
	}
	double GetTotSynBeta()	{
		return totsynbeta;
	}
	double GetTotNonSynBeta()	{
		return totnonsynbeta;
	}

	protected:

	virtual void ResetSufficientStatistic()	{
		PathConjugate::ResetSufficientStatistic();
		nucflag = false;
	}

	void ComputeNucSuffStat()	{

		for (int i=0; i<Nnuc; i++)	{
			nucrootcount[i] = 0;
			for (int j=0; j<Nnuc; j++)	{
				synpaircount[i][j] = 0;
				nonsynpaircount[i][j] = 0;
				synpairbeta[i][j] = 0;
				nonsynpairbeta[i][j] = 0;
			}
		}

		CodonStateSpace* cod = codonmatrix->GetCodonStateSpace();

		// root part
		if (! GetTime())	{
			for (map<int,int>::iterator i = rootcount.begin(); i!= rootcount.end(); i++)	{
				int codon = i->first;
				nucrootcount[cod->GetCodonPosition(0,codon)] += i->second;
				nucrootcount[cod->GetCodonPosition(1,codon)] += i->second;
				nucrootcount[cod->GetCodonPosition(2,codon)] += i->second;
			}
		}

		else	{
			for (map<int,double>::iterator i = waitingtime.begin(); i!= waitingtime.end(); i++)	{
				int codon = i->first;
				for (int c2 = 0; c2 < cod->GetNstate(); c2++)	{
					if (c2 != codon)	{
						int pos = cod->GetDifferingPosition(codon,c2);
						if (pos < 3)	{
							int n1 = cod->GetCodonPosition(pos,codon);
							int n2 = cod->GetCodonPosition(pos,c2);
							if (cod->Synonymous(codon,c2))	{
								synpairbeta[n1][n2] += i->second;
							}
							else	{
								nonsynpairbeta[n1][n2] += i->second;
							}
						}
					}
				}
			}

			for (map<pair<int,int>, int>::iterator i = paircount.begin(); i!= paircount.end(); i++)	{
				int cod1 = i->first.first;
				int cod2 = i->first.second;
				int pos = cod->GetDifferingPosition(cod1,cod2);
				if (pos == 3)	{
					cerr << "error in codon conj path suffstat\n";
					exit(1);
				}
				int n1 = cod->GetCodonPosition(pos,cod1);
				int n2 = cod->GetCodonPosition(pos,cod2);
				if (cod->Synonymous(cod1,cod2))	{
					synpaircount[n1][n2] += i->second;
				}
				else	{
					nonsynpaircount[n1][n2] += i->second;
				}
			}
		}
	}

	double SuffStatLogProb()	{


		if (! nucflag)	{
			ComputeNucSuffStat();
			nucflag = true;
		}
		double total = 0;

		// root part
		if (! GetTime())	{
			int nroot = 0;
			const double* rootstat = codonmatrix->GetNucMatrix()->GetStationary();
			for (int i=0; i<Nnuc; i++)	{
				total += nucrootcount[i] * log(rootstat[i]);
				nroot += nucrootcount[i];
			}
			total -= nroot / 3 * log(codonmatrix->GetNormStat());
		}

		// non root part
		else	{
			int totnsub = 0;
			double totscale = 0;
			double** synq = codonmatrix->GetSynNucArray();
			double** nonsynq = codonmatrix->GetNonSynNucArray();
			for (int i=0; i<Nnuc; i++)	{
				for (int j=0; j<Nnuc; j++)	{
					if (i != j)	{
						total += synpaircount[i][j] * log(synq[i][j]);
						totnsub += synpaircount[i][j];
						total += nonsynpaircount[i][j] * log(nonsynq[i][j]);
						totnsub += nonsynpaircount[i][j];

						totscale += synpairbeta[i][j] * synq[i][j];
						totscale += nonsynpairbeta[i][j] * nonsynq[i][j];
					}
				}
			}

			total -= GetRate() * GetTime() * totscale;
			if (totnsub)	{
				total += totnsub * log(GetRate() * GetTime());
			}
		}

		return total;
	}

	int* nucrootcount;
	int** synpaircount;
	int** nonsynpaircount;
	double** synpairbeta;
	double** nonsynpairbeta;

	int totsyncount;
	int totnonsyncount;
	double totsynbeta;
	double totnonsynbeta;

	bool nucflag;

	MGCodonSubMatrix* codonmatrix;

};



class MGCodonBranchMatrixPathConjugateTree : public PathConjugateTree	{

	public:

	MGCodonBranchMatrixPathConjugateTree(LengthTree* inlengthtree, BranchValPtrTree<RandomSubMatrix>* inmatrixtree,  SequenceAlignment* indata) : PathConjugateTree(inlengthtree, indata) {
		matrixtree = inmatrixtree;
		SetWithRoot(true);
		if (! matrixtree->WithRoot())	{
			cerr << "error in PathConjugateTree: matrixtree does not have root value\n";
			exit(1);
		}
		RecursiveCreate(GetRoot());
	}

	~MGCodonBranchMatrixPathConjugateTree()	{
		RecursiveDelete(GetRoot());
	}

	RandomMGCodonSubMatrix* GetRandomMGCodonSubMatrix(const Branch* branch)	{
		RandomMGCodonSubMatrix* tmp = dynamic_cast<RandomMGCodonSubMatrix*>(matrixtree->GetBranchVal(branch));
		if (! tmp)	{
			cerr << "error in mg codon branch matrix path conjugate tree\n";
			exit(1);
		}
		return tmp;
	}

	void ComputeTotSuffStat()	{
		/*
		totsyncount = 0;
		totnonsyncount = 0;
		totsynbeta = 0;
		totnonsynbeta = 0;
		*/
		RecursiveComputeTotSuffStat(GetRoot());
	}

	MGCodonPathConjugate* GetBranchMGCodonPathConjugate(const Branch* branch)	{
		return dynamic_cast<MGCodonPathConjugate*> (GetBranchVal(branch));
	}

	protected:

	void RecursiveComputeTotSuffStat(const Link* from)	{
		if (! from->isRoot())	{
			/*
			totsyncount += GetBranchVal(from->GetBranch())->GetTotSynCount();
			totnonsyncount += GetBranchVal(from->GetBranch())->GetTotNonSynCount();
			totsynbeta+= GetBranchVal(from->GetBranch())->GetTotSynBeta();
			totnonsynbeta+= GetBranchVal(from->GetBranch())->GetTotNonSynBeta();
			*/
			GetBranchMGCodonPathConjugate(from->GetBranch())->ComputeTotSuffStat();
		}

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveComputeTotSuffStat(link->Out());
		}
	}

	PathConjugate* CreateBranchVal(const Link* link)	{
		return new MGCodonPathConjugate(lengthtree->GetBranchVal(link->GetBranch()), 0, matrixtree->GetBranchVal(link->GetBranch()), 0);
	}

	BranchValPtrTree<RandomSubMatrix>* matrixtree;

	/*
	int totsyncount;
	int totnonsyncout;
	int totsynbeta;
	int totnonsynbeta;
	*/

};

class MGCodonSiteMatrixPathConjugateProcess : public PathConjugateProcess	{

	public:

	MGCodonSiteMatrixPathConjugateProcess(LengthTree* inlengthtree, RandomSubMatrix** inmatrixarray, SequenceAlignment* indata) : PathConjugateProcess(inlengthtree, indata) {
		matrixarray = inmatrixarray;
		RecursiveCreate(GetRoot());
		CreateTotSuffStat();
	}

	~MGCodonSiteMatrixPathConjugateProcess()	{
		DeleteTotSuffStat();
		RecursiveDelete(GetRoot());
	}

	RandomMGCodonSubMatrix* GetRandomMGCodonSubMatrix(int site)	{
		RandomMGCodonSubMatrix* tmp = dynamic_cast<RandomMGCodonSubMatrix*>(matrixarray[site]);
		if (! tmp)	{
			cerr << "error in mg codon branch matrix path conjugate tree\n";
			exit(1);
		}
		return tmp;
	}

	MGCodonPathConjugate* GetMGCodonPathConjugate(const Branch* branch, int site)	{
		return dynamic_cast<MGCodonPathConjugate*> (GetPathConjugate(branch,site));
	}

	void CreateTotSuffStat()	{
		totsyncount = new int[GetNsite()];
		totnonsyncount = new int[GetNsite()];
		totsynbeta = new double[GetNsite()];
		totnonsynbeta = new double[GetNsite()];
	}

	void DeleteTotSuffStat()	{
		delete[] totsyncount;
		delete[] totnonsyncount;
		delete[] totsynbeta;
		delete[] totnonsynbeta;
	}

	void ResetTotSuffStat()	{
		for (int i=0; i<GetNsite(); i++)	{
			totsyncount[i] = 0;
			totnonsyncount[i] = 0;
			totsynbeta[i] = 0;
			totnonsynbeta[i] = 0;
		}
	}
		
	void ComputeTotSuffStat()	{
		ResetTotSuffStat();
		RecursiveComputeTotSuffStat(GetRoot());
	}

	int GetTotSynCount(int i)	{
		return totsyncount[i];
	}
	int GetTotNonSynCount(int i)	{
		return totnonsyncount[i];
	}
	double GetTotSynBeta(int i)	{
		return totsynbeta[i];
	}
	double GetTotNonSynBeta(int i)	{
		return totnonsynbeta[i];
	}

	protected:

	PathConjugate* CreatePathConjugate(const Link* link, int site)	{
		return new MGCodonPathConjugate(lengthtree->GetBranchVal(link->GetBranch()), 0, matrixarray[site], 0);
	}

	void RecursiveComputeTotSuffStat(const Link* from)	{
		if (! from->isRoot())	{
			for (int i=0; i<GetNsite(); i++)	{
				GetMGCodonPathConjugate(from->GetBranch(),i)->ComputeTotSuffStat();
				totsyncount[i] += GetMGCodonPathConjugate(from->GetBranch(),i)->GetTotSynCount();
				totnonsyncount[i] += GetMGCodonPathConjugate(from->GetBranch(),i)->GetTotNonSynCount();
				totsynbeta[i] += GetMGCodonPathConjugate(from->GetBranch(),i)->GetTotSynBeta();
				totnonsynbeta[i] += GetMGCodonPathConjugate(from->GetBranch(),i)->GetTotNonSynBeta();
			}
		}

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveComputeTotSuffStat(link->Out());
		}
	}

	RandomSubMatrix** matrixarray;

	int* totsyncount;
	int* totnonsyncount;
	double* totsynbeta;
	double* totnonsynbeta;

};

#endif
