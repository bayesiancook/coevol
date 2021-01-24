
#ifndef BRANCHSITESUBPROCESS_H
#define BRANCHSITESUBPROCESS_H

#include "BranchSitePath.h"
#include "Var.h"
#include "BranchSiteSubstitutionProcess.h"
#include "RandomSubMatrix.h"

#include "EmpiricalSubMatrix.h"

// RandomBranchSitePath
// - represents a branch- and site-specific substitution process (the BranchSiteSubstitutionProcess superclass)
// - maintains a realisation of the substitution history along a branch/site (a path) conditional on the data (the BranchSitePath super class)
//
// - it is an abstract class (BranchSiteSubstitutionProcess has many pure virtual methods that need to be implemented in subclasses)
// - it is a random node of the model's graph (the Rnode superclass)
// - subclasses will be responsible for making the required connections between this Rnode and the relevant parameters of the process
//   on a model-specific basis

class PhyloProcess;

class RandomBranchSitePath : public virtual Rnode, public BranchSitePath, public BranchSiteSubstitutionProcess	{


	public:

	RandomBranchSitePath(PhyloProcess* inprocess) : rate(0), length(0), matrix(0), stationary(0), myprocess(inprocess), propmatrix(0), matswap(false) {
		SetName("path");
	}

	RandomBranchSitePath(PhyloProcess* inprocess, Var<PosReal>* inlength, Var<PosReal>* inrate, RandomSubMatrix* inmatrix, Var<Profile>* instationary) : myprocess(inprocess), propmatrix(0), matswap(false)	{
		SetName("path");
		length = inlength;
		rate = inrate;
		matrix = inmatrix;
		stationary = instationary;
		active_flag = true;
		if (stationary)	{
			cerr << "in random branch site path: stat is non null\n";
			exit(1);
		}

		if ((! matrix) && (! stationary))	{
			cerr << "error in RandomBranchSitePath: should specify a matrix or a stationary\n";
			exit(1);
		}

		Register(length);
		Register(rate);
		Register(matrix);
		Register(stationary);
		// Sample();
	}

	RandomBranchSitePath(PhyloProcess* inprocess, RandomTransitionMatrix* inmatrix, Var<Profile>* instationary, Var<PosReal>* inlength) : myprocess(inprocess){
		SetName("path");
		propmatrix = 0;
		matrix = 0;
		length = inlength;
		rate = 0;
		matswap = false;
		transitionmatrix = inmatrix;
		stationary = instationary;
		Register(transitionmatrix);
		if (stationary)	{
			cerr << "non null stat in random branch site path\n";
			exit(1);
		}
		// Register(stationary);
		// Register(length);
	}

	PhyloProcess*			GetPhyloProcess() 	{return myprocess;}

	bool SampleBranchMapping();

	virtual AbstractTransitionMatrix*	GetTransitionMatrix()	{
		return transitionmatrix;
	}

	virtual SubMatrix* 		GetSubMatrix()		{
		if (matswap)	{
			if (! propmatrix)	{
				cerr << "error in random branch site path: null prop matrix\n";
				exit(1);
			}
			return propmatrix;
		}
		if (! matrix)	{
			cerr << "error in random branch site path: getsubmatrix\n";
			exit(1);
		}
		return matrix;
	}

	virtual void			SwapMatrices()	{
		if ((! matswap) && (! propmatrix))	{
			cerr << "error in random branch site path swap matrix: null prop matrix\n";
			exit(1);
		}
		matswap = ! matswap;
		/*
		RandomSubMatrix* tmp = matrix;
		matrix = propmatrix;
		propmatrix = tmp;
		*/
	}

	// void	SetProposalMatrix(EmpiricalSubMatrix* inmatrix)	{
	void	SetProposalMatrix(SubMatrix* inmatrix)	{
		propmatrix = inmatrix;
	}

	void BackwardPropagate(const double* down, double* up)	{
		if (SampleBranchMapping())	{
			BranchSiteSubstitutionProcess::BackwardPropagate(down,up);
		}
		else	{
			transitionmatrix->BackwardPropagate(down,up,0);
		}
	}

	void ForwardPropagate(const double* up, double* down)	{
		if (SampleBranchMapping())	{
			BranchSiteSubstitutionProcess::ForwardPropagate(up,down);
		}
		else	{
			transitionmatrix->ForwardPropagate(up,down,0);
		}
	}

	void GetFiniteTimeTransitionProb(int state, double* aux)	{
		if (SampleBranchMapping())	{
			BranchSiteSubstitutionProcess::GetFiniteTimeTransitionProb(state,aux);
		}
		else	{
			TransitionGetFiniteTimeTransitionProb(state,aux);
		}
	}

	int DrawStationary()	{
		if (SampleBranchMapping())	{
			return BranchSiteSubstitutionProcess::DrawStationary();
		}
		return TransitionDrawStationary();
	}

	int DrawFiniteTime(int state)	{
		if (SampleBranchMapping())	{
			return BranchSiteSubstitutionProcess::DrawFiniteTime(state);
		}
		return TransitionDrawFiniteTime(state);
	}

	void TransitionGetFiniteTimeTransitionProb(int state, double* aux)	{
		const double* p = transitionmatrix->GetRow(state);
		for (int i=0; i<GetNstate(); i++)	{
			aux[i] = p[i];
		}
	}

	int TransitionDrawStationary()	{
		const double* p = transitionmatrix->GetStationary();
		int newstate = Random::DrawFromDiscreteDistribution(p,GetNstate());
		return newstate;
	}

	int TransitionDrawFiniteTime(int state)	{
		const double* p = transitionmatrix->GetRow(state);
		int newstate = Random::DrawFromDiscreteDistribution(p,GetNstate());
		return newstate;
	}

	virtual int GetNstate() {
		if (matrix)	{
			return GetSubMatrix()->GetNstate();
		}
		else if (transitionmatrix)	{
			return transitionmatrix->GetNstate();
		}
		else	{
			cerr << "error in RandomBranchSitePath: GetNstate\n";
			exit(1);
		}
	}

	virtual double 			GetRate() 		{return rate ?((double) rate->val()) : 1;}
	virtual const double* 		GetStationary()		{
			if (stationary)	{
				cerr << "error : non null stationary in log prob path\n";
				exit(1);
				return stationary->GetArray();
			}
			if (SampleBranchMapping())	{
				return GetSubMatrix()->GetStationary();
			}
			return transitionmatrix->GetStationary();
		}
			// return stationary ? stationary->GetArray() : matrix->GetStationary();}

	virtual double 			GetTime()		{
					// 		return return length ? ((double) length->val()) : 0;}
						if (length)	{
							if (std::isnan(((double) (length->val()))))	{
								cerr << "length is nan\n";
							}
							return length->val();
						}
						return 0;
					}
		

	bool 				isRoot() 		{return (length == 0);}

	bool				isActivated()		{return active_flag;}

	void				SetActivated(bool inflag)	{active_flag = inflag;}

	StateSpace* GetStateSpace();

	void 			SetUp(RandomBranchSitePath* inup);

	double 			GetTotalTime() {return GetTime();}
	void			SetTotalTime(double intime) {
		// what about corruption ?
		if (! length)	{
			cerr << "error in RandomBranchSitePath::SetTotalTime\n";
			exit(1);
		}
		length->setval(intime);
	}

	int			GetMaxTrial();

	virtual void 			Resample()	{
					if (! SampleBranchMapping())	{
						cerr << "error in random branch site path : resample map called\n";
						exit(1);
					}
					if (isRoot())	{
						cerr << "error in RandomBranchSitePath::Resample : called on root\n";
						exit(1);
					}
					else	{
						if (GetTime() == 0)	{
							cerr << "error in resample : null bl\n";
							exit(1);
						}
						bool ok = ResampleAcceptReject(GetMaxTrial());
						if (! ok)	{
							ResampleUniformized();
						}
						/*
						if (propmatrix)	{
							ResampleUniformized();
						}
						else	{
							bool ok = ResampleAcceptReject(GetMaxTrial());
							if (! ok)	{
								ResampleUniformized();
							}
						}
						*/
					}
				}

	double 			Move(double tuning)	{
					if (!propmatrix)	{
						Resample();
						return 1;
					}
					return Rnode::Move(tuning);
				}

	virtual void	localRestore()	{
		if (propmatrix)	{
			stateup = bkstateup;
			statedown = bkstatedown;
			RestorePath();
			logprob = bklogprob;
			/*
			if (logprob != logProb())	{
				cerr << "error in local restore\n";
				cerr << logprob - logProb() << '\n';
				exit(1);
			}
			*/
			flag = true;
		}
		else	{
			Rnode::localRestore();
		}
	}

	virtual void	localCorrupt(bool bk)	{
		if (propmatrix)	{
			if (bk)	{
				bkstateup = stateup;
				bkstatedown = statedown;
				BackupPath();
				bklogprob = logprob;
			}
			flag = false;
		}
		else	{
			Rnode::localCorrupt(bk);
		}
	}

	void 			drawSample()	{
					cerr << "in random branch site path drawsample\n";
					exit(1);
					Resample();
				}

	double 			logProb();
	double			PathLogProb();
	double			NoPathLogProb();

	void 			ResampleRoot(int state)	{
					stateup = statedown = state;
					Reset(state);
					localUpdate();
				}

	void 			Resample(int instateup, int instatedown)	{
					stateup = instateup;
					statedown = instatedown;
					if (SampleBranchMapping())	{
						if (propmatrix)	{
							SwapMatrices();
							Resample();
							SwapMatrices();
						}
						else	{
							Resample();
						}
					}
					else	{
						Reset(stateup);
					}
					localUpdate();
				}

	void 			SetMatrix(RandomSubMatrix* inmatrix) {matrix = inmatrix;}

	double ProposeMove(double tuning)	{
		cerr << "error : in random branch site path propose move(tuning)\n";
		exit(1);
	}

	double ProposeMove(int instateup, int instatedown)	{
		double logbefore = logProb();
		stateup = instateup;
		statedown = instatedown;
		Resample();
		double logafter = logProb();
		return logbefore-logafter;
	}

	double ProposeMoveRoot(int instate)	{
		double logbefore = logProb();
		stateup = statedown = instate;
		Reset(instate);
		double logafter = logProb();
		return logbefore-logafter;
	}

	protected:

	bool 			ResampleAcceptReject(int maxtrial);
	void			ResampleUniformized();
	// double			RecordResampleUniformizedLogProb();

	// RandomBranchSitePath* pathup;
	Var<PosReal>* rate;
	Var<PosReal>* length;

	RandomSubMatrix* matrix;
	RandomTransitionMatrix* transitionmatrix;

	protected:

	Var<Profile>* stationary;

	int stateup;
	int statedown;
	int bkstateup;
	int bkstatedown;
	PhyloProcess* myprocess;
	SubMatrix* propmatrix;
	// EmpiricalSubMatrix* propmatrix;
	bool active_flag;

	bool matswap;
	/*
	int uninsub;
	vector<int> unistate;
	vector<double> unitime;
	*/

};


#endif // RANDOMSITEPATH_H
