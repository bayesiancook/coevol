#ifndef RANDOMSEGMENTEDBRANCHSITEPATH_H
#define	RANDOMSEGMENTEDBRANCHSITEPATH_H

#include "RandomBranchSitePath.h"
#include "EvolutionSegmentedProcess.h"
#include "Chrono.h"

class RandomSegmentedBranchSitePath : public RandomBranchSitePath {

	protected :

		EvolutionSegmentedBranch* evol;
		int* interstate;				//The intermedial state at each node
		int* bkinterstate;			  //backup of the interstate


	public :



		RandomSegmentedBranchSitePath(PhyloProcess *phyloprocess) : RandomBranchSitePath(phyloprocess) {
		}

		RandomSegmentedBranchSitePath(PhyloProcess *phyloprocess, EvolutionSegmentedBranch* inevol) :
		RandomBranchSitePath(phyloprocess) {
			evol = inevol;
			if(isRoot())
				  SetName("Root Segmented Branch Site Path");
			else
				  SetName("Segmented Branch Site Path");

			Register(evol);
			Register(evol->GetGlobalTransMatrix());
		
			interstate = new int[GetNsegment()+1];
			bkinterstate = new int[GetNsegment()+1];
			transitionmatrix = evol->GetGlobalTransMatrix();
					  
		}


		virtual AbstractTransitionMatrix*	GetTransitionMatrix()	{
			return evol->GetGlobalTransMatrix();			  
	}

		virtual int GetNstate() {
			return evol->GetNstate();
		}

		virtual int GetNsegment() {
			return evol->GetNsegment();
		}

		virtual bool isRoot() {
			return evol->isRoot();
		}

		virtual double logProb() {
			if (isRoot())	{
		return log(evol->GetStationary()[stateup]);
			}

			double p = 1;

			for(int k=0; k<GetNsegment(); k++) {
				p *= evol->GetP(k)->GetR()[interstate[k]][interstate[k+1]];
			}
			double x = log(p);	 
		  
			return x;
		}


		virtual void Resample(int instateup, int instatedown) {
			stateup = instateup;
			statedown = instatedown;

			interstate[0] = stateup;
			interstate[GetNsegment()] = statedown;
			double *prob = new double[GetNstate()];
			for(int k=GetNsegment()-1; k>0; k--) {
				int statenext = interstate[k+1];
				for(int i=0; i<GetNstate(); i++)
					prob[i] = evol->GetP(k)->GetR()[i][statenext] * evol->GetCumul(k)[stateup][i] / evol->GetCumul(k+1)[stateup][statenext];
				interstate[k] = Random::DrawFromDiscreteDistribution(prob, GetNstate());

			}
			localUpdate();
			delete[] prob;

		}

		virtual void localRestore()	{
				for(int i = 0; i<GetNsegment()+1; i++)
					interstate[i] = bkinterstate[i];
				stateup = bkstateup;
		statedown = bkstatedown;
				RandomBranchSitePath::localRestore();
	}

	virtual void localCorrupt(bool bk)	{
				if (bk)	{
					for(int i = 0; i<GetNsegment()+1; i++)
						bkinterstate[i] = interstate[i];
					bkstateup = stateup;
					bkstatedown = statedown;
				}

				RandomBranchSitePath::localCorrupt(bk);
	}

		virtual int countA() {
			int count = 0;
			for(int i=0; i<GetNsegment(); i++) {
				if(interstate[i] == 0)
					count++;
			}
			return count;
		}


};



class BrownianPhyloProcess : public PhyloProcess {

	EvolutionSegmentedProcess *process;

public :
	BrownianPhyloProcess(LengthTree* intree, SequenceAlignment* indata, EvolutionSegmentedProcess* inprocess) :
		PhyloProcess(intree, indata, false)
	{
		process = inprocess;
	}

	RandomBranchSitePath* CreateRandomBranchSitePath(const Link* from, int site) {
		EvolutionSegmentedBranch* mat = process->GetBranchVal(from->GetBranch());
  
		return new RandomSegmentedBranchSitePath(this, mat);
	}  
  


};


#endif	/* RANDOMSEGMENTEDBRANCHSITEPATH_H */

