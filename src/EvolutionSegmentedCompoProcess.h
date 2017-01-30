/*
 * File:   EvolutionSegmentedCompoProcess.h
 * Author: bhorvill
 *
 * Created on May 29, 2013, 4:17 PM
 */

#ifndef EVOLUTIONSEGMENTEDCOMPOPROCESS_H
#define	EVOLUTIONSEGMENTEDCOMPOPROCESS_H

#include "EvolutionSegmentedProcess.h"
#include "Chrono.h"
#include "BrownianMove.h"
#include "SumConstrained.h"

class CompoMatrix : virtual public EvolutionSegmentedMatrix {
	public :


		CompoMatrix(RandomBrownianPath* bridge, Var<RealVector>* up, Var<RealVector>* down, Profile *inrelrate, SumConstrainedMapping* inmapping, int k, int inoffset, int Nstate)  :
			TransitionMatrix(Nstate), EvolutionSegmentedMatrix(bridge, up, down, Nstate, k) {
			relrate = inrelrate;
			rStationary = 0;
			offset = inoffset;
			mapping = inmapping;
		  
			ComputeArray();
		}

			//Implements TransitionMatrix::ComputeArray()
			virtual void ComputeR(double lengthLeft, double l) {

			  
				double length = brownianBridge->getLength();
				double dt = brownianBridge->getSegmentLength(k);

				//Compute the values at the sub-branch extremities, without accounting the brownian bridge (using Thales theorem)					 
				double r1 = up[POS_R] + lengthLeft/length * (down[POS_R]-up[POS_R]);
				double r2 = up[POS_R] + (lengthLeft + dt)/length * (down[POS_R]-up[POS_R]);

				//Compute the values of the integrated rate
				double rm = 0.5 * ( exp(r1 + brownianBridge->getValue(k, POS_R)) + exp(r2 + brownianBridge->getValue(k+1, POS_R)) );

				// vector of eq. freq.
				double** b = mapping->base;

				double tempup[GetNstate()];
				double tempdown[GetNstate()];

				double stat[GetNstate()];

				tempup[0] = 0;
				tempdown[0] = 0;

				for (int i=0; i<GetNstate()-1; i++)	{
					tempup[i+1] = up[offset] + lengthLeft/length * (down[offset]-up[offset]) + brownianBridge->getValue(k,offset);
					tempdown[i+1] = up[offset] + (lengthLeft + dt)/length * (down[offset]-up[offset]) + brownianBridge->getValue(k,offset);
				}

				double total = 0;
				double totup = 0;
				double totdown = 0;
				for (int i=0; i<GetNstate(); i++)	{
					double tmpup = 0;
					double tmpdown = 0;
					for (int j=0; j<GetNstate(); j++)	{
						tmpup += b[j][i] * tempup[j];
						tmpdown += b[j][i] * tempdown[j];
					}
					totup += tmpup;
					totdown += tmpdown;
					stat[i] = exp(tmpup) + exp(tmpdown);
					total += stat[i];
				}
				if (fabs(totup) > 1e-6)	{
					cerr << "error : non matching sum\n";
					cerr << "up : " << totup << '\n';
				}
				if (fabs(totdown) > 1e-6)	{
					cerr << "error : non matching sum\n";
					cerr << "down : " << totdown << '\n';
				}
				for (int i=0; i<GetNstate(); i++)	{
					stat[i] /= total;
				}

				for(int i=0; i<GetNstate(); i++) {
					double tot = 0;
					for(int j=0; j<GetNstate(); j++) {
						if(i!=j) {
							R[i][j] = rm*l*RelativeRate(i,j) * stat[j];
							tot+=R[i][j];
						}
					}
					R[i][i]=-tot;
				}
			}

			double RelativeRate(int i, int j) {return (*relrate)[rrindex(i,j,GetNstate())];}

			static int rrindex(int i, int j, int nstate)	{
					return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
			}

			//Compute the stationary from the beginning
			virtual void ComputeStationary() {
				for(int i=0; i<GetNstate(); i++) {
					stationary[i] = 1.0 / GetNstate();
				}		  
			}

	private :
	  
		Profile *relrate;
		Profile* rStationary;
		SumConstrainedMapping* mapping;
		int offset;
};


class RootCompoMatrix : virtual public TransitionMatrix, virtual public RandomTransitionMatrix {
	public :
	  RootCompoMatrix(Var<Profile>* instationary, int Nstate) : TransitionMatrix(Nstate, 0, instationary->GetArray()) {
		  Register(instationary);

	  }

	  virtual void ComputeArray() {
	  }
	  virtual void ComputeStationary() {
	  }
	  virtual void SetParameters() {
	  }

	private :
};


class EvolutionSegmentedCompoBranch : public EvolutionSegmentedBranch {
  
	public :

	EvolutionSegmentedCompoBranch(RandomBrownianPath* inbridge, Var<RealVector>* inup, Var<RealVector>* indown, Var<Profile>* inrelrate, SumConstrainedMapping* inmapping, int inoffset, int nState, bool paral) :
		EvolutionSegmentedBranch(inbridge, inup, indown, nState, paral) {
		relrate = inrelrate;
		mapping = inmapping;
		offset = inoffset;
		Register(relrate);
		CreateAllP();
	}

	EvolutionSegmentedCompoBranch(Var<RealVector>* inup, Var<Profile>* instationary) :
		EvolutionSegmentedBranch(inup, instationary->GetDim()) {
		rStationary = instationary;
		Register(rStationary);
		globalTransMatrix = CreateRootMatrix();
	}

	virtual EvolutionSegmentedMatrix* CreateP(int i) {
		CompoMatrix* tmp = new CompoMatrix(bridge, up, down, relrate, mapping, i, offset, GetNstate());
		return tmp;
	}

	virtual RandomTransitionMatrix* CreateRootMatrix() { 
		return new RootCompoMatrix(rStationary, rStationary->GetDim());
	}

	Var<Profile> *relrate;
	Var<Profile> *rStationary;
	SumConstrainedMapping* mapping;
	int offset;
};

class RootCompoStationary : public Dvar<Profile> {
	public :
		RootCompoStationary(Var<UnitReal> *ingamma) {
			dim = 4;
			profile = new double[dim];
			gamma = ingamma;
			Register(gamma);
			specialUpdate();

		}

		void specialUpdate() {
			for(int i=0; i<dim; i++)
				profile[i] = ((i == 1 || i == 2) ? gamma->val()/2 : (1-gamma->val())/2);
		}

		Var<UnitReal> *gamma;
};


class EvolutionSegmentedCompoProcess : public EvolutionSegmentedProcess {

public :

	EvolutionSegmentedCompoProcess(BrownianProcess *inprocess, Var<Profile>* inrelrate, Var<Profile>* rootstat, SumConstrainedMapping* inmapping, int inoffset, bool paral) :
		EvolutionSegmentedProcess(inprocess, paral) {
		relrate = inrelrate;
		rStationary = rootstat;
		mapping = inmapping;
		offset = inoffset;
		RecursiveCreate(GetRoot());
	}

	virtual EvolutionSegmentedBranch* CreateBranchVal(const Link* from) {
	   // grep the instant value associated to this node
		MultiNormal* vup = process->GetInstantProcess()->GetMultiNormal(from);

		if(from->isRoot())	{
			EvolutionSegmentedBranch* tmp = new EvolutionSegmentedCompoBranch(vup, rStationary);
			return tmp;
		}
		// grep the instant value associated to the node immediately upstream
		MultiNormal* vdown = process->GetInstantProcess()->GetMultiNormal(from->Out());

		// grep the brownian path of the branch
		RandomBrownianPath* brownianPath = process->GetPureBrownianProcess()->GetRandomBrownianPath(from);

		EvolutionSegmentedBranch* tmp = new EvolutionSegmentedCompoBranch(brownianPath, vup, vdown, relrate, mapping, offset, rStationary->GetDim(), paral);
		return tmp;
	}



	Var<Profile> *relrate;
	Var<Profile> *rStationary;
	SumConstrainedMapping* mapping;
	int offset;


};


#endif	/* EVOLUTIONSEGMENTEDCOMPOPROCESS_H */

