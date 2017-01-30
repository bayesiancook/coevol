#ifndef BROWNIANMOVE_H
#define	BROWNIANMOVE_H

#include "Move.h"
#include "BrownianProcess.h"
#include "Chrono.h"

class BrownianHorizontalMove : public MCUpdate {

private :
	Chronogram *chrono;
	BrownianProcess *process;
	double horizontalTuning;

	  

public:
	BrownianHorizontalMove(Chronogram *inchrono, BrownianProcess *inprocess, double inhorizontalTuning) {
		chrono = inchrono;
		process = inprocess;
		horizontalTuning = inhorizontalTuning;
	}

	 ~BrownianHorizontalMove() {
	 }


	double Move(double tuning_modulator) {
		int n = 0;
		double total = RecursiveMove(tuning_modulator, chrono->GetRoot(), n);
		return total/n;
	}

	double RecursiveMove(double tuning_modulator, const Link* from, int& n) {
 
		double total = 0;
		for(const Link* link = from->Next(); link != from; link = link->Next())
			total += RecursiveMove(tuning_modulator, link->Out(), n);

		if(!(from->isLeaf() || from->isRoot())) {
			total += LocalMove(tuning_modulator, from);
			n++;
		}

		return total;
	}

	 
	double LocalMove(double tuning_modulator, const Link* from) {
		//Creation of mnode
		Mnode* mnode = new Mnode(false);
		Rvar<PosReal>* nodeAge = chrono->GetNodeVal(from->GetNode());

		RandomBrownianPath* up = process->GetPureBrownianProcess()->GetRandomBrownianPath(from);
		RandomBrownianPath* left = process->GetPureBrownianProcess()->GetRandomBrownianPath(from->Next());
		RandomBrownianPath* right = process->GetPureBrownianProcess()->GetRandomBrownianPath(from->Next()->Next());
		MultiNormal *nodeVal = process->GetInstantProcess()->GetMultiNormal(from);
		int dim = nodeVal->GetDim();

		//Registering of moving parameters
		nodeAge->Register(mnode);
		up->Register(mnode);
		left->Register(mnode);
		right->Register(mnode);
		nodeVal->Register(mnode);

		double verif = 0;
		verif-=up->GetLogProb();
		verif-=left->GetLogProb();
		verif-=right->GetLogProb();
		verif-=nodeVal->GetLogProb();
		verif-=process->GetInstantProcess()->GetMultiNormal(from->Next()->Out())->GetLogProb();
		verif-=process->GetInstantProcess()->GetMultiNormal(from->Next()->Next()->Out())->GetLogProb();


		 //Morphing bridges in absolute bridges
		mnode->Corrupt(true);
		up->addSlope(process->GetInstantProcess()->GetNodeVal(from->Out()->GetNode())->GetArray(), nodeVal->GetArray());
		left->addSlope(nodeVal->GetArray(), process->GetInstantProcess()->GetNodeVal(from->Next()->Out()->GetNode())->GetArray());
		right->addSlope(nodeVal->GetArray(), process->GetInstantProcess()->GetNodeVal(from->Next()->Next()->Out()->GetNode())->GetArray());

		//Moving the ages
		double oldAge = nodeAge->val();
		chrono->SetMinMax(from);
		double logHastings = nodeAge->ProposeMove(horizontalTuning*tuning_modulator);
		double newAge = nodeAge->val();


		//Moving the instant value		  
		double t0, t1, t2;	  //Time at the window extremities
		double *x0,*x1,*x2;	   //Pure brownian values at the window extremities

		up->findPointBefore(oldAge>newAge ? oldAge : newAge, t0, x0);
		left->findPointAfter(oldAge<newAge ? oldAge : newAge, t1, x1);
		right->findPointAfter(oldAge<newAge ? oldAge : newAge, t2, x2);

		double r0 = 1.0/(t0-oldAge), r1 = 1.0/(oldAge-t1), r2 = 1.0/(oldAge-t2);   //Segments inverse lengths
		double r = r0+r1+r2;
		double* mean = new double[dim];
		for(int i=0; i<dim; i++)
			mean[i] = (r0*x0[i]+r1*x1[i]+r2*x2[i])/r;


		logHastings+=up->logProbNormal(nodeVal->GetArray(), mean, 1.0/r, up->getUnitVariance().GetInvMatrix(), 1.0/up->getUnitVariance().GetDeterminant());

		r0 = 1.0/(t0-newAge); r1 = 1.0/(newAge-t1); r2 = 1.0/(newAge-t2);   //Segments inverse lengths

		r = r0+r1+r2;
		for(int i=0; i<dim; i++)
			mean[i] = (r0*x0[i]+r1*x1[i]+r2*x2[i])/r;

		up->getUnitVariance().drawVal(nodeVal->GetArray());

		nodeVal->ScalarMultiplication(sqrt(1.0/r));
		for(int i=0; i<up->getDim(); i++) {
			nodeVal->GetArray()[i] += (r0*x0[i]+r1*x1[i]+r2*x2[i])/r;
		}

		logHastings-=up->logProbNormal(nodeVal->GetArray(), mean, 1.0/r, up->getUnitVariance().GetInvMatrix(), 1.0/up->getUnitVariance().GetDeterminant());

		delete[] mean;


	  

		//Moving the bridges
		double s = up->finalAgeMove(newAge, nodeVal->GetArray());
		logHastings +=s;
		up->removeSlope();

		s = left->initAgeMove(newAge, nodeVal->GetArray());
		logHastings +=s;
		left->removeSlope();

		s = right->initAgeMove(newAge, nodeVal->GetArray());
		logHastings +=s;
		right->removeSlope();
	  
	  
		verif+=up->GetLogProb();
		verif+=left->GetLogProb();
		verif+=right->GetLogProb();
		verif+=nodeVal->GetLogProb();
		verif+=process->GetInstantProcess()->GetMultiNormal(from->Next()->Out())->GetLogProb();
		verif+=process->GetInstantProcess()->GetMultiNormal(from->Next()->Next()->Out())->GetLogProb();

 
		//Validating the move
		double logratio = mnode->Update();
		bool accepted = (log(Random::Uniform()) < logratio+logHastings);
		if (! accepted)	{
			mnode->Corrupt(false);
			mnode->Restore();
		}
			  
		delete mnode;


		return (double) accepted;

	}

	void ToStream(ostream& os) {
	 
	}

};

#endif	/* BROWNIANMOVE_H */

