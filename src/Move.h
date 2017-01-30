

#ifndef MOVE_H
#define MOVE_H

#include "MCMC.h"
#include <vector>

#include "DAGnode.h"
#include "RandomTypes.h"

// #include "Normal.h"

class ProbModel;


// MCMC Update mechanisms (Metroplois Hastings, or Gibbs)
// can be done in 2 different ways
//
// 1. the simplest is to make a Move() method in the Model
// which calls the Move() methods of all the random variables of the model, each in turn
//
// 2. a more elaborate but more powerful method consists in
// making objects deriving from MCUpdate
// each object is responsible for making the update of one random node, or a group of nodes
//
// these objects are REGISTERED onto a MCScheduler object
// MCScheduler will then call the updates each one by one
// while maintaining statistics such as time spent in each update, success rates, etc
//
// All this is created at the level of the model:
// the model has a MCScheduler (called 'scheduler')
// the model has a method called MakeScheduler(), whose job is to create and register all required MCUpdate objects
// When the Move() method of the ProbModel object is called
// this calls the Cycle() function of the MCScheduler object
// which in turns calls the Move() function of each MCUpdate objects
// which themselves call the relevant function (e.g. the Move() function) of the random variables they are associated to

// a MCUpdate object has a pure virtual Move function
//
class MCUpdate	{

	public:

	virtual ~MCUpdate() {}

	virtual double Move(double tuning_modulator = 1) = 0;

	virtual void ToStream(ostream& os) {}
};



// The MCScheduler class
// a ProbModel has a MCScheduler object called scheduler
//
class MCScheduler : public MCUpdate {

	public:

	MCScheduler(ProbModel* inmodel) : model(inmodel), closed(false) , random(false) {}

	virtual void Cycle(double tuning_modulator, int nrep, bool verbose, bool check);
	void RandomCycle(double tuning_modulator, int nrep, bool verbose, bool check);
	void Move(double tuning_modulator, int i, bool verbose, bool check, int nrep);

	double Move(double tuning_modulator = 1)	{
		if (random)	{
			RandomCycle(tuning_modulator,1,false,false);
		}
		else	{
			Cycle(tuning_modulator,1,false,false);
		}
		return 1;
	}

	void SetRandom(bool inrand)	{
		random = inrand;
	}
	void Register(MCUpdate* inupdate, int inweight = 1, string inname = "");
	void Reset();

	double GetTotalTime() {return totaltime;}
	double GetMeanTimePerCycle() {return ncycle ? totaltime / ncycle : 0;}
	double GetTotalCycleNumber() {return ncycle;}

	void ToStream(ostream& os, ostream& osdetail);

		void OpenLoop(int n) {
			ostringstream oss;
			oss << '(' << n << ":";
			command += oss.str();
		}
		void CloseLoop() {
			command += ')';
		}
		vector<int> ReadCommand(unsigned int &n);


	protected:

	vector<MCUpdate*> update;

	vector<int> weight;
	vector<double> time;
	vector<double> success;
	vector<string> name;
	vector<int> ncall;

	double totalweight;
	double totaltime;
	int ncycle;
	int size;

		string command;

	ProbModel* model;

	bool closed;
	bool random;
};


// The simplest MCUpdate object:
// when its Move() function is called
// it calls the Move() function of the random variable it is associated to
//
class SimpleMove : public MCUpdate	{

	public:

	SimpleMove(MCMC* invar, double intuning) : var(invar), tuning(intuning) {}

	double Move(double tuning_modulator = 1)	{
		return var->Move(tuning * tuning_modulator);
	}

	protected:

	MCMC* var;
	double tuning;
};


class JointSimpleMove : public MCUpdate, public Mnode {

	public:

	JointSimpleMove(Rnode* ina1, Rnode* ina2, double intuning);

	double Move(double tuning_modulator = 1);

	private:

	Rnode* a1;
	Rnode* a2;
	double tuning;

};

template <class P, class S> class SemiConjugateMove : public MCUpdate	{

	public:

	SemiConjugateMove(P* inprior, S* insampling, double intuningprior, int innprior, double intuningsampling, int innsampling) : prior(inprior), sampling(insampling), tuningprior(intuningprior), nprior(innprior), tuningsampling(intuningsampling), nsampling(innsampling) {}

	double Move(double tuning_modulator = 1)	{
		prior->ActivateSufficientStatistic();
		double total = 0;
		for (int i=0; i<nprior; i++)	{
			total += prior->Move(tuningprior * tuning_modulator);
		}
		for (int i=0; i<nsampling; i++)	{
			total += sampling->Move(tuningsampling * tuning_modulator);
		}
		total /= nsampling + nprior;
		prior->InactivateSufficientStatistic();
		return total;
	}

	protected:

	P* prior;
	S* sampling;
	double tuningprior;
	int nprior;
	double tuningsampling;
	int nsampling;
};

template <class P, class S> class ConjugateMove : public MCUpdate	{

	public:

	ConjugateMove(P* inprior, S* insampling, double intuning, int inn) : prior(inprior), sampling(insampling), tuning(intuning), n(inn) {}

	double Move(double tuning_modulator = 1)	{
		prior->Integrate();
		double total = 0;
		for (int i=0; i<n; i++)	{
			total += sampling->Move(tuning * tuning_modulator);
		}
		total /= n;
		prior->Resample();
		return total;
	}

	protected:

	P* prior;
	S* sampling;
	double tuning;
	int n;

};

// Compensatory compensatory move functions
class MultiplicativeCompensatoryMove : public MCUpdate, public Mnode {

	public:

	MultiplicativeCompensatoryMove(Multiplicative* inm1, Multiplicative* inm2, double intuning);

	double Move(double tuning_modulator = 1);

	private:

	Multiplicative* m1;
	Multiplicative* m2;
	double tuning;

};

// Compensatory compensatory move functions
class AdditiveCompensatoryMove : public MCUpdate, public Mnode {

	public:

	AdditiveCompensatoryMove(Additive* ina1, Additive* ina2, double intuning);

	double Move(double tuning_modulator = 1);

	private:

	Additive* a1;
	Additive* a2;
	double tuning;

};

// Compensatory compensatory move functions
class AdditiveAntiCompensatoryMove : public MCUpdate, public Mnode {

	public:

	AdditiveAntiCompensatoryMove(Additive* ina1, Additive* ina2, double intuning);

	double Move(double tuning_modulator = 1);

	private:

	Additive* a1;
	Additive* a2;
	double tuning;

};

class RealVectorMove : public MCUpdate	{

	public:

	RealVectorMove(Rvar<RealVector>* invar, double intuning, int inm);

	double Move(double tuning_modulator = 1);

	private:

	Rvar<RealVector>* var;
	double tuning;
	int m;
};

class RealVectorTranslationMove : public MCUpdate	{

	public:

	RealVectorTranslationMove(Rvar<RealVector>* invar, double intuning);

	double Move(double tuning_modulator = 1);

	private:

	Rvar<RealVector>* var;
	double tuning;
};

class RealVectorComponentwiseCompensatoryMove : public MCUpdate, public Mnode {

	public:

	RealVectorComponentwiseCompensatoryMove(Rvar<RealVector>* ina1, Rvar<RealVector>* ina2, double intuning);

	double Move(double tuning_modulator = 1);

	private :

	Rvar<RealVector>* a1;
	Rvar<RealVector>* a2;
	double tuning;

};

/*
class OneToManyRealVectorComponentwiseCompensatoryMove : public MCUpdate, public Mnode {

	public:

	OneToManyRealVectorComponentwiseCompensatoryMove(Rvar<RealVector>* ina1, IIDNormal** ina2, int inK, double intuning);

	double Move(double tuning_modulator = 1);

	private :

	Rvar<RealVector>* a1;
	IIDNormal** a2;
	int K;
	double tuning;

};
*/

class ProfileMove : public MCUpdate	{

	public:

	ProfileMove(Rvar<Profile>* invar, double intuning, int inn);

	double Move(double tuning_modulator = 1);

	private:

	Rvar<Profile>* var;
	double tuning;
	int n;
};


#endif
