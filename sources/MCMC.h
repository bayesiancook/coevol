#ifndef MCMC_H
#define MCMC_H

#include "Exception.h"

/// An interface for Monte-carlo behavior
/**
 * is implemented by various components of the model such as
 * elemtary random variables (e.g. Rnode, Rvar<>)
 * entire modules made of several variables (e.g. PhyloProcess, ProbModel)
 */

#include <iostream>
using namespace std;


class MCMC {

	public:

	static const double MAXDIFF;

	MCMC() : clamp(false) {}

	virtual ~MCMC() {}

	/// resample current state, conditional on parents and children.
	// monte carlo sampling of a model proceeds by repeatedly calling the Move function of all its components taken in turn
	virtual double Move(double) = 0;

	// resample current state directly at equilibrium
	// random draws from the prior (i.e. simulations under the model)
	// proceed by recursively calling Sample, from the roots to the tips of the model
	virtual void	Sample()	{
		if (! isClamped())	{
			drawSample();
		}
	}

	virtual void 	drawSample() = 0;

	bool		isClamped() {return clamp;}
	void		Clamp() {clamp = true;}
	void		Unclamp() {clamp = false;}

	// returns the probability of the component given the current state of its parents
	virtual double GetLogProb() = 0;


	private:
	bool clamp;

};


// A more specific interface, meant for Metropolis Hastings moves

class MH : public MCMC	{

	public:
	virtual ~MH() {}

	virtual double ProposeMove(double tuning = 1) = 0;
	protected:
	virtual void Restore() = 0;

};

#endif // MCMC_H

