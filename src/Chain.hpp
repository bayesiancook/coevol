#ifndef CHAIN_H
#define CHAIN_H

#include <fstream>
#include <cstdlib>
using namespace std;

#include "ProbModel.h"

/// Chain is a Monte Carlo Markov Chain
//  it is responsible for creating a model, applying it to data
//  running a MCMC, to obtain a sample approximately from the posterior distribution
//  saving the sample onto a file, restarting a run
//
// file nomenclature:
// <chainname>.param   : current state
// <chainname>.chain   : list of all points since the beginning of the Monte Carlo (burnin included)
// <chainname>.trace   : trace file, each row corresponding to one point of the .chain file
// <chainname>.monitor : monitoring the success rate, time spent in each move, numerical errors, etc
// <chainname>.run     : put 0 in this file to stop the chain


class Chain	{

	public:

			Chain();

	virtual		~Chain() {};

	virtual void 	MakeFiles(int force = 0);
			// overwrites files if force == 1

	virtual void 	Monitor();
			// write the .trace and .monitor files

	virtual void 	SavePoint();
			// save one point in the .chain file

	virtual void 	Reset(int force = 0);
			// initialise model and make the files (overwrite if force == 1)

	virtual void 	Move();
			// perform one cycle of Monte Carlo "moves" (updates)

	virtual void 	Start();
			// start Monte Carlo

	virtual int 	GetRunningStatus();
			// returns 0 (means STOP) if one the following conditions holds true
			// 	.run file contains a 0 ("echo 0 > <chainname>.run" is the proper way to stop a chain from a shell)
			// 	size >= until (and until != -1)

	virtual void 	Run();
			// Move, Monitor amd Save while running status == 1


	virtual string 	GetModelType() = 0;
			// string meant as a check when opening files (model name)
			// you should give different names to chains based on different models

	virtual void 	New(int force = 0) = 0;
			// new chain (force = 1 : this will overwrite files)

	virtual void 	Open() = 0;
			// open a chain from files

	virtual void 	Save() = 0;
			// save a chain to files

	string		GetName() {return name;}

	void 		SetEvery(int inevery)	{every = inevery;}
	void		SetUntil(int inuntil) {until = inuntil;}

	ProbModel* 	GetModel() {return model;}
	int		GetSize() {return size;}


	protected:

	int every; 		// saving frequency
	int until;		// intended size of the run (number of saved points)
	int size;		// current size
	ProbModel* model;	// the model
	string name;		// the name of the chain in the filesystem
				// all files for this chain will be of the form : <name>.<ext>

};

#endif // CHAIN_H
