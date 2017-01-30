
#ifndef STATESPACE_H
#define STATESPACE_H

#include "BiologicalSequences.h"

// pure interface
//
class StateSpace	{

	public:

	virtual ~StateSpace() {}

	virtual string GetState(int state) = 0;
	virtual int GetNstate() = 0;

	virtual int GetState(string from) = 0;

	virtual bool isCompatible(int state1, int state2)	{
			return ((state1 == unknown) || (state2 == unknown) || (state1 == state2));
	}

};

// simple state space: assumes that states are referred to using a one-letter code
//
class SimpleStateSpace : public StateSpace	{

	public:


	int GetState(string from);

	int GetNstate() {
		return Nstate;
	}

	string GetState(int state);

	protected:
	int Nstate;
	char* Alphabet;
	int NAlphabetSet;
	char* AlphabetSet;
};

class DNAStateSpace : public SimpleStateSpace	{

	public:

	DNAStateSpace();
	~DNAStateSpace();
};

class RNAStateSpace : public SimpleStateSpace	{

	public:

	RNAStateSpace();
	~RNAStateSpace();
};

class ProteinStateSpace : public SimpleStateSpace	{

	public:

	ProteinStateSpace();
	~ProteinStateSpace();
};

class RYStateSpace : public SimpleStateSpace	{

	public:

	RYStateSpace();
	~RYStateSpace();

	int GetRYCoding(int from);
};

class GenericStateSpace : public SimpleStateSpace	{

	public:

	GenericStateSpace(int inNstate, char* inAlphabet, int inNAlphabetSet, char* inAlphabetSet)	{
		Nstate = inNstate;
		Alphabet = new char[Nstate];
		for (int i=0; i<Nstate; i++)	{
			Alphabet[i] = inAlphabet[i];
		}
		NAlphabetSet = inNAlphabetSet;
		AlphabetSet = new char[NAlphabetSet];
		for (int i=0; i<NAlphabetSet; i++)	{
			AlphabetSet[i] = inAlphabetSet[i];
		}
	}

	~GenericStateSpace()	{
		delete[] Alphabet;
		delete[] AlphabetSet;
	}

};

#endif // STATESPACE_H

