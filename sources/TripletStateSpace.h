
#ifndef TRIPLET_H
#define TRIPLET_H

#include "StateSpace.h"

const char dna5[] = {'A','C','G','T','?'};

class TripletStateSpace : public StateSpace	{

	// 0 1 2 3 4 : A C G T ?
	// 0 : AAA
	// 1 : AAC
	// 2 : AAG
	// 3 : AAT
	// 4 : AA?
	// 5 : ACA
	// etc

	public	:


	int GetNstate()	{
		return 125;
	}

	bool isCompatible(int state125, int state64)	{
		int s1 = state125 / 25;
		int s2 = (state125 -25*s1) / 5;
		int s3 = state125 - 25*s1 - 5*s2;

		int t1 = state64 / 16;
		int t2 = (state64 -16*t1) / 4;
		int t3 = state64 - 16*t1 - 4*t2;

		return (((s1 == 4) || (s1 == t1)) && ((s2 == 4) || (s2 == t2)) && ((s3 == 4) || (s3 == t3)));
	}

	string GetState(int state)	{
		int s1 = state / 25;
		int s2 = (state -25*s1) / 5;
		int s3 = state - 25*s1 - 5*s2;
		ostringstream s;
		s << dna5[s1] << dna5[s2] << dna5[s3];
		return s.str();
	}

	int GetState(int s1, int s2, int s3)	{
		if (s1 == -1)	{
			s1 = 4;
		}
		if (s2 == -1)	{
			s2 = 4;
		}
		if (s3 == -1)	{
			s3 = 4;
		}
		return s1 * 25 + s2 * 5 + s3;
	}

	int GetState(string from)	{
		cerr << "get state in triplet\n";
		exit(1);
		return 0;
	}

};


#endif

