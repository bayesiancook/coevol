#ifndef TRIPLETALI_H
#define TRIPLETALI_H

#include "SequenceAlignment.h"
#include "TripletStateSpace.h"

class TripletSequenceAlignment : public SequenceAlignment	{

	public:

	TripletSequenceAlignment(SequenceAlignment* from)	{

		if (from->GetNstate() != Nnuc)	{
			cerr << "error : triplets only on dna alignments\n";
			exit(1);
		}
		Ntaxa = from->Ntaxa;
		Nsite = from->Nsite / 3;
		taxset = from->taxset;
		tripstatespace = new TripletStateSpace();
		statespace = tripstatespace;

		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				int tmp = 0;
				if (3*j+3 < from->Nsite)	{
					tmp = from->Data[i][3*j+3];
				}
				else	{
					tmp = -1;
				}
				Data[i][j] = tripstatespace->GetState(from->Data[i][3*j+1], from->Data[i][3*j+2], tmp);
			}
		}
	}

	protected:

	TripletStateSpace* tripstatespace;
};


#endif

