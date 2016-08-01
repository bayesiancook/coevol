
#ifndef RYSEQUENCEALIGNMENT_H
#define RYSEQUENCEALIGNMENT_H


#include "ContinuousData.h"
#include "SequenceAlignment.h"
#include <cmath>

class RYSequenceAlignment : public SequenceAlignment	{

	public:

	RYSequenceAlignment(RYSequenceAlignment* from) : SequenceAlignment((SequenceAlignment*) from) {}

	RYSequenceAlignment(SequenceAlignment* from)	{

		try	{
			DNAsource = from;

			Nsite = from->Nsite;
			Ntaxa = from->Ntaxa;
			statespace = new RYStateSpace();

			taxset = DNAsource->GetTaxonSet();

			// make my own arrays
			// make translation
			Data = new int*[Ntaxa];
			for (int i=0; i<Ntaxa; i++)	{
				Data[i] = new int[Nsite];
				for (int j=0; j<Nsite; j++)	{
					try {
						Data[i][j] = GetRYStateSpace()->GetRYCoding(DNAsource->GetState(i,j));
					}
					catch(...)	{
					// catch(Exception e)	{
						cerr << "in RYSequenceAlignment: taxon " << i << " and codon " << j << "\n";
						cerr << "taxon : " << taxset->GetTaxon(i) << '\n';
						throw;
					}
				}
			}

		}
		catch(Exception)	{
			cerr << "RY Sequence Alignment: failed to read the datafile\n";
			exit(1);
		}
	}

	RYStateSpace* GetRYStateSpace()	{
		return (RYStateSpace*) (statespace);
	}

	private:

	SequenceAlignment* DNAsource;

};

#endif
