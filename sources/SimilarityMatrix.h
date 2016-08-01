#ifndef SIMILARITYMATRIX_H
#define SIMILARITYMATRIX_H

#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;

	enum AminoAcidModelType {POLARITY=0, VOLUME=1, CHARGE=2, POLARITYVOLUME=3, NONE = 4};

	const int NSmallSet = 9;
	const int SmallSet [] = {0,1,2,5,11,12,15,16,17};   
	const int NNonpolarSet = 8;
	const int NonpolarSet [] = {0,4,5,7,9,10,12,17};
	const int NPositivelyChargedSet = 3;
	const int PositivelyChargedSet [] = {6,8,14};
	const int NNegativelyChargedSet = 2;
	const int NegativelyChargedSet [] = {2,3};
	const int NLargePolarSet = 7;
	const int LargePolarSet [] = {3,6,8,13,14,18,19};
	const int NSmallPolarSet = 5;
	const int SmallPolarSet [] = {1,2,11,15,16};
	const int NSmallNonpolarSet = 4;
	const int SmallNonpolarSet [] = {0,5,12,17};
	const int NLargeNonpolarSet = 4;
	const int LargeNonpolarSet [] = {4,7,9,10};

class SimilarityMatrix{
	//create a 20*20 matrix of booleen 
	public :
	SimilarityMatrix();

	bool isRadical(int i, int j){ return omegaMatrix[i][j];}

	int GetDim()	{
		return 20;
	}

	//Fonction qui affiche la matrice OmegaMatrix
	void Affiche();

	protected :
	/*
	Depeding on the Model, each Amino Acid is assigned to a group
	0--> group1 1--> group2, 2--> group3, 3-->group4
	exemple: for the polarity-Based Model: Polar AA --> 0, Nonpolar AA -->1
	*/
	int* Groups;
	void initGroups(const int n1,const int group1[], const int n2, const int group2[], const int n3, const int group3[], const int n4, const int goup4[]);

	//  return true if i and j are in the sameGroup
	bool isConservative(int i,int j);

	/*
	a 20*20 matrix
	if omegaMatrix[i][j]= true --> going from i to j is radical, otherwise going from i to j is conservative
	*/
	bool** omegaMatrix;
	void initOmegaMatrix();
};

class PolarityBasedMatrix : public SimilarityMatrix{
	public:
	PolarityBasedMatrix();
};
class VolumeBasedMatrix : public SimilarityMatrix{
	public :
	VolumeBasedMatrix();
};
class ChargeBasedMatrix : public SimilarityMatrix {
	public:
	ChargeBasedMatrix ();
};
class PolarityAndVolumeBasedMatrix : public SimilarityMatrix{
	public:
	PolarityAndVolumeBasedMatrix();
};


inline istream& operator>>(istream& is, AminoAcidModelType & type)	{
	string t;
	is >> t;
	if (t == "NONE")	{
		type = NONE;
	}
	else if (t == "POLARITY")	{
		type = POLARITY;
	}
	else if (t == "CHARGE")	{
		type = CHARGE;
	}
	else if (t == "VOLUME")	{
		type = VOLUME;
	}
	else if (t == "POLARITYVOLUME")	{
		type = POLARITYVOLUME;
	}
	else	{
		cerr << "error in istream amino acid model type\n";
		cerr << type << '\n';
		exit(1);
	}
	return is;
}
inline ostream& operator<<(ostream& os, AminoAcidModelType type)	{
	if (type == NONE)	{
		os << "NONE\n";
	}
	else if (type == POLARITY)	{
		os << "POLARITY\n";
	}
	else if (type == CHARGE)	{
		os << "CHARGE\n";
	}
	else if (type == VOLUME)	{
		os << "VOLUME\n";
	}
	else if (type == POLARITYVOLUME)	{
		os << "POLARITYVOLUME\n";
	}
	else	{
		cerr << "error in ostream genetic code type\n";
		cerr << (int) type << '\n';
		exit(1);
	}
	return os;
}
#endif
