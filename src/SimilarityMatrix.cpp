
#include "SimilarityMatrix.h"
#include "BiologicalSequences.h"
#include <iostream>
using namespace std;

//*******************************************//
//	Similarity Matrix
//******************************************//
SimilarityMatrix :: SimilarityMatrix(){

	omegaMatrix = new bool*[Naa];
	Groups= new int[Naa];
	for(int i=0; i< Naa; i++){
		omegaMatrix[i] = new bool[Naa];
		Groups[i] = 0;
	}
	for(int i=0 ; i<Naa; i++){
		for (int j=0; j<Naa; j++){
			omegaMatrix[i][j] = false;
		}
	}
}
void SimilarityMatrix :: Affiche (){
		for (int i=0; i<Naa; i++){ 
			cout << AminoAcids[i] << " ";
			for(int j=0; j<Naa; j++){
				cout  << "|" << omegaMatrix[i][j] <<"|";
			}
			cout<< "\n";
		}
	}

void SimilarityMatrix :: initGroups(const int N1,const int* group1, const int N2,const int* group2, const int N3, const int* group3, const int N4, const int* group4){
	if(group1){
		for(int i=0; i< N1; i++)
			Groups[group1[i]] = 1;
	}
	if(group2){
		for(int i=0; i< N2; i++)
			Groups[group2[i]] = 2;
	}
	if(group3){
		for(int i=0; i< N3; i++)
			Groups[group3[i]] = 3;
	}
	if(group4){
		for(int i=0; i< N4; i++)
			Groups[group1[i]] = 4;
	}
}

void SimilarityMatrix :: initOmegaMatrix(){
	for(int i=0; i<Naa; i++){
		for(int j=0; j<Naa; j++){
			if(!isConservative(i,j))
				omegaMatrix[i][j] = true;
		}
	}
}

bool SimilarityMatrix :: isConservative (int i, int j){
	return (Groups[i]==Groups[j]);
}


//******************************************************//
//		Polarity-Based Model
//*****************************************************//
// 0--> polar amino acids 1--> nonpolar amino acids

PolarityBasedMatrix :: PolarityBasedMatrix(){
	initGroups(NNonpolarSet, NonpolarSet, 0,0,0,0,0,0);
	initOmegaMatrix();
}

//******************************************************//
//		Volume-Based Model
//*****************************************************//
// 0--> Large amino acids 1--> Small amino acids

VolumeBasedMatrix :: VolumeBasedMatrix(){
	initGroups(NSmallSet, SmallSet,0,0,0,0,0,0);
	initOmegaMatrix();
}


//****************************************************************************************************//
//					Charge-Based Model
//****************************************************************************************************//
// 0--> No Charged amino acids  1--> Positively charged amino acids 2--> Negatively charged amino acids

ChargeBasedMatrix :: ChargeBasedMatrix(){
	initGroups(NPositivelyChargedSet,PositivelyChargedSet,NNegativelyChargedSet,NegativelyChargedSet,0,0,0,0);
	initOmegaMatrix();
}

//************************************************************************//
//		Polarity And Volume-Based Model
//***********************************************************************//
// 0--> Large and polar amino acids  1--> Small and polar amino acids
// 2--> Small and nonpolar amino acids 3--> Large and nonpolar amino acids

PolarityAndVolumeBasedMatrix :: PolarityAndVolumeBasedMatrix(){
	initGroups(NSmallPolarSet,SmallPolarSet,NSmallNonpolarSet,SmallNonpolarSet,NLargeNonpolarSet,LargeNonpolarSet, 0,0);
	initOmegaMatrix();
}
