#ifndef MSCODONSUBMATRIX_H
#define MSCODONSUBMATRIX_H

#include "CodonSubMatrix.hpp"

// The Muse and Gaut codon substitution process
// with an fitness = exp(selection profile) parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGFitnessCodonSubMatrix : public MGCodonSubMatrix	{

public:

  MGFitnessCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* infitness, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
    fitness(infitness) {}


  double	GetFitness(int aastate) {
    return fitness[aastate] + 1e-10;
  }

protected:

  // look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
  void      ComputeArray(int state);
  void      ComputeStationary();
  void			SetFitnessProfile(double* infitness) {fitness = infitness;}

  // data members

  double* fitness;
};

class RandomMGFitnessCodonSubMatrix : public MGFitnessCodonSubMatrix, public RandomCodonSubMatrix	{

public:

  RandomMGFitnessCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inRandomFitnessProfile, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGFitnessCodonSubMatrix(instatespace,inmatrix,inRandomFitnessProfile->GetArray(), innormalise) ,
    RandomSubMatrix(instatespace->GetNstate(), innormalise),
    RandomCodonSubMatrix(instatespace, innormalise),
    matrix(inmatrix),
    RandomFitnessProfile(inRandomFitnessProfile) {

    Register(matrix);
    Register(RandomFitnessProfile);
    specialUpdate();

  }

  // before updating the matrix instant rates and stationary probabilities
  // we need to check that we are pointing on the right nucleotide mutation process
  // and we need to update the value of fitness profile stored by the MGFitnessCodonSubMatrix object,
  // based on the value stored by the RandomFitnessprofile parent
  void SetParameters()	{
    SetNucMatrix(matrix);
    SetFitnessProfile(RandomFitnessProfile->GetArray());
  }

protected:

  RandomSubMatrix* matrix;
  Var<Profile>* RandomFitnessProfile;
};


//if selection profiles and codonusageselection is of type Dirichlet
class MGFitnessCodonUsageSubMatrix : public MGCodonSubMatrix	{

public:

  MGFitnessCodonUsageSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* infitness, double* incodonusageselection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
    fitness(infitness),
    codonusageselection(incodonusageselection)	{}


  double	GetFitness(int aastate) {
    return fitness[aastate] + 1e-10;
  }

  double	GetCodonUsageSelection(int codonstate) {
    return codonusageselection[codonstate];
  }


protected:

  // look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
  void		ComputeArray(int state);
  void		ComputeStationary();

  void		SetFitnessProfile(double* infitness) {fitness = infitness;}
  void		SetCodonUsageSelection(double* incodonusageselection) {codonusageselection = incodonusageselection;}

  // data members
  double* fitness;
  double* codonusageselection;
};




class RandomMGFitnessCodonUsageSubMatrix : public MGFitnessCodonUsageSubMatrix, public RandomCodonSubMatrix	{

public:

  RandomMGFitnessCodonUsageSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inRandomFitnessProfile, Var<Profile>* inRandomCodonUsageSelection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGFitnessCodonUsageSubMatrix(instatespace,inmatrix,inRandomFitnessProfile->GetArray(), inRandomCodonUsageSelection->GetArray(), innormalise),
    RandomSubMatrix(instatespace->GetNstate(), innormalise),
    RandomCodonSubMatrix(instatespace, innormalise),
    matrix(inmatrix),
    RandomFitnessProfile(inRandomFitnessProfile),
    RandomCodonUsageSelection(inRandomCodonUsageSelection)  {

    Register(matrix);
    Register(RandomFitnessProfile);
    Register(RandomCodonUsageSelection);

    specialUpdate();

  }

  CodonStateSpace*	GetCodonStateSpace() {return RandomCodonSubMatrix::GetCodonStateSpace();}

  // before updating the matrix instant rates and stationary probabilities
  // we need to check that we are pointing on the right nucleotide mutation process
  // and we need to update the value of fitness profile stored by the MGFitnessCodonSubMatrix object,
  // based on the value stored by the RandomFitnessprofile parent
  void SetParameters()	{
    SetNucMatrix(matrix);
    SetFitnessProfile(RandomFitnessProfile->GetArray());
    SetCodonUsageSelection(RandomCodonUsageSelection->GetArray());
  }

protected:

  RandomSubMatrix* matrix;
  Var<Profile>* RandomFitnessProfile;
  Var<Profile>* RandomCodonUsageSelection;

};


//if selection profiles and codonusageselection is of type Normal or gamma
class MGSRFitnessNormalCodonUsageSubMatrix : public MGCodonSubMatrix	{

public:

  MGSRFitnessNormalCodonUsageSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* infitness, double* incodonusageselection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
    fitness(infitness),
    codonusageselection(incodonusageselection) {}


  double	GetFitness(int aastate) {
    return fitness[aastate] + 1e-10;
  }

  double	GetCodonUsageSelection(int codonstate) {
    return codonusageselection[codonstate];
  }


protected:

  // look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
  void			ComputeArray(int state);
  void			ComputeStationary();
  void			SetFitnessProfile(double* infitness) {fitness = infitness;}
  void			SetCodonUsageSelection(double* incodonusageselection) {codonusageselection = incodonusageselection;}

  // data members
  double* fitness;
  double* codonusageselection;
};


class RandomMGSRFitnessNormalCodonUsageSubMatrix : public MGSRFitnessNormalCodonUsageSubMatrix, public RandomCodonSubMatrix	{

public:

  RandomMGSRFitnessNormalCodonUsageSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<RealVector>* inRandomFitnessProfile, Var<RealVector>* inRandomCodonUsageSelection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGSRFitnessNormalCodonUsageSubMatrix(instatespace,inmatrix,inRandomFitnessProfile->GetArray(), inRandomCodonUsageSelection->GetArray(), innormalise),
    RandomSubMatrix(instatespace->GetNstate(), innormalise),
    RandomCodonSubMatrix(instatespace, innormalise),
    matrix(inmatrix),
    RandomFitnessProfile(inRandomFitnessProfile),
    RandomCodonUsageSelection(inRandomCodonUsageSelection)  {

    Register(matrix);
    Register(RandomFitnessProfile);
    Register(RandomCodonUsageSelection);
    specialUpdate();



    //cerr << "nucmatrix= "  << *matrix << '\n';
    //cerr << "stationary= " << *(RandomFitnessProfile->GetArray()) << '\n';
    //cerr << "codonusage= " << *(RandomCodonUsageSelection->GetArray()) <<'\n';

  }

  CodonStateSpace*	GetCodonStateSpace() {return RandomCodonSubMatrix::GetCodonStateSpace();}

  // before updating the matrix instant rates and stationary probabilities
  // we need to check that we are pointing on the right nucleotide mutation process
  // and we need to update the value of fitness profile stored by the MGFitnessCodonSubMatrix object,
  // based on the value stored by the RandomFitnessprofile parent
  void SetParameters()	{
    SetNucMatrix(matrix);
    SetFitnessProfile(RandomFitnessProfile->GetArray());
    SetCodonUsageSelection(RandomCodonUsageSelection->GetArray());
  }

protected:

  RandomSubMatrix* matrix;
  Var<RealVector>* RandomFitnessProfile;
  Var<RealVector>* RandomCodonUsageSelection;

};




//if selection profiles and codonusageselection is of type Normal or gamma
class MGMSFitnessNormalCodonUsageSubMatrix : public MGCodonSubMatrix	{

public:

  MGMSFitnessNormalCodonUsageSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* infitness, double* incodonusageselection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
    fitness(infitness),
    codonusageselection(incodonusageselection) {}


  double	GetFitness(int aastate) {
    return fitness[aastate];
  }

  double	GetCodonUsageSelection(int codonstate) {
    return codonusageselection[codonstate];
  }

protected:

  // look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
  void			ComputeArray(int state);
  void			ComputeStationary();
  void			SetFitnessProfile(double* infitness) {fitness = infitness;}
  void			SetCodonUsageSelection(double* incodonusageselection) {codonusageselection = incodonusageselection;}

  // data members
  double* fitness;
  double* codonusageselection;
};


class RandomMGMSFitnessNormalCodonUsageSubMatrix : public MGMSFitnessNormalCodonUsageSubMatrix, public RandomCodonSubMatrix	{

public:

  RandomMGMSFitnessNormalCodonUsageSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<RealVector>* inRandomFitnessProfile, Var<RealVector>* inRandomCodonUsageSelection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGMSFitnessNormalCodonUsageSubMatrix(instatespace,inmatrix,inRandomFitnessProfile->GetArray(), inRandomCodonUsageSelection->GetArray(), innormalise),
    RandomSubMatrix(instatespace->GetNstate(), innormalise),
    RandomCodonSubMatrix(instatespace, innormalise),
    matrix(inmatrix),
    RandomFitnessProfile(inRandomFitnessProfile),
    RandomCodonUsageSelection(inRandomCodonUsageSelection)  {

    Register(matrix);
    Register(RandomFitnessProfile);
    Register(RandomCodonUsageSelection);
    specialUpdate();



    //cerr << "nucmatrix= "  << *matrix << '\n';
    //cerr << "stationary= " << *(RandomFitnessProfile->GetArray()) << '\n';
    //cerr << "codonusage= " << *(RandomCodonUsageSelection->GetArray()) <<'\n';

  }

  CodonStateSpace*	GetCodonStateSpace() {return RandomCodonSubMatrix::GetCodonStateSpace();}

  // before updating the matrix instant rates and stationary probabilities
  // we need to check that we are pointing on the right nucleotide mutation process
  // and we need to update the value of fitness profile stored by the MGFitnessCodonSubMatrix object,
  // based on the value stored by the RandomFitnessprofile parent
  void SetParameters()	{
    SetNucMatrix(matrix);
    SetFitnessProfile(RandomFitnessProfile->GetArray());
    SetCodonUsageSelection(RandomCodonUsageSelection->GetArray());
  }

protected:

  RandomSubMatrix* matrix;
  Var<RealVector>* RandomFitnessProfile;
  Var<RealVector>* RandomCodonUsageSelection;

};



//square root
class MGSRFitnessCodonUsageSubMatrix : public MGCodonSubMatrix	{

public:

  MGSRFitnessCodonUsageSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* infitness, double* incodonusageselection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
    fitness(infitness),
    codonusageselection(incodonusageselection) {}


  double	GetFitness(int aastate) {
    return fitness[aastate] + 1e-6;
  }

  double	GetCodonUsageSelection(int codonstate) {
    return codonusageselection[codonstate];
  }


protected:

  // look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp

  void		ComputeArray(int state);
  void		ComputeStationary();
  void		SetFitnessProfile(double* infitness) {fitness = infitness;}
  void		SetCodonUsageSelection(double* incodonusageselection) {codonusageselection = incodonusageselection;}

  // data members
  double* fitness;
  double* codonusageselection;
};



//if selection profiles and codonusageselection is of type Dirichlet
class RandomMGSRFitnessCodonUsageSubMatrix : public MGSRFitnessCodonUsageSubMatrix, public RandomCodonSubMatrix	{

public:

  RandomMGSRFitnessCodonUsageSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inRandomFitnessProfile, Var<Profile>* inRandomCodonUsageSelection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGSRFitnessCodonUsageSubMatrix(instatespace,inmatrix,inRandomFitnessProfile->GetArray(), inRandomCodonUsageSelection->GetArray(), innormalise),
    RandomSubMatrix(instatespace->GetNstate(), innormalise),
    RandomCodonSubMatrix(instatespace, innormalise),
    matrix(inmatrix),
    RandomFitnessProfile(inRandomFitnessProfile),
    RandomCodonUsageSelection(inRandomCodonUsageSelection)  {

    Register(matrix);
    Register(RandomFitnessProfile);
    Register(RandomCodonUsageSelection);

    specialUpdate();

  }

  CodonStateSpace*	GetCodonStateSpace() {return RandomCodonSubMatrix::GetCodonStateSpace();}

  // before updating the matrix instant rates and stationary probabilities
  // we need to check that we are pointing on the right nucleotide mutation process
  // and we need to update the value of fitness profile stored by the MGFitnessCodonSubMatrix object,
  // based on the value stored by the RandomFitnessprofile parent
  void SetParameters()	{
    SetNucMatrix(matrix);
    SetFitnessProfile(RandomFitnessProfile->GetArray());
    SetCodonUsageSelection(RandomCodonUsageSelection->GetArray());
  }

protected:

  RandomSubMatrix* matrix;
  Var<Profile>* RandomFitnessProfile;
  Var<Profile>* RandomCodonUsageSelection;

};



//mutation selection
class MGMSFitnessCodonUsageSubMatrix : public MGCodonSubMatrix	{

public:

  MGMSFitnessCodonUsageSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double* infitness, double* incodonusageselection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
    fitness(infitness),
    codonusageselection(incodonusageselection) {}


  double	GetFitness(int aastate) {
    return fitness[aastate] + 1e-6;
  }

  double	GetCodonUsageSelection(int codonstate) {
    return codonusageselection[codonstate];
  }


protected:

  // look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp

  void		ComputeArray(int state);
  void		ComputeStationary();
  void		SetFitnessProfile(double* infitness) {fitness = infitness;}
  void		SetCodonUsageSelection(double* incodonusageselection) {codonusageselection = incodonusageselection;}

  void ToStream(ostream& os)	{
    cerr << "nucmatrix : \n";
    NucMatrix->CheckReversibility();
    cerr << '\n';
    NucMatrix->CorruptMatrix();
    NucMatrix->UpdateMatrix();
    NucMatrix->ToStream(cerr);
    cerr << '\n';
    NucMatrix->CheckReversibility();
    cerr << '\n';

    for (int i=0; i<Naa; i++)	{
      cerr << fitness[i] << '\t';
    }
    cerr << '\n';
    cerr << '\n';
    double total = 0;
    for (int i=0; i<GetNstate(); i++)	{
      total += NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i)) * GetFitness(GetCodonStateSpace()->Translation(i)) * GetCodonUsageSelection(i);
    }
    for (int i=0; i<GetNstate(); i++)	{
      cerr << i << '\t' << NucMatrix->Stationary(GetCodonPosition(0,i));
      cerr << '\t' << NucMatrix->Stationary(GetCodonPosition(1,i));
      cerr << '\t' << NucMatrix->Stationary(GetCodonPosition(2,i));
      cerr << '\t' << GetFitness(GetCodonStateSpace()->Translation(i));
      cerr << '\t' << GetCodonUsageSelection(i);
      cerr << '\t' << NucMatrix->Stationary(GetCodonPosition(0,i)) * NucMatrix->Stationary(GetCodonPosition(1,i)) * NucMatrix->Stationary(GetCodonPosition(2,i)) * GetFitness(GetCodonStateSpace()->Translation(i)) * GetCodonUsageSelection(i) / total;
      cerr << '\n';
    }
    cerr << '\n';
    for (int i=0; i<GetNstate(); i++)	{
      for (int j=0; j<GetNstate(); j++)	{
        if (i!=j)	{
          int pos = GetDifferingPosition(i,j);
          if ((pos != -1) && (pos != 3))	{
            int a = GetCodonPosition(pos,i);
            int b = GetCodonPosition(pos,j);

            double mut = (*NucMatrix)(a,b);

            double deltaS;
            if (! Synonymous(i,j))  {
              deltaS = ( log(GetFitness(GetCodonStateSpace()->Translation(j))) - log(GetFitness(GetCodonStateSpace()->Translation(i))) ) + ( log(GetCodonUsageSelection(j)) - log(GetCodonUsageSelection(i)) );
            }
            else	{
              deltaS = log(GetCodonUsageSelection(j)) - log(GetCodonUsageSelection(i));
            }
            double fix = 1;
            if (deltaS != 0)	{
              fix = deltaS/(1.0 - exp(-deltaS));
            }
            cerr << i << '\t' << j << '\t' << mut << '\t' << deltaS << '\t' << fix << '\t' << mut*fix << '\t' << Q[i][j] << '\n';
          }
        }
      }
    }
  }

  // data members
  double* fitness;
  double* codonusageselection;
};



//if selection profiles and codonusageselection is of type Dirichlet
class RandomMGMSFitnessCodonUsageSubMatrix : public MGMSFitnessCodonUsageSubMatrix, public RandomCodonSubMatrix	{

public:

  RandomMGMSFitnessCodonUsageSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inRandomFitnessProfile, Var<Profile>* inRandomCodonUsageSelection, bool innormalise = false) :
    SubMatrix(instatespace->GetNstate(), innormalise),
    CodonSubMatrix(instatespace, innormalise),
    MGMSFitnessCodonUsageSubMatrix(instatespace,inmatrix,inRandomFitnessProfile->GetArray(), inRandomCodonUsageSelection->GetArray(), innormalise),
    RandomSubMatrix(instatespace->GetNstate(), innormalise),
    RandomCodonSubMatrix(instatespace, innormalise),
    matrix(inmatrix),
    RandomFitnessProfile(inRandomFitnessProfile),
    RandomCodonUsageSelection(inRandomCodonUsageSelection)  {

    Register(matrix);
    Register(RandomFitnessProfile);
    Register(RandomCodonUsageSelection);

    specialUpdate();

  }

  CodonStateSpace*	GetCodonStateSpace() {return RandomCodonSubMatrix::GetCodonStateSpace();}

  // before updating the matrix instant rates and stationary probabilities
  // we need to check that we are pointing on the right nucleotide mutation process
  // and we need to update the value of fitness profile stored by the MGFitnessCodonSubMatrix object,
  // based on the value stored by the RandomFitnessprofile parent
  void SetParameters()	{
    SetNucMatrix(matrix);
    SetFitnessProfile(RandomFitnessProfile->GetArray());
    SetCodonUsageSelection(RandomCodonUsageSelection->GetArray());
  }

protected:

  RandomSubMatrix* matrix;
  Var<Profile>* RandomFitnessProfile;
  Var<Profile>* RandomCodonUsageSelection;

};


#endif
