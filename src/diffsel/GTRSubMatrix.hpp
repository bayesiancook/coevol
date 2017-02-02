#ifndef GTRSUBMATRIX_H
#define GTRSUBMATRIX_H

#include "RandomSubMatrix.hpp"
#include "BiologicalSequences.hpp"

class GTRSubMatrix : public virtual SubMatrix	{

public:

  GTRSubMatrix(int nstate, const double* rr, const double* stat, bool innormalise = false);
  ~GTRSubMatrix() {};

  int			GetNRelativeRate() {return Nrr;}
  double      RelativeRate(int i, int j) {return mRelativeRate[rrindex(i,j,GetNstate())];}

  static int rrindex(int i, int j, int nstate)	{
    return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
  }

protected:

  // make a copy of the entries (not of the pointer)
  void      CopyStationary(const double* instat);

  // copy the pointer
  void			SetRelativeRate(const double* inrelrate) {mRelativeRate = inrelrate;}

protected:

  void      ComputeArray(int state);
  void      ComputeStationary() {}

  // data members
  const double* mRelativeRate;
  int Nrr;
};

class GTRRandomSubMatrix : public RandomSubMatrix, public GTRSubMatrix  {

public:
  GTRRandomSubMatrix(Var<PosRealVector>* inrelrate, Var<Profile>* instat, bool innormalise = false) : SubMatrix(instat->GetDim(), innormalise), RandomSubMatrix(instat->GetDim(), innormalise) , GTRSubMatrix(instat->GetDim(),inrelrate->val().GetArray(),instat->val().GetArray(), innormalise)	{
    relrate = inrelrate;
    stat = instat;
    Register(relrate);
    Register(stat);
    specialUpdate();
  }

  Var<Profile>* GetRandomStationaries() {return stat;}
  Var<PosRealVector>* GetRandomRelativeRates() {return relrate;}

protected:

  void SetParameters()	{
    SetRelativeRate(relrate->val().GetArray());
    CopyStationary(stat->val().GetArray());
  }

private:
  Var<PosRealVector>* relrate;
  Var<Profile>* stat;
};

class GTRRandomSubMatrixWithNormRates : public RandomSubMatrix, public GTRSubMatrix  {

public:
  GTRRandomSubMatrixWithNormRates(Var<Profile>* inrelrate, Var<Profile>* instat, bool innormalise = false) : SubMatrix(instat->GetDim(), innormalise), RandomSubMatrix(instat->GetDim(), innormalise) , GTRSubMatrix(instat->GetDim(),inrelrate->val().GetArray(),instat->val().GetArray(), innormalise)	{
    relrate = inrelrate;
    stat = instat;
    rescaledrelrate = new double[relrate->GetDim()];
    Register(relrate);
    Register(stat);
    specialUpdate();
  }

  ~GTRRandomSubMatrixWithNormRates()	{
    delete[] rescaledrelrate;
  }

  /*
    Var<Profile>* GetRandomStationaries() {return stat;}
    Var<Profile>* GetRandomRelativeRates() {return relrate;}
  */

protected:

  void SetParameters()	{
    for (int i=0; i<relrate->GetDim(); i++)	{
      //	rescaledrelrate[i] = (*relrate)[i] * relrate->GetDim();
      rescaledrelrate[i] = (*relrate)[i];
    }
    SetRelativeRate(rescaledrelrate);
    CopyStationary(stat->val().GetArray());
  }

private:
  Var<Profile>* relrate;
  Var<Profile>* stat;
  double* rescaledrelrate;
};

const double LG_RR[] = {2.43501,0.386559,1.01598,0.248188,2.02114,0.35106,0.146574,0.524859,0.386747,1.09961,0.270803,1.15206,0.94882,0.415855,4.62446,2.09301,2.49251,0.176789,0.214201,0.0611966,0.00342322,1.08123,0.556892,0.626623,0.313658,0.0129776,0.581098,0.874258,0.517278,0.0737435,0.0829654,0.522935,2.72396,1.11863,1.91672,0.655562,1.1402,5.12993,0.0170374,0.826565,0.906972,0.010458,0.27681,0.014748,0.0249932,4.96585,0.385885,0.512014,0.121261,1.21333,0.416606,0.0371423,0.0292405,0.132172,0.0184023,0.341267,0.414673,0.0433031,1.7679,0.0681587,0.16996,0.529941,0.410296,4.03888,0.356061,0.598676,0.591407,0.23971,0.0761602,0.117429,0.0876389,0.667316,1.08855,0.0233983,2.53635,1.75976,0.0875798,0.0924114,0.035076,0.0515759,0.353957,0.161416,0.640457,2.40372,7.63435,0.304716,0.00851629,0.290191,0.0432987,0.136505,1.40641,0.192681,0.262136,0.381714,1.70218,0.127015,0.0750342,0.262657,0.0534908,0.106516,0.682111,0.358356,0.432853,4.41125,0.497794,4.70891,2.37387,0.968501,0.571566,0.116427,0.584078,5.19151,0.155611,4.055,4.18072,0.187342,0.0765801,0.0712711,0.124231,0.0627125,1.01127,10.4177,0.109234,0.22747,0.134512,0.642335,2.09847,0.381839,3.16401,6.1886,0.732413,1.11216,0.181177,0.048821,0.129065,6.17517,0.0669398,0.243648,0.5698,0.295289,0.178326,0.296353,1.66574,0.606165,0.293136,0.362942,0.0976794,1.63622,0.473613,0.339422,1.97647,1.85746,0.681046,0.470848,0.158271,1.6589,0.73554,3.92125,1.95721,0.081869,0.0443899,0.598725,0.610728,0.325308,1.30906,0.559051,0.290058,0.0930628,0.0876663,2.74689,1.19724,1.05666,0.205761,0.231066,0.251745,0.839504,0.566405,0.167174,0.580707,0.307606,6.33164,0.096231,0.243453,0.391844,2.14061,0.137765,0.240498,0.185392,0.243896,3.08335};


class LGSubMatrix : public GTRSubMatrix	{

public:

  LGSubMatrix(const double* stat, bool innormalise = false) : SubMatrix(Naa,innormalise), GTRSubMatrix(Naa,0,stat,innormalise)	{
    mRelativeRate = LG_RR;
  }
};

class LGRandomSubMatrix : public RandomSubMatrix, public LGSubMatrix  {

public:
  LGRandomSubMatrix(Var<Profile>* instat, bool innormalise = false) : SubMatrix(instat->GetDim(), innormalise), RandomSubMatrix(instat->GetDim(), innormalise) , LGSubMatrix(instat->val().GetArray(), innormalise)	{
    stat = instat;
    Register(stat);
    specialUpdate();
  }

protected:

  void SetParameters()	{
    CopyStationary(stat->val().GetArray());
  }

private:
  Var<Profile>* stat;

};

#endif
