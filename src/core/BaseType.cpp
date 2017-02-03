#include "core/BaseType.hpp"

const double Profile::MIN = 1e-20;

double	Profile::ProposeMove(double tuning, int n)	{ // n==0dirichlet resampling, otherwise, vase communiquants
  double ret = 0;
  if (!n)	{ // dirichlet
    double* oldprofile = new double[dim];
    for (int i=0; i<dim; i++)	{
      oldprofile[i] = profile[i];
    }
    double total = 0;
    for (int i=0; i<dim; i++)	{
      profile[i] = Random::sGamma(tuning*oldprofile[i]);
      if (profile[i] == 0)	{
        std::cerr << "error in dirichlet resampling : 0 \n";
        exit(1);
      }
      total += profile[i];
    }

    double logHastings = 0;
    for (int i=0; i<dim; i++)	{
      profile[i] /= total;

      logHastings += - Random::logGamma(tuning*oldprofile[i]) + Random::logGamma(tuning*profile[i])
        -  (tuning*profile[i] -1.0) * log(oldprofile[i]) + (tuning * oldprofile[i] -1.0) * log(profile[i]);
    }

    delete[] oldprofile;
    return logHastings;
  }
  else	{
    if (2*n > dim)	{
      n = dim / 2;
    }
    int* indices = new int[2*n];
    Random::DrawFromUrn(indices,2*n,dim);
    for (int i=0; i<n; i++)	{
      int i1 = indices[2*i];
      int i2 = indices[2*i+1];
      double tot = profile[i1] + profile[i2];
      double x = profile[i1];

      // double h = tuning * (Random::Uniform() - 0.5);
      double h = tot * tuning * (Random::Uniform() - 0.5);
      /*
        int c = (int) (h / (2 * tot));
        h -= c*2*tot;
      */
      x += h;
      while ((x<0) || (x>tot))	{
        if (x<0)	{
          x = -x;
        }
        if (x>tot)	{
          x = 2*tot - x;
        }
      }
      profile[i1] = x;
      profile[i2] = tot - x;
    }
    delete[] indices;
  }
  return ret;
}

void Additive::Register(DAGnode*) {std::cerr << "error in Additive::Register\n"; throw;}

void Multiplicative::Register(DAGnode*)  {std::cerr << "error in Multiplicative::Register\n"; throw;}

Real& Real::operator=(const Real& from) {
  value = from.value;
  return *this;
}

Real& Real::operator=(double from) {
  value = from;
  return *this;
}

Real& Real::operator+=(const Real& from) {
  value += from.value;
  return *this;
}

Real& Real::operator/=(double from) {
  value /= from;
  return *this;
}

int Real::ScalarAddition(double d) {
  value += d;
  return 1;
}

double Real::ProposeMove(double tuning) {
  // simple additive move
  double m = tuning*(Random::Uniform() - 0.5);
  value += m;
  return 0;
}

int Real::Check() {return 1;}

UnitReal& UnitReal::operator=(const UnitReal& from) {
  value = from.value;
  return *this;
}

UnitReal& UnitReal::operator=(double from) {
  value = from;
  return *this;
}

UnitReal& UnitReal::operator+=(const UnitReal from) {
  value += from.value;
  return *this;
}

double UnitReal::ProposeMove(double tuning) {
  // simple additive move
  double m = tuning*(Random::Uniform() - 0.5);
  value += m;
  while ((value<0) || (value>1)) {
    if (value<0) {
      value = -value;
    }
    if (value>1) {
      value = 2 - value;
    }
  }
  return 0.0;
}

int UnitReal::Check() {return 1;}

PosReal& PosReal::operator=(const PosReal& from) {
  value = from.value;
  return *this;
}

PosReal& PosReal::operator=(const double& from) {
  value = from;
  return *this;
}

PosReal& PosReal::operator+=(const PosReal from) {
  value += from.value;
  return *this;
}

PosReal& PosReal::operator/=(double from) {
  value /= from;
  return *this;
}

int PosReal::ScalarMultiplication(double d) {
  value *= d;
  return 1;
}

double PosReal::ProposeMove(double tuning) {
  // simple multiplicative move
  double m = tuning*(Random::Uniform() - 0.5);
  value *= exp(m);
  return m;
}

int PosReal::Check() {
  if (value<=0) {
    std::cerr << "error : positive double is not positive : " << value << '\n';
    return 0;
  }
  return 1;
}
