#include "core/BaseType.hpp"
#include <cmath>

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

std::istream& operator>>(std::istream& is, Real& r) {
  is >> r.value;
  return is;
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

std::istream& operator>>(std::istream& is, UnitReal& r) {
  is >> r.value;
  return is;
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

std::istream& operator>>(std::istream& is, PosReal& r) {
  is >> r.value;
  return is;
}

Profile::Profile(int indim, double* v) {
  dim = indim;
  profile = new double[dim];
  if (v) {
    double total = 0;
    for (int k=0; k<dim; k++) {
      if (!(v[k]>0)) {
        std::cerr << "error : profiles should be strictly positive\n";
        exit(1);
      }
      profile[k] = v[k];
      total += profile[k];
    }
    for (int k=0; k<dim; k++) {
      profile[k] /= total;
    }
  }
}

Profile::Profile(const Profile& from) {
  dim = from.dim;
  profile = new double[dim];
  for (int k=0; k<dim; k++) {
    profile[k] = from.profile[k];
  }
}

Profile::~Profile() {
  delete[] profile;
}

Profile& Profile::operator=(const Profile& from) {
  if (!dim) {
    dim = from.dim;
    profile = new double[dim];
  }
  if (dim != from.dim) {
    std::cerr << "error : non matching dimenstion for profiles\n";
    std::cerr << dim << '\t' << from.dim << '\n';
    exit(1);
    dim = from.dim;
    delete[] profile;
    profile = new double[dim];
  }
  for (int k=0; k<dim; k++) {
    profile[k] = from.profile[k];
  }
  return *this;
}

void Profile::setuniform() {
  for (int i=0; i<dim; i++) {
    profile[i] = 1.0 / dim;
  }
}

void Profile::setarray(double* in) {
  for (int i=0; i<dim; i++) {
    profile[i] = in[i];
  }
}

const double* Profile::GetArray() const {return profile;}

double* Profile::GetArray() {return profile;}

double& Profile::operator[](int i) {
  return profile[i];
}

double& Profile::operator[](int i) const  {
  return profile[i];
}

int Profile::GetDim() const {return dim;}

void Profile::SetAtZero() {
  for (int i=0; i<dim; i++) {
    profile[i] = 0;
  }
}

void Profile::ScalarMultiplication(double d) {
  for (int i=0; i<dim; i++) {
    profile[i] *= d;
  }
}

void Profile::Add(const Profile& in) {
  for (int i=0; i<dim; i++) {
    profile[i] += in[i];
  }
}

int Profile::Check() {return 1;}

double Profile::GetEntropy() const {
  double total = 0;
  for (int i=0; i<dim; i++) {
    total += (profile[i]>1e-8) ? -profile[i]*log(profile[i]) : 0;
  }
  return total;
}

double Profile::ProposeMove(double tuning, int dim);

double Profile::ProposeMove(double tuning) {
  return ProposeMove(tuning,dim);
}

double PosRealVector::ProposeMove(double tuning, int n) {
  if ((n<=0) || (n > dim)) {
    n = dim;
  }
  int* indices = new int[n];
  Random::DrawFromUrn(indices,n,dim);
  double ret = 0;
  for (int i=0; i<n; i++) {
    double m = tuning * (Random::Uniform() - 0.5);
    vec[indices[i]] *= exp(m);
    ret += m;
  }
  delete[] indices;
  return ret;
}

double PosRealVector::GetEntropy() const {
  double total = 0;
  for (int i=0; i<dim; i++) {
    total += vec[i];
  }
  double ent = 0;
  for (int i=0; i<dim; i++) {
    double tmp = vec[i]/total;
    ent += (tmp>1e-8) ? -tmp*log(tmp) : 0;
  }
  return ent;
}
