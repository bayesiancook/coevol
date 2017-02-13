#ifndef RANDOMTYPES_H
#define RANDOMTYPES_H

#include "BaseType.hpp"
#include "Var.hpp"


// if variance is null (either as a pointer, or as a value), then this is an improper uniform
// distribution
// drawSample is then an ad-hoc sampling procedure, just to draw something when first instantiating
// the variable
class Normal : public virtual Rvar<Real> {
  public:
    Normal(Var<Real>* inmean, Var<PosReal>* invariance);
    Normal(Var<RealVector>* inmeanvec, Var<PosReal>* invariance, int inindex);

    double logProb();

  protected:
    void drawSample();

  private:
    int index;
    Var<Real>* mean;
    Var<RealVector>* meanvec;
    Var<PosReal>* variance;
};


class Exponential : public Rvar<PosReal> {
  public:
    enum ParentType { MEAN, RATE };

    Exponential(Var<PosReal>* inscale, ParentType intype);
    Exponential(const Exponential& from);

    virtual ~Exponential() {}

    void drawSample();

  private:
    virtual double logProb();

    Var<PosReal>* scale;
    ParentType type;
};


class Gamma : public virtual Rvar<PosReal> {
  public:
    static const double GAMMAMIN;

    Gamma(Var<PosReal>* inshape, Var<PosReal>* inscale, bool inmeanvar = false);
    Gamma(const Gamma& from);

    virtual ~Gamma() {}

    void drawSample();

  protected:
    virtual double logProb();

    Var<PosReal>* scale;
    Var<PosReal>* shape;
    bool meanvar;
};


class Beta : public Rvar<UnitReal> {
  public:
    Beta(Var<PosReal>* inalpha, Var<PosReal>* inbeta);
    Beta(const Beta& from);

    virtual ~Beta() {}

    void drawSample();

    double GetAlpha() { return alpha->val(); }
    double GetBeta() { return beta->val(); }

  private:
    virtual double logProb();

    Var<PosReal>* alpha;
    Var<PosReal>* beta;
};


class PosUniform : public Rvar<PosReal> {
  public:
    PosUniform(Var<PosReal>* inroot, double inmax);
    PosUniform(const PosUniform& from);

    virtual ~PosUniform() {}

    void drawSample();

  private:
    virtual double logProb();

    Var<PosReal>* root;
    double max;
};


class Binomial : public Rvar<Int> {
  public:
    Binomial(int inN, Var<UnitReal>* intheta);
    virtual ~Binomial() {}

    void drawSample();
    double ProposeMove(double tuning);
    // int Check();

  private:
    virtual double logProb();

    int N;
    Var<UnitReal>* theta;
};


class Poisson : public virtual Rvar<Int> {
  public:
    Poisson(Var<PosReal>* inmu);
    virtual ~Poisson() {}

    void drawSample();

  protected:
    virtual double logProb();

    Var<PosReal>* mu;
};


class Dirichlet : public virtual Rvar<Profile> {
  public:
    Dirichlet(int dimension);
    Dirichlet(Var<Profile>* incenter, Var<PosReal>* inconcentration);

    virtual ~Dirichlet() {}

    void drawSample();

    virtual double Move(double tuning, int n);
    virtual double Move(double tuning);
    virtual double ProposeMove(double tuning);

  protected:
    Var<Profile>* center;
    Var<PosReal>* concentration;

    virtual double logProb();
};


class Multinomial : public virtual Rvar<IntVector> {
  public:
    Multinomial(Var<Profile>* inprobarray, int inN);
    virtual ~Multinomial() {}

    void drawSample();

  private:
    Var<Profile>* probarray;
    int N;

    virtual double logProb();
};


class FiniteDiscrete : public virtual Rvar<Int> {
  public:
    FiniteDiscrete(Var<Profile>* inprobarray);
    virtual ~FiniteDiscrete() {}

    void drawSample();

  private:
    Var<Profile>* probarray;

    virtual double logProb();
};


class IIDExp : public Rvar<PosRealVector> {
  public:
    IIDExp(int indim, Var<PosReal>* inmean);
    IIDExp(int indim);

    virtual ~IIDExp() {}

    void drawSample();

    void setall(double in);

  private:
    Var<PosReal>* mean;
    int dim;

    virtual double logProb();
};


class IIDGamma : public virtual Rvar<PosRealVector> {
  public:
    IIDGamma(int indim, Var<PosReal>* inalpha, Var<PosReal>* inbeta);
    IIDGamma(int indim);

    virtual ~IIDGamma() {}

    void drawSample();

    void setall(double in);

  protected:
    Var<PosReal>* alpha;
    Var<PosReal>* beta;
    int dim;

    virtual double logProb();
};


class Product : public Dvar<PosReal> {
  public:
    Product(Var<PosReal>* ina, Var<PosReal>* inb) {
        SetName("product");
        a = ina;
        b = inb;
        Register(a);
        Register(b);
    }

    void specialUpdate() { setval(a->val() * b->val()); }

  protected:
    Var<PosReal>* a;
    Var<PosReal>* b;
};

#endif  // RANDOMTYPES_H
