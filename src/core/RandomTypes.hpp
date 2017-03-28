#ifndef RANDOMTYPES_H
#define RANDOMTYPES_H

#include "BaseType.hpp"
#include "Var.hpp"

// if variance is null (either as a pointer, or as a value), then this is an
// improper uniform
// distribution
// drawSample is then an ad-hoc sampling procedure, just to draw something when
// first instantiating
// the variable
class Normal : public virtual Rvar<Real> {
  public:
    Normal(Var<Real> *inmean, Var<PosReal> *invariance);
    Normal(Var<RealVector> *inmeanvec, Var<PosReal> *invariance, int inindex);

    double logProb() override;

    void drawSample() override;

  private:
    int index;
    Var<Real> *mean;
    Var<RealVector> *meanvec;
    Var<PosReal> *variance;
};

class Exponential : public Rvar<PosReal> {
  public:
    enum ParentType { MEAN, RATE };

    Exponential(Var<PosReal> *inscale, ParentType intype);
    Exponential(const Exponential &from);

    ~Exponential() override = default;

    void drawSample() override;

  private:
    double logProb() override;

    Var<PosReal> *scale;
    ParentType type;
};

class Gamma : public virtual Rvar<PosReal> {
  public:
    static const double GAMMAMIN;

    Gamma(Var<PosReal> *inshape, Var<PosReal> *inscale, bool inmeanvar = false);
    Gamma(const Gamma &from);

    ~Gamma() override = default;

    void drawSample() override;

  protected:
    double logProb() override;

    Var<PosReal> *scale;
    Var<PosReal> *shape;
    bool meanvar;
};

class Beta : public Rvar<UnitReal> {
  public:
    Beta(Var<PosReal> *inalpha, Var<PosReal> *inbeta);
    Beta(const Beta &from);

    ~Beta() override = default;

    void drawSample() override;

    double GetAlpha() { return alpha->val(); }
    double GetBeta() { return beta->val(); }

  private:
    double logProb() override;

    Var<PosReal> *alpha;
    Var<PosReal> *beta;
};

class PosUniform : public Rvar<PosReal> {
  public:
    PosUniform(Var<PosReal> *inroot, double inmax);
    PosUniform(const PosUniform &from);

    ~PosUniform() override = default;

    void drawSample() override;

  private:
    double logProb() override;

    Var<PosReal> *root;
    double max;
};

class Binomial : public Rvar<Int> {
  public:
    Binomial(int inN, Var<UnitReal> *intheta);
    ~Binomial() override = default;

    void drawSample() override;
    double ProposeMove(double tuning) override;
    // int Check();

  private:
    double logProb() override;

    int N;
    Var<UnitReal> *theta;
};

class Poisson : public virtual Rvar<Int> {
  public:
    Poisson(Var<PosReal> *inmu);
    ~Poisson() override = default;

    void drawSample() override;

  protected:
    double logProb() override;

    Var<PosReal> *mu;
};

class Dirichlet : public virtual Rvar<Profile> {
  public:
    Dirichlet(int dimension);
    Dirichlet(Var<Profile> *incenter, Var<PosReal> *inconcentration);

    ~Dirichlet() override = default;

    void drawSample() override;

    virtual double Move(double tuning, int n);
    double Move(double tuning) override;
    double ProposeMove(double tuning) override;

  protected:
    Var<Profile> *center;
    Var<PosReal> *concentration;

    double logProb() override;
};

class Multinomial : public virtual Rvar<IntVector> {
  public:
    Multinomial(Var<Profile> *inprobarray, int inN);
    ~Multinomial() override = default;

    void drawSample() override;

  private:
    Var<Profile> *probarray;
    int N;

    double logProb() override;
};

class FiniteDiscrete : public virtual Rvar<Int> {
  public:
    FiniteDiscrete(Var<Profile> *inprobarray);
    ~FiniteDiscrete() override = default;

    void drawSample() override;

  private:
    Var<Profile> *probarray;

    double logProb() override;
};

class IIDExp : public Rvar<PosRealVector> {
  public:
    IIDExp(int dimension, Var<PosReal> *inmean);
    IIDExp(int dimension);

    ~IIDExp() override = default;

    void drawSample() override;

    void setall(double in);

  private:
    Var<PosReal> *mean;
    int dim;

    double logProb() override;
};

class IIDGamma : public virtual Rvar<PosRealVector> {
  public:
    IIDGamma(int dimension, Var<PosReal> *inalpha, Var<PosReal> *inbeta);
    IIDGamma(int dimension);

    ~IIDGamma() override = default;

    void drawSample() override;

    void setall(double in);

  protected:
    Var<PosReal> *alpha;
    Var<PosReal> *beta;
    int dim;

    double logProb() override;
};

class Product : public Dvar<PosReal> {
  public:
    Product(Var<PosReal> *ina, Var<PosReal> *inb) {
        SetName("product");
        a = ina;
        b = inb;
        Register(a);
        Register(b);
    }

    void specialUpdate() override { setval(a->val() * b->val()); }

  protected:
    Var<PosReal> *a;
    Var<PosReal> *b;
};

#endif  // RANDOMTYPES_H
