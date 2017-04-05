#ifndef IIDNORMALIIDARRAY_H
#define IIDNORMALIIDARRAY_H

#include "IID.hpp"
#include "Normal.hpp"
#include "utils/Chrono.hpp"

class IIDNormalIIDArray : public IIDArray<RealVector> {
  public:
    IIDNormalIIDArray(int insize, int indim, Var<Real> *inmean, Var<PosReal> *invariance)
        : IIDArray<RealVector>(insize) {
        mean = inmean;
        variance = invariance;
        dim = indim;
        Create();
    }

    void ClampAtZero() {
        for (int i = 0; i < GetSize(); i++) {
            GetNormalVal(i)->ClampAtZero();
        }
    }

    // make a MoveN(double tuning, int m)  calls over all sites
    virtual double MoveN(double tuning, int m) {
        double tot = 0;
        for (int i = 0; i < GetSize(); i++) {
            tot += this->GetNormalVal(i)->Move(tuning, m);
        }

        return tot / GetSize();
    }

    IIDNormal *GetNormalVal(int site) { return dynamic_cast<IIDNormal *>(GetVal(site)); }

    double GetMeanVar() {
        double mean = 0;
        for (int i = 0; i < GetSize(); i++) {
            mean += GetNormalVal(i)->GetVar();
        }
        mean /= GetSize();
        return mean;
    }

  protected:
    Rvar<RealVector> *CreateVal(int /*site*/) override {
        return new IIDNormal(dim, mean, variance);
    }

    Var<Real> *mean;
    Var<PosReal> *variance;
    int dim;
};

class IIDNormalIIDArrayMove : public MCUpdate {
  public:
    IIDNormalIIDArrayMove(IIDNormalIIDArray *inselectarray, double intuning, int inm) {
        selectarray = inselectarray;
        tuning = intuning;
        m = inm;
    }

    double Move(double tuning_modulator = 1) override {
        return selectarray->MoveN(tuning_modulator * tuning, m);
    }

  private:
    IIDNormalIIDArray *selectarray;
    double tuning;
    int m;
};

#endif  // IIDNORMALIIDARRAY_H