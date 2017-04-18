// SubMatrix:
// this class implements
// the instant rate matrix of a substitution process
// but only the mathematical aspects of it
//
// RandomSubMatrix:
// this class derives from SubMatrix and from Dnode
// therefore, it knows everything about substitution processes (as a SubMatrix)
// and at the same time, it can be inserted as a deterministic node in a
// probabilistic model (as a
// Dnode)
//
// If you need to define a new substitution process
// - derive a new class deriving (directly or indirectly) from SubMatrix
// in this class, implements the ComputeArray and ComputeStationary functions
// (these are the functions that construct the instant rates and the stationary
// probabilities of the
// matrix)
//
// - derive a new class deriving both from RandomSubMatrix, and from your new
// class
// in this class, implement SetParameters
// (this function is responsible for updating the internal parameters that the
// SubMatrix uses in
// ComputeArray And ComputeStationary,
// based on the values stored by the parent nodes)
//
//

#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <cmath>
#include <cstdlib>
#include <iostream>

class AbstractTransitionMatrix {
  public:
    virtual ~AbstractTransitionMatrix() = default;

    virtual void BackwardPropagate(const double *down, double *up, double length) = 0;
    virtual void ForwardPropagate(const double *up, double *down, double length) = 0;
    virtual const double *GetStationary() = 0;
    virtual double Stationary(int i) = 0;

    virtual int GetNstate() = 0;
    virtual void CorruptMatrix() = 0;
    virtual double operator()(int, int) = 0;
    virtual const double *GetRow(int i) = 0;

    virtual bool check() { return true; }
};


class SubMatrix : public virtual AbstractTransitionMatrix {
  protected:
    // these 2 pure virtual functions are the most essential component of the
    // SubMatrix class
    // see GTRSubMatrix.cpp and CodonSubMatrix.cpp for examples

    // ComputeArray(int state) is in charge of computing the row of the rate
    // matrix
    // corresponding to all possible rates of substitution AWAY from state
    //
    virtual void ComputeArray(int state) = 0;

    // ComputeStationary() is in charge of computing the vector of stationary
    // probabilities
    // (equilibirum frequencies)
    // of the substitution process
    virtual void ComputeStationary() = 0;

  public:
    static const int UniSubNmax = 500;

    static int nuni;
    static int nunimax;

    static double GetMeanUni() { return ((double)nunimax) / nuni; }

    SubMatrix(int inNstate, bool innormalise = false);
    ~SubMatrix() override;

    void Create();

    double operator()(int /*i*/, int /*j*/) override;
    const double *GetRow(int i) override;

    const double *GetStationary() override;
    double Stationary(int i) override;

    int GetNstate() override { return Nstate; }

    double GetRate();
    void ScalarMul(double e);

    bool isNormalised() { return normalise; }
    void Normalise();

    void CorruptMatrix() override;
    void UpdateMatrix();

    void ActivatePowers();
    void InactivatePowers();
    double Power(int n, int i, int j);
    double GetUniformizationMu();

    double *GetEigenVal();
    double **GetEigenVect();
    double **GetInvEigenVect();

    virtual void ToStream(std::ostream &os);
    void CheckReversibility();

    int GetDiagStat() { return ndiagfailed; }

    void BackwardPropagate(const double *up, double *down, double length) override;
    void ForwardPropagate(const double *down, double *up, double length) override;
    // virtual void     FiniteTime(int i0, double* down, double length);

    double **GetQ() { return Q; }
    void ComputeExponential(double range, double **expo);
    void ApproachExponential(double range, double **expo, int prec = 1024);
    void PowerOf2(double **y, int z);

    static double meanz;
    static double maxz;
    static double nz;

  protected:
    void UpdateRow(int state);
    void UpdateStationary();

    void ComputePowers(int N);
    void CreatePowers(int n);

    bool ArrayUpdated();

    int Diagonalise();

    // data members

    bool powflag;
    bool diagflag;
    bool statflag;
    bool *flagarray;

    int Nstate;
    int npow;
    double UniMu;

    double ***mPow;

    // Q : the infinitesimal generator matrix
    double **Q;

    // the stationary probabilities of the matrix
    double *mStationary;

    bool normalise;

    // an auxiliary matrix
    double **aux;

  protected:
    // v : eigenvalues
    // vi : imaginary part
    // u : the matrix of eigen vectors
    // invu : the inverse of u

    double **u;
    double **invu;
    double *v;
    double *vi;

    int ndiagfailed;
};

//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------

inline double SubMatrix::operator()(int i, int j) {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q[i][j];
}

inline const double *SubMatrix::GetRow(int i) {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q[i];
}

inline const double *SubMatrix::GetStationary() {
    if (!statflag) {
        UpdateStationary();
    }
    return mStationary;
}

inline double SubMatrix::Stationary(int i) {
    if (!statflag) {
        UpdateStationary();
    }
    return mStationary[i];
}

inline void SubMatrix::CorruptMatrix() {
    diagflag = false;
    statflag = false;
    for (int k = 0; k < Nstate; k++) {
        flagarray[k] = false;
    }
    InactivatePowers();
}

inline bool SubMatrix::ArrayUpdated() {
    bool qflag = true;
    for (int k = 0; k < Nstate; k++) {
        qflag &= static_cast<int>(flagarray[k]);
    }
    return qflag;
}

inline void SubMatrix::UpdateStationary() {
    ComputeStationary();
    statflag = true;
}

inline void SubMatrix::UpdateRow(int state) {
    if (isNormalised()) {
        UpdateMatrix();
    } else {
        if (!statflag) {
            UpdateStationary();
        }
        ComputeArray(state);
        flagarray[state] = true;
    }
}

inline void SubMatrix::BackwardPropagate(const double *up, double *down, double length) {
    double **eigenvect = GetEigenVect();
    double **inveigenvect = GetInvEigenVect();
    double *eigenval = GetEigenVal();

    int matSize = GetNstate();

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] = 0;
    }
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += inveigenvect[i][j] * up[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] *= exp(length * eigenval[i]);
    }

    for (int i = 0; i < GetNstate(); i++) {
        down[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            down[i] += eigenvect[i][j] * aux[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        if (std::isnan(down[i])) {
            std::cerr << "error in back prop\n";
            for (int j = 0; j < GetNstate(); j++) {
                std::cerr << up[j] << '\t' << down[j] << '\t' << Stationary(j) << '\n';
            }
            exit(1);
        }
    }
    double maxup = 0;
    for (int k = 0; k < matSize; k++) {
        if (up[k] < 0) {
            std::cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
            // down[k] = 0;
        }
        if (maxup < up[k]) {
            maxup = up[k];
        }
    }
    double max = 0;
    for (int k = 0; k < matSize; k++) {
        if (down[k] < 0) {
            down[k] = 0;
        }
        if (max < down[k]) {
            max = down[k];
        }
    }
    if (maxup == 0) {
        std::cerr << "error in backward propagate: null up array\n";
        exit(1);
    }
    if (max == 0) {
        std::cerr << "error in backward propagate: null array\n";
        for (int k = 0; k < matSize; k++) {
            std::cerr << up[k] << '\t' << down[k] << '\n';
        }
        std::cerr << '\n';
        exit(1);
    }
    down[matSize] = up[matSize];

    delete[] aux;
}

inline void SubMatrix::ForwardPropagate(const double *down, double *up, double length) {
    double **eigenvect = GetEigenVect();
    double **inveigenvect = GetInvEigenVect();
    double *eigenval = GetEigenVal();

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += down[j] * eigenvect[j][i];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] *= exp(length * eigenval[i]);
    }

    for (int i = 0; i < GetNstate(); i++) {
        up[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            up[i] += aux[j] * inveigenvect[j][i];
        }
    }

    delete[] aux;
}

#endif  // SUBMATRIX_H
