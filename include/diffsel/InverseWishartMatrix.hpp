#ifndef WISHART_H
#define WISHART_H

#include <fstream>
#include "CovMatrix.hpp"
#include "GenericTimeLine.hpp"
#include "GlobalScalingFunction.hpp"
#include "Tree.hpp"
#include "ValArray.hpp"

class SigmaZero : public Dnode {
  public:
    SigmaZero(Var<PosReal>** inarray, int indim) {
        diag = nullptr;
        dim = indim;
        array = inarray;
        for (int i = 0; i < GetDim(); i++) {
            Register(array[i]);
        }
        specialUpdate();
    }

    SigmaZero(VarArray<PosReal>* indiag) {
        array = nullptr;
        diag = indiag;
        dim = diag->GetSize();
        for (int i = 0; i < GetDim(); i++) {
            Register(diag->GetVal(i));
        }
        specialUpdate();
    }

    ~SigmaZero() override = default;

    double val(int index) {
        if (diag != nullptr) {
            return diag->GetVal(index)->val();
        }
        return array[index]->val();
    }


    // this actually draws a value based on the INVERSE of the matrix
    void drawValInv(double* x) {
        for (int i = 0; i < GetDim(); i++) {
            x[i] = Random::sNormal() / sqrt(val(i));
        }
    }

    int GetDim() { return dim; }

    double GetDeterminant() {
        double p = 1;
        for (int i = 0; i < GetDim(); i++) {
            p *= val(i);
        }
        return p;
    }

    double GetLogDeterminant() {
        double logp = 0;
        for (int i = 0; i < GetDim(); i++) {
            logp += log(val(i));
        }
        return logp;
    }

    void specialUpdate() override {}

  private:
    int dim;
    VarArray<PosReal>* diag;
    Var<PosReal>** array;
};

class InverseWishartMatrix : public virtual Rvar<CovMatrix> {
  protected:
    SigmaZero* diagonalMatrix;
    int P;

  public:
    SigmaZero* GetDiagonalMatrix() { return diagonalMatrix; }

    int GetP() { return P; }

    InverseWishartMatrix(SigmaZero* inDiagonalMatrix, int inP)
        : Rvar<CovMatrix>(CovMatrix(inDiagonalMatrix->GetDim())) {
        diagonalMatrix = inDiagonalMatrix;
        P = inP;
        Register(diagonalMatrix);
        Sample();
    }


    void drawSample() override {
        int cont = 1;
        auto echantillon = new double*[P];
        for (int i = 0; i < P; i++) {
            echantillon[i] = new double[GetDim()];
        }

        while (cont != 0) {
            // draws a sample of P vectors from a normal distribution of variance covariance equal
            // to the inverse of SigmaZero
            for (int i = 0; i < P; i++) {
                diagonalMatrix->drawValInv(echantillon[i]);
            }

            // compute the scatter matrix, invert it, and this gives the sample
            ScatterizeOnZero(echantillon, P);

            cont = Invert();
            /*
              if (Invert())	{
              std::cerr << "in regular draw sample\n";

              std::cerr << '\n';
              std::cerr << "diag matrix\n";
              for (int i=0; i<GetDim(); i++)	{
              std::cerr << diagonalMatrix->val(i) << '\n';
              }
              std::cerr << '\n';

              std::cerr << "P : " << P << '\n';
              std::cerr << '\n';

              std::cerr << "sample\n";
              for (int i=0; i< P ; i++) {
              for (int j=0; j<GetDim(); j++)	{
              std::cerr << echantillon[i][j] << '\t';
              }
              std::cerr << '\n';
              }

              exit(1);
              }
            */
        }

        for (int i = 0; i < P; i++) {
            delete[] echantillon[i];
        }
        delete[] echantillon;
    }

    void drawSample(CovMatrix* A) {
        // algorithm of Odell and Feiveson, 1966
        auto v = new double[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            v[i] = Random::Gamma(0.5 * (P - i), 0.5);
        }
        auto n = new double*[GetDim()];
        auto b = new double*[GetDim()];
        auto a = new double*[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            n[i] = new double[GetDim()];
            b[i] = new double[GetDim()];
            a[i] = new double[GetDim()];
        }

        for (int i = 0; i < GetDim(); i++) {
            for (int j = i + 1; j < GetDim(); j++) {
                n[i][j] = Random::sNormal();
            }
        }
        for (int i = 0; i < GetDim(); i++) {
            b[i][i] = v[i];
            for (int k = 0; k < i; k++) {
                b[i][i] += n[k][i] * n[k][i];
            }
        }
        for (int i = 0; i < GetDim(); i++) {
            for (int j = i + 1; j < GetDim(); j++) {
                b[i][j] = n[i][j] * sqrt(v[i]);
                for (int k = 0; k < i; k++) {
                    b[i][j] += n[k][i] * n[k][j];
                }
                b[j][i] = b[i][j];
            }
        }

        A->CorruptDiag();
        A->Diagonalise();

        double** p = A->GetEigenVect();
        // double** invp = A->GetInvEigenVect();
        double* d = A->GetEigenVal();


        for (int i = 0; i < GetDim(); i++) {
            for (int j = 0; j < GetDim(); j++) {
                a[i][j] = p[i][j] / sqrt(d[j]);
            }
        }

        for (int i = 0; i < GetDim(); i++) {
            for (int j = 0; j < GetDim(); j++) {
                double tmp = 0;
                for (int k = 0; k < GetDim(); k++) {
                    tmp += b[i][k] * a[j][k];
                }
                n[i][j] = tmp;
            }
        }

        for (int i = 0; i < GetDim(); i++) {
            for (int j = 0; j < GetDim(); j++) {
                double tmp = 0;
                for (int k = 0; k < GetDim(); k++) {
                    tmp += a[i][k] * n[k][j];
                }
                (*this)[i][j] = tmp;
            }
        }

        Invert();

        for (int i = 0; i < GetDim(); i++) {
            delete[] n[i];
            delete[] b[i];
            delete[] a[i];
        }
        delete[] n;
        delete[] b;
        delete[] a;
        delete[] v;
    }

    // old naive versions: draw P iid normal and then calculate scattermatrix
    void drawSampleOld(CovMatrix* A) {
        auto echantillon = new double*[P];
        for (int i = 0; i < P; i++) {
            echantillon[i] = new double[GetDim()];
            A->drawValInv(echantillon[i]);
        }
        ScatterizeOnZero(echantillon, P);
        for (int i = 0; i < P; i++) {
            delete[] echantillon[i];
        }
        delete[] echantillon;
        if (Invert() != 0) {
            std::cerr << "in draw sample cov matrix A\n";
            std::cerr << '\n' << *A << '\n';
            std::cerr << "DF : " << P << '\n';
            exit(1);
        }
    }

    void SetAndClampAtSigmaZero() {
        for (int i = 0; i < GetDim(); i++) {
            value[i][i] = diagonalMatrix->val(i);
            for (int j = 0; j < i; j++) {
                value[i][j] = 0;
                value[j][i] = 0;
            }
        }
        Clamp();
    }

    void SetAtSigmaZero() {
        for (int i = 0; i < GetDim(); i++) {
            value[i][i] = diagonalMatrix->val(i);
            for (int j = 0; j < i; j++) {
                value[i][j] = 0;
                // value[i][j] = Random::Uniform()*1e-3 * sqrt(diagonalMatrix->val(i) *
                // diagonalMatrix->val(j));
                value[j][i] = value[i][j];
            }
        }
    }

    void localRestore() override {
        Rvar<CovMatrix>::localRestore();
        diagflag = false;
    }


    void localCorrupt(bool bk) override {
        Rvar<CovMatrix>::localCorrupt(bk);
        diagflag = false;
    }


    double logProb() override {
        if (isPositiveDefine()) {
            double sum = 0;
            for (int i = 0; i < GetDim(); i++) {
                sum += GetInvMatrix()[i][i] * diagonalMatrix->val(i);
            }
            double d = -((GetLogDeterminant() * (GetDim() + P + 1)) + sum) * 0.5;
            d += diagonalMatrix->GetLogDeterminant() * P * 0.5;
            return d;
        }
        // std::cerr << "singular cov matrix\n";
        return -1000;
    }

    double MarginalLogProb(CovMatrix& scalestat, int shapestat) {
        for (int i = 0; i < dim; i++) {
            scalestat[i][i] += diagonalMatrix->val(i);
        }
        auto postscale = new CovMatrix(scalestat);
        int postshape = P + shapestat;
        double l = diagonalMatrix->GetLogDeterminant() * P * 0.5 -
                   postscale->GetLogDeterminant() * postshape * 0.5;
        delete postscale;
        for (int i = 0; i < dim; i++) {
            scalestat[i][i] -= diagonalMatrix->val(i);
        }
        return l;
    }

    void GibbsResample(CovMatrix& scalestat, int shapestat) {
        for (int i = 0; i < dim; i++) {
            scalestat[i][i] += diagonalMatrix->val(i);
        }
        P += shapestat;
        scalestat.Diagonalise();
        drawSample(&scalestat);
        for (int i = 0; i < GetDim(); i++) {
            scalestat[i][i] -= diagonalMatrix->val(i);
        }
        P -= shapestat;
    }

    void LeftRightMultiply(double** P) {
        auto aux = new double*[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            aux[i] = new double[GetDim()];
        }

        for (int i = 0; i < GetDim(); i++) {
            for (int j = 0; j < GetDim(); j++) {
                double tmp = 0;
                for (int k = 0; k < GetDim(); k++) {
                    tmp += (*this)[i][k] * P[j][k];
                }
                aux[i][j] = tmp;
            }
        }

        for (int i = 0; i < GetDim(); i++) {
            for (int j = 0; j < GetDim(); j++) {
                double tmp = 0;
                for (int k = 0; k < GetDim(); k++) {
                    tmp += P[i][k] * aux[k][j];
                }
                (*this)[i][j] = tmp;
            }
        }

        for (int i = 0; i < GetDim(); i++) {
            delete[] aux[i];
        }
        delete[] aux;
    }
};

class DiagonalCovMatrix : public virtual Rvar<CovMatrix> {
  protected:
    SigmaZero* diagonalMatrix;
    int P;

  public:
    DiagonalCovMatrix(SigmaZero* inDiagonalMatrix, int inP)
        : Rvar<CovMatrix>(CovMatrix(inDiagonalMatrix->GetDim())) {
        diagonalMatrix = inDiagonalMatrix;
        P = inP;
        Register(diagonalMatrix);
        // SetAtSigmaZero();
        Sample();
    }

    void drawSample() override {
        // draw each entry along the diagonal from an inverse gamma of parameters alpha and beta
        for (int i = 0; i < GetDim(); i++) {
            double alpha = 0.5 * (P + GetDim() - 1);
            // double alpha = 0.5 * (P - GetDim() + 1);
            double beta = 0.5 * diagonalMatrix->val(i);
            double tmp = Random::Gamma(alpha, beta);
            value[i][i] = 1.0 / tmp;
        }
        diagflag = false;
    }

    void SetAndClampAtSigmaZero() {
        SetAtSigmaZero();
        Clamp();
    }

    void SetAtSigmaZero() {
        for (int i = 0; i < GetDim(); i++) {
            value[i][i] = diagonalMatrix->val(i);
            bkvalue[i][i] = diagonalMatrix->val(i);
            for (int j = 0; j < i; j++) {
                value[i][j] = 0;
                bkvalue[i][j] = 0;
                value[j][i] = 0;
                bkvalue[j][i] = 0;
            }
        }
    }

    double GetLogDeterminant() {
        double total = 0;
        for (int i = 0; i < GetDim(); i++) {
            total += log(diagonalMatrix->val(i));
        }
        return total;
    }

    double GetDeterminant() {
        double total = 1;
        for (int i = 0; i < GetDim(); i++) {
            total *= diagonalMatrix->val(i);
        }
        return total;
    }

    void localRestore() override {
        Rvar<CovMatrix>::localRestore();
        diagflag = false;
    }


    void localCorrupt(bool bk) override {
        Rvar<CovMatrix>::localCorrupt(bk);
        diagflag = false;
    }

    double ProposeMove(double tuning) override {
        // multiplicative move on one entry chosen at random along the diagonal
        int k = (int)(GetDim() * Random::Uniform());
        double h = tuning * (Random::Uniform() - 0.5);
        double e = exp(h);
        value[k][k] *= e;
        return h;
    }

    double logProb() override {
        double sum = 0;
        for (int i = 0; i < GetDim(); i++) {
            double alpha = 0.5 * (P + GetDim() - 1);
            // double alpha = 0.5 * (P - GetDim() + 1);
            double beta = 0.5 * diagonalMatrix->val(i);
            sum += alpha * log(beta) - (alpha + 1) * log(value[i][i]) - beta / value[i][i] -
                   Random::logGamma(alpha);
        }
        return sum;
    }
};


class BoundForMultiNormal {
  public:
    BoundForMultiNormal() : tax1(""), tax2(""), lowerbound(-1), upperbound(-1), index(-1) {}

    BoundForMultiNormal(std::string intax1, std::string intax2, int inindex, double inlowerbound,
                        double inupperbound)
        : tax1(std::move(intax1)),
          tax2(std::move(intax2)),
          lowerbound(inlowerbound),
          upperbound(inupperbound),
          index(inindex) {}

    ~BoundForMultiNormal() = default;

    double GetUpperBound() const { return upperbound; }
    double GetLowerBound() const { return lowerbound; }

    int GetIndex() const { return index; }

    std::string GetTaxon1() const { return tax1; }
    std::string GetTaxon2() const { return tax2; }

    void ToStream(std::ostream& os) const {
        os << tax1 << '\t' << tax2 << '\t' << index << '\t' << lowerbound << '\t' << upperbound
           << '\n';
    }

    void FromStream(std::istream& is) {
        is >> tax1 >> tax2 >> index >> lowerbound >> upperbound;
        std::cerr << index << '\n';
    }

    friend std::ostream& operator<<(std::ostream& os, const BoundForMultiNormal& cal) {
        cal.ToStream(os);
        return os;
    }

    friend std::istream& operator>>(std::istream& is, BoundForMultiNormal& cal) {
        cal.FromStream(is);
        return is;
    }

  private:
    std::string tax1;
    std::string tax2;
    double lowerbound;
    double upperbound;
    int index;
};

class BoundSet {
  public:
    BoundSet(Tree* intree) : tree(intree) {}

    void ToStream(std::ostream& os) const {
        int Ncalib = 0;
        for (auto i = nodemap.begin(); i != nodemap.end(); i++) {
            Ncalib++;
        }
        os << Ncalib << '\n';
        for (const auto& i : nodemap) {
            os << i.second << '\n';
        }
    }

    void FromStream(std::istream& is) {
        int Ncalib;
        is >> Ncalib;
        for (int i = 0; i < Ncalib; i++) {
            BoundForMultiNormal bound;
            is >> bound;
            const Link* link = tree->GetLCA(bound.GetTaxon1(), bound.GetTaxon2());
            if (link == nullptr) {
                std::cerr << "error in calibration set: did not find common ancestor of " << bound
                     << '\n';
                exit(1);
            }
            nodemap[link->GetNode()] = bound;
        }
    }

    const std::map<const Node*, BoundForMultiNormal>& GetNodeMap() const { return nodemap; }

  private:
    std::map<const Node*, BoundForMultiNormal> nodemap;
    Tree* tree;
};

class FileBoundSet : public BoundSet {
  public:
    FileBoundSet(std::string filename, Tree* intree) : BoundSet(intree) {
        std::ifstream is(filename.c_str());
        FromStream(is);
    }
};


// Vector of reals draw from a Wishart distribution.
class MultiNormal : public virtual Rvar<RealVector> {
  public:
    MultiNormal(Var<CovMatrix>* insigma, Var<RealVector>* inup, Var<PosReal>* intime,
                Var<PosReal>* indate, Var<PosReal>* inagescale,
                GlobalScalingFunction* inscalefunction) {
        scale = nullptr;
        drift = nullptr;
        driftphi = nullptr;
        drift2 = nullptr;
        driftphi2 = nullptr;
        timeline = nullptr;
        alpha = nullptr;
        rootmean = nullptr;
        rootvar = nullptr;

        sigma = insigma;
        up = inup;
        time = intime;
        date = indate;
        agescale = inagescale;
        scalefunction = inscalefunction;

        Register(sigma);
        Register(up);
        Register(time);
        Register(date);
        Register(agescale);
        Register(scalefunction);

        setval(RealVector(insigma->GetDim()));
        ClampVector = new bool[GetDim()];
        hasupperbound = new bool[GetDim()];
        haslowerbound = new bool[GetDim()];
        upperbound = new double[GetDim()];
        lowerbound = new double[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            hasupperbound[i] = false;
            haslowerbound[i] = false;
        }
        hasbounds = false;
        rootmax = 1000;
        for (int i = 0; i < GetDim(); i++) {
            UnClamp(i);
        }
        Sample();
    }

    MultiNormal(int /*unused*/, Var<CovMatrix>* insigma, Var<RealVector>* inrootmean,
                Var<PosRealVector>* inrootvar) {
        scalefunction = nullptr;

        sigma = insigma;
        up = nullptr;
        time = nullptr;
        scale = nullptr;
        drift = nullptr;
        driftphi = nullptr;
        drift2 = nullptr;
        driftphi2 = nullptr;
        date = nullptr;
        timeline = nullptr;
        alpha = nullptr;

        rootmean = inrootmean;
        rootvar = inrootvar;
        Register(sigma);
        if (rootmean != nullptr) {
            Register(rootmean);
        }
        if (rootvar != nullptr) {
            Register(rootvar);
        }

        setval(RealVector(insigma->GetDim()));
        ClampVector = new bool[GetDim()];
        hasupperbound = new bool[GetDim()];
        haslowerbound = new bool[GetDim()];
        upperbound = new double[GetDim()];
        lowerbound = new double[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            hasupperbound[i] = false;
            haslowerbound[i] = false;
        }
        hasbounds = false;
        rootmax = 1000;
        for (int i = 0; i < GetDim(); i++) {
            UnClamp(i);
        }
        Sample();
    }

    MultiNormal(Var<CovMatrix>* insigma, Var<RealVector>* inup = nullptr,
                Var<PosReal>* intime = nullptr, Var<PosReal>* inscale = nullptr,
                Var<RealVector>* indrift = nullptr, Var<PosReal>* indriftphi = nullptr,
                Var<PosReal>* indate = nullptr, Var<RealVector>* indrift2 = nullptr,
                Var<PosReal>* indriftphi2 = nullptr, Var<PosReal>* inagescale = nullptr,
                double inkt = 0, GenericTimeLine* intimeline = nullptr,
                Var<Real>* inalpha = nullptr) {
        scalefunction = nullptr;

        rootmean = nullptr;
        rootvar = nullptr;

        sigma = insigma;
        up = inup;
        time = intime;
        scale = inscale;
        drift = indrift;
        driftphi = indriftphi;
        drift2 = indrift2;
        driftphi2 = indriftphi2;
        date = indate;
        agescale = inagescale;
        kt = inkt;
        timeline = intimeline;
        alpha = inalpha;
        setval(RealVector(insigma->GetDim()));
        ClampVector = new bool[GetDim()];
        hasupperbound = new bool[GetDim()];
        haslowerbound = new bool[GetDim()];
        upperbound = new double[GetDim()];
        lowerbound = new double[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            hasupperbound[i] = false;
            haslowerbound[i] = false;
        }
        hasbounds = false;
        rootmax = 1000;
        if (isRoot()) {
            // Clamp(GetDim());
        } else {
            Register(up);
            Register(time);
        }
        Register(sigma);
        if (alpha != nullptr) {
            Register(alpha);
        }
        if (scale != nullptr) {
            Register(scale);
        }
        if (drift != nullptr) {
            Register(drift);
        }
        if (driftphi != nullptr) {
            Register(driftphi);
            if (date == nullptr) {
                std::cerr << "error in MultiNormal::MultiNormal: null date\n";
                exit(1);
            }
            Register(date);
        }
        if (drift2 != nullptr) {
            Register(drift2);
        }
        if (driftphi2 != nullptr) {
            Register(driftphi2);
        }
        if (agescale != nullptr) {
            Register(agescale);
        }
        for (int i = 0; i < GetDim(); i++) {
            UnClamp(i);
        }
        Sample();
    }

    ~MultiNormal() override {
        delete[] ClampVector;
        delete[] hasupperbound;
        delete[] haslowerbound;
        delete[] upperbound;
        delete[] lowerbound;
    }

    bool HasBounds() { return hasbounds; }

    void SetBound(const BoundForMultiNormal& bound, int offset) {
        hasbounds = true;
        int index = bound.GetIndex() + offset;
        double upper = bound.GetUpperBound();
        std::cerr << bound.GetIndex() << '\t' << offset << '\t' << index << '\t' << upper << '\n';
        if (upper != -1) {
            if (upper <= 0) {
                std::cerr << "error : negative upper bound : " << upper << '\n';
                exit(1);
            }
            SetUpperBound(index, log(upper));
        }
        double lower = bound.GetLowerBound();
        if (lower != -1) {
            if (lower <= 0) {
                std::cerr << "error : negative lower bound : " << lower << '\n';
                exit(1);
            }
            SetLowerBound(index, log(lower));
        }
    }

    void SetUpperBound(int i, double in) {
        hasupperbound[i] = true;
        upperbound[i] = in;
    }

    void SetLowerBound(int i, double in) {
        haslowerbound[i] = true;
        lowerbound[i] = in;
    }

    void CheckBounds() {
        for (int i = 0; i < GetDim(); i++) {
            if (haslowerbound[i] && hasupperbound[i]) {
                while ((val()[i] > static_cast<double>(hasupperbound[i])) &&
                       (val()[i] < static_cast<double>(haslowerbound[i]))) {
                    if (val()[i] > upperbound[i]) {
                        val()[i] = 2 * upperbound[i] - val()[i];
                    }
                    if (val()[i] < lowerbound[i]) {
                        val()[i] = 2 * lowerbound[i] - val()[i];
                    }
                }
            } else if (haslowerbound[i]) {
                if (val()[i] < lowerbound[i]) {
                    val()[i] = 2 * lowerbound[i] - val()[i];
                }
            } else if (hasupperbound[i]) {
                if (val()[i] > upperbound[i]) {
                    val()[i] = 2 * upperbound[i] - val()[i];
                }
            }
        }
    }

    double GetRootMax() { return rootmax; }

    void SetRootMax(double d) { rootmax = d; }

    void drawSampleUnclamped() {
        if (HasBounds()) {
            std::cerr << "error: resampling a multinormal with bounds\n";
            exit(1);
        }
        if (isRoot()) {
            if (rootmean != nullptr) {
                for (int i = 0; i < GetDim(); i++) {
                    if ((*rootvar)[i] != 0.0) {
                        val()[i] = Random::sNormal() * sqrt((*rootvar)[i]) + (*rootmean)[i];
                    }
                }
            } else {
                for (int i = 0; i < GetDim(); i++) {
                    val()[i] = 0;
                    // val()[i] = 2*rootmax * (Random::Uniform() - 0.5);
                }
            }
        } else {
            auto t = new double[GetDim()];
            sigma->drawVal(t);

            for (int i = 0; i < GetDim(); i++) {
                double tt = time->val();
                if (scalefunction != nullptr) {
                    double t2 = date->val() * agescale->val();
                    double t1 = (date->val() + time->val()) * agescale->val();
                    double scalefactor = scalefunction->GetScalingFactor(t1, t2);
                    tt *= scalefactor;
                } else if (scale != nullptr) {
                    tt *= scale->val();
                }
                if (drift != nullptr) {
                    if ((driftphi != nullptr) && (driftphi2 != nullptr)) {
                        //  not yet done
                        val()[i] = t[i] * sqrt(tt) + up->val()[i];
                    } else if (driftphi != nullptr) {
                        double u = time->val();
                        if (driftphi->val() > 1e-10) {
                            u = exp(-driftphi->val() * date->val()) *
                                (1 - exp(-driftphi->val() * time->val())) / driftphi->val();
                        }
                        val()[i] = t[i] * sqrt(tt) + (*drift)[i] * u + up->val()[i];
                    } else {
                        val()[i] = t[i] * sqrt(tt) + (*drift)[i] * tt + up->val()[i];
                    }
                } else {
                    val()[i] = t[i] * sqrt(tt) + up->val()[i];
                }
            }
            delete[] t;
        }
    }

    void drawSample() override {
        if (HasBounds()) {
            std::cerr << "error: resampling a multinormal with bounds\n";
            exit(1);
        }
        if (isRoot()) {
            if (rootmean != nullptr) {
                for (int i = 0; i < GetDim(); i++) {
                    if (!ClampVector[i]) {
                        if ((*rootvar)[i] != 0.0) {
                            val()[i] = Random::sNormal() * sqrt((*rootvar)[i]) + (*rootmean)[i];
                        }
                    }
                }
            } else {
                for (int i = 0; i < GetDim(); i++) {
                    if (!ClampVector[i]) {
                        val()[i] = 0;
                        // val()[i] = 2*rootmax * (Random::Uniform() - 0.5);
                    }
                }
            }
        } else {
            auto t = new double[GetDim()];
            sigma->drawVal(t);
            for (int i = 0; i < GetDim(); i++) {
                if (!ClampVector[i]) {
                    double tt = time->val();
                    if (scalefunction != nullptr) {
                        double t2 = date->val() * agescale->val();
                        double t1 = (date->val() + time->val()) * agescale->val();
                        double scalefactor = scalefunction->GetScalingFactor(t1, t2);
                        tt *= scalefactor;
                    } else if (scale != nullptr) {
                        tt *= scale->val();
                    }
                    if (drift != nullptr) {
                        if ((driftphi != nullptr) && (driftphi2 != nullptr)) {
                            //  not yet done
                            val()[i] = t[i] * sqrt(tt) + up->val()[i];
                        } else if (driftphi != nullptr) {
                            double u = time->val();
                            if (driftphi->val() > 1e-10) {
                                u = exp(-driftphi->val() * date->val()) *
                                    (1 - exp(-driftphi->val() * time->val())) / driftphi->val();
                            }
                            val()[i] = t[i] * sqrt(tt) + (*drift)[i] * u + up->val()[i];
                        } else {
                            val()[i] = t[i] * sqrt(tt) + (*drift)[i] * tt + up->val()[i];
                        }
                    } else {
                        val()[i] = t[i] * sqrt(tt) + up->val()[i];
                    }
                }
            }
            delete[] t;
        }
    }

    void DrawNormalFromPrecision(double** mean, CovMatrix* cov, int offset = 0) {
        int k = cov->GetDim();
        /*
          if (k > 1)	{
          std::cerr << "dim error\n";
          exit(1);
          }
        */
        auto t = new double[k];
        cov->drawValInv(t);
        for (int i = 0; i < k; i++) {
            if (ClampVector[i + offset]) {
                std::cerr << "error in MultiNormal::DrawNormal: clamped\n";
                exit(1);
            }
            val()[i + offset] = t[i] + mean[i][0];
        }
        delete[] t;
    }

    void ConditionalDrawNormalFromPrecision(double** mean, CovMatrix* cov, int* index, int k) {
        auto t = new double[k];
        cov->drawValInv(t);
        for (int i = 0; i < k; i++) {
            if (ClampVector[index[i]]) {
                std::cerr << "error in MultiNormal::DrawNormal: clamped\n";
                exit(1);
            }
            val()[index[i]] = t[i] + mean[i][0];
        }
        delete[] t;
    }

    void DrawNormalFromPrecisionWithConstraint(double** mean, CovMatrix* cov) {
        int K = 0;
        for (int i = 0; i < GetDim(); i++) {
            if (!ClampVector[i]) {
                K++;
            }
        }
        int L = GetDim() - K;
        auto index1 = new int[K];
        auto index2 = new int[L];
        auto mean1 = new double[K];
        auto mean2 = new double[L];
        int k = 0;
        int l = 0;

        for (int i = 0; i < GetDim(); i++) {
            if (ClampVector[i]) {
                mean2[l] = mean[i][0] - (*this)[i];
                index2[l] = i;
                l++;
            } else {
                mean1[k] = mean[i][0];
                index1[k] = i;
                k++;
            }
        }

        auto Cov11 = new CovMatrix(K);
        double** cov11 = Cov11->GetMatrix();
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                cov11[i][j] = (*cov)[index1[i]][index1[j]];
            }
        }
        Cov11->CorruptDiag();
        Cov11->Diagonalise();
        double** invcov11 = Cov11->GetInvMatrix();

        auto cov12 = new double*[K];
        for (int i = 0; i < K; i++) {
            cov12[i] = new double[L];
            for (int j = 0; j < L; j++) {
                cov12[i][j] = (*cov)[index1[i]][index2[j]];
            }
        }
        auto tmp = new double[K];
        for (int i = 0; i < K; i++) {
            tmp[i] = 0;
            for (int j = 0; j < L; j++) {
                tmp[i] += cov12[i][j] * mean2[j];
            }
        }
        for (int i = 0; i < K; i++) {
            double temp = 0;
            for (int j = 0; j < K; j++) {
                temp += invcov11[i][j] * tmp[j];
            }
            mean1[i] += temp;
        }
        Cov11->drawValInv(tmp);
        for (int i = 0; i < K; i++) {
            (*this)[index1[i]] = mean1[i] + tmp[i];
        }

        for (int i = 0; i < K; i++) {
            delete[] cov12[i];
        }
        delete[] cov12;
        delete cov11;
        delete[] tmp;
        delete[] index1;
        delete[] index2;
        delete[] mean1;
        delete[] mean2;
    }

    void LeftMultiply(double** P) {
        auto aux = new double[GetDim()];
        for (int i = 0; i < GetDim(); i++) {
            double tmp = 0;
            for (int j = 0; j < GetDim(); j++) {
                tmp += P[i][j] * (*this)[j];
            }
            aux[i] = tmp;
        }
        for (int i = 0; i < GetDim(); i++) {
            (*this)[i] = aux[i];
        }
        delete[] aux;
    }

    bool CheckIdentity() {
        if (GetDim() != bkvalue.GetDim()) {
            std::cerr << "non matching dimension : " << GetDim() << '\t' << bkvalue.GetDim() << '\n';
            exit(1);
        }
        bool ok = true;
        for (int i = 0; i < GetDim(); i++) {
            ok &= static_cast<int>(fabs((*this)[i] - (bkvalue)[i]) < 1e-8);
        }
        if (!ok) {
            for (int i = 0; i < GetDim(); i++) {
                std::cerr << (*this)[i] << '\t' << (bkvalue)[i] << '\n';
            }
        }
        return ok;
    }

    double SegmentMove(double tuning, int imin, int size) {
        Corrupt(true);
        double logHastings = ProposeSegmentMove(tuning, imin, size);
        double deltaLogProb = Update();
        double logRatio = deltaLogProb + logHastings;
        bool accepted = (log(Random::Uniform()) < logRatio);
        if (!accepted) {
            Corrupt(false);
            Restore();
            return 0;
        }
        return 1;
    }

    double PiecewiseTranslationMove(double tuning, int index, int n) {
        if (!isClamped()) {
            Corrupt(true);
            double logHastings = ProposePiecewiseTranslationMove(tuning, index, n);
            double deltaLogProb = Update();
            double logRatio = deltaLogProb + logHastings;
            bool accepted = (log(Random::Uniform()) < logRatio);
            if (!accepted) {
                Corrupt(false);
                Restore();
                return 0;
            }
        }
        return 1;
    }

    double ProposePiecewiseTranslationMove(double tuning, int index, int n) {
        double u = tuning * (Random::Uniform() - 0.5);
        return PiecewiseTranslation(u, index, n);
    }

    double PiecewiseTranslation(double u, int index, int n) {
        for (int i = 0; i < n; i++) {
            if (i + index >= GetDim()) {
                std::cerr << "error in MultiNormal::PiecewiseTranslationMove : " << index << '\t' << n
                     << '\t' << GetDim() << "\n";
                exit(1);
            }
            if ((!isClamped()) && (!ClampVector[index + i])) {
                (*this)[index + i] += u;
            }
        }
        /*
          bool noclamped = true;
          for ( int i=0; i <n; i++){
          if (i+index >= GetDim())	{
          std::cerr << "error in MultiNormal::PiecewiseTranslationMove : " << index << '\t' << n << '\t'
          << GetDim() << "\n";
          exit(1);
          }
          if (ClampVector[index+i])	{
          noclamped = false;
          }
          }

          if (noclamped)	{
          for ( int i=0; i <n; i++){
          (*this)[index + i] += u;
          }
          }
        */
        return 0;
    }

    void Shift(double* delta, double f = 1) {
        for (int i = 0; i < GetDim(); i++) {
            if (!ClampVector[i]) {
                val()[i] += f * delta[i];
            }
        }
    }

    double ProposeMove(double tuning) override {
        if (HasBounds()) {
            return SimpleProposeMove(tuning);
        }

        if ((!isRoot()) && Random::Uniform() < 0.5) {
            return MatrixBasedProposeMove(tuning);
        }
        return SimpleProposeMove(tuning);

        return 0;
    }

    double MatrixBasedProposeMove(double tuning) {
        auto t = new double[GetDim()];
        sigma->drawVal(t);
        for (int i = 0; i < GetDim(); i++) {
            if (!ClampVector[i]) {
                double tt = time->val();
                if (scale != nullptr) {
                    tt *= scale->val();
                }
                val()[i] += tuning * t[i] * sqrt(tt);
            }
        }
        delete[] t;
        return 0;
    }

    double ProposeSegmentMove(double tuning, int imin, int size) {
        if ((imin < 0) || (imin + size >= GetDim())) {
            std::cerr << "error in propose segment move: overflow\n";
            std::cerr << imin << '\t' << size << '\t' << GetDim() << '\n';
            exit(1);
        }

        int choose = (int)(3 * Random::Uniform());
        if (choose == 0) {
            for (int i = imin; i < imin + size; i++) {
                if (!ClampVector[i]) {
                    val()[i] += tuning * (Random::Uniform() - 0.5);
                }
            }
        } else if (choose == 1) {
            double u = tuning * (Random::Uniform() - 0.5);
            for (int i = imin; i < imin + size; i++) {
                if (!ClampVector[i]) {
                    val()[i] += u;
                }
            }
        } else {
            int i = (int)(size * Random::Uniform()) + imin;
            if (!ClampVector[i]) {
                val()[i] += tuning * (Random::Uniform() - 0.5);
            }
        }
        CheckBounds();
        return 0;
    }

    double SimpleProposeMove(double tuning) {
        int choose = (int)(3 * Random::Uniform());
        if (choose == 0) {
            for (int i = 0; i < GetDim(); i++) {
                if (!ClampVector[i]) {
                    val()[i] += tuning * (Random::Uniform() - 0.5);
                }
            }
        } else if (choose == 1) {
            double u = tuning * (Random::Uniform() - 0.5);
            for (int i = 0; i < GetDim(); i++) {
                if (!ClampVector[i]) {
                    val()[i] += u;
                }
            }
        } else {
            int i = (int)(GetDim() * Random::Uniform());
            if (!ClampVector[i]) {
                val()[i] += tuning * (Random::Uniform() - 0.5);
            }
        }
        CheckBounds();
        return 0;
    }


    bool isRoot() { return (this->up == nullptr); }

    /*
      bool isclamped(int index)	{
      return ClampVector[index];
      }
    */

    void Clamp(int index) { ClampVector[index] = true; }

    void ClampAt(double inval, int index) {
        val()[index] = inval;
        ClampVector[index] = true;
    }

    void ClampAtZero() {
        for (int i = 0; i < GetDim(); i++) {
            val()[i] = 0;
            ClampVector[i] = true;
        }
    }

    void UnClamp(int index) { ClampVector[index] = false; }

    double GetTrend(double t, int index) {
        if (drift2 != nullptr) {
            if (t > kt) {
                return 0;
            }
            return (*drift)[index] * exp(-driftphi->val() * t) +
                   (*drift2)[index] * (1 - exp(-driftphi2->val() * (kt - t)));
        }
        return (*drift)[index] * exp(-driftphi->val() * t);
    }

    double logProb() override {
        if (isRoot()) {
            double total = 0;
            if (rootmean != nullptr) {
                for (int i = 0; i < GetDim(); i++) {
                    if ((*rootvar)[i] != 0.0) {
                        double tmp = val()[i] - (*rootmean)[i];
                        total += -0.5 * (log((*rootvar)[i]) + tmp * tmp / (*rootvar)[i]);
                    } else {
                        if (fabs(val()[i]) < rootmax) {
                            total -= log(2 * rootmax);
                        } else {
                            total -= 100;
                        }
                    }
                }
            } else {
                for (int i = 0; i < GetDim(); i++) {
                    if (fabs(val()[i]) < rootmax) {
                        total -= log(2 * rootmax);
                    } else {
                        total -= 100;
                    }
                }
            }
            return total;
        }
        // calcul de transposé de X * sigma * X
        double tXSX = 0;
        auto dval = new double[GetDim()];
        double tt = time->val();
        if (scalefunction != nullptr) {
            double t2 = date->val() * agescale->val();
            double t1 = (date->val() + time->val()) * agescale->val();
            double scalefactor = scalefunction->GetScalingFactor(t1, t2);
            tt *= scalefactor;
        } else if (scale != nullptr) {
            tt *= scale->val();
        }
        double roott = sqrt(tt);
        if (drift != nullptr) {
            if (timeline != nullptr) {
                double t2 = date->val() * agescale->val();
                double t1 = (date->val() + time->val()) * agescale->val();
                for (int i = 0; i < GetDim(); i++) {
                    double f1 = timeline->GetVal(t1, i);
                    double f2 = timeline->GetVal(t2, i);
                    dval[i] = ((val()[i] - up->val()[i]) - alpha->val() * (f2 - f1)) / roott;
                }
            } else if (driftphi2 != nullptr) {
                double t2 = date->val() * agescale->val();
                double t1 = (date->val() + time->val()) * agescale->val();
                if ((t1 < 0) || (t1 > 2000) || (t2 < 0) || (t2 > 2000)) {
                    std::cerr << "error in timeline2\n";
                    std::cerr << t1 << '\t' << t2 << '\n';
                    exit(1);
                }
                for (int i = 0; i < GetDim(); i++) {
                    double f1 = GetTrend(t1, i);
                    double f2 = GetTrend(t2, i);
                    dval[i] = ((val()[i] - up->val()[i]) - (f2 - f1)) / roott;
                }
            } else if (driftphi != nullptr) {
                double t2 = date->val();
                double t1 = (date->val() + time->val());
                for (int i = 0; i < GetDim(); i++) {
                    double f1 = GetTrend(t1, i);
                    double f2 = GetTrend(t2, i);
                    dval[i] = ((val()[i] - up->val()[i]) - (f2 - f1)) / roott;
                }
            } else {
                for (int i = 0; i < GetDim(); i++) {
                    dval[i] = (val()[i] - up->val()[i] - tt * (*drift)[i]) / roott;
                }
            }
        } else {
            for (int i = 0; i < GetDim(); i++) {
                dval[i] = (val()[i] - up->val()[i]) / roott;
            }
        }
        for (int i = 0; i < GetDim(); i++) {
            tXSX += sigma->GetInvMatrix()[i][i] * dval[i] * dval[i];
            for (int j = 0; j < i; j++) {
                tXSX += 2 * sigma->GetInvMatrix()[i][j] * dval[j] * dval[i];
            }
        }
        double d = -0.5 * (sigma->GetLogDeterminant() + tXSX + GetDim() * 2 * log(roott));
        delete[] dval;
        return d;
    }


    bool* ClampVector;

  protected:
    bool hasbounds;
    bool* hasupperbound;
    bool* haslowerbound;
    double* upperbound;
    double* lowerbound;
    Var<PosReal>* scale;
    Var<PosReal>* time;
    Var<RealVector>* up;
    Var<CovMatrix>* sigma;
    Var<RealVector>* drift;
    Var<PosReal>* driftphi;
    Var<RealVector>* drift2;
    Var<PosReal>* driftphi2;
    Var<PosReal>* date;
    Var<PosReal>* agescale;
    double kt;

    Var<RealVector>* rootmean;
    Var<PosRealVector>* rootvar;

    GenericTimeLine* timeline;
    Var<Real>* alpha;

    GlobalScalingFunction* scalefunction;

    double rootmax;
};

// Vector of reals draw from a Wishart distribution.
class MultivariateNormal : public virtual Rvar<RealVector> {
  public:
    MultivariateNormal(Var<CovMatrix>* insigma) {
        sigma = insigma;
        setval(RealVector(insigma->GetDim()));
        Register(sigma);
        Sample();
    }

    ~MultivariateNormal() override = default;

    void drawSample() override { sigma->drawVal(val().GetArray()); }

    double ProposeMove(double tuning) override {
        if (Random::Uniform() < 0.3) {
            for (int i = 0; i < GetDim(); i++) {
                val()[i] += tuning * (Random::Uniform() - 0.5);
            }
        } else if (Random::Uniform() < 0.6) {
            double u = tuning * (Random::Uniform() - 0.5);
            for (int i = 0; i < GetDim(); i++) {
                val()[i] += u;
            }
        } else {
            int i = (int)(GetDim() * Random::Uniform());
            val()[i] += tuning * (Random::Uniform() - 0.5);
        }
        return 0;
    }


    double logProb() override {
        double tXSX = 0;
        for (int i = 0; i < GetDim(); i++) {
            tXSX += sigma->GetInvMatrix()[i][i] * val()[i] * val()[i];
            for (int j = 0; j < i; j++) {
                tXSX += 2 * sigma->GetInvMatrix()[i][j] * val()[i] * val()[j];
            }
        }
        double d = -0.5 * (sigma->GetLogDeterminant() + tXSX);
        return d;
    }

  protected:
    Var<CovMatrix>* sigma;
};


#endif
