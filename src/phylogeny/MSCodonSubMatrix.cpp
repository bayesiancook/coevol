#include "MSCodonSubMatrix.hpp"
using namespace std;

void MGFitnessCodonSubMatrix::ComputeStationary() {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGFitnessCodonSubMatrix::ComputeArray(int i) {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                if (a == b) {
                    cerr << GetCodonStateSpace()->GetState(i) << '\t'
                         << GetCodonStateSpace()->GetState(j) << '\n';
                    cerr << pos << '\n';
                    exit(1);
                }
                Q[i][j] = (*NucMatrix)(a, b);
                if (!Synonymous(i, j)) {
                    Q[i][j] *= sqrt((GetFitness(GetCodonStateSpace()->Translation(j))) /
                                    (GetFitness(GetCodonStateSpace()->Translation(i))));
                }
            } else {
                Q[i][j] = 0;
            }
            total += Q[i][j];
        }
    }
    Q[i][i] = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

// when the fitness profile and codonusageselection are from Dirichlet
void MGFitnessCodonUsageSubMatrix::ComputeStationary() {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i)) *
                         GetCodonUsageSelection(i);
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGFitnessCodonUsageSubMatrix::ComputeArray(int i) {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                if (a == b) {
                    cerr << GetCodonStateSpace()->GetState(i) << '\t'
                         << GetCodonStateSpace()->GetState(j) << '\n';
                    cerr << pos << '\n';
                    exit(1);
                }
                Q[i][j] = (*NucMatrix)(a, b) *
                          sqrt(GetCodonUsageSelection(j) / GetCodonUsageSelection(i));
                if (!Synonymous(i, j)) {
                    Q[i][j] *= sqrt(((GetFitness(GetCodonStateSpace()->Translation(j))) /
                                     (GetFitness(GetCodonStateSpace()->Translation(i)))));
                }
            } else {
                Q[i][j] = 0;
            }
            total += Q[i][j];
        }
    }
    Q[i][i] = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

// when the fitness profile and codonusageselection are from Normal or gamma //
// phenomenological
// square root model
void MGSRFitnessNormalCodonUsageSubMatrix::ComputeStationary() {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i)) *
                         exp(GetCodonUsageSelection(i));

        total += mStationary[i];
    }
    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGSRFitnessNormalCodonUsageSubMatrix::ComputeArray(int i) {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                if (a == b) {
                    cerr << GetCodonStateSpace()->GetState(i) << '\t'
                         << GetCodonStateSpace()->GetState(j) << '\n';
                    cerr << pos << '\n';
                    exit(1);
                }
                Q[i][j] = (*NucMatrix)(a, b) *
                          sqrt(exp(GetCodonUsageSelection(j) - GetCodonUsageSelection(i)));
                if (!Synonymous(i, j)) {
                    Q[i][j] *= sqrt(((GetFitness(GetCodonStateSpace()->Translation(j))) /
                                     (GetFitness(GetCodonStateSpace()->Translation(i)))));
                }
            } else {
                Q[i][j] = 0;
            }
            total += Q[i][j];
        }
    }

    Q[i][i] = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

// mechanical mutation-selection model
void MGMSFitnessNormalCodonUsageSubMatrix::ComputeStationary() {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i)) *
                         GetCodonUsageSelection(i);

        total += mStationary[i];
    }
    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGMSFitnessNormalCodonUsageSubMatrix::ComputeArray(int i) {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                if (a == b) {
                    cerr << GetCodonStateSpace()->GetState(i) << '\t'
                         << GetCodonStateSpace()->GetState(j) << '\n';
                    cerr << pos << '\n';
                    exit(1);
                }
                Q[i][j] = (*NucMatrix)(a, b);
                // Q[i][j] = (*NucMatrix)(a,b) *
                // sqrt(exp(GetCodonUsageSelection(j)-GetCodonUsageSelection(i))) ;

                double deltaS;
                if (!Synonymous(i, j)) {
                    deltaS = log((GetFitness(GetCodonStateSpace()->Translation(j))) /
                                 (GetFitness(GetCodonStateSpace()->Translation(i)))) +
                             log(GetCodonUsageSelection(j) / GetCodonUsageSelection(i));
                } else {
                    deltaS = log(GetCodonUsageSelection(j) / GetCodonUsageSelection(i));
                }

                if ((fabs(deltaS)) < 1e-10) {
                    Q[i][j] *= 1 / (1 - (deltaS / 2));
                }
                /*
                  else if (deltaS > 50)	{
                  Q[i][j] *= deltaS;
                  }
                  else if (deltaS < -50)	{
                  Q[i][j] = 0;
                  }
                */

                else {
                    Q[i][j] *= (deltaS) / (1.0 - exp(-deltaS));
                }
                if (isinf(Q[i][j])) {
                    cerr << "Q matrix infinite: " << Q[i][j] << '\n';
                    cerr << deltaS << '\t' << GetFitness(GetCodonStateSpace()->Translation(i))
                         << '\t' << GetFitness(GetCodonStateSpace()->Translation(j)) << '\n';
                    exit(1);
                }

            } else {
                Q[i][j] = 0;
            }
            total += Q[i][j];

            if (isinf(Q[i][j])) {
                cerr << "Q matrix infinite: " << Q[i][j] << '\n';
                exit(1);
            }

            if (Q[i][j] < 0) {
                cerr << "Q matrix negative: " << Q[i][j] << '\n';
                exit(1);
            }
        }
    }

    Q[i][i] = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

void MGSRFitnessCodonUsageSubMatrix::ComputeStationary() {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i)) *
                         GetCodonUsageSelection(i);
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGSRFitnessCodonUsageSubMatrix::ComputeArray(int i) {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                if (a == b) {
                    cerr << GetCodonStateSpace()->GetState(i) << '\t'
                         << GetCodonStateSpace()->GetState(j) << '\n';
                    cerr << pos << '\n';
                    exit(1);
                }
                Q[i][j] = (*NucMatrix)(a, b) *
                          sqrt(GetCodonUsageSelection(j) / GetCodonUsageSelection(i));
                if (!Synonymous(i, j)) {
                    Q[i][j] *= sqrt(((GetFitness(GetCodonStateSpace()->Translation(j))) /
                                     (GetFitness(GetCodonStateSpace()->Translation(i)))));
                }
            } else {
                Q[i][j] = 0;
            }
            total += Q[i][j];
        }
    }

    Q[i][i] = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

void MGMSFitnessCodonUsageSubMatrix::ComputeStationary() {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i)) *
                         GetCodonUsageSelection(i);
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    // double min = 1;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
        /*
          if (min > mStationary[i])	{
          min = mStationary[i];
          }
        */
    }
    /*
      if (min < 1e-10)	{
      cerr << "small stat\n";
      }
    */
}

void MGMSFitnessCodonUsageSubMatrix::ComputeArray(int i) {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);

                Q[i][j] = (*NucMatrix)(a, b);

                double deltaS;
                if (!Synonymous(i, j)) {
                    deltaS = (log(GetFitness(GetCodonStateSpace()->Translation(j))) -
                              log(GetFitness(GetCodonStateSpace()->Translation(i)))) +
                             (log(GetCodonUsageSelection(j)) - log(GetCodonUsageSelection(i)));
                } else {
                    deltaS = log(GetCodonUsageSelection(j)) - log(GetCodonUsageSelection(i));
                }
                if ((fabs(deltaS)) < 1e-30) {
                    Q[i][j] *= 1 + deltaS / 2;
                } else if (deltaS > 50) {
                    Q[i][j] *= deltaS;
                } else if (deltaS < -50) {
                    Q[i][j] = 0;
                }
                if (deltaS != 0) {
                    // else	{
                    Q[i][j] *= deltaS / (1.0 - exp(-deltaS));
                }
                if (isinf(Q[i][j])) {
                    cerr << "Q matrix infinite: " << Q[i][j] << '\n';
                    cerr << deltaS << '\t' << GetFitness(GetCodonStateSpace()->Translation(i))
                         << '\t' << GetFitness(GetCodonStateSpace()->Translation(j)) << '\n';
                    cerr << deltaS << '\t' << log(GetFitness(GetCodonStateSpace()->Translation(i)))
                         << '\t' << log(GetFitness(GetCodonStateSpace()->Translation(j))) << '\n';
                    cerr << log(GetFitness(GetCodonStateSpace()->Translation(i))) -
                                log(GetFitness(GetCodonStateSpace()->Translation(j)))
                         << '\n';
                    cerr << GetCodonUsageSelection(i) << '\t' << GetCodonUsageSelection(j) << '\n';
                    exit(1);
                }
            } else {
                Q[i][j] = 0;
            }
            total += Q[i][j];

            if (isinf(Q[i][j])) {
                cerr << "Q matrix infinite: " << Q[i][j] << '\n';
                exit(1);
            }

            if (Q[i][j] < 0) {
                cerr << "Q matrix negative: " << Q[i][j] << '\n';
                exit(1);
            }
        }
    }

    Q[i][i] = -total;

    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}
