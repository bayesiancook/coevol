#ifndef STABLEDISTRIBUTION_H
#define STABLEDISTRIBUTION_H

class StableDistribution {
 public:

  StableDistribution(double alpha = 2, double scale = 1);

  double alpha() const {return a;}
  double scale() const {return c;}

  double logPDF(double x) const;
  double logPDF(double X, double t) const;

  void setAlpha(double alpha);
  void setScale(double scale);

  double cubicInterpolate (double p[4], double x) const;

  double bicubicInterpolate (double p[4][4], double x, double y) const;

 private:

  static const double grid[40][10000];
  double c;
  double a;
  double s_a;
  double cpowa;
};

#endif // STABLEDISTRIBUTION_H
