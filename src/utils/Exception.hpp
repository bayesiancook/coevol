#ifndef EXCEPTION_H
#define EXCEPTION_H

class Exception {};

class CheckSumException : public Exception {
  public:
    CheckSumException(double in) : checksum(in) {}

    inline double GetCheckSum() { return checksum; }

  private:
    double checksum;
};

#endif  // EXCEPTION_H
