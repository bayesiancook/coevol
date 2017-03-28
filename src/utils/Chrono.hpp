#ifndef CHRONO_H
#define CHRONO_H

#include <sys/time.h>
#include <ctime>

class Chrono {
  public:
    void Reset();
    void Start();
    void Stop();

    inline int operator++() { return N++; }

    inline double GetTime() { return TotalTime; }

    inline double GetTimePerCount() { return TotalTime / N; }

    inline int GetCount() { return N; }

  private:
    // this is in milli seconds
    double sec1;
    double sec2;
    double milli1;
    double milli2;
    double TotalTime;
    int N;
};

#endif  // CHRONO_H
