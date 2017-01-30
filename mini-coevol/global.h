#ifndef GLOBAL_H
#define GLOBAL_H

#include "fmath.hpp"

#define MULTITHREADED // UNDEFINE FOR GET A SINGLE THREADED VERSION

#define SQRT2PI 2.5066282746310002
#define SQRT2 1.4142135623730951

inline double fastPow(double a, double b) {
  union {
	double d;
	int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

inline double fastLog(double a) {
	return fmath::log(a);
}
	
// undefine LINEAR_INTERPOLATION TO GET CUBIC_INTERPOLATION, SLOWER BUT MORE ACCURATE
#define LINEAR_INTERPOLATION

// define POW(a,b) as pow(a,b) for slow accurate(?) and fastPow(a,b) for fast
//#define POW(a,b) pow(a,b)
#define POW(a,b) fastPow(a,b)

// define LOG(a) as log(a) for slow accurate(?) and fastLog(a) for fast
//#define LOG(a) log(a)
#define LOG(a) fastLog(a)

#endif // GLOBAL_H
