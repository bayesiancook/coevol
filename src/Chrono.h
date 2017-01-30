#ifndef CHRONO_H
#define CHRONO_H

#include <ctime>

#include <sys/time.h>

class Chrono	{

	public:

	Chrono(){};
	~Chrono() {};
	void Reset();
	void Start();
	void Stop();
	int operator++();

	double GetTime();
	double GetTimePerCount();
	int GetCount();

	private:

	// this is in milli seconds
	double sec1;
	double sec2;
	double milli1;
	double milli2;

	double TotalTime;
	int N;

}
;



#endif

