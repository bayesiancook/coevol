#include "Chrono.h"

void Chrono::Reset()	{

	TotalTime = 0;
	N = 0;
}

void Chrono::Start()	{
	// clock_gettime(CLOCK_REALTIME, &nano1);

	struct timeval tod;
	gettimeofday(&tod,NULL);
	sec1 = ((double) tod.tv_sec);
	milli1 = ((double) tod.tv_usec) / 1000;
}

void Chrono::Stop()	{
	/*
	clock_gettime(CLOCK_REALTIME, &nano2);
	double t1 = ((double) (nano2.tv_sec))  - ((double) (nano1.tv_sec));
	double t2 = (nano2.tv_nsec - nano1.tv_nsec) * (1e-9);
	double duration = t1 + t2;
	*/

	struct timeval tod;
	gettimeofday(&tod,NULL);
	sec2 = ((double) tod.tv_sec);
	milli2 = ((double) tod.tv_usec) / 1000;
	double duration = 1000*(sec2 - sec1) + milli2 - milli1;

	TotalTime += duration;
}

int Chrono::operator++()	{
	return N++;
}

double Chrono::GetTime()	{
	return TotalTime;
}

double Chrono::GetTimePerCount()	{
	return TotalTime / N;
}

int Chrono::GetCount()	{
	return N;
}

