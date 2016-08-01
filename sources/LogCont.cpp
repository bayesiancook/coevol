
#include "Random.h"
#include "ContinuousData.h"

#include <fstream>
int main(int argc, char* argv[])	{

	FileContinuousData* contdata = new FileContinuousData(argv[1]);
	ofstream os (argv[2]);
	contdata->ToStreamLog(os);
}

