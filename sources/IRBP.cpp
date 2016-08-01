
#include <fstream>
using namespace std;
#include "SequenceAlignment.h"


int main(int argc, char* argv[])	{

	FileSequenceAlignment* irbp = new FileSequenceAlignment(argv[1]);
	ofstream os(argv[2]);
	irbp->ToStream(os);
}
