#include "SequenceAlignment.h"

int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);
	cout << "alphabet size : " << data->GetNstate() << '\n';
	cout << "mean diversity: " << data->GetMeanDiversity() << '\n';

}
