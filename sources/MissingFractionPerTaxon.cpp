#include "SequenceAlignment.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	FileSequenceAlignment* data = new FileSequenceAlignment(argv[1]);
    data->MissingFractionPerTaxon(cout);
}

