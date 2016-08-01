
#include "SequenceAlignment.h"


class JCSequenceAlignment : public SequenceAlignment	{


	public:
	
	JCSequenceAlignment(string filename, bool pir = false);

	private:

	void ReadPIR(string filename);
	void ReadDataFromFile(string filename);
	
	void Convert();

	map<string,string> ali;
};

/*
class Concatenation : public SequenceAlignment	{

	public:

	Concatenation(SequenceAlignment** genelist, int Ngene)	{

		map<string,string> ali;
	}

};
*/
