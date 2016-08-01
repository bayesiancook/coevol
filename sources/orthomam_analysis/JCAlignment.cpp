#include "JCAlignment.h"
#include "StringStreamUtils.h"
#include "BiologicalSequences.h"

#include <fstream>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 SequenceAlignment
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

JCSequenceAlignment::JCSequenceAlignment(string filename, bool pir)	{

	if (pir)	{
		ReadPIR(filename);
	}
	else	{
		ReadDataFromFile(filename);
	}
	Convert();

}

void JCSequenceAlignment::ReadPIR(string filename)	{

	ifstream is(filename.c_str());
	string name = "begin";
	while (! is.eof())	{

		if (name == "begin")	{
			name = ReadLine(is);
		}
		string finalname = name.substr(1,name.length()-1);
		ostringstream seq;
		do {
			name = ReadLine(is);
			if (name[0] != '>')	{
				for (unsigned int k=0; k<name.length(); k++)	{
					if ((name[k] != '\n') && (name[k] != ' '))	{
						seq << name[k];
					}
				}
			}
		} while ((! is.eof()) && (name[0] != '>'));
		ali[finalname] = seq.str();
		// cerr << "seq : " << seq.str() << '\n';
		// cerr << seq.str().length() << '\n';
	}
}

void JCSequenceAlignment::ReadDataFromFile(string filename)	{

	ifstream is(filename.c_str());
	string name = "begin";
	while (! is.eof())	{

		if (name == "begin")	{
			name = ReadLine(is);
		}
		string name2 = StringReplace('_'," ",name);
		string name3 = StringReplace('>'," ",name2);
		istringstream s(name3);
		string genus;
		string species;
		s >> genus >> species;
		string finalname = genus + "_" + species;
		ostringstream seq;
		do {
			name = ReadLine(is);
			if (name[0] != '>')	{
				for (unsigned int k=0; k<name.length(); k++)	{
					if (name[k] != '\n')	{
						seq << name[k];
					}
				}
			}
		} while ((! is.eof()) && (name[0] != '>'));
		ali[finalname] = seq.str();
		// cerr << "seq : " << seq.str() << '\n';
		// cerr << seq.str().length() << '\n';
	}
}

void JCSequenceAlignment::Convert()	{

	statespace = new DNAStateSpace();
	Ntaxa = ali.size();
	string* names = new string[Ntaxa];
	string* seq = new string[Ntaxa];
	int k = 0;
	Nsite = -1;
	for (map<string,string>::const_iterator i=ali.begin(); i!= ali.end(); i++)	{
		names[k] = i->first;
		seq[k] = i->second;
		if (Nsite == -1)	{
			Nsite = seq[k].length();
		}
		else	{
			if (((unsigned int) Nsite) != (seq[k].length()))	{
				cerr << "error: non matching sequence lengths\n";
				cerr << Nsite << '\t' << seq[k].length() << '\n';
				cerr << seq[k] << '\n';
				exit(1);
			}
		}
		k++;
	}

	taxset = new TaxonSet(names,Ntaxa);
	Data = new int*[Ntaxa]; 
	for (int i=0; i<Ntaxa; i++)	{
		Data[i] = new int[Nsite];
		for (int j=0; j<Nsite; j++)	{
			ostringstream s;
			s << seq[i][j];
			Data[i][j] = statespace->GetState(s.str());
		}
	}
}

int main(int argc, char* argv[])	{

	string name = argv[1];
	string out = argv[2];

	JCSequenceAlignment* ali = new JCSequenceAlignment(name,true);
	ofstream os(out.c_str());
	ali->ToStream(os);
}

