
#include "Tree.h"


int main(int argc, char* argv[])	{

	string file = "";
	string out = "";

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-o")	{
				i++;
				out = argv[i];
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				file = argv[i];
			}
			i++;
		}
		if (file == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "drawtree [-x <xsize> -y <ysize> -o <outfile> -ns <node_val_scale> -mm] <inputfile> \n";
		cerr << '\n';
		exit(1);
	}

	if (out == "")	{
		out = file;
	}

	SimpleDrawTree* tree = new SimpleDrawTree(file);
	tree->Draw(out);
}

