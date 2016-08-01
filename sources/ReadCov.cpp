
#include "Random.h"

#include "MeanCovMatrix.h"

int main(int argc, char* argv[])	{

	int burnin = 1;
	string name = "";

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-b")	{
				i++;
				burnin = atoi(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readcoevol [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	ifstream is(name.c_str());

	int size = 0;

	int dim = -1;

	ReadLine(is);


	string tmp;
	is >> tmp >> tmp >> tmp >> tmp;
	is >> tmp;
	if (tmp != "[")	{
		cerr << "error when reading file\n";
		exit(1);
	}
	is >> tmp;
	if (tmp != "[")	{
		cerr << "error when reading file\n";
		exit(1);
	}
	is >> tmp;
	while ((tmp != "]") && (is))	{
		is >> tmp;
		dim++;
	}
	if (! is)	{
		cerr << "error: end of file\n";
		exit(1);
	}
	cerr << "dimension: " << dim << '\n';
	
	ReadLine(is);

	size ++;

	while (is && (size < burnin))	{
		ReadLine(is);
		size++;
	}

	if (! is)	{
		cerr << "error: end of file before reaching burnin\n";
		exit(1);
	}
			
	size = 0;

	MeanCovMatrix* meancov = new MeanCovMatrix(dim);

	CovMatrix* cov = new CovMatrix(dim);

	while (is)	{

		string tmp;
		is >> tmp >> tmp >> tmp >> tmp;
		is >> tmp;
		if (tmp != "[")	{
			cerr << "error when reading file\n";
			exit(1);
		}

		for (int i=0; i<dim; i++)	{
			is >> tmp;
			if (tmp != "[")	{
				cerr << "error when reading file\n";
				exit(1);
			}
			for (int j=0; j<dim; j++)	{
				is >> tmp;
				(*cov)[i][j] = atof(tmp.c_str());
			}
			is >> tmp;
			if (tmp != "]")	{
				cerr << "error when reading file\n";
				exit(1);
			}
			if (i < dim-1)	{
				is >> tmp;
				if (tmp != ",")	{
					cerr << "error when reading file\n";
					exit(1);
				}
			}
		}
		
		size++;
	}


		
}

