
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N;
	is >> N;
	double cov[N][N];
	double p[N][N];
	int permut[N];
	for (int i=0; i<N; i++)	{
		is >> permut[i];
	}
	for (int i=0; i<N; i++)	{
		for (int j=0; j<N; j++)	{
			is >> cov[i][j];
		}
	}

	for (int i=0; i<N; i++)	{
		for (int j=0; j<N; j++)	{
			if (i != j)	{
				is >> p[i][j];
			}
			else	{
				string c;
				is >> c;
				p[i][j] = 0;
			}
		}
	}

	string target = argv[2];

	double size = atof(argv[3]);
	double scale = atof(argv[4]);

	bool header = true;
	if (argc == 6)	{
		header = false;
	}

	string texfile = target + ".tex";
	if (header)	{
		string auxpath = "~/coevol1.1/aux_ps/";
		string header = "header_tikz.tex";
		// tex output
		string appel = "cp " + auxpath + header + " " + texfile;
		system(appel.c_str());
	}
	else	{
		string appel = "echo \"\" >  " + texfile;
		system(appel.c_str());
	}
	ofstream os(texfile.c_str(), ios_base::app);

	if (header)	{
		os << "\\begin{document}\n";
		os << "\\begin{frame}\n";
	}

	os << "\\begin{tikzpicture}\n";

	for (int i=0; i<N; i++)	{
		for (int j=i; j<N; j++)	{
			double tmp = cov[permut[i]][permut[j]];
			double op = 0.8;
			if ((p[permut[i]][permut[j]] > 0.025) && (p[permut[i]][permut[j]] < 0.975))	{
				op = 0.4;
			}
			if (tmp > 0)	{
				os << "\\path [fill=red,opacity=" << op << ",anchor=center] (" << size*j << ',' << size*(N-i) << ") circle (" << scale * sqrt(tmp) << "pt);\n";
			}
			else	{
				os << "\\path [fill=blue,opacity=" << op << ",anchor=center] (" << size*j << ',' << size*(N-i) << ") circle (" << scale * sqrt(-tmp) << "pt);\n";
			}
		}
	}
	os << '\n';
	os << "\\end{tikzpicture}\n";

	if (header)	{
		os << "\\end{frame}\n";
		os << "\\end{document}\n";
	}
	os.close();
}

