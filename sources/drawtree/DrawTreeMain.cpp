
// #include "CompoTree.h"
#include "DrawTree.h"

#include "BiologicalSequences.h"

int main(int argc, char* argv[])	{

	string file = "";
	string out = "";
	double nodescale = 5;
	// double x = 6;
	// double y = 8;
	double x = 10;
	double y = 17;
	bool sizeset = false;
	bool minmax = false;
	int withnodetext = 0;
	int prec = 2;
	bool withheader = true;
	bool unibranch = false;

	double fontsize = 8;
	// double fontsize = 4;
	string style = "article";
	double xoffset = 0;

	int Nstate = Nnuc;
	double thickness = 0.06;

	bool bubbletree = false;
	bool chronotree = false;
	bool heattree = false;


	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-fs")	{
				i++;
				fontsize = atof(argv[i]);
			}
			else if (s == "-aa")	{
				Nstate = Naa;
			}
			else if (s == "-branch")	{
				unibranch = true;
			}
			else if (s == "-th")	{
				i++;
				thickness = atof(argv[i]);
			}
			else if (s == "-xoff")	{
				i++;
				xoffset = atof(argv[i]);
			}
			else if (s == "-x")	{
				i++;
				x = atof(argv[i]);
				sizeset = true;
			}
			else if (s == "-y")	{
				i++;
				y = atof(argv[i]);
				sizeset = true;
			}
			else if (s == "-beamer")	{
				style = "beamer";
				if (! sizeset)	{
					x = 6;
					y = 8;
				}
			}
			else if (s == "-article")	{
				style = "article";
				if (! sizeset)	{
					x = 15;
					y = 25;
				}
			}
			else if (s == "-o")	{
				i++;
				out = argv[i];
			}
			else if (s == "-h")	{
				withheader = false;
			}
			else if (s == "-mm")	{
				minmax = true;
			}
			else if (s == "-nt")	{
				withnodetext = 1;
			}
			else if (s == "-lt")	{
				withnodetext = 2;
			}
			else if (s == "-bubble")	{
				bubbletree = true;
			}
			else if (s == "-chrono")	{
				chronotree = true;
			}
			else if (s == "-color")	{
				heattree = true;
			}
			else if (s == "-p")	{
				i++;
				prec = atoi(argv[i]);
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

	if (heattree)	{
		FileHeatTree* tree = new FileHeatTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetFontSize(fontsize);
		tree->SetUniBranch(unibranch);
		tree->SetThickness(thickness);
		tree->Draw(out);
	}
	else if (bubbletree)	{
		FileBubbleTree* tree = new FileBubbleTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetFontSize(fontsize);
		tree->SetNodeScale(nodescale);
		tree->Draw(out);
	}
	else if (chronotree)	{
		FileChronoDrawTree* tree = new FileChronoDrawTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetFontSize(fontsize);
		tree->SetTimeScale(false);
		tree->Draw(out);
	}
	else	{
		FileDrawTree* tree = new FileDrawTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetFontSize(fontsize);
		tree->Draw(out);
	}
}

