
// #include "CompoTree.h"
#include "FileDrawTree.h"

int main(int argc, char* argv[])	{

	string file = "";
	string out = "";
	double nodescale = 5;
	double nodepower = 1;
	double x = 6;
	double y = 8;
	double z = 0;
	// double x = 10;
	// double y = 17;
	bool sizeset = false;
	double min = -1;
	double max = -1;
	int withnodetext = 0;
	bool withleafnames = true;
	int prec = 2;
	bool withheader = true;
	bool unibranch = false;

	double barwidth = 0.04;

	// double fontsize = 8;
	double fontsize = 6;
	double groupfontsize = 8;
	// string style = "article";
	string style = "beamer";
	double xoffset = 0;
	bool fontset = false;

	int Nstate = Nnuc;
	double thickness = 0.06;

	bool bubbletree = false;
	bool chronotree = false;
	bool heattree = false;
	bool justtab = false;

	double maxtime = 0;

	double texscale = 1;

	string branchvalfile = "";
	string nodevalfile = "";
	string nodecifile = "";

	string groupfile = "";
	bool onlygroups = false;

	double rescale = 0;

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
			else if (s == "-np")	{
				i++;
				nodepower = atof(argv[i]);
			}
			else if (s == "-fs")	{
				i++;
				fontsize = atof(argv[i]);
				fontset = true;
			}
			else if (s == "-rescale")	{
				i++;
				rescale = atof(argv[i]);
			}
			else if (s == "-aa")	{
				Nstate = Naa;
			}
			else if (s == "-branch")	{
				unibranch = true;
			}
			else if (s == "-leafnames")	{
				withleafnames = false;
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
				if (! fontset)	{
					fontsize = 4;
				}
			}
			else if (s == "-article")	{
				style = "article";
				if (! sizeset)	{
					x = 15;
					y = 25;
				}
				if (! fontset)	{
					fontsize = 8;
				}
			}

			else if (s == "-scale")	{
				i++;
				texscale = atof(argv[i]);
			}
			else if (s == "-o")	{
				i++;
				out = argv[i];
			}
			else if (s == "-h")	{
				withheader = false;
			}
			else if (s == "-mm")	{
				i++;
				min = atof(argv[i]);
				i++;
				max = atof(argv[i]);
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
			else if (s == "-maxtime")	{
				i++;
				maxtime = atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
			}
			else if (s == "-color")	{
				heattree = true;
			}
			else if (s == "-branchval")	{
				i++;
				branchvalfile = argv[i];
			}
			else if (s == "-nodeval")	{
				i++;
				nodevalfile = argv[i];
			}
			else if (s == "-nodeci")	{
				i++;
				nodecifile = argv[i];
			}
			else if (s == "-groups")	{
				i++;
				groupfile = argv[i];
			}
			else if (s == "-groupmargin")	{
				i++;
				z = atof(argv[i]);
			}
			else if (s == "-groupfs")	{
				i++;
				groupfontsize = atof(argv[i]);
			}
			else if (s == "-onlygroups")	{
				onlygroups = true;
			}
			else if (s == "-tab")	{
				justtab = true;
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
		cerr << "drawtree [-x <xsize> -y <ysize> -o <outfile> -fs <font_size (default 8)> -color] <inputfile> \n";
        cerr << "\t-x <xsize> : tree depth (horizontal, default: 6)\n";
        cerr << "\t-y <xsize> : tree height (vertical, default: 8)\n";
        cerr << "\t-fs <font_size> : font size (default: 8)\n";
        cerr << "\t-color : tree with colored branches (heat map from yellow to red): input file should have point estimates at nodes\n";
        cerr << "\t-bubble : tree with bubbles attached to nodes: input file should have min_max interval estimates at nodes\n";
        // cerr << "\t-chrono : chronogram: input file should have min_max date estimates at nodes\n";
        cerr << "\t-nt : with node text (node values or names printed full text next to each node)\n";
        // cerr << "\t-scale <scale>: set scale in tex picture (default: 1)\n";
        cerr << "\t-leafnames : tip names not printed out\n";
		cerr << '\n';
		exit(1);
	}

	if (out == "")	{
		out = file;
	}

	ifstream testis(file.c_str());
	if (testis)	{
	if (justtab)	{
		FileTree* tree = new FileTree(file);
		ofstream os((out + ".tab").c_str());
		tree->Tabulate(os);
	}
	else if (heattree)	{
		FileHeatTree* tree = new FileHeatTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetTexScale(texscale);
		tree->SetWithHeader(withheader);
		tree->SetWithLeafNames(withleafnames);
		tree->SetFontSize(fontsize);
		tree->SetUniBranch(unibranch);
		tree->SetThickness(thickness);
		tree->SetMinMax(min,max);
		tree->SetNodeText(withnodetext);
		tree->SetRescalingFactor(rescale);
		if (branchvalfile != "")	{
			tree->SetBranchVal(branchvalfile);
		}
		if (nodevalfile != "")	{
			tree->SetExternalNodeVal(nodevalfile);
		}
		if (nodecifile != "")	{
			tree->SetExternalNodeCI(nodecifile);
		}
		if (groupfile != "")	{
			tree->SetGroupOffset(z);
			tree->SetGroups(groupfile);
			tree->SetGroupFontSize(groupfontsize);
		}
		tree->Draw(out);
	}
	else if (bubbletree)	{
		cerr << "bubble\n";
		FileBubbleTree* tree = new FileBubbleTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetTexScale(texscale);
		tree->SetWithHeader(withheader);
		tree->SetWithLeafNames(withleafnames);
		tree->SetFontSize(fontsize);
		tree->SetNodeScale(nodescale);
		tree->SetNodePower(nodepower);
		tree->Draw(out);
	}
	else if (chronotree)	{
		FileChronoDrawTree* tree = new FileChronoDrawTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetTexScale(texscale);
		tree->SetWithHeader(withheader);
		tree->SetWithLeafNames(withleafnames);
		tree->SetFontSize(fontsize);
		tree->SetTimeScale(false);
		tree->SetBarWidth(barwidth);
		if (groupfile != "")	{
			tree->SetGroupOffset(z);
			tree->SetGroups(groupfile);
			tree->SetGroupFontSize(groupfontsize);
			tree->SetOnlyGroups(onlygroups);
			tree->SetMaxTime(maxtime);
		}
		tree->Draw(out);
	}
	else	{
		FileDrawTree* tree = new FileDrawTree(file);
		tree->SetDocumentStyle(style);
		tree->SetScale(x,y);
		tree->SetTexScale(texscale);
		tree->SetWithHeader(withheader);
		tree->SetWithLeafNames(withleafnames);
		tree->SetNodeText(withnodetext);
		tree->SetFontSize(fontsize);
		tree->Draw(out);
	}
	}

	else	{
		cerr << "file " << file << " does not exist\n";
	}
}

