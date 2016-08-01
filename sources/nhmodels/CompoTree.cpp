#include "CompoTree.h"

void CompoTree::ToTikz(string target, double sizeX, double sizeY, double fontsize)	{

	string texfile = target + ".tex";
	if (header)	{
		string auxpath = "~/mcmc/aux_ps/";
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
		os << "\\begin{frame}[plain,fragile]\n";
	}

	ComputeSize(root);
	ComputeDepth(root);
	if (minmax)	{
		double max = ComputeMaxNodeVal(root);
		double min = ComputeMinNodeVal(root);
		cerr << min << '\t' << max << '\n';
		minnodeval = min;
		nodescale = maxnodeval / sqrt(max - min);
	}
	else	{
		double max = ComputeMaxNodeVal(root);
		nodescale = maxnodeval / sqrt(max);
	}

	double scaleX = sizeX / GetNode(root)->GetDepth();
	double scaleY = sizeY / (GetNode(root)->GetSize() - 1);
	
	os << "\\begin{tikzpicture}\n";
	os << "[anchor=west,font=\\fontsize{" << fontsize << "}{" << fontsize << "}\\selectfont,aa/.style={font=\\fontsize{" << alphabetfontsize << "}{7}\\selectfont,minimum height=4mm}]\n";
	os << '\n';

	DrawTikz(root,os,0,sizeY/2,scaleX,scaleY);

	os << '\n';
	os << "\\end{tikzpicture}\n";

	if (header)	{
		os << "\\end{frame}\n";
		os << "\\end{document}\n";
	}
	os.close();
}


inline double texapprox(double f)	{
	return ((double) ((int) (f * 1000))) / 1000;
}

void CompoTree::DrawProfile(ostream& os)	{

	os << "\\matrix[above=-" << 2 + 10 * maxminus * composcale << "mm,column sep=" << colsep << "mm,row sep =" << rowsep << "mm] at (node" << nodecount << ") {\n";
	for (int i=0; i<nstate; i++)	{
		if (compo[i] > 0)	{
			os << "\\draw[fill=blue] (-" << barwidth<< ",0) rectangle (" << barwidth << "," << texapprox(compo[i] * composcale) << ");";
		}
		if (i < nstate-1)	{
			os << "&\n";
		}
		else	{
			os << "\\\\\n";
		}
	}
	for (int i=0; i<nstate; i++)	{
		os << "\\node[aa]{\\verb!" << alphabet[i] << "!};";
		if (i < nstate-1)	{
			os << "&\n";
		}
		else	{
			os << "\\\\\n";
		}
	}
	for (int i=0; i<nstate; i++)	{
		if (compo[i] < 0)	{
			os << "\\draw[fill=blue] (-" << barwidth << "," <<  texapprox((maxminus + compo[i]) * composcale) << ") rectangle (" << barwidth << "," << texapprox(maxminus * composcale) << ");";
		}
		if (i < nstate-1)	{
			os << "&\n";
		}
		else	{
			os << "\\\\\n";
		}
	}
	os << "};\n";
}

double CompoTree::DrawTikz(Link* from, ostream& os, double X, double Y, double scaleX, double scaleY)	{

	// string format = "";
	// string format = "\\tiny ";
	// double thickness = 0.5;
//	os << "\\linethickness{" << thickness << "mm}\n";

	double ret = 0;

	if (from->isLeaf())	{
		// write species name
		os << "\\path (" << texapprox(X) + 0.1 << "," << texapprox(Y) << ") node {" <<  GetNodeName(from) << " };\n";
		ret = Y;
	}

	else	{

		double y = Y - 0.5 * GetNode(from)->GetSize() * scaleY;

		Link* link = from->Next();
		double yfirst = 0;
		double ylast = 0;

		while (link != from)	{

			// compute horizontal offset for child node
			double xoffset = GetLength(link) * scaleX;

			// compute vertical span of child node 
			double yspan = GetNode(link->Out())->GetSize() * scaleY;

			// shift half the vertical offset
			y += 0.5 * yspan;

			// draw child node
			double x = X + xoffset;
			double downret = DrawTikz(link->Out(), os, x, y, scaleX, scaleY);

			// draw horizontal line
			// here : prepare a special node in the middle of the path
			os << "\\path [draw] (" << texapprox(X) << ',' << texapprox(downret) << ") -- node[above](node" << nodecount << "){} +(" << texapprox(xoffset) << ",0);\n";

			// at this special node,
			// draw the frequency skyline
			GetComposition(link->Out());
			DrawProfile(os);
			nodecount++;

			/*
			if (withnodecircles)	{
				double d = GetCircleNodeVal(link);
				if (d > 0)	{
					os << "\\path [fill=red,opacity=0.8,anchor=center] (" << texapprox(x) << ',' << texapprox(downret) << ") circle (" << d << "pt);\n";
				}
				else	{
					os << "\\path [fill=blue,opacity=0.8,anchor=center] (" << texapprox(x) << ',' << texapprox(downret) << ") circle (" << -d << "pt);\n";
				}
			}
			*/
			if ((withnodetext == 1) && (! link->Out()->isLeaf()))	{
				os << "\\path [anchor=west] (" << texapprox(x) << "," << texapprox(downret) << ") node {" <<  GetNodeName(link) << " };\n";
			}

			if (link == from->Next())	{
				yfirst = downret;
			}

			if (link->Next() == from)	{
				ylast = downret;
			}

			// shift half the vertical offset again
			y += 0.5 * yspan;

			link = link->Next();
		}

		// draw vertical line
		os << "\\path [draw] (" << texapprox(X) << ',' << texapprox(yfirst) << ") -- +(0," << texapprox(ylast-yfirst) << ");\n";

		ret = (yfirst + ylast) / 2;
		
	}
	return ret;
}
