
#include "DrawTree.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/regex.hpp>  // Boost.Regex lib
#include "StringStreamUtils.h"

using namespace boost;

void FileTree::Check(const Link* from)	{

	string name = from->GetNode()->GetName();
	regex re;
	int minmax = 2;
	string s = from->isLeaf() ? ".+_" : "";
	s += "(\\d*(\\.\\d+)?)_(\\d*(\\.\\d+)?)";
	re.assign(s,regex_constants::icase);
	if (!regex_match(name,re))	{
		minmax = 1;
		string s = from->isLeaf() ? ".+_" : "";
		s += "(\\d*(\\.\\d+)?)";
		re.assign(s,regex_constants::icase);
		if (!regex_match(name,re))	{
			minmax = 0;
		}
	}
	if (from->isLeaf())	{
		if (leafval == -1)	{
			leafval = minmax;
		}
		else if (leafval != minmax)	{
			cerr << "error : not all leaves have same format\n";
			cerr << leafval << '\t' << minmax << '\n';
			cerr << name << '\n';
			exit(1);
		}
		else	{
			leafval = minmax;
		}
	}
	else	{
		if (internalval == -1)	{
			internalval = minmax;
		}
		else if (internalval != minmax)	{
			cerr << "error : not all internal nodes have same format\n";
			cerr << name << '\n';
			exit(1);
		}
		else	{
			internalval = minmax;
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Check(link->Out());
	}
}

string FileTree::GetLeafNodeName(const Link* from)	const {

	if (!from->isLeaf())	{
		cerr << "error :get leaf node name\n";
		exit(1);
	}
	string name = from->GetNode()->GetName();
	string s = "(.+)";
	regex re;
	if (HasLeafMinMax())	{
		s += "_(\\d*(\\.\\d+)?)_(\\d*(\\.\\d+)?)";
	}
	else if (HasLeafVal())	{
		s += "_(\\d*(\\.\\d+)?)";
	}
	else	{
	}
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		cerr << "get leaf node name\n";
		cerr << name << '\n';
		exit(1);
	}
	string rets(c[1].first,c[1].second);
	return StringReplace('_'," ",rets);
}

double FileTree::GetNodeVal(const Link* from)	{

	string name = from->GetNode()->GetName();
	string s = from->isLeaf() ? ".+_" : "";
	regex re;
	s += "(\\d*(\\.\\d+)?)";
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		exit(1);
	}
	string rets(c[1].first,c[1].second);
	double ret = atof(rets.c_str());
	return ret;
}


double FileTree::GetNodeMin(const Link* from)	{

	string name = from->GetNode()->GetName();
	string s = from->isLeaf() ? ".+_" : "";
	regex re;
	s += "(\\d*(\\.\\d+)?)_(\\d*(\\.\\d+)?)";
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		exit(1);
	}
	string rets(c[1].first,c[1].second);
	double ret = atof(rets.c_str());
	return ret;
}

double FileTree::GetNodeMax(const Link* from)	{

	string name = from->GetNode()->GetName();
	string s = from->isLeaf() ? ".+_" : "";
	regex re;
	s += "(\\d*(\\.\\d+)?)_(\\d*(\\.\\d+)?)";
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		exit(1);
	}
	string rets(c[3].first,c[3].second);
	double ret = atof(rets.c_str());
	return ret;
}


int DrawTree::RecursiveComputeSize(const Link* from)	{
	if (from->isLeaf())	{
		size[from] = 1;
	}
	else	{
		int tot = 0;
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			tot += RecursiveComputeSize(link->Out());
		}
		size[from] = tot;
	}
	return size[from];
}
		
double DrawTree::RecursiveComputeDepth(const Link* from)	{
	if (from->isLeaf())	{
		depth[from] = GetLength(from);
	}
	else	{
		double max=0;
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = RecursiveComputeDepth(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		max += GetLength(from);
		depth[from] = max;
	}
	return depth[from];
}

void DrawTree::Draw(string target)	{

	string texfile = target + ".tex";
	ofstream os(texfile.c_str());

	if (header)	{
		os << "\\documentclass{";
		if (beamer)	{
			os << "beamer";
		}
		else	{
			os << "amsart";
		}
		os << "}\n";

		os << "\\usepackage{times}\n";
		os << "\\usepackage{pgf}\n";
		os << "\\usepackage[english]{babel}\n";
		os << "\\usepackage{tikz}\n";
		os << "\\usepackage{nopageno}\n";
		os << '\n';

		os << "\\begin{document}\n";
		if (beamer)	{
			os << "\\begin{frame}[plain]\n";
		}
		os << "\\begin{center}\n";
	}
	else	{
		os << "\\begin{center}\n";
	}

	PrepareDrawing();

	double scaleX = sizeX / GetDepth(GetRoot());
	double scaleY = sizeY / (GetSize(GetRoot()) - 1);
	// os << "\\setlength{\\unitlength}{0.6cm}\n";
	os << "\\begin{tikzpicture}\n";
	os << "[anchor=west,font=\\fontsize{" << fontsize << "}{" << fontsize << "}\\selectfont,scale=1]\n";
	// os << "\\begin{tikzpicture}(" << sizeX << ',' << sizeY << ")\n";
	os << '\n';

	RecursiveDraw(GetRoot(),os,0,sizeY/2,scaleX,scaleY);

	FinishDrawing(os,scaleX,scaleY);
	// scale 
	
	os << '\n';
	os << "\\end{tikzpicture}\n";

	if (header)	{
		os << "\\end{center}\n";
		if (beamer)	{
			os << "\\end{frame}\n";
		}
		os << "\\end{document}\n";
	}
	else	{
		os << "\\end{center}\n";
	}
	os.close();
}

void DrawTree::PrepareDrawing()	{

	ComputeSize();
	ComputeDepth();
}



double DrawTree::RecursiveDraw(const Link* from, ostream& os, double X, double Y, double scaleX, double scaleY)	{

	// string format = "";
	// string format = "\\tiny ";
	// double thickness = 0.5;
//	os << "\\linethickness{" << thickness << "mm}\n";

	double ret = 0;

	if (from->isLeaf())	{
		// write species name
		os << "\\path (" << texapprox(X) + shiftname << "," << texapprox(Y) << ") node { \\it " <<  GetLeafNodeName(from) << " };\n";
		// os << "\\path (" << texapprox(X) + 10 * shiftname << "," << texapprox(Y) << ") node { \\it " <<  GetLeafNodeName(from) << " };\n";
		ret = Y;
	}

	else	{

		double y = Y - 0.5 * GetSize(from) * scaleY;

		const Link* link = from->Next();
		double yfirst = 0;
		double ylast = 0;

		while (link != from)	{

			// compute horizontal offset for child node
			double xoffset = GetLength(link) * scaleX;

			// compute vertical span of child node 
			double yspan = GetSize(link->Out()) * scaleY;

			// shift half the vertical offset
			y += 0.5 * yspan;

			// draw child node
			double x = X + xoffset;
			double downret = RecursiveDraw(link->Out(), os, x, y, scaleX, scaleY);

			// draw horizontal line
			LocalDrawBranch(link,os,X,downret,xoffset);

			LocalDraw(link,os,x,downret,scaleX,scaleY);

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
		LocalDrawVerticalTrait(from, os, X, yfirst, ylast);

		ret = (yfirst + ylast) / 2;

		if (from->isRoot())	{
			LocalDraw(from,os,X,ret,scaleX,scaleY);
		}
	}
	return ret;
}

void DrawTree::LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx)	{
	os << "\\path [draw] (" << texapprox(x) << ',' << texapprox(y) << ") -- +(" << texapprox(dx) << ",0);\n";
}

void DrawTree::LocalDrawVerticalTrait(const Link* link, ostream& os, double x, double y1, double y2)	{
	os << "\\path [draw] (" << texapprox(x) << ',' << texapprox(y1) << ") -- +(0," << texapprox(y2-y1) << ");\n";
}

void DrawTree::DrawNodeText(const Link* link, ostream& os, double x, double y)	{
	if (! link->Out()->isLeaf())	{
		// os << "\\path [anchor=west] (" << texapprox(x) << "," << texapprox(y) << ") node {" <<  GetNodeName(link) << " };\n";
	}
}

		
double BubbleTree::ComputeMaxNodeVal(const Link* from)	{
	double max = 0;
	if (from->isLeaf())	{
		max = fabs(GetNodeMax(from));
	}
	else	{
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = ComputeMaxNodeVal(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		// double tmp = GetNodeVal(from);
		double tmp = GetNodeMax(from);
		if (max < tmp)	{
			max = tmp;
		}
	}
	return max;
}
		
		
double BubbleTree::ComputeMinNodeVal(const Link* from)	{
	double min = 0;
	if (from->isLeaf())	{
		min = fabs(GetNodeMin(from));
	}
	else	{
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = ComputeMinNodeVal(link->Out());
			if (((! from->isLeaf()) && (link == from->Next())) || (min > tmp))	{
				min = tmp;
			}
		}
		double tmp = GetNodeMin(from);
		if (min > tmp)	{
			min = tmp;
		}
	}
	return min;
}

void BubbleTree::PrepareDrawing()	{

	/*
	if (minmax)	{
		double max = ComputeMaxNodeVal(root);
		double min = ComputeMinNodeVal(root);
		cerr << min << '\t' << max << '\n';
		// exit(1);
		minnodeval = min;
		nodescale = maxnodeval / sqrt(max - min);
	}
	else	{
	*/
		double max = ComputeMaxNodeVal(GetRoot());
		nodescale = maxnodeval / exp(0.5 * nodepower * log(max));
	// }
}


void BubbleTree::WriteBubbleText(const Link* link, ostream& os, double x, double y)	{
	os << "\\path [anchor=west] (" << texapprox(x) << "," << texapprox(y) << ") node {" <<  "(" << GetNodeMin(link) << "," << GetNodeMax(link) << ")" << " };\n";
}

double BubbleTree::GetCircleNodeMin(const Link* link)	{
	return approx(nodescale * exp(0.5 * nodepower * log(GetNodeMin(link->Out()))));
}

double BubbleTree::GetCircleNodeMax(const Link* link)	{
	return approx(nodescale * exp(0.5 * nodepower * log(GetNodeMax(link->Out()))));
}

double BubbleTree::GetCircleNodeVal(const Link* link)	{
	double d = 0;
	/*
	if (minmax)	{
		double tmp = GetNodeVal(link->Out()) - minnodeval;
		if (tmp < 0)	{
			cerr << "??? " << GetNodeVal(link->Out()) << '\t' << minnodeval << '\n';
			exit(1);
		}
		d = approx(nodescale * sqrt(GetNodeVal(link->Out()) - minnodeval));
	}
	else	{
	*/
		d = approx(nodescale * exp(0.5 * nodepower * log(GetNodeVal(link->Out()))));
	// }
	return d;
}

void BubbleTree::DrawNodeBubble(const Link* link, ostream& os, double x, double y)	{

	if (1)	{
	// if (! link->Out()->isLeaf())	{
		double dmin = GetCircleNodeMin(link);
		double dmax = GetCircleNodeMax(link);
		os << "\\path [fill=red,opacity=0.7,anchor=center] (" << texapprox(x) << ',' << texapprox(y) << ") circle (" << dmin << "pt);\n";
		os << "\\path [fill=red,opacity=0.3,anchor=center] (" << texapprox(x) << ',' << texapprox(y) << ") circle (" << dmax << "pt);\n";
	}
	else	{
		double d = GetCircleNodeVal(link);
		if (d > 0)	{
			os << "\\path [fill=red,opacity=0.8,anchor=center] (" << texapprox(x) << ',' << texapprox(y) << ") circle (" << d << "pt);\n";
		}
		else	{
			os << "\\path [fill=blue,opacity=0.8,anchor=center] (" << texapprox(x) << ',' << texapprox(y) << ") circle (" << -d << "pt);\n";
		}
	}
}

void BubbleTree::FinishDrawing(ostream& os, double xscale, double yscale)	{

	int dmin = (int) (log(ComputeMinNodeVal(GetRoot())) / log(10.0));
	int dmax = (int) (log(ComputeMaxNodeVal(GetRoot())) / log(10.0));
	cerr << dmin << '\t' << dmax << '\n';
	DrawBubbleScale(os, -sizeX / 4, 0.7, dmin, dmax, 1);
	// DrawTree::FinishDrawing();
}

void BubbleTree::DrawBubbleScale(ostream& os, double x, double y, double logmin, double logmax, double logstep)	{

	double delta = nodescale / 15 * exp(0.5 * nodepower * logmax * log(10.0));
	double yy = y;
	for (double d = logmin; d<=logmax; d+= logstep)	{
		double m = nodescale * exp(0.5 * nodepower * d * log(10.0));
		os << "\\path [fill=red,opacity=0.8,anchor=center] (" << x << ',' << sizeY - yy << ") circle (" << m << "pt);\n";
		os << "\\path (" << x << ',' << sizeY - yy << ") node { " <<  ((int) exp(d * log(10.0)))  << " };\n";
		yy += delta;
	}
}

void ChronoDrawTree::DrawTimeInterval(const Link* link, ostream& os, double x, double y, double xscale, double yscale)	{
	if (! link->Out()->isLeaf())	{
		if (link->isRoot())	{
			// cerr << GetDepth(GetRoot()) << '\t' << GetMinTime(GetRoot()) << '\t' << GetMaxTime(GetRoot()) << '\n';
		}
		double tmin = xscale * (GetDepth(GetRoot()) - GetMinTime(link->Out()));
		double tmax = xscale * (GetDepth(GetRoot()) - GetMaxTime(link->Out()));
		double offset = 0.5 * barwidth;
		os << "\\path [fill=blue,opacity=0.6,anchor=center] (" << texapprox(tmin) << ',' << texapprox(y-offset) << ") rectangle (" << texapprox(tmax) << ',' << texapprox(y+offset)  << ");\n";
	}
}

void ChronoDrawTree::DrawTimeScale(ostream& os, double xscale, double yscale)	{

	// x = 0: root
	// sizeX = age of root
	// double rootage = GetDepth(GetRoot());
	double rootmaxage = GetMaxTime(GetRoot());
	// cerr << sizeX << '\t' << xscale * rootage << '\n';
	// int ntics = (int) (rootage / 10);

	double step = 10;

	double y = -0.04 * sizeY;
	double dy= 0.006 * sizeY;
	os << "\\path [draw] (" << sizeX - xscale * (rootmaxage+step)  << "," << y << ") -- (" << sizeX << "," << y << ");\n";
	int tic = 0;
	for (double age=0; age<rootmaxage+step; age += step)	{
		double x = sizeX - xscale * age; 
	// for (int tic=0; tic<ntics; tic++)	{
		// double x = sizeX * (1 - ((double) tic) / ntics);
		os << "\\path [draw] (" << x << "," << y + dy << ") -- +(0," << -2 *dy << ");\n";

		if (! (tic % 10))	{
			os << "\\path (" << x << "," << -8 *dy << ") node[below] {" << tic * 10 << "};\n";
		}
		tic ++;
	}
	os << "\\path (" << sizeX + 0.3 << "," << -8 *dy << ") node[below] {Myrs};\n";
	/*
	os << "\\path[draw,dashed,thin] (" << sizeX * (rootage - 65.0) / rootage << "," << y << ")-- +(0," << sizeY - y << ");\n";
	os << "\\path (" << sizeX * (rootage - 65.0) / rootage << "," <<  - 8*dy << ") node[below] {\\bf KT};\n";
	*/
}

void HeatTree::LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx)	{
	double z2 = GetNodeVal(link->Out());
	double z1 = unibranch ? GetNodeVal(link->Out()) : GetNodeVal(link);
	/*
	double z1 = GetNodeVal(link);
	double z2 = unibranch ? GetNodeVal(link) : GetNodeVal(link->Out());
	*/
	double zmid = 0.5 * (z1 + z2);
	string s1 = GetColorCode(z1);
	string s2 = GetColorCode(z2);
	string smid = GetColorCode(zmid);
	os << "\\shade[left color=";
	os << s1;
	os << ", right color =";
	os << s2;
	os << ", middle color = ";
	os << smid;
	os << "]"; 

	os << " (" << texapprox(x) << "," << texapprox(y-0.5*thickness) << ") rectangle (" << texapprox(x+dx) << ',' << texapprox(y+0.5*thickness) << ");\n";

	/*
	if ((!from->Out()->isLeaf()) && withnodetext)	{
		WriteNodeText(from,os,x,y);
	}
	*/
}

string HeatTree::GetColorCode(double z)	{

	double x = (z - minnodeval) / (maxnodeval - minnodeval);
	if ((x<0) || (x>1))	{
		cerr << "error in color code : " << x << '\n';
		exit(1);
	}
	// cerr << x << '\t';
	ostringstream s;
	s << "red!" << ((int) (100*x)) << "!yellow";
	return s.str();
}
		
double HeatTree::ComputeMaxNodeVal(const Link* from)	{
	double max = 0;
	if (from->isLeaf())	{
		max = fabs(GetNodeVal(from));
	}
	else	{
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = ComputeMaxNodeVal(link->Out());
			if (max < tmp)	{
				max = tmp;
			}
		}
		// double tmp = GetNodeVal(from);
		double tmp = GetNodeVal(from);
		if (max < tmp)	{
			max = tmp;
		}
	}
	return max;
}
		
		
double HeatTree::ComputeMinNodeVal(const Link* from)	{
	double min = 0;
	if (from->isLeaf())	{
		min = fabs(GetNodeVal(from));
	}
	else	{
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = ComputeMinNodeVal(link->Out());
			if (((! from->isLeaf()) && (link == from->Next())) || (min > tmp))	{
				min = tmp;
			}
		}
		double tmp = GetNodeVal(from);
		if (min > tmp)	{
			min = tmp;
		}
	}
	return min;
}

void HeatTree::PrepareDrawing()	{
	DrawTree::PrepareDrawing();
	maxnodeval = ComputeMaxNodeVal(GetRoot());
	minnodeval = ComputeMinNodeVal(GetRoot());
	// cerr << minnodeval << '\t' << maxnodeval << '\n';
}


void HeatTree::WriteNodeText(const Link* link, ostream& os, double x, double y)	{
	os << "\\path [anchor=west] (" << texapprox(x) << "," << texapprox(y) << ") node {" << GetNodeVal(link) << " };\n";
}
