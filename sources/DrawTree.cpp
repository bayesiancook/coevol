
#include "DrawTree.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include "StringStreamUtils.h"



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

	double scaleX = sizeX / GetDepth();
	double startX = (GetDepth() - GetDepth(GetRoot())) * scaleX;
	double scaleY = sizeY / (GetSize(GetRoot()) - 1);
	// os << "\\setlength{\\unitlength}{0.6cm}\n";
	os << "\\begin{tikzpicture}\n";
	os << "[anchor=west,font=\\fontsize{" << fontsize << "}{" << fontsize << "}\\selectfont,scale=" << texscale<< "]\n";
	// os << "\\begin{tikzpicture}(" << sizeX << ',' << sizeY << ")\n";
	os << '\n';

	RecursiveDraw(GetRoot(),os,startX,sizeY/2,scaleX,scaleY);

	RecursiveAfterDraw(GetRoot(),os,startX,sizeY/2,scaleX,scaleY);

	if (Z)	{
		cerr << "write group names\n";
		WriteGroupNames(os);
	}

	FinishDrawing(os,scaleX,scaleY);
	// texscalecale

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

void DrawTree::WriteGroupNames(ostream& os)	{

	for (map<const Node*,string>::iterator i=groupname.begin(); i!=groupname.end(); i++)	{

		const Node* node = i->first;
		string name = i->second;
		if (grouptext[node])	{
			cerr << name << '\n';
			double x = sizeX + Z;
			double y = groupy[node];
			cerr << x << '\t' << y << '\n';
			os << "\\path [anchor=west] (" << texapprox(x) << "," << texapprox(y) << ") node[font=\\fontsize{" << groupfontsize << "}{" << groupfontsize << "}\\selectfont] {" <<  name << " };\n";
		}
	}
}

double DrawTree::RecursiveDraw(const Link* from, ostream& os, double X, double Y, double scaleX, double scaleY)	{

	// string format = "";
	// string format = "\\tiny ";
	// double thickness = 0.5;
//	os << "\\linethickness{" << thickness << "mm}\n";

	double ret = 0;

	if (from->isLeaf())	{
		if (withleafnames)	{
			// write species name
			os << "\\path (" << texapprox(X) + shiftname << "," << texapprox(Y) << ") node { " << GetPreLeafNodeName(from) << " \\,  " << "\\it " <<  GetLeafNodeName(from) << " };\n";
			// os << "\\path (" << texapprox(X) + 10 * shiftname << "," << texapprox(Y) << ") node { \\it " <<  GetLeafNodeName(from) << " };\n";
		}
		ret = Y;
	}

	else	{

		double y = Y - 0.5 * GetSize(from) * scaleY;

		if (groupname[from->GetNode()] != "")	{
			groupy[from->GetNode()] = Y;
		}

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
			LocalDrawBranch(from,os,X,ret,0);
		}
	}
	return ret;
}

double DrawTree::RecursiveAfterDraw(const Link* from, ostream& os, double X, double Y, double scaleX, double scaleY)	{

	// string format = "";
	// string format = "\\tiny ";
	// double thickness = 0.5;
//	os << "\\linethickness{" << thickness << "mm}\n";

	double ret = 0;

	if (from->isLeaf())	{
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
			double downret = RecursiveAfterDraw(link->Out(), os, x, y, scaleX, scaleY);

			LocalAfterDraw(link->Out(),os,x,downret,scaleX,scaleY);

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

		ret = (yfirst + ylast) / 2;

		if (from->isRoot())	{
			LocalAfterDraw(from,os,X,ret,scaleX,scaleY);
		}
	}
	return ret;
}

void DrawTree::LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx)	{
	if (! link->isRoot())	{
	os << "\\path [draw] (" << texapprox(x) << ',' << texapprox(y) << ") -- +(" << texapprox(dx) << ",0);\n";
	}
	if ((!link->Out()->isLeaf()) && withnodetext)	{
		WriteNodeText(link->Out(),os,x+dx,y);
	}
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

	DrawTree::PrepareDrawing();
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
		cerr << dmin << '\t' << dmax << '\n';
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
	cerr << "min max\n";
	cerr << ComputeMinNodeVal(GetRoot()) << '\n';
	cerr << ComputeMaxNodeVal(GetRoot()) << '\n';
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
		if ((! onlygroups) || (groupname[link->Out()->GetNode()] != ""))	{
			if (link->isRoot())	{
				// cerr << GetDepth(GetRoot()) << '\t' << GetMinTime(GetRoot()) << '\t' << GetMaxTime(GetRoot()) << '\n';
			}
			double tmin = xscale * (GetDepth() - GetMinTime(link->Out()));
			double tmax = xscale * (GetDepth() - GetMaxTime(link->Out()));
			/*
			double tmin = xscale * (GetDepth(GetRoot()) - GetMinTime(link->Out()));
			double tmax = xscale * (GetDepth(GetRoot()) - GetMaxTime(link->Out()));
			*/
			double offset = 0.5 * barwidth;
			/*
			int color = 0;
			if (onlygroups)	{
				color = groupcolor[link->Out()->GetNode()];
			}
			*/
			int color = groupcolor[link->Out()->GetNode()];
			if (color == 0)	{
				os << "\\path [fill=blue,opacity=0.5,anchor=center] (" << texapprox(tmin) << ',' << texapprox(y-offset) << ") rectangle (" << texapprox(tmax) << ',' << texapprox(y+offset)  << ");\n";
			}
			else if (color == 1)	{
				offset *= 1.5;
				os << "\\path [fill=blue,opacity=0.6,anchor=center] (" << texapprox(tmin) << ',' << texapprox(y-offset) << ") rectangle (" << texapprox(tmax) << ',' << texapprox(y+offset)  << ");\n";
			}
			else if (color == 2)	{
				offset *= 1.5;
				os << "\\path [fill=orange,opacity=0.6,anchor=center] (" << texapprox(tmin) << ',' << texapprox(y-offset) << ") rectangle (" << texapprox(tmax) << ',' << texapprox(y+offset)  << ");\n";
			}
			else if (color == 3)	{
				offset *= 1.5;
				os << "\\path [fill=red,opacity=0.6,anchor=center] (" << texapprox(tmin) << ',' << texapprox(y-offset) << ") rectangle (" << texapprox(tmax) << ',' << texapprox(y+offset)  << ");\n";
			}
		}
	}
}

void ChronoDrawTree::DrawTimeScale(ostream& os, double xscale, double yscale)	{

	// x = 0: root
	// sizeX = age of root
	double rootmaxage = GetMaxTime();
	// double rootmaxage = GetMaxTime(GetRoot());
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

		if (! (tic % 5))	{
			os << "\\path (" << x << "," << -8 *dy << ") node[below,font=\\fontsize{" << groupfontsize << "}{" << groupfontsize << "}\\selectfont] {" << tic * 10 << "};\n";
		}
		tic ++;
	}
	os << "\\path (" << sizeX + 0.4 << "," << -8 *dy << ") node[below,font=\\fontsize{" << groupfontsize << "}{" << groupfontsize << "}\\selectfont] {Myr};\n";

	// KT
	// double rootage = GetDepth();
	// double rootage = GetDepth(GetRoot());
	// os << "\\path[draw,dashed,thin] (" << sizeX * (rootage - 65.0) / rootage << "," << y << ")-- +(0," << sizeY - y << ");\n";
	// os << "\\path (" << sizeX * (rootage - 65.0) / rootage << "," <<  - 8*dy << ") node[below,font=\\fontsize{" << groupfontsize << "}{" << groupfontsize << "}\\selectfont] {\\bf KPg};\n";
}


void HeatTree::DrawCI(const Link* link, ostream& os, double x, double y)	{

	map<const Node*, double>::iterator imin = extnodemin.find(link->GetNode());
	map<const Node*, double>::iterator imax = extnodemax.find(link->GetNode());
	ostringstream s;

	s.precision(2);

	double z = GetNodeVal(link);
	s << fixed << z << " ";

	s << "(";
	if (imin == extnodemin.end())	{
		s << "?";
	}
	else	{
		double tmp = imin->second;
		s << fixed << tmp;
	}
	s <<  " , ";
	if (imin == extnodemax.end())	{
		s << "?";
	}
	else	{
		double tmp = imax->second;
		s << fixed << tmp;
	}
	s << ")";

	double branchfontsize = 3;
	os << "\\path (" << texapprox(x) << "," << texapprox(y) << ") node[above=-2,font=\\fontsize{" << branchfontsize << "}{" << branchfontsize << "}\\selectfont] {" << s.str() << "};\n";

}

void HeatTree::LocalDrawBranch(const Link* link, ostream& os, double x, double y, double dx)	{
	if (! link->isRoot())	{
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

	if (withbranchval)	{
		double xmid = x + 0.5 * dx;
		map<const Branch*,string>::iterator i = branchval.find(link->GetBranch());
		// map<const Branch*,double>::iterator i = branchval.find(link->GetBranch());
		ostringstream s;
		if (i == branchval.end())	{
			// cerr << "did not find " << GetLeftMost(link) << '\t' << GetRightMost(link) << '\n';
			s << "?";
			// exit(1);
		}
		else	{
			/*
			double v = ((int) (100 * (i->second)));
			if (v)	{
				s << v;
			}
			else	{
				s << "";
				// s << "$<1$";
			}
			*/
			string tmp = i->second;
			for (unsigned int k = 0; k<tmp.length(); k++)	{
				char c = tmp[k];
				if (c == '/')	{
					s << " / ";
				}
				else	{
					s << c;
				}
			}
			// s << i->second;
		}
		double branchfontsize = 4;
		os << "\\path (" << texapprox(xmid) << "," << texapprox(y) << ") node[above=-2,font=\\fontsize{" << branchfontsize << "}{" << branchfontsize << "}\\selectfont] {" << s.str() << "};\n";
	}
	}

	if ((!link->Out()->isLeaf()) && withnodetext)	{
		WriteNodeText(link->Out(),os,x+dx,y);
	}
}

string HeatTree::GetColorCode(double z)	{

	double x = (z - minnodeval) / (maxnodeval - minnodeval);
	if (x > 1)	{
		x = 1;
	}
	if (x < 0)	{
		x = 0;
	}
	if ((x<0) || (x>1))	{
		cerr << "error in color code : " << x << '\n';
		exit(1);
	}
	// cerr << x << '\t';
	ostringstream s;
	s << "red!" << ((int) (100*x)) << "!yellow";
	return s.str();

	/*
	double y = (z - minnodeval) / (maxnodeval - minnodeval);
	double d2 = 1.0 / 3;
	double d1 = 2.0 / 3;
	double x = (y > d2) ? 0.5  + (y-d2)/2/d1 : (y /2/d2);
	ostringstream s;
	s << "red!" << ((int) (100*x)) << "!yellow";
	return s.str();
	*/

	/*
	x *= 2;
	ostringstream s;
	if (x > 1)	{
		// s << "black!" << ((int) (100*(x-1))) << "!red";
		s << "red!" << ((int) (100*(x-1))) << "!yellow";
	}
	else	{
		// s << "red!" << ((int) (100*x)) << "!yellow";
		s << "yellow!" << ((int) (100*x)) << "!green";
	}
	return s.str();
	*/
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
	cerr << minnodeval << '\t' << maxnodeval << '\n';
	if (maxnodeval == -1)	{
		maxnodeval = ComputeMaxNodeVal(GetRoot());
	}
	if (minnodeval == -1)	{
		minnodeval = ComputeMinNodeVal(GetRoot());
	}
	cerr << minnodeval << '\t' << maxnodeval << '\n';
}

void HeatTree::FinishDrawing(ostream& os, double xscale, double yscale)	{
	MakeScale(os);
}

void HeatTree::MakeScale(ostream& os)	{
// void HeatTree::MakeScale(double x, double y, double dx, double min, double max, int ngrad)	{
	double y = -0.04 * sizeY;
	double dy= 0.006 * sizeY;
	double x = 0.04 * sizeX;
	double dx = sizeX / 5;
	double min = minnodeval;
	double max = maxnodeval;
	double z1 = min;
	double z2 = max;
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

	os << " (" << texapprox(x) << "," << texapprox(y-0.5*thickness) << ") rectangle (" << texapprox(x+dx-0.04) << ',' << texapprox(y+0.5*thickness) << ");\n";

	int ngrad = 2;
	for (int i=0; i <= ngrad; i++)	{
		double xx = texapprox(x + dx * ((double) i) / ngrad);
		os << "\\path [draw] (" << xx << "," << y + dy << ") -- +(0," << -2 *dy << ");\n";
		os.precision(3);
		double x = ((double) ((int) (10 * (min + (max-min) * ((double) i) / ngrad)))) / 10;
		os << "\\path (" << xx << "," << -8 *dy << ") node[below,font=\\fontsize{" << groupfontsize << "}{" << groupfontsize << "}\\selectfont] {" << x << "};\n";
	// 	os << "\\path (" << xx << "," << -8 *dy << ") node[below] {" << min + (max - min) * ((double) i) / ngrad  << "};\n";
	}
}


void DrawTree::WriteNodeText(const Link* link, ostream& os, double x, double y)	{
	os << "\\path [anchor=south east,font=\\fontsize{" << fontsize << "}{" << fontsize << "}\\selectfont,scale=" << texscale<< "]";
	// os << "\\path [anchor=west,font=\\fontsize{" << fontsize - 1 << "}{" << fontsize -1 << "}\\selectfont,scale=" << texscale<< "]";
	os << "(" << texapprox(x-0.1) << "," << texapprox(y) << ") node {" << GetNodeName(link) << " };\n";
	// os << "(" << texapprox(x-0.1) << "," << texapprox(y) << ") node {" << ((double) ((int) (10 * GetNodeVal(link)))) / 10 << " };\n";
}

void HeatTree::SetBranchVal(string infile)	{
	ifstream is(infile.c_str());
	if (!is)	{
		cerr << "error in heattree: did not find " << infile << '\n';
		exit(1);
	}

	// int n;
	// is >> n;
	// for (int i=0; i<n; i++)	{
	bool cont = true;
	while (cont)	{
		string name1, name2;
		string tmp;
		// double tmp;
		is >> name1 >> name2 >> tmp;
		cerr << name1 << '\t' << name2 << '\t' << tmp << '\n';
		if (name1 != "END")	{
		const Link* link = GetLCA(name1,name2);
		if (! link)	{
			cerr << "error in heattree: did not find common ancestor of " << name1 << " and " << name2 << '\n';
			exit(1);
		}
		if (! link->isRoot())	{
			// cerr << link << '\t' << name1 << '\t' << name2 << '\t' << tmp << '\n';
			branchval[link->GetBranch()] = tmp;
		}
		}
		else	{
			cerr << "end\n";
			cont = false;
		}
	}
	withbranchval = true;
}

void HeatTree::SetExternalNodeVal(string infile)	{
	ifstream is(infile.c_str());
	if (!is)	{
		cerr << "error in heattree: did not find " << infile << '\n';
		exit(1);
	}
	int n;
	is >> n;
	for (int i=0; i<n; i++)	{
		string name1, name2;
		double tmp;
		is >> name1 >> name2 >> tmp;
		const Link* link = GetLCA(name1,name2);
		if (! link)	{
			cerr << "error in heattree: did not find common ancestor of " << name1 << " and " << name2 << '\n';
			exit(1);
		}
		cerr << link << '\t' << name1 << '\t' << name2 << '\n';
		extnodeval[link->GetNode()] = tmp;
	}
	withexternalnodeval = true;
}

void HeatTree::SetExternalNodeCI(string infile)	{
	ifstream is(infile.c_str());
	if (!is)	{
		cerr << "error in heattree: did not find " << infile << '\n';
		exit(1);
	}
	int n;
	is >> n;
	for (int i=0; i<n; i++)	{
		string name1, name2;
		double tmpmin, tmpmax;
		is >> name1 >> name2 >> tmpmin >> tmpmax;
		const Link* link = GetLCA(name1,name2);
		if (! link)	{
			cerr << "error in heattree: did not find common ancestor of " << name1 << " and " << name2 << '\n';
			exit(1);
		}
		cerr << link << '\t' << name1 << '\t' << name2 << '\n';
		extnodemin[link->GetNode()] = tmpmin;
		extnodemax[link->GetNode()] = tmpmax;
	}
	withexternalnodeci = true;
}

