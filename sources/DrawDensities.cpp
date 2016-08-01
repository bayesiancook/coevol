
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>

using namespace std;

int main(int argc, char* argv[])	{

	ifstream taxis(argv[1]);
	ifstream taxis2(argv[2]);
	string tabcov = argv[3];
	string tabdiag = argv[4];

	string basename = argv[5];

	double z = atof(argv[6]);

	ofstream os((basename + ".R").c_str());

	os << "cov <- read.table(\"" << tabcov << "\")\n";
	os << "diag <- read.table(\"" << tabdiag << "\")\n";
	os << "tcov <- t(cov)\n";
	os << "tdiag <- t(diag)\n";
	os << '\n';


	int ntaxa;
	taxis >> ntaxa;
	// string* leftname = new string[ntaxa];
	// string* rightname = new string[ntaxa];
	string* name = new string[ntaxa];
	map<pair<string,string>,string> c;
	for (int i=0; i<ntaxa; i++)	{
		string name1, name2, name3;
		taxis >> name1 >> name2 >> name3; // leftname[i] >> rightname[i] >> name[i];
		c[pair<string,string>(name1,name2)] = name3;
	}
	for (int i=0; i<ntaxa; i++)	{
		string name1, name2;
		taxis2 >> name1 >> name2;
		name[i] = c[pair<string,string>(name1,name2)];
		cerr << name[i] << '\n';
	}

	int n = (int) (sqrt(ntaxa));
	if (n*n != ntaxa)	{
		cerr << "error : not a perfect square\n";
		exit(1);
	}

	for (int i=0; i<n; i++)	{
		for (int j=0; j<n; j++)	{
			int k = i*n + j;
			os << "pdf(\"" << basename << "_" << k+1 << ".pdf\")\n";
			os << "mcov <- mean(tcov[," << k+1 << "])\n";
			os << "mdiag <- mean(tdiag[," << k+1 << "])\n";
			os << "scov <- sqrt(var(tcov[," << k+1 << "]))\n";
			os << "sdiag <- sqrt(var(tdiag[," << k+1 << "]))\n";
			os << "mmin <- min(mcov,mdiag)\n";
			os << "mmax <- max(mcov,mdiag)\n";
			os << "xmin <- mmin - " << z << "*0.5*(scov+sdiag)\n";
			os << "xmax <- mmax + " << z << "*0.5*(scov+sdiag)\n";
			os << "ymax <- 0.9 / (scov + sdiag)\n";
			os << "plot(density(tcov[," << k+1 << "], adjust=2), ylim=range(0,ymax), xlim=range(xmin,xmax), col=\"black\", cex.axis = 2, cex.lab=2.5, cex.main=3, lwd=3,";
			// os << "plot(density(tcov[," << k+1 << "], adjust=2), ylim=range(0,ymax), xlim=range(xmin,xmax), col=\"red\", cex.axis = 2, cex.lab=2.5, cex.main=3, lwd=3,";
			if (i==n-1)	{
				os << " xlab=expression(paste(log[10], \" Mass (g)\")),";
			}
			else	{
				os << " xlab = \"\",";
			}
			if (!j)	{
				os << " ylab = \"\",";
				// os << " ylab=\"density\",";
			}
			else	{
				os << " ylab = \"\",";
			}
			// os << "plot(density(tcov[," << k+1 << "], adjust=2), xlim=range(" << xmin << "," << xmax << "), ylim=range(" << ymin << "," << ymax << "), col=\"red\", cex.lab=1.5, cex.main=2, lwd=3,";
			os << "main=\"" << name[k] << "\")\n";
			os << "box()\n";
			os << "lines(density(tdiag[," << k+1 << "],adjust=2) ,col=\"black\", lwd=3, lt=2)\n";
			// os << "lines(density(tdiag[," << k+1 << "],adjust=2) ,col=\"blue\", lwd=3)\n";
			os << "dev.off()\n";
			os << '\n';
		}
	}
}

