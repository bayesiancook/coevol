
#include "Random.h"
#include "Ontology.h"

using namespace std;

int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cerr << "redont datafile ontologyfile mingene minconcept outfile\n";
		exit(1);
	}

	string datafile = argv[1];
	string ontologyfile = argv[2];
	int mingene = atoi(argv[3]);
	int minconcept = atoi(argv[4]);
	string outfile = argv[5];

	cerr << "datafile : " << datafile << '\n';
	ifstream is(datafile.c_str());
	int tmpNgene;
	is >> tmpNgene;
	string* tmpgenename = new string[tmpNgene];
	for (int gene=0; gene<tmpNgene; gene++)	{
		is >> tmpgenename[gene];
	}

	int* include = new int[tmpNgene];

	Ontology* tmpontology = new Ontology(ontologyfile);
	Ontology* ontology = new Ontology(tmpontology,tmpNgene,tmpgenename,mingene,minconcept,include);
	int Ngene = 0;
	for (int g=0; g<tmpNgene; g++)	{
		Ngene += include[g];
	}
	if (Ngene != ontology->GetNgene())	{
		cerr << "error : non matching reduxed gene sets\n";
		exit(1);
	}
	string* genename = new string[Ngene];
	int gene = 0;
	for (int g=0; g<tmpNgene; g++)	{
		if (include[g])	{
			genename[gene] = tmpgenename[g];
			gene++;
		}
	}

	cerr << "after reduction: " << ontology->GetNgene() << '\t' << ontology->GetNconcept() << '\n';

	ofstream os((outfile + ".list").c_str());

	os << Ngene << '\n';
	for (int gene=0; gene<Ngene; gene++)	{
		os << genename[gene] << '\n';
	}

	ofstream oos((outfile + ".mat").c_str());
	ontology->ToStream(oos);
	oos.close();
}

