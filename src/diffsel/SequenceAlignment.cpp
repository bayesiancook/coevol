#include "diffsel/SequenceAlignment.hpp"
#include <fstream>
#include "diffsel/BiologicalSequences.hpp"
#include "utils/Random.hpp"
#include "utils/StringStreamUtils.hpp"


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     SequenceAlignment
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
int Int(string s) { return atoi(s.c_str()); }

double Double(string s) { return atof(s.c_str()); }

void SequenceAlignment::GetEmpiricalFreq(double* in) {
    int n = GetNstate();
    for (int i = 0; i < GetNstate(); i++) {
        in[i] = 1;
    }
    for (int i = 0; i < GetNtaxa(); i++) {
        for (int j = 0; j < GetNsite(); j++) {
            if (GetState(i, j) != unknown) {
                in[GetState(i, j)]++;
                n++;
            }
        }
    }
    for (int i = 0; i < GetNstate(); i++) {
        in[i] /= n;
    }
}


void SequenceAlignment::GetSiteEmpiricalFreq(double** in, double pseudocount) {
    if (pseudocount < 0) {
        for (int j = 0; j < GetNsite(); j++) {
            for (int i = 0; i < GetNstate(); i++) {
                in[j][i] = -pseudocount * Random::sGamma(1.0);
            }
        }
    } else {
        for (int j = 0; j < GetNsite(); j++) {
            for (int i = 0; i < GetNstate(); i++) {
                in[j][i] = pseudocount;
            }
        }
    }

    for (int i = 0; i < GetNtaxa(); i++) {
        for (int j = 0; j < GetNsite(); j++) {
            if (GetState(i, j) != unknown) {
                in[j][GetState(i, j)]++;
            }
        }
    }
    for (int j = 0; j < GetNsite(); j++) {
        double total = 0;
        for (int i = 0; i < GetNstate(); i++) {
            total += in[j][i];
        }
        for (int i = 0; i < GetNstate(); i++) {
            in[j][i] /= total;
        }
    }
}

void SequenceAlignment::ToFasta(ostream& os) {
    for (int i = 0; i < Ntaxa; i++) {
        os << '>' << taxset->GetTaxon(i) << '\n';
        for (int j = 0; j < Nsite; j++) {
            os << statespace->GetState(GetState(i, j));
        }
        os << '\n';
    }
}
/*
  void SequenceAlignment::ToStream(ostream& os)	{

  os << Ntaxa << '\t' << Nsite << '\n';
  int max = 0;
  for (int i=0; i<Ntaxa; i++)	{
  int l = taxset->GetTaxon(i).length();
  if (max < l)	{
  max = l;
  }
  }

  for (int i=0; i<Ntaxa; i++)	{
  os << taxset->GetTaxon(i);
  for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
  os << ' ';
  }
  for (int j=0; j<Nsite; j++)	{
  os << statespace->GetState(GetState(i,j));
  }
  os << '\n';
  }
  os << '\n';
  }
*/

void SequenceAlignment::ToStream(ostream& os) {
    // os << Ntaxa << '\t' << 876<< '\n';
    os << Ntaxa << '\t' << Nsite << '\n';
    int max = 0;
    for (int i = 0; i < Ntaxa; i++) {
        int l = taxset->GetTaxon(i).length();
        if (max < l) {
            max = l;
        }
    }

    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i);
        for (unsigned int j = 0; j < 5 + max - taxset->GetTaxon(i).length(); j++) {
            os << ' ';
        }
        for (int j = 0; j < Nsite; j++) {
            os << statespace->GetState(GetState(i, j));
        }
        os << '\n';
    }
    os << '\n';
}

int SequenceAlignment::GetNonMissingTriplet() {
    int nsite = 0;
    for (int j = 2; j < Nsite - 3; j += 3) {
        if (NoMissingColumn(j - 1) && NoMissingColumn(j) && NoMissingColumn(j + 1)) {
            nsite++;
        }
    }
    return nsite;
}

void SequenceAlignment::ToStreamTriplet(ostream& os) {
    int nsite = 0;
    for (int j = 2; j < Nsite - 3; j += 3) {
        if (NoMissingColumn(j - 1) && NoMissingColumn(j) && NoMissingColumn(j + 1)) {
            nsite++;
        }
    }

    os << Ntaxa << '\t' << 3 * nsite << '\n';
    int max = 0;
    for (int i = 0; i < Ntaxa; i++) {
        int l = taxset->GetTaxon(i).length();
        if (max < l) {
            max = l;
        }
    }

    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i);
        for (unsigned int j = 0; j < 5 + max - taxset->GetTaxon(i).length(); j++) {
            os << ' ';
        }
        for (int j = 2; j < Nsite - 3; j += 3) {
            if (NoMissingColumn(j - 1) && NoMissingColumn(j) && NoMissingColumn(j + 1)) {
                os << statespace->GetState(GetState(i, j - 1));
                os << statespace->GetState(GetState(i, j));
                os << statespace->GetState(GetState(i, j + 1));
            }
        }
        os << '\n';
    }
    os << '\n';
}

FileSequenceAlignment::FileSequenceAlignment(string filename, int) {  // FIXME unused parameter

    SpeciesNames = 0;
    // cerr << "read data from file : " << filename << "\n";
    ReadDataFromFile(filename, 0);
    taxset = new TaxonSet(SpeciesNames, Ntaxa);
    delete[] SpeciesNames;
    /*
      cerr << "number of taxa  : " << GetNtaxa() << '\n';
      cerr << "number of sites : " << GetNsite() << '\n';
      cerr << "number of states: " << GetNstate() << '\n';
      cerr << "after delete\n";
    */
}

int FileSequenceAlignment::ReadDataFromFile(string filespec, int forceinterleaved) {
    string tmp;
    ifstream is((Path + filespec).c_str());
    if (!is) {
        cerr << "error : cannot find data file " << filespec << '\n';
        cerr << "\n";
        exit(1);
    }
    is >> tmp;
    try {
        if (tmp == "#NEXUS") {
            ReadNexus(filespec);
            return 1;
        } else if (tmp == "#SPECIALALPHABET") {
            ReadSpecial(filespec);
            return 1;
        } else {
            cerr << "phylip\n";
            if (!forceinterleaved) {
                int returnvalue = TestPhylipSequential(filespec);
                if (returnvalue) {
                    cerr << "sequential\n";
                    ReadPhylipSequential(filespec);
                    return 1;
                }
            }
            cerr << "interleaved\n";
            int returnvalue = TestPhylip(filespec, 1);
            if (returnvalue) {
                cerr << "interleaved, taxon names repeated\n";
                ReadPhylip(filespec, 1);
                return 1;
            }
            cerr << "interleaved, taxon names not repeated\n";
            TestPhylip(filespec, 0);
            ReadPhylip(filespec, 0);
            return 1;
        }
    } catch (...) {
        exit(1);
        return 0;
    }
    return 1;
}

int FileSequenceAlignment::ReadNexus(string filespec) {
    ifstream theStream((Path + filespec).c_str());
    try {
        GoPastNextWord(theStream, "dimensions");
        GoPastNext(theStream, '=');
        theStream >> Ntaxa;
        GoPastNext(theStream, '=');
        theStream >> Nsite;
        GoPastNextWord(theStream, "format");
        GoPastNextWord(theStream, "datatype");
        GoPastNext(theStream, '=');
        string type;
        theStream >> type;


        if (EquivalentStrings(type, "protein")) {
            statespace = new ProteinStateSpace();
        } else if (EquivalentStrings(type, "dna")) {
            statespace = new DNAStateSpace();
        } else if (EquivalentStrings(type, "rna")) {
            statespace = new RNAStateSpace();
        } else {
            cerr << "error cannot recognise data type\n";
            cerr << type << "\n";
            exit(1);
        }

        if (Data) {
            for (int i = 0; i < Ntaxa; i++) {
                delete Data[i];
            }
            delete[] Data;
        }
        Data = new int*[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }

        SpeciesNames = new string[Ntaxa];

        GoPastNextWord(theStream, "Matrix");

        int l = 0;
        while (l < Nsite) {
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                string temp;
                theStream >> temp;
                while (temp == "[") {
                    unsigned char c;
                    c = 'i';
                    while (c != ']') c = theStream.get();
                    theStream >> temp;
                }

                if (!l) {
                    SpeciesNames[i] = temp;
                } else {
                    if (temp != SpeciesNames[i]) {
                        cerr << "error when reading tree base: " << temp << '\t' << SpeciesNames[i]
                             << '\n';
                        exit(1);
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if (c == '[') {
                        while (c != ']') c = theStream.get();
                        c = theStream.get();
                    }
                    if ((c != ' ') && (c != '\t') && (c != '\n') && (c != 13)) {
                        if (c == '(') {
                            Data[i][k] = unknown;
                            while (c != ')') {
                                theStream >> c;
                            }
                        } else if (c == '{') {
                            Data[i][k] = unknown;
                            while (c != '}') {
                                theStream >> c;
                            }
                        } else {
                            ostringstream s;
                            s << c;
                            Data[i][k] = statespace->GetState(s.str());
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                if (!m) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "error when reading nexus : " << m << '\t' << k << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        if (m > k) {
                            while (k != m) {
                                Data[i][k] = unknown;
                                k++;
                            }
                        }
                    }
                }
            }
            l = m;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return 1;
}

// ---------------------------------------------------------------------------
//     ReadSpecial()
// ---------------------------------------------------------------------------


int FileSequenceAlignment::ReadSpecial(string filename) {
    int returnvalue = 0;
    try {
        ifstream theStream(filename.c_str());

        string tmp;
        theStream >> tmp;
        theStream >> Ntaxa;
        theStream >> Nsite;
        theStream >> tmp;
        cerr << tmp << '\n';
        int Nstate = tmp.length();
        char* Alphabet = new char[Nstate];
        int NAlphabetSet = Nstate + 5;
        char* AlphabetSet = new char[NAlphabetSet];
        cerr << "alphabet size : " << Nstate << '\n';
        cerr << "alphabet : ";
        for (int i = 0; i < Nstate; i++) {
            Alphabet[i] = tmp[i];
            AlphabetSet[i] = tmp[i];
            cerr << Alphabet[i] << ' ';
        }
        cerr << '\n';

        AlphabetSet[Nstate] = '?';
        AlphabetSet[Nstate + 1] = '-';
        AlphabetSet[Nstate + 2] = '*';
        AlphabetSet[Nstate + 3] = 'X';
        AlphabetSet[Nstate + 4] = 'x';

        statespace = new GenericStateSpace(Nstate, Alphabet, NAlphabetSet, AlphabetSet);

        if (Data) {
            for (int i = 0; i < Ntaxa; i++) {
                delete Data[i];
            }
            delete[] Data;
        }
        Data = new int*[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }

        delete[] SpeciesNames;
        SpeciesNames = new string[Ntaxa];

        int ntaxa = 0;
        string temp;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            char c;
            do {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != ')') {
                            theStream >> c;
                        }
                    } else if (c == '{') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != '}') {
                            theStream >> c;
                        }
                    } else {
                        int p = 0;
                        while ((p < NAlphabetSet) && (c != AlphabetSet[p])) p++;
                        if (p == NAlphabetSet) {
                            cout << "error: does not recognise character. taxon " << ntaxa << '\t'
                                 << SpeciesNames[ntaxa] << "  site  " << nsite << '\t' << c << '\n';
                            exit(1);
                        }
                        if (p >= Nstate) {
                            Data[ntaxa][nsite] = unknown;
                        } else {
                            for (int l = 0; l < Nstate; l++) {
                                if (c == Alphabet[l]) {
                                    Data[ntaxa][nsite] = l;
                                }
                            }
                        }
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            ntaxa++;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return returnvalue;
}

// ---------------------------------------------------------------------------
//     ReadPhylip()
// ---------------------------------------------------------------------------

int FileSequenceAlignment::TestPhylipSequential(string filespec) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = atoi(temp.c_str());
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = atoi(temp.c_str());

        delete[] SpeciesNames;
        SpeciesNames = new string[Ntaxa];

        int AAcomp = 1;
        int DNAcomp = 1;
        int RNAcomp = 1;

        int ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            char c = ' ';
            do {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        while (c != ')') {
                            theStream >> c;
                        }
                    } else if (c == '{') {
                        while (c != '}') {
                            theStream >> c;
                        }
                    } else {
                        int p = 0;
                        if (DNAcomp) {
                            while ((p < DNAN) && (c != DNAset[p])) p++;
                            if (p == DNAN) {
                                DNAcomp = 0;
                            }
                        }
                        p = 0;
                        if (RNAcomp) {
                            while ((p < RNAN) && (c != RNAset[p])) p++;
                            if (p == RNAN) RNAcomp = 0;
                        }
                        p = 0;
                        if (AAcomp) {
                            while ((p < AAN) && (c != AAset[p])) p++;
                            if (p == AAN) {
                                AAcomp = 0;
                            }
                        }
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            if (theStream.eof()) {
                if (nsite < Nsite) {
                    return 0;
                }
            }
            ntaxa++;
        }
        if (theStream.eof()) {
            if (ntaxa < Ntaxa) {
                return 0;
            }
        }
        if (DNAcomp) {
            statespace = new DNAStateSpace;
        } else if (RNAcomp) {
            statespace = new RNAStateSpace;
        } else if (AAcomp) {
            statespace = new ProteinStateSpace;
        } else {
            return 0;
        }
    } catch (...) {
        return 0;
    }
    return 1;
}


void FileSequenceAlignment::ReadPhylipSequential(string filespec) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = Int(temp);

        Data = new int*[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }
        SpeciesNames = new string[Ntaxa];

        int ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            char c;
            do {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != ')') {
                            theStream >> c;
                        }
                    } else if (c == '{') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != '}') {
                            theStream >> c;
                        }
                    } else {
                        ostringstream s;
                        s << c;
                        Data[ntaxa][nsite] = statespace->GetState(s.str());
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            ntaxa++;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        exit(1);
    }
}

int FileSequenceAlignment::TestPhylip(string filespec, int repeattaxa) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = Int(temp);

        delete[] SpeciesNames;
        SpeciesNames = new string[Ntaxa];

        int AAcomp = 1;
        int DNAcomp = 1;
        int RNAcomp = 1;

        int l = 0;
        int block = 0;
        while (l < Nsite) {
            block++;
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                if ((!l) || repeattaxa) {
                    string temp;
                    theStream >> temp;
                    if (!l) {
                        SpeciesNames[i] = temp;
                    } else {
                        if (temp != SpeciesNames[i]) {
                            return 0;
                        }
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') &&
                        (c != 13)) {
                        if (c == '(') {
                            while (c != ')') {
                                theStream >> c;
                            }
                        } else if (c == '{') {
                            while (c != '}') {
                                theStream >> c;
                            }
                        } else {
                            int p = 0;
                            if (DNAcomp) {
                                while ((p < DNAN) && (c != DNAset[p])) p++;
                                if (p == DNAN) DNAcomp = 0;
                            }
                            p = 0;
                            if (RNAcomp) {
                                while ((p < RNAN) && (c != RNAset[p])) p++;
                                if (p == RNAN) RNAcomp = 0;
                            }
                            p = 0;
                            if (AAcomp) {
                                while ((p < AAN) && (c != AAset[p])) p++;
                                if (p == AAN) {
                                    AAcomp = 0;
                                }
                            }
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13) && (c != 10));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13))) {
                    c = theStream.get();
                    c = theStream.peek();
                }
                if (!m) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "in test phylip\n";
                        cerr << "error when reading data non matching number of sequences in block "
                                "number "
                             << block << " for taxon " << i + 1 << " " << SpeciesNames[i] << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l = m;
        }
        if (l < Nsite) {
            cerr << "error : reached end of stream \n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        if (DNAcomp) {
            statespace = new DNAStateSpace;
            // cerr << "dna sequences\n";
        } else if (RNAcomp) {
            statespace = new DNAStateSpace;
            // cerr << "rna sequences\n";
        } else if (AAcomp) {
            statespace = new ProteinStateSpace;
            // cerr << "protein sequences\n";
        } else {
            // cerr << "format not recognised\n";
            return 0;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return 1;
}

void FileSequenceAlignment::ReadPhylip(string filespec, int repeattaxa) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (!IsInt(temp)) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = Int(temp);
        // cerr << Ntaxa << '\t' << Nsite << '\n';

        Data = new int*[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }
        SpeciesNames = new string[Ntaxa];

        int l = 0;
        int block = 0;
        while (l < Nsite) {
            block++;
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                if ((!l) || repeattaxa) {
                    string temp;
                    theStream >> temp;
                    if (!l) {
                        SpeciesNames[i] = temp;
                    } else {
                        if (temp != SpeciesNames[i]) {
                            cerr << "error when reading data: read " << temp << " instead of "
                                 << SpeciesNames[i] << '\n';
                            exit(1);
                        }
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') &&
                        (c != 13)) {
                        if (c == '(') {
                            Data[i][k] = unknown;
                            while (c != ')') {
                                theStream >> c;
                            }
                        } else if (c == '{') {
                            Data[i][k] = unknown;
                            while (c != '}') {
                                theStream >> c;
                            }
                        } else {
                            ostringstream s;
                            s << c;
                            Data[i][k] = statespace->GetState(s.str());
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13))) {
                    c = theStream.get();
                    c = theStream.peek();
                }

                if (!m) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "error when reading data non matching number of sequences in block "
                                "number "
                             << block << " for taxon " << i << " " << SpeciesNames[i] << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l = m;
        }
        if (l < Nsite) {
            cerr << "error : reached end of stream \n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
    } catch (...) {
        cerr << "error while reading data file\n";
    }
}

// ---------------------------------------------------------------------------
//     EliminateUnnownColumns
// ---------------------------------------------------------------------------


/*
  void FileSequenceAlignment::EliminateUnknownColumns()	{

  cerr << "eliminate unknown columns\n";

  // delete all the sites that are '-' for every species

  int i=0;
  int j=0;
  int Eliminated = 0;
  while (i<Nsite)	{
  int k=0;
  Boolean test = true;
  while (test &&  k<Ntaxa)	{
  test &= (Data[k][i] == unknown);
  k++;
  }
  if (! test)	{
  for (int k=0; k<Ntaxa; k++)	{
  Data[k][j] = Data[k][i];
  }
  j++;
  }
  else	{
  GeneSize[Gene[i]]--;
  Eliminated ++;
  }
  i++;
  }
  Nsite -= Eliminated;
  if (Eliminated)	{
  cerr << Eliminated << " columns completely undetermined (full of gaps, or unknown): eliminated\n";
  }
  }
*/

/*
// ---------------------------------------------------------------------------
//     EstimateEmpiricalFrequencies()
// ---------------------------------------------------------------------------

void FileSequenceAlignment::EstimateEmpiricalFrequencies()	{

cerr << "estimate emp freq\n";
for (int k=0; k<Nstate; k++)	{
EmpiricalFreq[k] = 1;
}
for (int i=0; i<Ntaxa; i++)	{
for (int j=0; j<Nsite; j++)	{
if (Data[i][j] != unknown)	{
if ( (Data[i][j] < 0 ) || (Data[i][j] >= Nstate) )	{
cerr << "error in data matrix : " << Data[i][j] << '\n';
exit(1);
}
EmpiricalFreq[Data[i][j]] += 1.0;
}
}
}
double total = 0;
for (int k = 0; k<Nstate; k++)	{
total += EmpiricalFreq[k];
}
for (int k=0; k<Nstate; k++)	{
EmpiricalFreq[k] /= total;
}
}

// ---------------------------------------------------------------------------
//     ConstantColumns()
// ---------------------------------------------------------------------------

int FileSequenceAlignment::ConstantColumns()	{

int n=0;
for (int i=0; i<Nsite; i++)	{
int j = 0;
int test = 1;
while ( (j<Ntaxa) && (Data[j][i] == unknown))	{
j++;
}
if (j != Ntaxa)	{
int k = j+1;
while (k < Ntaxa)	{
if (Data[k][i] != unknown)	{
if (Data[j][i] != Data[k][i])	{
test = 0;
}
}
k++;
}
}
if (test)	{
n++;
}
}
return n;
}


// ---------------------------------------------------------------------------
//     EliminateConstantPositions()
// ---------------------------------------------------------------------------

void FileSequenceAlignment::EliminateConstantPositions()	{

cerr << "eliminate constant positions\n";
int i=0;
int j=0;
int Eliminated = 0;
while (i<Nsite)	{
int k = 0;
while ((k<Ntaxa) && (Data[k][i] == unknown)) k++;
if (k<Ntaxa)	{
int a = Data[k][i];
k++;
while ((k<Ntaxa) && ((Data[k][i] == unknown) || (Data[k][i] == a))) k++;
if (k==Ntaxa)	{
Eliminated ++;
GeneSize[Gene[i]]--;
}
else	{
for (int k=0; k<Ntaxa; k++)	{
Data[k][j] = Data[k][i];
}
j++;
}
}
i++;
}

Nsite -= Eliminated;
cout << "number of positions eliminated : " << Eliminated << '\n';
WriteDataToFile("woconst.puz");
}


// ---------------------------------------------------------------------------
//     RegisterWithData()
// ---------------------------------------------------------------------------

void FileSequenceAlignment::RegisterWithData()	{

// test that no taxon name appears twice
for (int i=0; i<Ntaxa; i++)     {
for (int j=i+1; j<Ntaxa; j++)   {
if (SpeciesNames[i] == SpeciesNames[j]) {
cerr << "error: taxa " << i << " and " << j << " in datafile have same name\n";
exit(1);
}
}
}

if (Ngene == 0)	{
Ngene = 1;
Gene = new int[Nsite];
for (int i=0; i<Nsite; i++)	{
Gene[i] =0;
}
GeneSize = new int[Ngene];
GeneSize[0] = Nsite;
GeneFirstSite = new int[Ngene];
GeneFirstSite[0] = 0;
}

if (Recoding)	{
RecodeData();
}
EliminateUnknownColumns();
int total = 0;
for (int i=0; i<Ngene; i++)	{
GeneFirstSite[i] = total;
int temp = GeneSize[i];
for (int j=0; j<temp; j++)	{
Gene[j+total] = i;
}
total += temp;
}
if (total != Nsite)	{
cerr << "total number of sites does not match\n";
exit(1);
}
if (DeleteConstant)	{
EliminateConstantPositions();
int total = 0;
for (int i=0; i<Ngene; i++)	{
GeneFirstSite[i] = total;
int temp = GeneSize[i];
for (int j=0; j<temp; j++)	{
Gene[j+total] = i;
}
total += temp;
}
if (total != Nsite)	{
cerr << "total number of sites does not match\n";
exit(1);
}
}
ComputeZipArrays();
EstimateEmpiricalFrequencies();
int nconst = 0;
for (int i=0; i<Nsite; i++)	{
if (OrbitSize[i] == 1)	{
nconst ++;
}
}
// cerr << "empirical proportion of invariable sites : " << ((double) nconst) / Nsite << '\n';
}

// ---------------------------------------------------------------------------
//     ComputeZipArrays()
// ---------------------------------------------------------------------------

void FileSequenceAlignment::ComputeZipArrays()	{

EmpiricalFreq = new double[Nstate];
for (int i=0; i<Nstate; i++)	{
EmpiricalFreq[i] = 1.0 / Nstate;
}
SiteEmpiricalCount = new int*[Nsite];
for (int i=0; i<Nsite; i++)	{
SiteEmpiricalCount[i] = new int[Nstate];
for (int j=0; j<Nstate; j++)	{
SiteEmpiricalCount[i][j] = 0;
}
for (int j=0; j<Ntaxa; j++)	{
if (Data[j][i] != unknown)	{
SiteEmpiricalCount[i][Data[j][i]]++;
}
}
}

Indices = new int*[Nsite];
ZipIndices = new int*[Nsite];
for (int i=0; i<Nsite; i++)	{
Indices[i] = new int[Nstate];
ZipIndices[i] = new int[Nstate];
}

ZipSize = new int[Nsite];

OrbitSize = new int[Nsite];
Orbit = new int *[Nsite];
for (int i=0; i<Nsite; i++)	{
Orbit[i] = new int[Nstate];
}

ZipData = new int*[Ntaxa];
for (int i=0; i<Ntaxa; i++)	{
ZipData[i] = new int[Nsite];
}

for (int i=0; i<Nsite; i++)	{

OrbitSize[i] = 0;
for (int k=0 ; k< Nstate; k++)	{
Orbit[i][k]= false;
}

for (int j=0; j<Ntaxa; j++)	{
int d = Data[j][i];
if (d != unknown)	{
if (! Orbit[i][d])	{
Orbit[i][d] = true;
Indices[i][OrbitSize[i]] = d;
OrbitSize[i] ++;
}
}
}

// sort Indices[i]
for (int j=0; j<OrbitSize[i]; j++)	{
for (int k=OrbitSize[i]-1; k>j; k--)	{
if (Indices[i][j] > Indices[i][k])	{
int tmp = Indices[i][j];
Indices[i][j] = Indices[i][k];
Indices[i][k] = tmp;
}
}
}

// reverse translation table
for (int j=0; j<OrbitSize[i]; j++)	{
ZipIndices[i][Indices[i][j]] = j;
}
for (int j=0; j<Nstate; j++)	{
if (! Orbit[i][j])	{
ZipIndices[i][j] = OrbitSize[i];
}
}

if (OrbitSize[i] == 0)	{
cerr << "PhyloParameters::RegisterWithData : missing column";
// exit(1);
}

if (OrbitSize[i] < Nstate)	{
ZipSize[i] = OrbitSize[i] + 1;
}
else	{
ZipSize[i] = OrbitSize[i];
}

for (int j=0; j<Ntaxa; j++)	{
ZipData[j][i] = -2;
int d = (Data)[j][i];
if (d == unknown)	{
ZipData[j][i] = unknown;
}
else	{
for (int k=0; k<OrbitSize[i]; k++)	{
if (Indices[i][k] == d)	{
ZipData[j][i] = k;
}
}

// Zip identity:
// Indices[i][ZipData[j][i]] == mParam->Data[j][i] , for every i and j
// within their respective range

}

// here, may be check that ZipData != -2

if (ZipData[j][i] == -2)	{
cerr << "PhyloParameters::RegisterWithData : error in zip data making\n";
}
}

int observed[Nstate];
int orbitsize= OrbitSize[i];
for (int k=0; k<Nstate; k++)	{
observed[k] = 0;
}
for (int k=0; k<orbitsize; k++)	{
observed[Indices[i][k]] = 1;
}
for (int k=0; k<Nstate; k++)	{
if (! observed[k])	{
Indices[i][orbitsize++] = k;
}
}
if (orbitsize != Nstate)	{
cerr << "error in FileSequenceAlignment::RegisterWithData\n";
cerr << "site : " << i << '\n';
cerr << '\n';
exit(1);
}
}

double temp = 0;
for (int i=0; i<Nsite; i++)	{
temp += ((double) ZipSize[i] * ZipSize[i]) / Nstate / Nstate;
}
SpeedFactor =  temp / Nsite ;
cerr << "speed factor : " << SpeedFactor << '\n';
int nconst = 0;
for (int i=0; i<Nsite; i++)	{
if (OrbitSize[i] == 1)	{
nconst++;
}
}
cerr << "number of constant columns: " << nconst << " (" << ((double) nconst) / Nsite * 100 << ")"
<< '\n';

}




*/
