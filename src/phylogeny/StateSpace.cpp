#include <cstdlib>
#include <iostream>
using namespace std;

#include <sstream>
#include "BiologicalSequences.hpp"
#include "StateSpace.hpp"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     StateSpace
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// static inline int EquivalentStrings(string a, string b) {
//     if (a.length() != b.length()) {
//         return 0;
//     }
//     unsigned int k = 0;
//     int cont = 1;
//     while ((k < a.length()) && ((cont) != 0)) {
//         char ca = a[k];
//         char cb = b[k];
//         if ((ca >= 65) && (ca <= 90)) {
//             ca += 32;
//         }
//         if ((cb >= 65) && (cb <= 90)) {
//             cb += 32;
//         }
//         if (ca != cb) {
//             cont = 0;
//         }
//         k++;
//     }
//     return cont;
// }

string SimpleStateSpace::GetState(int state) {
    ostringstream s;
    if (state == unknown) {
        s << "-";
    } else {
        s << Alphabet[state];
    }
    return s.str();
}

int SimpleStateSpace::GetState(string from) {
    if (from.length() != 1) {
        cerr << "error in SingleLetterStateSpace\n";
        exit(1);
    }
    char c = from[0];
    int p = 0;
    while ((p < NAlphabetSet) && (c != AlphabetSet[p])) {
        p++;
    }
    if (p == NAlphabetSet) {
        cout << "error: does not recognise character " << c << '\n';
        exit(1);
    }
    if (p >= 2 * Nstate) {
        return unknown;
    }
    int k = 0;
    for (int l = 0; l < Nstate; l++) {
        if ((c == Alphabet[l]) || (c == Alphabet[l] + 32)) {
            k = l;
        }
    }
    return k;

    return 0;
}

DNAStateSpace::DNAStateSpace() {
    Nstate = 4;
    Alphabet = new char[Nstate];
    for (int i = 0; i < Nstate; i++) {
        Alphabet[i] = DNAletters[i];
    }
    NAlphabetSet = DNAN;
    AlphabetSet = new char[NAlphabetSet];
    for (int i = 0; i < NAlphabetSet; i++) {
        AlphabetSet[i] = DNAset[i];
    }
}

DNAStateSpace::~DNAStateSpace() {
    delete[] Alphabet;
    delete[] AlphabetSet;
}

RYStateSpace::RYStateSpace() {
    Nstate = 2;
    Alphabet = new char[Nstate];
    for (int i = 0; i < Nstate; i++) {
        Alphabet[i] = RYletters[i];
    }
    NAlphabetSet = Nstate;
    AlphabetSet = new char[NAlphabetSet];
    for (int i = 0; i < NAlphabetSet; i++) {
        AlphabetSet[i] = RYletters[i];
    }
}

RYStateSpace::~RYStateSpace() {
    delete[] Alphabet;
    delete[] AlphabetSet;
}

int RYStateSpace::GetRYCoding(int from) {
    if (from == -1) {
        return -1;
    }
    if ((from == 0) || (from == 2)) {
        return 0;
    }
    return 1;
}

RNAStateSpace::RNAStateSpace() {
    Nstate = 4;
    Alphabet = new char[Nstate];
    for (int i = 0; i < Nstate; i++) {
        Alphabet[i] = RNAletters[i];
    }
    NAlphabetSet = RNAN;
    AlphabetSet = new char[NAlphabetSet];
    for (int i = 0; i < NAlphabetSet; i++) {
        AlphabetSet[i] = RNAset[i];
    }
}

RNAStateSpace::~RNAStateSpace() {
    delete[] Alphabet;
    delete[] AlphabetSet;
}

ProteinStateSpace::ProteinStateSpace() {
    Nstate = 20;
    Alphabet = new char[Nstate];
    for (int i = 0; i < Nstate; i++) {
        Alphabet[i] = AminoAcids[i];
    }
    NAlphabetSet = AAN;
    AlphabetSet = new char[NAlphabetSet];
    for (int i = 0; i < NAlphabetSet; i++) {
        AlphabetSet[i] = AAset[i];
    }
}

ProteinStateSpace::~ProteinStateSpace() {
    delete[] Alphabet;
    delete[] AlphabetSet;
}
