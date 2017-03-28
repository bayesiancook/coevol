#ifndef STATESPACE_H
#define STATESPACE_H

#include <string>
#include "BiologicalSequences.hpp"  //FIXME only here because constant unknown

// pure interface
//
class StateSpace {
  public:
    virtual ~StateSpace() = default;

    virtual std::string GetState(int state) = 0;
    virtual int GetNstate() = 0;

    virtual int GetState(std::string from) = 0;

    virtual bool isCompatible(int state1, int state2) {
        return ((state1 == unknown) || (state2 == unknown) || (state1 == state2));
    }
};

// simple state space: assumes that states are referred to using a one-letter
// code
//
class SimpleStateSpace : public StateSpace {
  public:
    int GetState(std::string from) override;

    int GetNstate() override { return Nstate; }

    std::string GetState(int state) override;

  protected:
    int Nstate;
    char *Alphabet;
    int NAlphabetSet;
    char *AlphabetSet;
};

class DNAStateSpace : public SimpleStateSpace {
  public:
    DNAStateSpace();
    ~DNAStateSpace() override;
};

class RNAStateSpace : public SimpleStateSpace {
  public:
    RNAStateSpace();
    ~RNAStateSpace() override;
};

class ProteinStateSpace : public SimpleStateSpace {
  public:
    ProteinStateSpace();
    ~ProteinStateSpace() override;
};

class RYStateSpace : public SimpleStateSpace {
  public:
    RYStateSpace();
    ~RYStateSpace() override;

    int GetRYCoding(int from);
};

class GenericStateSpace : public SimpleStateSpace {
  public:
    GenericStateSpace(int inNstate, char *inAlphabet, int inNAlphabetSet, char *inAlphabetSet) {
        Nstate = inNstate;
        Alphabet = new char[Nstate];
        for (int i = 0; i < Nstate; i++) {
            Alphabet[i] = inAlphabet[i];
        }
        NAlphabetSet = inNAlphabetSet;
        AlphabetSet = new char[NAlphabetSet];
        for (int i = 0; i < NAlphabetSet; i++) {
            AlphabetSet[i] = inAlphabetSet[i];
        }
    }

    ~GenericStateSpace() override {
        delete[] Alphabet;
        delete[] AlphabetSet;
    }
};

#endif  // STATESPACE_H
