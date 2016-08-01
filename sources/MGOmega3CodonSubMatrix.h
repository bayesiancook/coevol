
#ifndef OMEGA3CODON
#define OMEGA3CODON

#include "CodonSubMatrix.h"


// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGOmega2CodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGOmega2CodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double inomegats, double inomegatv, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise) ,
			omegats(inomegats), omegatv(inomegatv) {}


	double			GetOmegaTs() {return omegats;}
	double			GetOmegaTv() {return omegatv;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void			SetOmegaTs(double inomegats) {omegats = inomegats;}
	void			SetOmegaTv(double inomegatv) {omegatv = inomegatv;}

	// data members

	double omegats;
	double omegatv;
};

class RandomMGOmega2CodonSubMatrix : public MGOmega2CodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGOmega2CodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inRandomOmegaTs, Var<PosReal>* inRandomOmegaTv, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MGOmega2CodonSubMatrix(instatespace,inmatrix,inRandomOmegaTs->val(), inRandomOmegaTv->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix),
			RandomOmegaTs(inRandomOmegaTs),
			RandomOmegaTv(inRandomOmegaTv) {

		Register(matrix);
		Register(RandomOmegaTs);
		Register(RandomOmegaTv);
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetOmegaTs(RandomOmegaTs->val());
		SetOmegaTv(RandomOmegaTv->val());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* RandomOmegaTs;
	Var<PosReal>* RandomOmegaTv;
};

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGOmega3CodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGOmega3CodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double inomegats, double inomegatv0, double inomegatvgc, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise) ,
			omegats(inomegats), omegatv0(inomegatv0), omegatvgc(inomegatvgc)  {}


	double			GetOmegaTs() {return omegats;}
	double			GetOmegaTv0() {return omegatv0;}
	double			GetOmegaTvGC() {return omegatvgc;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void			SetOmegaTs(double inomegats) {omegats = inomegats;}
	void			SetOmegaTv0(double inomegatv0) {omegatv0 = inomegatv0;}
	void			SetOmegaTvGC(double inomegatvgc) {omegatvgc = inomegatvgc;}

	// data members

	double omegats;
	double omegatv0;
	double omegatvgc;
};

class RandomMGOmega3CodonSubMatrix : public MGOmega3CodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGOmega3CodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inRandomOmegaTs, Var<PosReal>* inRandomOmegaTv0, Var<PosReal>* inRandomOmegaTvGC, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MGOmega3CodonSubMatrix(instatespace,inmatrix,inRandomOmegaTs->val(), inRandomOmegaTv0->val(), inRandomOmegaTvGC->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix),
			RandomOmegaTs(inRandomOmegaTs),
			RandomOmegaTv0(inRandomOmegaTv0),
			RandomOmegaTvGC(inRandomOmegaTvGC) {

		Register(matrix);
		Register(RandomOmegaTs);
		Register(RandomOmegaTv0);
		Register(RandomOmegaTvGC);
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetOmegaTs(RandomOmegaTs->val());
		SetOmegaTv0(RandomOmegaTv0->val());
		SetOmegaTvGC(RandomOmegaTvGC->val());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* RandomOmegaTs;
	Var<PosReal>* RandomOmegaTv0;
	Var<PosReal>* RandomOmegaTvGC;
};

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGOmega3X3CodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGOmega3X3CodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double insynratets, double insynratetvgc, double inomegats, double inomegatv0, double inomegatvgc, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise) ,
			synratets(insynratets), synratetvgc(insynratetvgc),
			omegats(inomegats), omegatv0(inomegatv0), omegatvgc(inomegatvgc)  {}


	double			GetSynRateTs() {return synratets;}
	double			GetSynRateTvGC() {return synratetvgc;}

	double			GetOmegaTs() {return omegats;}
	double			GetOmegaTv0() {return omegatv0;}
	double			GetOmegaTvGC() {return omegatvgc;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void			SetOmegaTs(double inomegats) {omegats = inomegats;}
	void			SetOmegaTv0(double inomegatv0) {omegatv0 = inomegatv0;}
	void			SetOmegaTvGC(double inomegatvgc) {omegatvgc = inomegatvgc;}

	void			SetSynRateTs(double insynratets) {synratets = insynratets;}
	void			SetSynRateTvGC(double insynratetvgc) {synratetvgc = insynratetvgc;}

	// data members

	double synratets;
	double synratetvgc;

	double omegats;
	double omegatv0;
	double omegatvgc;
};

class RandomMGOmega3X3CodonSubMatrix : public MGOmega3X3CodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGOmega3X3CodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inRandomSynRateTs, Var<PosReal>* inRandomSynRateTvGC, Var<PosReal>* inRandomOmegaTs, Var<PosReal>* inRandomOmegaTv0, Var<PosReal>* inRandomOmegaTvGC, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MGOmega3X3CodonSubMatrix(instatespace,inmatrix, inRandomSynRateTs->val(), inRandomSynRateTvGC->val(), inRandomOmegaTs->val(), inRandomOmegaTv0->val(), inRandomOmegaTvGC->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix),
			RandomSynRateTs(inRandomSynRateTs),
			RandomSynRateTvGC(inRandomSynRateTvGC),
			RandomOmegaTs(inRandomOmegaTs),
			RandomOmegaTv0(inRandomOmegaTv0),
			RandomOmegaTvGC(inRandomOmegaTvGC) {

		Register(matrix);
		Register(RandomSynRateTs);
		Register(RandomSynRateTvGC);
		Register(RandomOmegaTs);
		Register(RandomOmegaTv0);
		Register(RandomOmegaTvGC);
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetSynRateTs(RandomSynRateTs->val());
		SetSynRateTvGC(RandomSynRateTvGC->val());
		SetOmegaTs(RandomOmegaTs->val());
		SetOmegaTv0(RandomOmegaTv0->val());
		SetOmegaTvGC(RandomOmegaTvGC->val());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* RandomSynRateTs;
	Var<PosReal>* RandomSynRateTvGC;
	Var<PosReal>* RandomOmegaTs;
	Var<PosReal>* RandomOmegaTv0;
	Var<PosReal>* RandomOmegaTvGC;
};

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
class MGOmega2X2CodonSubMatrix : public MGCodonSubMatrix	{

	public:

	MGOmega2X2CodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, double insynratets, double inomegats, double inomegatv0, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			MGCodonSubMatrix(instatespace, inNucMatrix, innormalise) ,
			synratets(insynratets),
			omegats(inomegats), omegatv0(inomegatv0) {}


	double			GetSynRateTs() {return synratets;}

	double			GetOmegaTs() {return omegats;}
	double			GetOmegaTv0() {return omegatv0;}

	protected:

	// look at how ComputeArray and ComputeStationary are implemented in CodonSubMatrix.cpp
	void 			ComputeArray(int state);
	void			SetOmegaTs(double inomegats) {omegats = inomegats;}
	void			SetOmegaTv0(double inomegatv0) {omegatv0 = inomegatv0;}

	void			SetSynRateTs(double insynratets) {synratets = insynratets;}

	// data members

	double synratets;

	double omegats;
	double omegatv0;
};

class RandomMGOmega2X2CodonSubMatrix : public MGOmega2X2CodonSubMatrix, public RandomCodonSubMatrix	{

	public:

	RandomMGOmega2X2CodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<PosReal>* inRandomSynRateTs, Var<PosReal>* inRandomOmegaTs, Var<PosReal>* inRandomOmegaTv0, bool innormalise = false) :
			SubMatrix(instatespace->GetNstate(), innormalise),
			CodonSubMatrix(instatespace,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
			MGOmega2X2CodonSubMatrix(instatespace,inmatrix, inRandomSynRateTs->val(), inRandomOmegaTs->val(), inRandomOmegaTv0->val(), innormalise) ,
			RandomCodonSubMatrix(instatespace, innormalise),
			matrix(inmatrix),
			RandomSynRateTs(inRandomSynRateTs),
			RandomOmegaTs(inRandomOmegaTs),
			RandomOmegaTv0(inRandomOmegaTv0)	{

		Register(matrix);
		Register(RandomSynRateTs);
		Register(RandomOmegaTs);
		Register(RandomOmegaTv0);
		specialUpdate();

	}

	// before updating the matrix instant rates and stationary probabilities
	// we need to check that we are pointing on the right nucleotide mutation process
	// and we need to update the value of omega stored by the MGOmegaCodonSubMatrix object,
	// based on the value stored by the RandomOmega parent
	void SetParameters()	{
		SetNucMatrix(matrix);
		SetSynRateTs(RandomSynRateTs->val());
		SetOmegaTs(RandomOmegaTs->val());
		SetOmegaTv0(RandomOmegaTv0->val());
	}

	protected:

	RandomSubMatrix* matrix;
	Var<PosReal>* RandomSynRateTs;
	Var<PosReal>* RandomOmegaTs;
	Var<PosReal>* RandomOmegaTv0;
};


void MGOmega2CodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)	{
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = (*NucMatrix)(a,b);
				if (! Synonymous(i,j))	{
					if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
						Q[i][j] *= GetOmegaTs();
					}
					else	{
						Q[i][j] *= GetOmegaTv();
					}
				}
			}
			else	{
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

void MGOmega3CodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)	{
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = (*NucMatrix)(a,b);
				if (! Synonymous(i,j))	{
					if (((a==0) && (b==2)) || ((a==2) && (b==0)) || ((a==1) && (b==3)) || ((a==3) && (b==1)))	{
						Q[i][j] *= GetOmegaTs();
					}
					else	{
						if (((a==1) && (b==2)) || ((a==2) && (b==1)) || ((a==0) && (b==3)) || ((a==3) && (b==0)))	{
							Q[i][j] *= GetOmegaTv0();
						}
						else	{
							Q[i][j] *= GetOmegaTvGC();
						}
					}
				}
			}
			else	{
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

void MGOmega3X3CodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)	{
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = (*NucMatrix)(a,b);
				if (((a==0) && (b==2)) || ((a==2) && (b==0)) || ((a==1) && (b==3)) || ((a==3) && (b==1)))	{
					Q[i][j] *= GetSynRateTs();
				}
				else	{
					if (((a==1) && (b==2)) || ((a==2) && (b==1)) || ((a==0) && (b==3)) || ((a==3) && (b==0)))	{
						// Q[i][j] *= GetSynRateTv0();
					}
					else	{
						Q[i][j] *= GetSynRateTvGC();
					}
				}
				if (! Synonymous(i,j))	{
					if (((a==0) && (b==2)) || ((a==2) && (b==0)) || ((a==1) && (b==3)) || ((a==3) && (b==1)))	{
						Q[i][j] *= GetOmegaTs();
					}
					else	{
						if (((a==1) && (b==2)) || ((a==2) && (b==1)) || ((a==0) && (b==3)) || ((a==3) && (b==0)))	{
							Q[i][j] *= GetOmegaTv0();
						}
						else	{
							Q[i][j] *= GetOmegaTvGC();
						}
					}
				}
			}
			else	{
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

void MGOmega2X2CodonSubMatrix::ComputeArray(int i)	{

	double total = 0;
	for (int j=0; j<GetNstate(); j++)	{
		if (i!=j)	{
			int pos = GetDifferingPosition(i,j);
			if ((pos != -1) && (pos != 3))	{
				int a = GetCodonPosition(pos,i);
				int b = GetCodonPosition(pos,j);
				if (a == b)	{
					cerr << GetCodonStateSpace()->GetState(i) << '\t' << GetCodonStateSpace()->GetState(j) << '\n';
					cerr << pos << '\n';
					exit(1);
				}
				Q[i][j] = (*NucMatrix)(a,b);
				if (((a==0) && (b==2)) || ((a==2) && (b==0)) || ((a==1) && (b==3)) || ((a==3) && (b==1)))	{
					Q[i][j] *= GetSynRateTs();
				}
				if (! Synonymous(i,j))	{
					if (((a==0) && (b==2)) || ((a==2) && (b==0)) || ((a==1) && (b==3)) || ((a==3) && (b==1)))	{
						Q[i][j] *= GetOmegaTs();
					}
					else	{
						Q[i][j] *= GetOmegaTv0();
					}
				}
			}
			else	{
				Q[i][j] = 0;
			}
			total += Q[i][j];
		}
	}
	Q[i][i] = -total;
	if (total <0)	{
		cerr << "negative rate away\n";
		exit(1);
	}
}

#endif

