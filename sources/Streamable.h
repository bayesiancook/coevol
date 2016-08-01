
#ifndef STREAMABLE_H
#define STREAMABLE_H

class Streamable	{

	public:

	virtual ~Streamable() {}

	virtual void ToStream(ostream& os) const = 0;
	virtual void FromStream(istream& is) = 0;

	friend operator<<(ostream& os, const Streamable& s)	{
		s.ToStream(os);
		return os;
	}

	friend operator>>(istream& is, Streamable& s)	{
		s.FromStream(is);
		return is;
	}

};


#endif

