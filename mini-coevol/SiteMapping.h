
#ifndef SITEMAPPING_H
#define SITEMAPPING_H

#include <iostream>
#include "BranchSitePath.h"

class SiteMapping	{

	public:

	virtual			~SiteMapping() {}

	virtual BranchSitePath*	GetPath(const Branch* branch) = 0;
	virtual Tree*		GetTree() = 0;
	Link*			GetRoot();

	virtual void		Print(ostream& os, bool redundant);
	void			Print(ostream& os, Link* from, bool redundant);

};

inline Link* SiteMapping::GetRoot() {return GetTree()->GetRoot();}

#endif // SITEMAPPING_H
