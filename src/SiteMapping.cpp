#include "SiteMapping.h"

//-------------------------------------------------------------------------
//	* SiteMapping
//-------------------------------------------------------------------------

void SiteMapping::Print(ostream& os, bool redundant)	{

	Print(os,GetRoot(), redundant);
	os << ";\n";
}


void SiteMapping::Print(ostream& os, Link* from, bool redundant)	{

	if (from->isLeaf())	{
		os << GetTree()->GetNodeName(from);
		if (redundant)	{
			os << '_';
		}
		else	{
			os << ':';
		}
	}
	else	{
		os << '(';
		for (Link* link=from->Next(); link!=from; link=link->Next())	{
			Print(os,link->Out(), redundant);
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	if (from->isRoot())	{
		os << GetPath(from->GetBranch())->GetCharInitState();
	}
	else	{
		os << GetPath(from->GetBranch())->ToString(redundant);
	}
}
