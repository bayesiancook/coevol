#ifndef SITEMAPPING_H
#define SITEMAPPING_H

#include <iostream>
#include "diffsel/ValTree.hpp"
// #include "diffsel/BranchSitePath.hpp"

class BranchSitePath;
class Branch;
class Tree;
class Link;

class SiteMapping {
  public:
    virtual ~SiteMapping() = default;

    virtual BranchSitePath* GetPath(const Branch* branch) = 0;
    virtual Tree* GetTree() = 0;
    Link* GetRoot();

    virtual void Print(std::ostream& os, bool redundant);
    void Print(std::ostream& os, Link* from, bool redundant);
};

inline Link* SiteMapping::GetRoot() { return GetTree()->GetRoot(); }

#endif  // SITEMAPPING_H
