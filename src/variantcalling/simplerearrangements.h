#ifndef SIMPLEREARRANGEMENTS_H
#define SIMPLEREARRANGEMENTS_H

#include <tclap/CmdLine.h>
#include "../common.h"
#include "../util.h"
#include "variantcaller.h"

namespace SyntenyFinder
{
    struct Genome
    {
        int n;
        std::vector<std::vector<int> > vertices;
        std::vector<std::vector<int> > synteny;

        Genome(int n_, const std::vector<int>&);
        std::string ShowArrangement();
    };

    std::vector<std::string> GetRearrangements(const std::vector<int>&, const std::vector<int>&);
    std::vector<std::string> GetRearrangements(const std::vector<BlockInstance>& blocks);
}
#endif // SIMPLEREARRANGEMENTS_H
