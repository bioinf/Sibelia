#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <set>
#include <sstream>
#include "simplerearrangements.h"

namespace SyntenyFinder
{
    Genome::Genome(int n_, const std::vector<int> &permutation) : n(n_)
    {
      //synteny = std::vector<std::vector<int>(2) >(n);
      for (size_t i = 0; i < n + 1; i++) synteny.push_back(std::vector<int>(2));
      int prev = -1;
      for (size_t i = 0; i < n; ++i) {
          vertices.push_back(std::vector<int>(2));                            //
          vertices[i][0] = 2*abs(permutation[i]) - 1 * (permutation[i] > 0);  //add one synteny block
          vertices[i][1] = prev;                                              //
          prev = 2*abs(permutation[i]) - 1 * (permutation[i] < 0);            //

                                                                              //
          synteny[abs(permutation[i])][1 * (permutation[i] < 0)] = i;              //
          synteny[abs(permutation[i])][1 * (permutation[i] > 0)] = i + 1;          //
      }                                                                       //
      vertices.push_back(std::vector<int>(2));                                //
      vertices[n][0] = -1;                                                    //
      vertices[n][1] = prev;                                                  //
    }

    std::string Genome::ShowArrangement()
    {
        int vAddr, next;
        std::stringstream ss;
        std::vector<bool> visited = std::vector<bool>(vertices.size());
        for (size_t i = 0; i < vertices.size(); i++) //show linear
        {
            if (!visited[i] && (vertices[i][0] == -1 || vertices[i][1] == -1))
            {
                visited[i] = true;
                next = vertices[i][1 * vertices[i][0] == -1];
                while (next != -1)
                {
                    if (next%2 != 0)
                    {
                        ss << "+" << (next/2 + 1) << " ";
                        vAddr = synteny[next/2 + 1][1];
                        next = vertices[vAddr][1 * (vertices[vAddr][0] == (next + 1))];
                    } else {
                        ss << -next/2  << " ";
                        vAddr = synteny[next/2][0];
                        next = vertices[vAddr][1 * (vertices[vAddr][0] == (next - 1))];
                    }
                    visited[vAddr] = true;
                }
                ss<< "$ ";
            }
        }
        for (size_t i = 0; i < vertices.size(); i++) //show circular
        {
            if (!visited[i])
            {
                visited[i] = true;
                next = vertices[i][0];//luboe voobshe?
                while (vAddr != i)
                {
                    if (next%2 != 0)
                    {
                        ss << "+" << (next/2 + 1) << " ";
                        vAddr = synteny[next/2 + 1][1];
                        next = vertices[vAddr][1 * (vertices[vAddr][0] == (next + 1))];
                    } else {
                        ss << -next/2 << " ";
                        vAddr = synteny[next/2][0];
                        next = vertices[vAddr][1 * (vertices[vAddr][0] == (next - 1))];
                    }
                    visited[vAddr] = true;
                }
                ss << (next/2 + 1 * (next%2 != 0)) << " ";
            }
        }
        return ss.str();
    }

std::vector<std::string> GetRearrangements(const std::vector<int> &pa, const std::vector<int> &pb)
{
    std::vector<std::string> result;
    assert(pa.size() == pb.size());
    int n = pa.size();
    Genome a = Genome(n, pa);
    Genome b = Genome(n, pb);
    int p, q, s, t, uAddr, vAddr;
    std::vector<int> u(2), v(2);
    for (size_t i = 0; i < n + 1; ++i)
    {
        if (b.vertices[i][0] != -1 && b.vertices[i][1] != -1) // is adjacency
        {
            p = b.vertices[i][0];
            q = b.vertices[i][1];
            uAddr = a.synteny[p/2 + 1 * (p%2 != 0)][1 * (p%2 == 0)];
            vAddr = a.synteny[q/2 + 1 * (q%2 != 0)][1 * (q%2 == 0)];
            u[0] = a.vertices[uAddr][0]; // contains edge p
            u[1] = a.vertices[uAddr][1];
            v[0] = a.vertices[vAddr][0]; // contains edge q
            v[1] = a.vertices[vAddr][1];
            if ((u[0] != v[0] && u[0] != v[1]) || (u[1] != v[1] && u[1] != v[1])) // u != v
            {
                s = u[1 * (u[0] == p)]; // not p edge from u
                t = v[1 * (v[0] == q)]; // not q edge from v
                a.vertices[uAddr][0] = p;
                a.vertices[uAddr][1] = q;
                a.synteny[q/2 +  1 * (q%2 != 0)][1 * (q%2 == 0)] = uAddr; // new vertex for edge q

                a.vertices[vAddr][0] = s;
                a.vertices[vAddr][1] = t;
                if (s != -1) a.synteny[s/2 + 1 * (s%2 != 0)][1 * (s%2 == 0)] = vAddr; // new vertex for edge s, if it exist
                result.push_back(a.ShowArrangement());
            }
        }
    }

    for (size_t i = 0; i < n + 1; ++i)
    {
        p = -1;
        if (b.vertices[i][0] == -1) p = b.vertices[i][1];
        if (b.vertices[i][1] == -1) p = b.vertices[i][0];
        if (p != -1) // is telomere
        {
            uAddr = a.synteny[p/2 + 1 * (p%2 != 0)][1 * (p%2 == 0)];
            u[0] = a.vertices[uAddr][0]; //contains edge p
            u[1] = a.vertices[uAddr][1];
            if (u[0] != -1 && u[1] != -1) //u is adjacency
            {
                s = u[1 * (u[0] == p)]; // not p edge from u
                a.vertices[uAddr][1 * (u[0] == p)] = -1;

                a.vertices.push_back(std::vector<int>(2)); // create new telomere
                a.vertices[a.vertices.size()][0] = s;
                a.synteny[s/2 + 1 * (s%2 != 0)][1 * (s%2 == 0)] = a.vertices.size(); // new vertex for edge s
                result.push_back(a.ShowArrangement());
            }
        }
    }
    return result;
}

    std::vector<std::string> GetRearrangements(const std::vector<BlockInstance> &blocks)
    {
        std::set<BlockInstance> ReferenceSet;
        std::set<BlockInstance> ContigsSet;

        for (size_t i; i < blocks.size(); ++i)
        {
            BlockInstance block = blocks[i];
            if (block.GetChrId() == 1) {
                ReferenceSet.insert(block);
            } else {
                ContigsSet.insert(block);
            }
        }
        std::vector<int> pa;
        std::vector<int> pb;

        for(std::set<BlockInstance>::iterator it = ReferenceSet.begin(); it != ReferenceSet.end(); it++)
        {
           BlockInstance block = *it;
           if (ContigsSet.count(block) != 0) pa.push_back(block.GetSignedBlockId());
        }
        for(std::set<BlockInstance>::iterator it = ContigsSet.begin(); it != ContigsSet.end(); it++)
        {
            BlockInstance block = *it;
            if (ReferenceSet.count(block) == 0) pb.push_back(block.GetSignedBlockId());
        }
        return GetRearrangements(pa, pb);
    }
}
