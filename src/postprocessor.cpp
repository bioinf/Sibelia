//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "postprocessor.h"

namespace SyntenyFinder
{
	namespace
	{
		struct Stripe
		{
			int firstBlock;
			int secondBlock;
			Stripe() {}
			Stripe(int firstBlock, int secondBlock): firstBlock(firstBlock), secondBlock(secondBlock) {}
			bool operator < (const Stripe & toCompare) const
			{
				return firstBlock < toCompare.firstBlock;
			}
		};		
	}

	void Postprocessor::GlueStripes(std::vector<BlockInstance> & block, const std::vector<FASTARecord> & chr)
	{
		std::vector<std::vector<BlockInstance> > perm(chr.size());
		for(size_t i = 0; i < block.size(); i++)
		{
			perm[block[i].GetChrId()].push_back(block[i]);
		}

		for(size_t i = 0; i < perm.size(); i++)
		{
			std::sort(perm[i].begin(), perm[i].end(), compareByStart);
		}

		int sentinel = INT_MAX >> 1;
		bool glue = false;
		do
		{
			std::vector<Stripe> stripe;
			for(size_t chr = 0; chr < perm.size(); chr++)
			{
				for(size_t i = 0; i < perm[chr].size(); i++)
				{
					int bid = perm[chr][i].GetSignedBlockId();
					if(bid > 0)
					{
						int nextBid = i < perm[chr].size() - 1 ? perm[chr][i + 1].GetSignedBlockId() : sentinel;
						stripe.push_back(Stripe(bid, nextBid));
					}
					else
					{
						int prevBid = i > 0 ? perm[chr][i - 1].GetSignedBlockId() : -sentinel;
						stripe.push_back(Stripe(-bid, -prevBid));
					}
				}
			}

			size_t now = 0;
			size_t next = 0;
			std::sort(stripe.begin(), stripe.end());
			for(; now < stripe.size(); now = next)
			{
				glue = true;
				for(; next < stripe.size() && stripe[next].firstBlock == stripe[now].firstBlock; next++)
				{
					if(stripe[next].secondBlock != stripe[now].secondBlock || stripe[next].secondBlock == sentinel || Abs(stripe[next].secondBlock) == stripe[next].firstBlock)
					{
						glue = false;
					}
				}

				if(glue)
				{
					typedef std::vector<Stripe>::iterator It;
					std::pair<It, It> range = std::equal_range(stripe.begin(), stripe.end(), Stripe(Abs(stripe[now].secondBlock), 0));
					if(range.second - range.first != next - now)
					{
						glue = false;
					}
					else
					{
						break;
					}
				}
			}

			if(glue)
			{
				assert(next - now > 1);
				int glueBid = stripe[now].firstBlock;
				for(size_t chr = 0; chr < perm.size(); chr++)
				{
					for(size_t i = 0; i < perm[chr].size(); i++)
					{
						int bid = perm[chr][i].GetBlockId();
						if(bid == glueBid)
						{
							bid = perm[chr][i].GetSignedBlockId();
							if(bid > 0)
							{
								BlockInstance & a = perm[chr][i];
								BlockInstance & b = perm[chr][i + 1];
								a = BlockInstance(a.GetSignedBlockId(), &a.GetChrInstance(), a.GetStart(), b.GetEnd());
								perm[chr].erase(perm[chr].begin() + i + 1);
							}
							else
							{
								BlockInstance & a = perm[chr][--i];
								BlockInstance & b = perm[chr][i + 1];
								a = BlockInstance(b.GetSignedBlockId(), &a.GetChrInstance(), a.GetStart(), b.GetEnd());
								perm[chr].erase(perm[chr].begin() + i + 1);
							}
						}
					}
				}
			}
		}
		while(glue);

		block.clear();
		std::vector<int> oldId;
		for(size_t chr = 0; chr < perm.size(); chr++)
		{
			for(size_t i = 0; i < perm[chr].size(); i++)
			{
				block.push_back(perm[chr][i]);
				oldId.push_back(perm[chr][i].GetBlockId());
			}
		}

		std::sort(oldId.begin(), oldId.end());
		oldId.erase(std::unique(oldId.begin(), oldId.end()), oldId.end());
		for(std::vector<BlockInstance>::iterator it = block.begin(); it != block.end(); ++it)
		{
			int sign = it->GetSignedBlockId() > 0 ? +1 : -1;
			size_t newId = std::lower_bound(oldId.begin(), oldId.end(), it->GetBlockId()) - oldId.begin() + 1;
			*it = BlockInstance(static_cast<int>(newId) * sign, &it->GetChrInstance(), it->GetStart(), it->GetEnd());
		}
	}
}