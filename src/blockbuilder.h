//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "platform.h"
#include "debruijnindex.h"

#ifndef _BLOCK_BUILDER_H_
#define _BLOCK_BUILDER_H_

namespace SyntenyFinder
{	
	class BlockBuilder
	{
	public:
		enum State
		{
			start,
			run,
			end
		};

		typedef boost::function<void(size_t, State)> ProgressCallBack;
		
		static const size_t PROGRESS_STRIDE;

		BlockBuilder(const std::vector<FastaRecord> * originalChr, const std::string & tempDir);
		void ConstructIndex(size_t k);
		void WriteIndexToDot(std::ostream & out) const;		
		size_t Simplify(size_t maxBranchSize, size_t maxIterations, ProgressCallBack callBack = ProgressCallBack());		
		void GenerateSyntenyBlocks(size_t k, size_t trimK, size_t minSize, std::vector<BlockInstance> & block) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(BlockBuilder);
		size_t lastK_;
		std::auto_ptr<DeBruijnIndex> index_;
		std::vector<size_t> virtualChrSize_;
		const std::vector<FastaRecord> * originalChr_;
		std::string tempDir_;

		struct VisitData
		{
			size_t kmerId;
			size_t distance;
			VisitData() {}
			VisitData(size_t kmerId, size_t distance): kmerId(kmerId), distance(distance) {}
		};

		struct BranchData
		{
			BranchData() {}
			BranchData(char ch, size_t maxBifMlp): endChar(ch), maxBifMlp(maxBifMlp) {}
			char endChar;
			size_t maxBifMlp;
			std::vector<VisitData> branch;
		};

		struct PotentialBlockInstance
		{
			size_t chr;
			size_t start;
			size_t end;
			FastaRecord::Direction dir;
			PotentialBlockInstance(size_t chr, size_t start, size_t end, FastaRecord::Direction dir): chr(chr), start(start), end(end), dir(dir)
			{
			}
		};

		typedef std::vector<PotentialBlockInstance> PotentialInstanceVector;

		typedef std::vector< std::vector<VisitData> > BulgedBranches;

		static bool BlockBuilder::CmpPotVector(const PotentialInstanceVector & a, const PotentialInstanceVector & b);
		static bool BlockBuilder::CmpBifIt(const DeBruijnIndex::BifurcationIterator & a, const DeBruijnIndex::BifurcationIterator & b);
		size_t RemoveBulges(size_t minBranchSize, size_t bifId);		
		bool AnyBulges(std::vector<DeBruijnIndex::BifurcationIterator> bif, BulgedBranches & bulges, size_t maxBranchSize) const;
		bool Overlap(const std::vector<DeBruijnIndex::BifurcationIterator> & bif, VisitData sourceData, VisitData targetData) const;
		void CollapseBulge(const std::vector<DeBruijnIndex::BifurcationIterator> & edge, VisitData sourceData, VisitData targetData);	
		void BlockBuilder::ResolveOverlap(PotentialInstanceVector & v, size_t minSize, std::vector<std::vector<bool> > & overlap) const;
		size_t EnumerateBifurcationsSArray(const std::vector<std::string> & data, size_t k_, const std::string & tempDir, std::vector<DeBruijnIndex::ChrBifVector> & ret) const;
	};
}

#endif