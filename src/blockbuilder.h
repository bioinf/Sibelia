//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
		void GenerateBlocks(std::vector<BlockInstance> & ret, size_t minBlockSize) const;
		size_t Simplify(size_t minBranchSize, size_t maxIterations, ProgressCallBack callBack = ProgressCallBack());
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

		size_t RemoveBulges(size_t minBranchSize, size_t bifId);
		bool Overlap(const std::vector<DeBruijnIndex::Edge> & edge, VisitData sourceData, VisitData targetData) const;
		void CollapseBulgeGreedily(std::vector<DeBruijnIndex::Edge> & edge, VisitData sourceData, VisitData targetData);

	#ifdef _DEBUG		
		std::vector<std::string> virtualChr_;
		std::map<std::string, size_t> debugIndex_;
		void PrintRaw(size_t chr0, size_t chr1, std::ostream & out) const;
		void PrintPath(DeBruijnIndex::Edge e, size_t distance, std::ostream & out) const;
	#endif
	};
}

#endif