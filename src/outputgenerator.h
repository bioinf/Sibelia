//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _OUTPUT_GENERATOR_H_
#define _OUTPUT_GENERATOR_H_

#include "resource.h"
#include "blockfinder.h"
#include "variantcalling/variant.h"

namespace SyntenyFinder
{	
	class OutputGenerator
	{
	public:
		typedef std::vector<BlockInstance> BlockList;
		OutputGenerator(const ChrList & chrList, const BlockList & blockList): chrList_(chrList), blockList_(blockList) {}
		void GenerateReport(const std::string & fileName) const;		
		void GenerateCircosOutput(const std::string & outFile, const std::string & outDir) const;
		void GenerateHierarchyCircosOutput(const std::string & outFile, const std::string & outDir, const std::vector<BlockList> & history) const;
		void GenerateD3Output(const std::string & outFile) const;
		void ListBlocksIndices(const std::string & fileName) const;
		void ListBlocksSequences(const std::string & fileName) const;		
		void ListChromosomesAsPermutations(const std::string & fileName) const;
		void OutputBuffer(const std::string & fileName, const std::string & buffer) const;		
	private:
		DISALLOW_COPY_AND_ASSIGN(OutputGenerator);
		const ChrList & chrList_;
		mutable BlockList blockList_;
		void ListChrs(std::ostream & out) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;
		void TryOpenResourceFile(const std::string & fileName, std::ifstream & stream) const;
	};
}

#endif
