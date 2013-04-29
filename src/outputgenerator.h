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
		OutputGenerator(const ChrList & chrList): chrList_(chrList) {}
		void GenerateReport(const BlockList & blockList, const std::string & fileName) const;		
		void GenerateCircosOutput(const BlockList & blockList, const std::string & outFile, const std::string & outDir) const;
		void GenerateHierarchyCircosOutput(const std::vector<BlockList> & history, const std::string & outFile, const std::string & outDir) const;
		void GenerateD3Output(const BlockList & blockList, const std::string & outFile) const;
		void GenerateVariantOutput(const std::vector <Variant> & variants, const std::string & assemblyFile, const std::string & outFile) const;
		void ListBlocksIndices(const BlockList & blockList, const std::string & fileName) const;
		void ListBlocksIndicesHeirarchy(const std::vector<BlockList> & history, const std::string & fileName) const;
		void OutputTree(const std::vector<BlockList> & history, const std::string & fileName) const;
		void ListBlocksSequences(const BlockList & blockList, const std::string & fileName) const;		
		void ListChromosomesAsPermutations(const BlockList & blockList, const std::string & fileName) const;
		void RearrangementScenario(const std::vector<std::string> & steps, const std::string & fileName) const;
		void OutputBuffer(const std::string & fileName, const std::string & buffer) const;		
        void OutputBlocksInSAM(const BlockList & block, const std::string & fileName) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(OutputGenerator);
		static const int CIRCOS_MAX_COLOR;
		static const int CIRCOS_DEFAULT_RADIUS;
		static const int CIRCOS_RESERVED_FOR_LABEL;
		static const int CIRCOS_HIGHLIGHT_THICKNESS;
		const ChrList & chrList_;		
		void ListChrs(std::ostream & out) const;
		void TryOpenFile(const std::string & fileName, std::ofstream & stream) const;
		void TryOpenResourceFile(const std::string & fileName, std::ifstream & stream) const;		
		void WriteCircosKaryoType(const std::string & outDir, const std::string & fileName) const;
		void WriteCircosImageConfig(const std::string & outDir, const std::string & fileName, int r) const;
		void WriteCircosLinks(const std::string & outDir, const std::string & fileName, const BlockList & block) const;
		void WriteCircosHighlight(const std::string & outDir, const std::string & fileName, const BlockList & block, int r0, int r1, bool ideogram, std::ofstream & config) const;		
	};
}

#endif
