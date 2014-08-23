//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _INDEXED_SEQUENCE_H_
#define _INDEXED_SEQUENCE_H_

#include "hashing.h"
#include "platform.h"

namespace SyntenyFinder
{
	class IndexedSequence
	{
	public:
		static const char SEPARATION_CHAR;
/*		void Test();
		DNASequence& Sequence();
		const DNASequence& Sequence() const;
		BifurcationStorage& BifStorage();
		const BifurcationStorage& BifStorage() const;
		void ConstructChrIndex();
		size_t GetChr(StrandIterator it) const;
		IndexedSequence(const std::vector<std::string> & record, size_t k, const std::string & tempDir);
		IndexedSequence(const std::vector<std::string> & record, std::vector<std::vector<Pos> > & original, size_t k, const std::string & tempDir, bool clear = false);
		static bool StrandIteratorPosGEqual(StrandIterator a, StrandIterator b);		
		static size_t StrandIteratorDistance(StrandIterator start, StrandIterator end);*/

		struct BifurcationInstance
		{
		public:
			
			BifurcationInstance() {}
			BifurcationInstance(size_t bifId, size_t pos): bifId_(bifId), pos_(pos) {}
			bool operator < (const BifurcationInstance & toCompare) const
			{
				return pos_ < toCompare.pos_;
			}

			size_t GetBifurcationId() const
			{
				return bifId_;
			}

			size_t GetBifurcationPostion() const
			{
				return pos_;
			}

		private:
			size_t bifId_;	
			size_t pos_;
		};

		typedef std::vector<BifurcationInstance> BifVector;
		typedef std::vector<BifVector> ChrBifVector;

		static size_t EnumerateBifurcationsSArray(const std::vector<std::string> & data, size_t k_, const std::string & tempDir, std::vector<ChrBifVector> & ret);

	private:
		DISALLOW_COPY_AND_ASSIGN(IndexedSequence);
	/*	typedef std::pair<StrandIterator, size_t> IteratorChrPair;
		size_t k_;
		std::auto_ptr<DNASequence> sequence_;
		std::auto_ptr<BifurcationStorage> bifStorage_;
		std::vector<IteratorChrPair> chrIndex_;		

		struct IteratorHash
		{
			size_t operator() (const StrandIterator & it) const
			{
				return it.GetElementId();
			}
		};

		typedef std::vector<Pos> PosVector;		
		typedef boost::unordered_map<std::string, size_t> KMerBifMap;
		typedef boost::unordered_map<StrandIterator, size_t, IteratorHash> IteratorIndexMap;

	#ifdef _DEBUG
		KMerBifMap idMap_;
	#endif		

		size_t GetMustBeBifurcation(StrandIterator it);
		void Init(std::vector<std::string> record, std::vector<std::vector<Pos> > & original, size_t k, const std::string & tempDir, bool clear);		*/
	};
}

#endif