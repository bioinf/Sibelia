//****************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "fasta.h"
#include "blockinstance.h"

#ifndef _DE_BRUIJN_GRAPH_H_
#define _DE_BRUIJN_GRAPH_H_

namespace SyntenyFinder
{	
	class DeBruijnIndex
	{	
	private:
		class BifurcationData
		{
		public:
			BifurcationData();						
			BifurcationData(size_t pos, size_t bifId, size_t projection, char inMark, char outMark);
			bool IsValid() const;
			char GetInMark() const;
			char GetOutMark() const;
			size_t GetPosition() const;
			size_t GetProjection() const;
			size_t GetBifurcationId() const;
			void Invalidate();
		private:
			uint32_t pos_;
			uint32_t bifId_;
			uint32_t projection_;
			char inMark_;
			char outMark_;
			bool valid_;			
		};

		typedef std::vector<BifurcationData> BifurcationVector;
	public:
		static const char SEPARATION_CHAR;
		static const size_t MAX_POSITION;
		static const size_t MAX_BIFURCATION_ID;		
		static const size_t MAX_SEQUENCE_NUMBER;		
		
		class BifurcationIterator
		{
		public:
			BifurcationIterator();
			bool AtEnd() const;
			bool IsValid() const;			
			char GetOutMark() const;		
			size_t GetPosition() const;
			size_t GetProjection() const;
			size_t GetChromosomeId() const;
			size_t GetBifurcationId() const;			
			size_t GetPositivePosition() const;
			size_t GetPositiveEndingPosition() const;
			BifurcationIterator& operator++();
			BifurcationIterator operator++(int);
			BifurcationIterator operator + (size_t shift) const;
			bool operator == (const BifurcationIterator & it) const;
			bool operator != (const BifurcationIterator & it) const;
		private:
			BifurcationIterator(const DeBruijnIndex * parent, size_t chrId, size_t index, FastaRecord::Direction dir);
			size_t chrId_;
			size_t index_;
			FastaRecord::Direction dir_;
			const DeBruijnIndex * parent_;
			friend class DeBruijnIndex;
		};

		struct BifurcationInstance
		{
		public:
			BifurcationInstance() {}
			BifurcationInstance(size_t bifId, size_t pos, size_t proj): bifId_(bifId), pos_(pos), proj_(proj) {}
			bool operator < (const BifurcationInstance & toCompare) const
			{
				return pos_ < toCompare.pos_;
			}

			size_t GetId() const
			{
				return bifId_;
			}

			size_t GetProjection() const
			{
				return proj_;
			}

			size_t GetPostion() const
			{
				return pos_;
			}

		private:
			size_t bifId_;	
			size_t pos_;
			size_t proj_;
		};

		typedef std::vector<BifurcationInstance> BifVector;
		typedef std::vector<BifVector> ChrBifVector;
		
		void ApplyChanges();		
		size_t GetBifurcationsNumber() const;
		size_t CountInstances(size_t bifId) const;
		size_t GetBifurcationInstances(size_t bifId, std::vector<BifurcationIterator> & ret) const;
		DeBruijnIndex(const std::vector<ChrBifVector> & bifurcation, const std::vector<std::string> & record, size_t k, size_t bifurcationNumber);
		void Replace(BifurcationIterator sourceStart, BifurcationIterator sourceEnd, BifurcationIterator targetStart, BifurcationIterator targetEnd);
		BifurcationIterator Begin(size_t chr, FastaRecord::Direction dir);
		BifurcationIterator End(size_t chr, FastaRecord::Direction dir);
	private:
		DISALLOW_COPY_AND_ASSIGN(DeBruijnIndex);
		
		class Location
		{
		public:
			Location() {}
			Location(size_t chrId, size_t index);
			size_t GetIndex() const;
			size_t GetChromosomeId() const;
		private:
			uint32_t chrId_;
			uint32_t index_;
		};
		
		static size_t GetStrand(FastaRecord::Direction dir);
		typedef std::vector<Location> LocationVector;		

		size_t k_;
		const std::vector<FastaRecord> * chr_;
		std::vector<size_t> revCompDictionary_;		
		std::vector<LocationVector> bifurcationPlace_;		
		std::vector<BifurcationVector> bifurcationData_;
		friend class BifurcationIterator;		
	};
}

#endif
