//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "blockfinder.h"

namespace SyntenyFinder
{
	BlockFinder::Edge::Edge(size_t chr, DNASequence::Direction direction, size_t startVertex, size_t endVertex, size_t actualPosition, size_t actualLength, size_t originalPosition, size_t originalLength, char firstChar):
		chr(chr), direction(direction), startVertex(startVertex), endVertex(endVertex), actualPosition(actualPosition), actualLength(actualLength), originalPosition(originalPosition), originalLength(originalLength), firstChar(firstChar) {}

	bool BlockFinder::CompareEdgesByDirection(const Edge & a, const Edge & b)
	{
		return a.direction < b.direction;
	}

	bool BlockFinder::Edge::Coincide(const Edge & edge) const
	{
		return startVertex == edge.startVertex && endVertex == edge.endVertex && firstChar == edge.firstChar;
	}

	bool BlockFinder::Edge::Overlap(const Edge & edge) const
	{
		size_t a1 = originalPosition;
		size_t b1 = originalPosition + originalLength;
		size_t a2 = edge.originalPosition;
		size_t b2 = edge.originalPosition + edge.originalLength;
		size_t overlap = 0;
		if(a1 >= a2 && a1 <= b2)
			overlap = std::min(b1, b2) - a1;
		if(a2 >= a1 && a2 <= b1)
			overlap = std::min(b1, b2) - a2;
		return edge.chr == chr && overlap > 0;
	}
}