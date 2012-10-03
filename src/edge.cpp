#include "blockfinder.h"

namespace SyntenyFinder
{
	BlockFinder::Edge::Edge(size_t chr, DNASequence::Direction direction, size_t startVertex, size_t endVertex, size_t actualPosition, size_t actualLength, size_t originalPosition, size_t originalLength, char firstChar):
		chr(chr), direction(direction), startVertex(startVertex), endVertex(endVertex), actualPosition(actualPosition), actualLength(actualLength), originalPosition(originalPosition), originalLength(originalLength), firstChar(firstChar) {}

	bool BlockFinder::Edge::Coincide(const Edge & edge) const
	{
		return startVertex == edge.startVertex && endVertex == edge.endVertex && firstChar == edge.firstChar;
	}
}