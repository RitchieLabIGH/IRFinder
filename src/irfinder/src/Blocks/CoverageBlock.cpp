#include "CoverageBlock.h"
// using namespace std;

CoverageBlock::CoverageBlock(uint start, uint end) {
	blockStart = start;
	blockEnd = end;
	firstDepth[0] = 0;
	firstDepth[1] = 0;
	blockExtents = NULL;
	blockExtentsL = NULL;
}

//direction -- 0=False/Neg, 1=True/Pos.
void CoverageBlock::RecordCover(uint readStart, uint readEnd, bool dir) {

	if (readStart <= blockStart && readEnd > blockStart) {
		firstDepth[dir]++;
	} else if (readStart < blockEnd) {
		// Need to increment the starts vector.
		uint inc_index = readStart - blockStart - 1;
		if (blockExtentsL) { //already an int vector
			blockExtentsL->at(inc_index).start[dir]++;
		} else if (!blockExtents) { //don't have a char vector either - create first.
			blockExtents = new std::vector<start_stops>(vectorLen());
			blockExtents->at(inc_index).start[dir]++;
		} else {
			if (blockExtents->at(inc_index).start[dir] == 254) {
				blockExtentsL = new std::vector<start_stopsL>(
						blockExtents->begin(), blockExtents->end());
				delete blockExtents;
				blockExtents = NULL;
				blockExtentsL->at(inc_index).start[dir]++;
			} else {
				blockExtents->at(inc_index).start[dir]++;
			}
		}
	} else {
		return;
	}

	if (readEnd >= blockEnd) {
		return;
	} else {
		// Need to increment the ends vector.
		uint inc_index = readEnd - blockStart - 1;

		if (blockExtentsL) { //already an int vector
			blockExtentsL->at(inc_index).end[dir]++;
		} else if (!blockExtents) { //don't have a char vector either - create first.
			blockExtents = new std::vector<start_stops>(vectorLen());
			blockExtents->at(inc_index).end[dir]++;
		} else {
			if (blockExtents->at(inc_index).end[dir] == 254) {
				blockExtentsL = new std::vector<start_stopsL>(
						blockExtents->begin(), blockExtents->end());
				delete blockExtents;
				blockExtents = NULL;
				blockExtentsL->at(inc_index).end[dir]++;
			} else {
				blockExtents->at(inc_index).end[dir]++;
			}
		}
	}
	// Can Throw: Out of range exception.
}

void CoverageBlock::updateCoverageHist(std::map<uint, uint> &hist, uint start,
		uint end) const {

	if (!blockExtentsL && !blockExtents) {
		// how many bases in this block?
		hist[firstDepth[0] + firstDepth[1]] += std::min(blockEnd, end)
				- std::max(blockStart, start);
	} else {
		// There are read starts and ends -- need to walk the positions from the start of this block
		//  even if not in the region of interest.

		//special handling for the first base -- the one before the vector starts.
		uint depth = firstDepth[0] + firstDepth[1];
		if (start <= blockStart) {
			// use the first depth, before commencing in the vector.
			hist[depth]++;
		}

		uint startindex = std::max(blockStart + 1, start) - blockStart - 1;
		uint endindex = std::min(blockEnd, end) - blockStart - 1;
		if (blockExtents) {
			for (uint i = 0; i < endindex; i++) {
				depth += -(*blockExtents)[i].end[0] - (*blockExtents)[i].end[1]
						+ (*blockExtents)[i].start[0]
						+ (*blockExtents)[i].start[1];
				if (i >= startindex) {
					hist[depth]++;
				}
			}
		} else {
			for (uint i = 0; i < endindex; i++) {
				depth += -(*blockExtentsL)[i].end[0]
						- (*blockExtentsL)[i].end[1]
						+ (*blockExtentsL)[i].start[0]
						+ (*blockExtentsL)[i].start[1];
				if (i >= startindex) {
					hist[depth]++;
				}
			}
		}
		//  When in the region of interest, update the hist each step.
	}
}

void CoverageBlock::updateCoverageHist(std::map<uint, uint> &hist, uint start,
		uint end, bool dir) const {
	if (!blockExtentsL && !blockExtents) {
		// how many bases in this block?
		hist[firstDepth[dir]] += std::min(blockEnd, end)
				- std::max(blockStart, start);
	} else {
		//special handling for the first base -- the one before the vector starts.
		uint depth = firstDepth[dir];
		if (start <= blockStart) {
			// use the first depth, before commencing in the vector.
			hist[depth]++;
		}

		uint startindex = std::max(blockStart + 1, start) - blockStart - 1;
		uint endindex = std::min(blockEnd, end) - blockStart - 1;
		if (blockExtents) {
			for (uint i = 0; i < endindex; i++) {
				depth += -(*blockExtents)[i].end[dir]
						+ (*blockExtents)[i].start[dir];
				if (i >= startindex) {
					hist[depth]++;
				}
			}
		} else {
			for (uint i = 0; i < endindex; i++) {
				depth += -(*blockExtentsL)[i].end[dir]
						+ (*blockExtentsL)[i].start[dir];
				if (i >= startindex) {
					hist[depth]++;
				}
			}
		}
	}
}

void CoverageBlock::updateCoverageArray(std::vector<uint> &arr,
		std::vector<bool> &covered, uint start, uint end) const {
	uint depth = firstDepth[0] + firstDepth[1],
			startindex = std::max( blockStart, start-1) - blockStart,
			endindex = std::min(blockEnd, end) - blockStart,
			startarray = std::max(blockStart+1, start ) - start ,
			endarray = std::min(blockEnd, end) - start ;

	if (!blockExtentsL && !blockExtents) {
		for (uint i = startindex; i < endindex && startarray < endarray;
				i++, startarray++) {
			arr[startarray] += depth;
			covered[startarray] = true;
		}
	} else {
		// There are read starts and ends -- need to walk the positions from the start of this block
		//  even if not in the region of interest.
		if (blockExtents) {
			for (uint i = 0; i < endindex && startarray < endarray; i++) {
				depth += -(*blockExtents)[i].end[0] - (*blockExtents)[i].end[1]
						+ (*blockExtents)[i].start[0]
						+ (*blockExtents)[i].start[1];
				if (i >= startindex) {
					arr[startarray] += depth;
					covered[startarray] = true;
					startarray++;
				}
			}
		} else {
			for (uint i = 0; i < endindex && startarray < endarray; i++) {
				depth += -(*blockExtentsL)[i].end[0]
						- (*blockExtentsL)[i].end[1]
						+ (*blockExtentsL)[i].start[0]
						+ (*blockExtentsL)[i].start[1];
				if (i >= startindex) {
					arr[startarray] += depth;
					covered[startarray] = true;
					startarray++;
				}
			}
		}
		//  When in the region of interest, update the hist each step.
	}
}

void CoverageBlock::updateCoverageArray(std::vector<uint> &arr,
		std::vector<bool> &covered, uint start, uint end, bool dir) const {

	uint depth = firstDepth[0] + firstDepth[1],
				startindex = std::max( blockStart, start-1) - blockStart,
				endindex = std::min(blockEnd, end) - blockStart,
				startarray = std::max(blockStart+1, start ) - start ,
				endarray = std::min(blockEnd, end) - start ;
	if (!blockExtentsL && !blockExtents) {
		for (uint i = startindex; i < endindex && startarray < endarray;
				i++, startarray++) {
			arr[startarray] += depth;
			covered[startarray] = true;
		}
	} else {
		// There are read starts and ends -- need to walk the positions from the start of this block
		//  even if not in the region of interest.
		if (blockExtents) {
			for (uint i = 0; i < endindex && startarray < endarray; i++) {
				depth += -(*blockExtents)[i].end[dir]
						+ (*blockExtents)[i].start[dir];
				if (i >= startindex) {
					arr[startarray] += depth;
					covered[startarray] = true;
					startarray++;
				}
			}
		} else {
			for (uint i = 0; i < endindex && startarray < endarray; i++) {
				depth += -(*blockExtentsL)[i].end[dir]
						+ (*blockExtentsL)[i].start[dir];
				if (i >= startindex) {
					arr[startarray] += depth;
					covered[startarray] = true;
					startarray++;
				}
			}
		}
		//  When in the region of interest, update the hist each step.
	}
}

void CoverageBlock::print(std::ostream &os) const {
	os << "Coverage block " << blockStart << " - " << blockEnd << "\n";
	os << "First depth 0 : " << firstDepth[0] << "\n";
	os << "First depth 1 : " << firstDepth[0] << "\n";
	uint i=0;
	if (blockExtents) {
		os << "BlockExtents: \n";
		for (auto &a : (*blockExtents)) {
			os << i+blockStart << " " << (uint) a.start[0] << ":" << (uint) a.start[1] << " - " << (uint) a.end[0] <<  ":"<< (uint) a.end[1] << "\n";
			i++;
		}
	}
	if (blockExtentsL) {
		os << "BlockExtentsL: \n";
		for (auto &a : (*blockExtentsL)) {
			os << i+blockStart << " " << (uint) a.start[0] << ":" << (uint) a.start[1] << " - " << (uint) a.end[0] <<  ":"<< (uint) a.end[1] << "\n";
			i++;
		}
	}
}

std::ostream& operator<<(std::ostream &os, const CoverageBlock &cb) {
	cb.print(os);
	return os;
}
