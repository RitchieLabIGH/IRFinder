#ifndef CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS
#define CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS

#include "../Blocks/CoverageBlock.h"
#include "ReadBlockProcessor.h"
#include "../Blocks/FragmentBlocks.h"

struct BEDrecord {
	std::string chrName;
	std::string name;
	uint start;
	uint end;
	bool direction;
	std::vector<std::pair<uint,uint>> blocks;
};


class CoverageBlocks : public ReadBlockProcessor {
	//Store the Blocked BED record for each ROI/intron. This won't be referred to again until the end.
	//XX Create the temporary vectors (per Chr) which simply list the blocks sequentially as read.
	//XX Sort the temporary vectors
	//XX Build the final vectors of "blocks of interest"
	//xx Delete the temporary vectors
	//xx Create the parallel vectors with counter objects. (do these as a batch at the end, once vector size is known - for best memory layout)
	//xx Process fragments against the counter structure. (have I already written a class/object for this?)
	
	//Produce summary statistical output for each Blocked BED record, using the counter structure.

	protected:

		// Coverage depth data-structures.
		std::map<std::string, std::vector<CoverageBlock>> chrName_CoverageBlocks;
		std::map<std::string, std::vector<CoverageBlock>> chrName_FlankCoverageBlocks;
		// Shortcut pointers to depth data-structures.
		std::vector<std::vector<CoverageBlock>*> chrID_CoverageBlocks;
		std::vector<std::vector<CoverageBlock>*> chrID_FlankCoverageBlocks;

		// TODO: what is optimal for speed & memory usage?
//		static const uint coverage_block_max_length = 5000;
		static const uint coverage_block_max_length = 500;


		std::vector<BEDrecord> BEDrecords;
		bool long_read=false;
		int jitter = 3;

	public:
		CoverageBlocks(std::string read_type) {
			long_read = read_type == "LR";
		}
		void setJitter(int j){jitter=j;};
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<std::string> &chrmap);
		void loadRef(std::istream &IN);
		int WriteOutput(std::ostream *os) const;
		
		void fillHist(std::map<uint,uint> &hist, const std::string &chrName, const std::vector<std::pair<uint,uint>> &blocks) const;
		void fillHist(std::map<uint,uint> &hist, const std::string &chrName, const std::vector<std::pair<uint,uint>> &blocks, bool direction) const;
		void getCoverageArray(std::vector<uint> &coverages,
				std::vector<bool> & covered,
				const std::string &chrName,
				const uint arr_start, const uint arr_end) const;
		void getCoverageArray(std::vector<uint> &coverages,
				std::vector<bool> & covered,
				const std::string &chrName,
				const uint arr_start, const uint arr_end,
				bool direction) const;
		double meanFromHist(const std::map<uint,uint> &hist) const;
		double coverageFromHist(const std::map<uint,uint> &hist) const;
		double percentileFromHist(const std::map<uint,uint> &hist, uint percentile) const;
		double trimmedMeanFromHist(const std::map<uint,uint> &hist, uint centerPercent) const;
};

class CoverageBlocksIRFinder : public CoverageBlocks {
	private:
		uint AI_warn=0;
		uint AI_intron=1;
		double AI_ratio=0.05;
	public:

	CoverageBlocksIRFinder(std::string read_type) : CoverageBlocks(read_type){
	}
	void setAI(uint AI_warning_level, uint AI_min_intron_coverage, double AI_IRratio){
		AI_warn=AI_warning_level;
		AI_intron=AI_min_intron_coverage;
		AI_ratio=AI_IRratio;
	}
		int WriteOutput(std::ostream *os, std::ostream *osAI, const JunctionCount &JC, const SpansPoint &SP, int directionality = 0) const;

};


#endif
