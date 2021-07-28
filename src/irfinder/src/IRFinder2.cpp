#include "Utils/includedefine.h"
#include "ReadBlock/ReadBlockProcessor.h"
#include "ReadBlock/CoverageBlocks.h"
#include "Blocks/BAM2blocks.h"
#include <regex>

inline bool file_exists(const std::string &name) {
	return ( access( name.c_str(), F_OK ) != -1 );
}

void checkArguments(std::string outputDir, std::string s_inCoverageBlocks,
		std::string s_inSJ, std::string s_inSpansPoint, std::string s_inROI,
		std::string read_type, std::string input_BAM, std::string AI_level, std::string jitter) {

	if (read_type != "LR" && read_type != "SR") {
		std::cerr << "Error! read type must be LR or SR\n";
		exit(1);
	}
	if (file_exists(outputDir + "/IRFinder-IR-nondir.txt")) {
		std::cerr
				<< "Error! The output directory already contains an IRFinder's result.\n";
		exit(1);
	}
	for (std::string fname : { s_inCoverageBlocks, s_inSJ, s_inSpansPoint,
			input_BAM }) {
		if (!file_exists(fname)) {
			std::cerr << "Error! File " << fname << " not found!\n";
			exit(1);
		}
	}
	if (s_inROI != "NULL") {
		if (!file_exists(s_inROI)) {
			std::cerr << "Error! File " << s_inROI
					<< " not found!\nROI file is optional and, in case of absence, NULL must be given.\n";
			exit(1);
		}
	}
	if ( ! std::regex_match(AI_level, std::regex("^[0-9]+:[0-9]+:0.[0-9]+$") )){
		std::cerr << "Error! AI level has to be a string of two integers and a double separated by a semicolumn.\n";
		std::cerr << "Example: 1:2  -> output all the introns with at most the warning level of 1 and with at least 2 as IntronDepth.\n";
		exit(1);
	}
	if ( ! std::regex_match(jitter, std::regex("^[0-9]+$") )){
		std::cerr << "Error! Jitter must be an integer number.\n";
		exit(1);
	}
}


int main(int argc, char *argv[]) {

	if (argc != 10) {
		std::cerr << argc << " "
				<< "Usage: IRFinder output-directory ref-coverage.bed ref-SJ.ref ref-spans-point.ref ROI-named.bed(or NULL) [LR|SR] AI_level jitter input.bam"<< std::endl;
		for ( int i=0 ; i< argc; i++){
			std::cerr << " - " << argv[i] << "\n";
		}
		exit(1);
	}
	std::string outputDir = argv[1];
	std::string s_inCoverageBlocks = argv[2];
	std::string s_inSJ = argv[3];
	std::string s_inSpansPoint = argv[4];
	std::string s_inROI = argv[5];
	std::string read_type = argv[6];
	std::string AI_level = argv[7];
	std::string jitter = argv[8];
	std::string input_BAM = argv[9];
	checkArguments(outputDir, s_inCoverageBlocks, s_inSJ, s_inSpansPoint,
			s_inROI, read_type, input_BAM, AI_level, jitter);

	std::smatch ai_match;
	std::regex_match(AI_level, ai_match, std::regex("([0-9]+):([0-9]+):(0.[0-9]+)"));
	uint AI_warn = std::stoi(ai_match[1]), AI_intron=std::stoi(ai_match[2]);
	double AI_ratio=std::stod(ai_match[3]);
	std::cout << "IRFinder run with options:\n"
			<< " - Output Dir:          \t" << outputDir << "\n"
			<< " - Main intron ref.:    \t" << s_inCoverageBlocks << "\n"
			<< " - Splice junction ref.:\t" << s_inSJ << "\n"
			<< " - Read spans ref.:     \t" << s_inSpansPoint << "\n"
			<< " - Optional ROI ref.:   \t" << s_inROI << "\n"
			<< " - Read type:           \t" << read_type << "\n";
	if ( read_type == "LR") {
		std::cout << " - Jitter:              \t" << jitter << "\n";
	} else {
		std::cout << " - AI levels:           \t" << AI_warn << ":" << AI_intron << ":" << AI_ratio << "\n";
	}

	std::cout 	<< " - Input BAM:           \t" << input_BAM << "\n\n";
	std::cout.flush();

	std::cout << "Preparing the reference:\n";
	std::cout << " - Junction count...";std::cout.flush();
	JunctionCount oJuncCount;
	std::ifstream inJuncCount;
	inJuncCount.open(s_inSJ, std::ifstream::in);
	oJuncCount.loadRef(inJuncCount);
	inJuncCount.close();
	std::cout << "done.\n - Span points...";std::cout.flush();
	SpansPoint oSpansPoint;
	oSpansPoint.setSpanLength(5, 4);
	std::ifstream inSpansPoint;
	inSpansPoint.open(s_inSpansPoint, std::ifstream::in);
	oSpansPoint.loadRef(inSpansPoint);
	inSpansPoint.close();
	std::cout << "done.\n - Coverage blocks...";std::cout.flush();
	CoverageBlocksIRFinder oCoverageBlocks(read_type);
	std::ifstream inCoverageBlocks;
	inCoverageBlocks.open(s_inCoverageBlocks, std::ifstream::in);
	oCoverageBlocks.loadRef(inCoverageBlocks);
	inCoverageBlocks.close();
	std::cout << "done.\n";std::cout.flush();
	BAM2blocks BB;

	BB.registerCallbackChrMappingChange(
			std::bind(&JunctionCount::ChrMapUpdate, &oJuncCount,
					std::placeholders::_1));
	BB.registerCallbackProcessBlocks(
			std::bind(&JunctionCount::ProcessBlocks, &oJuncCount,
					std::placeholders::_1));


	FragmentsInChr oFragmentsInChr;
	BB.registerCallbackChrMappingChange(
			std::bind(&FragmentsInChr::ChrMapUpdate, &oFragmentsInChr,
					std::placeholders::_1));
	BB.registerCallbackProcessBlocks(
			std::bind(&FragmentsInChr::ProcessBlocks, &oFragmentsInChr,
					std::placeholders::_1));

	BB.registerCallbackChrMappingChange(
			std::bind(&SpansPoint::ChrMapUpdate, &oSpansPoint,
					std::placeholders::_1));
	BB.registerCallbackProcessBlocks(
			std::bind(&SpansPoint::ProcessBlocks, &oSpansPoint,
					std::placeholders::_1));

	FragmentsInROI oFragmentsInROI;
	if (s_inROI != "NULL") {
		std::cout << " - ROI...";std::cout.flush();
		std::ifstream inFragmentsInROI;
		inFragmentsInROI.open(s_inROI, std::ifstream::in);
		oFragmentsInROI.loadRef(inFragmentsInROI);
		inFragmentsInROI.close();

		BB.registerCallbackChrMappingChange(
				std::bind(&FragmentsInROI::ChrMapUpdate, &oFragmentsInROI,
						std::placeholders::_1));
		BB.registerCallbackProcessBlocks(
				std::bind(&FragmentsInROI::ProcessBlocks, &oFragmentsInROI,
						std::placeholders::_1));
		std::cout << "done\n";std::cout.flush();
	}

	BB.registerCallbackChrMappingChange(
			std::bind(&CoverageBlocks::ChrMapUpdate, &oCoverageBlocks,
					std::placeholders::_1));
	BB.registerCallbackProcessBlocks(
			std::bind(&CoverageBlocks::ProcessBlocks, &oCoverageBlocks,
					std::placeholders::_1));
	BB.openFile(input_BAM);
	std::cout << "\nProcessing the BAM\n";std::cout.flush();
	BB.processAll();

	if (s_inROI != "NULL") {
		// Output computed statistics from data structures.
		//oFragmentsInROI -- this tells us if the data was directional or not -- if we need to know for other output modules.
		std::ofstream outFragmentsInROI;
		outFragmentsInROI.open(outputDir + "/IRFinder-ROI.txt",
				std::ifstream::out);
		oFragmentsInROI.WriteOutput(&outFragmentsInROI);
		outFragmentsInROI.flush();
		outFragmentsInROI.close();
	}

	std::ofstream outJuncCount;
	outJuncCount.open(outputDir + "/IRFinder-JuncCount.txt",
			std::ifstream::out);
	oJuncCount.WriteOutput(&outJuncCount);
	outJuncCount.flush();
	outJuncCount.close();

	int directionality = oJuncCount.Directional();
	std::cout << "RNA-Seq directionality -1/0/+1:\t" << directionality << "\n";

	std::ofstream outSpansPoint;
	outSpansPoint.open(outputDir + "/IRFinder-SpansPoint.txt",
			std::ifstream::out);
	oSpansPoint.WriteOutput(&outSpansPoint);
	outSpansPoint.flush();
	outSpansPoint.close();

	std::ofstream outFragmentsInChr;
	outFragmentsInChr.open(outputDir + "/IRFinder-ChrCoverage.txt",
			std::ifstream::out);
	oFragmentsInChr.WriteOutput(&outFragmentsInChr);
	outFragmentsInChr.flush();
	outFragmentsInChr.close();

	std::ofstream outCoverageBlocks, outCoverageBlocksAI;
	outCoverageBlocks.open(outputDir + "/IRFinder-IR-nondir.txt",
			std::ifstream::out);

	if ( AI_warn > 0 && read_type == "SR" ){
		outCoverageBlocksAI.open(outputDir + "/IRFinder-IR-nondir-AI.txt",
				std::ifstream::out);
	}
	if ( read_type == "LR"){
		oCoverageBlocks.setJitter(std::stoi(jitter));
	}
	oCoverageBlocks.setAI(AI_warn, AI_intron, AI_ratio);
	oCoverageBlocks.WriteOutput(&outCoverageBlocks, &outCoverageBlocksAI, oJuncCount, oSpansPoint);
	outCoverageBlocks.flush();
	outCoverageBlocks.close();
	if ( AI_warn > 0  && read_type == "SR" ){
		outCoverageBlocksAI.flush(); outCoverageBlocksAI.close();
	}
	if (directionality != 0) {
		outCoverageBlocks.open(outputDir + "/IRFinder-IR-dir.txt",
				std::ifstream::out);
		if ( AI_warn > 0 && read_type == "SR" ){
		outCoverageBlocksAI.open(outputDir + "/IRFinder-IR-dir-AI.txt",
			std::ifstream::out);
		}
		oCoverageBlocks.WriteOutput(&outCoverageBlocks,&outCoverageBlocksAI,  oJuncCount, oSpansPoint,
				directionality); // Directional.
		outCoverageBlocks.flush();
		outCoverageBlocks.close();
		if ( AI_warn > 0 && read_type == "SR" ){
			outCoverageBlocksAI.flush(); outCoverageBlocksAI.close();
		}
	}
}

