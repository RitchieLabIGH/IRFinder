#ifndef CODE_BAM2BLOCKS
#define CODE_BAM2BLOCKS

#include "FragmentBlocks.h"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/* Little Endian .. for big endian each group of 4 bytes needs to be reversed before individual members are accessed. */
// std c11 allows anonymous struct/union. -Wall may give a warning as non-portable to older c++ standards.



class BAM2blocks {

	// TODO -- are structs best hidden inside the class? Does doing so push them into namespace of the class only?
	struct bam_read_core {
		union {
		  char c[36];
		  struct {
			int32_t block_size;
			int32_t refID;
			int32_t pos;
			uint8_t l_read_name;
			uint8_t mapq;
			uint16_t bin;
			uint16_t n_cigar_op;
			uint16_t flag;
			int32_t l_seq;
			int32_t next_refID;
			int32_t next_pos;
			int32_t tlen;
		  }; // anonymous struct to allow easy access to members.
		};
		char read_name[256];
		union {
		  char cigar_buffer[20000];
		  int32_t cigar[5000];
		};
	};
 
	union bam_header {
		char c[8];
		struct {
		  char magic[4];
		  int32_t l_text;
		};
	};

	union stream_int32 {
		char c[4];
		int32_t i;
	};

	static const int BAM_HEADER_BYTES = 8;
	static const int BAM_READ_CORE_BYTES = 36;
	static const int BAM_READ_CORE_MAX_CIGAR = 20000;

	FragmentBlocks oBlocks;

	std::vector< std::function<void(const std::vector<std::string> &)> > callbacksChrMappingChange;
	std::vector< std::function<void(const FragmentBlocks &)> > callbacksProcessBlocks;

	// Statistics.
	ulong cShortPairs;
	ulong cIntersectPairs;
	ulong cLongPairs;
	ulong cSingleReads;
	ulong cPairedReads;
	ulong cErrorReads;
	ulong cSkippedReads;
	uint64_t totalNucleotides;
	std::map<uint16_t, uint> skippedReason;

	std::map<std::string, std::vector<char>> tmp_reads;
	bam_read_core tmp_read;
	bam_read_core tmp_mate;
	uint64_t current_read=0;

	bool getNextReadHead(bam_read_core &);
	void errorMessage();
	void getReadBody(bam_read_core &);
	void handlePairs(bam_read_core &, bam_read_core &);
	std::string getName(bam_read_core &);
	void setMate(std::vector<char> & mate);
	void saveMate();
	std::istream * IN;
	std::istream instream;
	void cigar2block(int32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len);

	unsigned int processPair(bam_read_core * read1, bam_read_core * read2);
	unsigned int processSingle(bam_read_core * read1);

	std::vector<unsigned char> stream_buffer;
	void fillBuffer();
	std::ifstream file;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	bool coord_sorted=false;
public:
  	BAM2blocks();
  	void openFile(std::istream * _IN);
  	void openFile(std::string in_file);
  	void readBamHeader();  // implied by openFile. So perhaps should be private.
  	int processAll();

	void registerCallbackChrMappingChange( std::function<void(const std::vector<std::string> &)> callback );
	void registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback );

	std::string samHeader;
	std::vector<std::string> chr_names;   //tab terminated chromosome names.
	std::vector<int32_t> chr_lens;	//length of each chromosome (not used when reading, used if optionally outputting an altered BAM file)
};


#endif
