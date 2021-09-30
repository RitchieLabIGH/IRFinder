// WARNING: code is little endian only!

#include "BAM2blocks.h"
// using namespace std;

//const char cigarChar[] = {'M','I','D','N','S','H','P','=','X'};

BAM2blocks::BAM2blocks() :
		instream(std::cin.rdbuf()) {
	oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?
	cShortPairs = 0;
	cIntersectPairs = 0;
	cLongPairs = 0;
	cSingleReads = 0;
	cPairedReads = 0;
	cErrorReads = 0;
	cSkippedReads = 0;
}

// OK.
void BAM2blocks::readBamHeader() {
	char buffer[1000];
	std::string chrName;
	//std::vector<std::string> chr_names;
	bam_header bamhead;

	IN->read(bamhead.c, BAM_HEADER_BYTES);

	//IN->ignore(bamhead.l_text);
	char headertext[bamhead.l_text + 1];
	IN->read(headertext, bamhead.l_text);
	samHeader = std::string(headertext, bamhead.l_text);
	coord_sorted = samHeader.find("SO:coordinate") != std::string::npos;
	stream_int32 i32;
	IN->read(i32.c, 4);
	uint n_chr = i32.i;

	for (uint i = 0; i < n_chr; i++) {
		IN->read(i32.c, 4);
		IN->read(buffer, i32.i);
		chrName = std::string(buffer, i32.i - 1);
		chr_names.push_back(chrName);

		IN->read(i32.c, 4);
		chr_lens.push_back(i32.i);
	}

	for (auto &callback : callbacksChrMappingChange) {
		callback(chr_names);
	}

}

void BAM2blocks::cigar2block(int32_t *cigar, uint16_t n_cigar_op,
		std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len) {
	bool inBlock = true;
	int relpos = 0;
	int curblock = 0;
	starts.resize(1); // Is this expensive or not -- does this call destroy on further items, or is it a single op, adjusting the end? If expensive we can revert to earlier behaviour where we keep track of how many blocks, just overwriting relevant parts of the vector.
	lens.resize(1);
	starts[curblock] = 0;
	lens[curblock] = 0;

	for (; n_cigar_op > 0; n_cigar_op--) {
		if (inBlock) {
			switch (*cigar & 15) {
			case 0:
			case 2:
			case 7:
			case 8:
				// increment len of last block
				lens[curblock] += (*cigar >> 4);
				relpos += (*cigar >> 4);
				break;
			case 3:
				curblock++;
				relpos += (*cigar >> 4);
				// extend arrays.
				starts.push_back(relpos);
				lens.push_back(0);
				inBlock = false;
				break;
			}
		} else {
			switch (*cigar & 15) {
			case 0:
			case 2:
			case 7:
			case 8:
				lens[curblock] = (*cigar >> 4);
				relpos += (*cigar >> 4);
				inBlock = true;
				break;
			case 3:
				// push start of next further out
				relpos += (*cigar >> 4);
				starts[curblock] = relpos;
				break;
			}
		}
		cigar++;
	}
	ret_genome_len = relpos;
//  *ret_blocks = curblock+1;  // Unnecessary if we are using vectors in the expected manner - changing their length as needed.
}

//OK - translated - doesn't call the callbacks yet though.
unsigned int BAM2blocks::processPair(bam_read_core *read1,
		bam_read_core *read2) {
	// R1 is to the left of R2 (or equal starts).
	int r1_genome_len;
	//int r1_blocks;
	int r2_genome_len;

	std::string debugstate;

	//int r2_blocks;
	//char dir;

	if (read1->flag & 0x40) {
		//this is first of pair.
		if (read1->flag & 0x10) {
			oBlocks.direction = 0;
		} else {
			oBlocks.direction = 1;
		}
	} else {
		if (read1->flag & 0x20) {
			oBlocks.direction = 0;
		} else {
			oBlocks.direction = 1;
		}
	}

	cigar2block(read1->cigar, read1->n_cigar_op, oBlocks.rStarts[0],
			oBlocks.rLens[0], r1_genome_len);
	cigar2block(read2->cigar, read2->n_cigar_op, oBlocks.rStarts[1],
			oBlocks.rLens[1], r2_genome_len);

	if (read1->pos + r1_genome_len < read2->pos) {
		cLongPairs++;
		//reads do not intersect
		oBlocks.readCount = 2;
		debugstate.append("-Long-");
	} else if (read1->pos + r1_genome_len >= read2->pos + r2_genome_len) {
		cShortPairs++;
		// Read 2 is a short read & read 1 fully contains it (or perhaps just a trimmed read with two exactly complementary reads remaining).
		oBlocks.readCount = 1;
		debugstate.append("-Short-");
	} else {
		debugstate.append("-Intersect-");
		cIntersectPairs++;
		bool goodPair = true;
		oBlocks.readCount = 1;
		// We have two reads that intersect - construct just one complete fragment.

// Guaranteed assumptions:
//   Read 1 starts to the left of Read 2.
//   Read 2 end extends beyond the end of Read 1 end.
		int r1pos = read1->pos;
		int r2pos = read2->pos;
		for (uint i = 0; i < oBlocks.rStarts[0].size(); i++) {
			if (r1pos + oBlocks.rStarts[0][i] + oBlocks.rLens[0][i] >= r2pos) {
				if (r1pos + oBlocks.rStarts[0][i] <= r2pos) {
					oBlocks.rLens[0][i] = r2pos - r1pos - oBlocks.rStarts[0][i]
							+ oBlocks.rLens[1][0];
					//r1_blocks = i + r2_blocks;
					oBlocks.rStarts[0].resize(i + oBlocks.rStarts[1].size());
					oBlocks.rLens[0].resize(i + oBlocks.rStarts[1].size());
					// Maybe this can be optimised by using push_back below instead of running resize.
					for (uint j = 1; j < oBlocks.rStarts[1].size(); j++) {
						i++;
						oBlocks.rLens[0][i] = oBlocks.rLens[1][j];
						oBlocks.rStarts[0][i] = oBlocks.rStarts[1][j] + r2pos
								- r1pos;
					}
					r1_genome_len = r2pos - r1pos + r2_genome_len;
					break;
				} else {
					//cerr << "Fault with this synthetic read, outputting each of the overlapping reads as singles: " << read1->read_name << endl;
					// This error is not worth reporting. The current version of STAR outputs a good number of these, concordance would be nice, but it is better to get at least one read illustrating the splice junction.
					goodPair = false;
					oBlocks.readCount = 2;
				}
			}
		}

		if (!goodPair) {
			oBlocks.readCount = 2;
		}
	}
	oBlocks.chr_id = read1->refID;
	oBlocks.readStart[0] = read1->pos;
	oBlocks.readEnd[0] = read1->pos + r1_genome_len;
	oBlocks.readName.resize(read1->l_read_name - 1);
	oBlocks.readName.replace(0, read1->l_read_name - 1, read1->read_name,
			read1->l_read_name - 1); // is this memory/speed efficient?

	uint totalBlockLen = 0;
	for (auto blockLen : oBlocks.rLens[0]) {
		totalBlockLen += blockLen;
	}
	if (oBlocks.readCount > 1) {
		oBlocks.readStart[1] = read2->pos;
		oBlocks.readEnd[1] = read2->pos + r2_genome_len;
		for (auto blockLen : oBlocks.rLens[1]) {
			totalBlockLen += blockLen;
		}
	}
	//DEBUG:
	oBlocks.readName.append(debugstate);
	oBlocks.readName.append(std::to_string(oBlocks.readCount));
// TODO - restructure -- we could instead do the manipulation from 2 reads-> 1 synthetic in a non-const callback.
//        not required until that future flexibility is needed if part of the framework is repurposed.
	for (auto &callback : callbacksProcessBlocks) {
		callback(oBlocks);
	}
	return totalBlockLen;
}

unsigned int BAM2blocks::processSingle(bam_read_core *read1) {
	int r1_genome_len;

	std::string debugstate;

	if (read1->flag & 0x10) {
		oBlocks.direction = 0;
	} else {
		oBlocks.direction = 1;
	}
	cigar2block(read1->cigar, read1->n_cigar_op, oBlocks.rStarts[0],
			oBlocks.rLens[0], r1_genome_len);

	oBlocks.readCount = 1;
	oBlocks.chr_id = read1->refID;
	oBlocks.readStart[0] = read1->pos;
	oBlocks.readEnd[0] = read1->pos + r1_genome_len;
	oBlocks.readName.resize(read1->l_read_name - 1);
	oBlocks.readName.replace(0, read1->l_read_name - 1, read1->read_name,
			read1->l_read_name - 1); // is this memory/speed efficient?
	//DEBUG:
	oBlocks.readName.append(debugstate);
	oBlocks.readName.append(std::to_string(oBlocks.readCount));
	//cout << "process pair - callbacks" << endl;  
	for (auto &callback : callbacksProcessBlocks) {
		callback(oBlocks);
	}
	uint totalBlockLen = 0;
	for (auto blockLen : oBlocks.rLens[0]) {
		totalBlockLen += blockLen;
	}
	return totalBlockLen;
}

void BAM2blocks::setMate(std::vector<char> & mate) {
	uint64_t i=0, j=0;
	while (i < BAM_READ_CORE_BYTES){
		tmp_mate.c[j++]=mate[i++];
	}

	j=0;
	while (i < BAM_READ_CORE_BYTES+tmp_mate.l_read_name ){
		tmp_mate.read_name[j++]=mate[i++];
	}
	j=0;
	while (i < mate.size()){
		tmp_mate.cigar_buffer[j++]=mate[i++];
	}
}
void BAM2blocks::saveMate() {
	std::vector<char> mate;
	for (int i=0; i< BAM_READ_CORE_BYTES; i++) {
		mate.push_back(tmp_read.c[i]);
	}
	for (int i = 0; i < tmp_read.l_read_name; i++) {
		mate.push_back(tmp_read.read_name[i]);
	}
	for (int i = 0; i < tmp_read.n_cigar_op * 4; i++) {
		mate.push_back(tmp_read.cigar_buffer[i]);
	}
	tmp_reads[getName(tmp_read)] = mate;
}

std::string BAM2blocks::getName(bam_read_core & read ){
	std::string name="";
	for (int i=0; i<read.l_read_name; i++) {
		name+=read.read_name[i];
	}
	return name;
}

int BAM2blocks::processAll() {

	totalNucleotides = 0;
	current_read = 0;
	//int bytesread = 0;
	bool running = getNextReadHead(tmp_read);
	while (running) {
		current_read++;
		if (IN->fail()) {
			errorMessage();
			return (1);
			//This is possibly also just about the end of the file (say an extra null byte).
			//IN->gcount() knows how many characters were actually read last time.
		}
		if (!(tmp_read.flag & 0x1) || !(tmp_read.flag & 0x2)
				|| (tmp_read.flag & 0x8)) {
			/* If is a single read ( or mate is unmapped or are not mapped in a proper way ) -- process it as a single -- then discard/overwrite */
			getReadBody(tmp_read);
			cSingleReads++;
			totalNucleotides += processSingle(&tmp_read);
		} else {
			/* If it is potentially a paired read, store it in our buffer, process the pair together when it is complete */
			getReadBody(tmp_read);
			if (coord_sorted) {
				std::string name=getName(tmp_read);
				if (tmp_reads.count(name) == 1) {
					setMate(tmp_reads[name]);
					handlePairs(tmp_mate, tmp_read);
					tmp_reads.erase(name);
				} else {
					saveMate();
				}
			} else {
				current_read++;
				running = getNextReadHead(tmp_mate);
				if (running) {
					getReadBody(tmp_mate);
					handlePairs(tmp_read, tmp_mate);
				}
			}
		}
		running = getNextReadHead(tmp_read);
	}
	if (coord_sorted && tmp_reads.size() > 0) {
		std::cerr << "WARNING! " << tmp_reads.size()
				<< " reads have no mates. Processed a single end." << std::endl;
		for (auto &pair : tmp_reads) {
			cSingleReads++;
			setMate(pair.second);
			totalNucleotides += processSingle(&tmp_mate);
		}
	}
	std::cout << "Total reads processed: " << current_read - 1 << std::endl;
	std::cout << "Total nucleotides: " << totalNucleotides << std::endl;
	std::cout << "Total singles processed: " << cSingleReads << std::endl;
	std::cout << "Total pairs processed: "
			<< cShortPairs + cIntersectPairs + cLongPairs << std::endl;
	std::cout << "Short pairs: " << cShortPairs << std::endl;
	std::cout << "Intersect pairs: " << cIntersectPairs << std::endl;
	std::cout << "Long pairs: " << cLongPairs << std::endl;
	std::cout << "Skipped reads: " << cSkippedReads << std::endl;
	for (auto reason : skippedReason) {
		std::cout << " - flag " << reason.first << ": " << reason.second
				<< std::endl;
	}
	std::cout << "Error reads: " << cErrorReads << std::endl;
	return (0);
}

void BAM2blocks::errorMessage() {
	std::cerr << "Input error at line:" << current_read << std::endl;
	std::cerr << "Characters read on last read call:" << IN->gcount()
			<< std::endl;
	std::cout << "ERR-Total reads processed: " << current_read - 1 << std::endl;
	std::cout << "ERR-Total nucleotides: " << totalNucleotides << std::endl;
	std::cout << "ERR-Total singles processed: " << cSingleReads << std::endl;
	std::cout << "ERR-Total pairs processed: "
			<< cShortPairs + cIntersectPairs + cLongPairs << std::endl;
	std::cout << "ERR-Short pairs: " << cShortPairs << std::endl;
	std::cout << "ERR-Intersect pairs: " << cIntersectPairs << std::endl;
	std::cout << "ERR-Long pairs: " << cLongPairs << std::endl;
	std::cout << "ERR-Skipped reads: " << cSkippedReads << std::endl;
	for (auto reason : skippedReason) {
		std::cout << " - flag " << reason.first << ": " << reason.second
				<< std::endl;
	}
	std::cout << "ERR-Error reads: " << cErrorReads << std::endl;
}

void BAM2blocks::handlePairs(bam_read_core &read, bam_read_core &mate) {
	if (read.refID == mate.refID && read.l_read_name == mate.l_read_name
			&& strncmp(read.read_name, mate.read_name, read.l_read_name) == 0) {
		if (read.pos <= mate.pos) {
			totalNucleotides += processPair(&read, &mate);
		} else {
			totalNucleotides += processPair(&mate, &read);
		}
		cPairedReads++;
	} else {
		std::cout <<" A: " << getName(read) << " B: " << getName(mate) << "\n";
		cErrorReads++;
		// Now this should never happend
	}
}

bool BAM2blocks::getNextReadHead(bam_read_core &read) {
	if (IN->eof()) {
		return false;
	}
	IN->read(read.c, BAM_READ_CORE_BYTES);
	if ((tmp_read.flag & 0x900) || (tmp_read.flag & 0x4)) {
		IN->ignore(read.block_size - BAM_READ_CORE_BYTES + 4);
		if (skippedReason.count(tmp_read.flag) == 0) {
			skippedReason[tmp_read.flag] = 1;
		} else {
			skippedReason[tmp_read.flag]++;
		}
		cSkippedReads++;
		return getNextReadHead(read);
	}
	return true;
}

void BAM2blocks::getReadBody(bam_read_core &read) {
	IN->read(read.read_name, read.l_read_name);
	IN->read(read.cigar_buffer, read.n_cigar_op * 4);
	IN->ignore(
			read.block_size - BAM_READ_CORE_BYTES + 4 - read.l_read_name
					- (read.n_cigar_op * 4));
}

void BAM2blocks::openFile(std::istream *_IN) {
	IN = _IN;
	readBamHeader(); // readBamHeader needs to call the ChrMappingChange callbacks.
}

void BAM2blocks::openFile(std::string in_file) {
	file.open(in_file, std::ios::binary);
	inbuf.empty();
	inbuf.push(boost::iostreams::gzip_decompressor());
	inbuf.push(file);
	//Convert streambuf to istream
	instream.rdbuf(&inbuf);
	IN = &instream;
	readBamHeader();
}

void BAM2blocks::registerCallbackChrMappingChange(
		std::function<void(const std::vector<std::string>&)> callback) {
	callbacksChrMappingChange.push_back(callback);
}

void BAM2blocks::registerCallbackProcessBlocks(
		std::function<void(const FragmentBlocks&)> callback) {
	callbacksProcessBlocks.push_back(callback);
}
