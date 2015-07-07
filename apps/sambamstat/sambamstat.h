// ==========================================================================
//                                 sambamstat
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================
// Version 0.1 = NH record
// Version 0.7 = replace BamSamheader instead of appending to it
// Version 0.15 = new roi format
// Version 0.16 remove position sorted options
// Version 0.17
// Version 0.18 based on problem in underlying libraries (problem with reading comments from sam file)
// Version 0.19 adapt to new command line parser
//              some clean-up (rm fasta support,...)
// Version 0.20 remove sam/bam header info
//              removed some variable in function declarations that were not used anymore
// Version 0.21 moved another header info to verbose >=2
// Version 0.22 moved header info (verbose >=2) to stderr

#ifndef SANDBOX_JAGLA_APPS_SAMBAMSTAT_SAMBAMSTAT_H_
#define SANDBOX_JAGLA_APPS_SAMBAMSTAT_SAMBAMSTAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sstream>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/bam_io/read_bam.h>
#include <seqan/bam_io/write_bam.h>
#include <seqan/bam_io/write_sam.h>


#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/arg_parse.h>

//#include "dataanalysis.h"

#define VERSION "0.22"
#define PROGNAME_ID "sambamstats-ID"
//#define PROGNAME_POS "sambamstats-POS"

//#if SEQAN_HAS_ZLIB

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

enum Format {
	FORMAT_AUTO, FORMAT_SAM, FORMAT_BAM
};

//enum SortOrder {
//	SORTORDER_NA, SORTORDER_ID, SORTORDER_POS
//};

struct Stats {
	__uint64 numRecords;
	__uint64 alignedRecords;
	String<double> avrgQuality;
	String<unsigned> nhHisto;
	String<unsigned> intronLHisto;
	String<unsigned> exonLHisto;
	String<unsigned> intronCHisto;
	String<unsigned> lengthsHisto;
	String<unsigned> cigarHisto;
	String<unsigned> flagsHisto;
	String<unsigned> rIDs;
	String<unsigned> insertHisto;
	String<unsigned> deletionHisto;
	String<unsigned> mismatchHisto;
	String<unsigned> editDistanceHisto;
	String<unsigned> roicountMax;
	String<unsigned> roiLenHisto;
	Stats() :
			numRecords(0), alignedRecords(0) {
	}
};

struct Options {
	bool showHelp;
	bool showVersion;
	unsigned int verbosity;
	CharString inFileName;
	CharString outFileName;
	CharString roiFileName;
	CharString tmpDir;
//	bool statistics;
	Format inFormat;
//	SortOrder sortOrder;
	bool forceRedo;
	CharString cmd;
	unsigned maxCoverage;
	bool strandSpecific;

	Options() {
		// Set defaults.
		showHelp = false;
		showVersion = false;
		verbosity = 0;
		//statistics = false;
		forceRedo = false;
		inFileName = "";
		outFileName = "";
		roiFileName = "";
		maxCoverage = 100000000;
		//statistics = false;
		inFormat = FORMAT_AUTO;
//		sortOrder = SORTORDER_NA;
		tmpDir = "";
		strandSpecific = false;
	}
};

double programStartTime=0.0;
int roiCount = 0;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void setupCommandLineParser(ArgumentParser & parser) {
	setVersion(parser, VERSION);
	setShortDescription(parser, "");
	setDate(parser, "March 2015");

	addUsageLine(parser, "[OPTIONS] -if inputFile -of outputFile");

	addDescription(parser, "Calculates basic statics of sam/bam formated files that are sorted by ID. " 
							"Adds attributes to sam/bam files: NH,IH,CC (see sam format specs); " 
							"il (intron length); cl (coverage length); ic (intron count). The program " 
							"also prints out extensive statistics on the distribution of these and similar values: " 
							"Basic statistics, NH histogram, Flags histogram, Cigar histogram, Chromosome coverage, " 
							"intron distribution, intron length distribution, exon length distribution, per position " 
							"quality stats, edit distance distribution. Output format is compliant with FASTQC " 
							"(http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)" 
							"The program makes use of the env variable TMPDIR [/tmp]");

	addSection(parser, "Main Options");
	addOption(parser,
			ArgParseOption("if", "input-file", "Input file.",
					ArgParseArgument::INPUT_FILE, "SAM/BAM" ));
	setValidValues(parser, "input-file", "sam bam");
    setRequired(parser, "input-file", true);
//	addOption(parser,
//			ArgParseOption("S", "input-sam",
//					"Input file is SAM (default: auto)."));
//	addOption(parser,
//			ArgParseOption("B", "input-bam",
//					"Input file is BAM (default: auto)."));
	addOption(parser,
			ArgParseOption("of", "output-file", "Output file.",
					ArgParseArgument::OUTPUT_FILE , "SAM/BAM"));
	setValidValues(parser, "output-file", "sam bam");

	addSection(parser, "optional parameters:");
	addOption(parser,
			ArgParseOption("ss", "strandspecific",
					"calculate strandspecific stats (only position sorted)"));
	addOption(parser,
			ArgParseOption("f", "forceRedo", "enforce recalculation"));
	addOption(parser,
			ArgParseOption("T", "tempdirectory",
					"use DIR for temporaries, not /tmp; multiple options specify multiple directories",
					ArgParseArgument::STRING, "PATH"));
	addSection(parser, "General Options");
	addOption(parser,
			ArgParseOption("vv", "very-verbose", "Very verbose output."));
	addOption(parser,
			ArgParseOption("v", "verbose", "verbose output."));

}

ArgumentParser::ParseResult
parseCommandLineAndCheck(Options & options, ArgumentParser & parser,
		int argc, char const ** argv) {
	ArgumentParser::ParseResult ret = parse(parser, argc, argv);
	if (ret != ArgumentParser::PARSE_OK){
		std::cerr << "error with parameter parsing " << std::endl;
		return  ret;
	}

	getOptionValue(options.inFileName, parser, "input-file" );
	getOptionValue(options.outFileName, parser, "output-file" );
	if (endsWith(options.outFileName,"bam"))
		options.inFormat = FORMAT_BAM;
	else
		options.inFormat = FORMAT_SAM;
	getOptionValue(options.tmpDir, parser, "tempdirectory");
	if (isSet(parser, "verbose"))
		options.verbosity = 1;
	if (isSet(parser, "very-verbose"))
		options.verbosity = 2;
	if (isSet(parser, "forceRedo"))
		options.forceRedo = true;
	if (isSet(parser, "strandspecific"))
		options.strandSpecific = true;
	for (int i = 0; i < argc; i++) {
		options.cmd += argv[i];
		options.cmd += " ";
	}

	return ArgumentParser::PARSE_OK;
}

// record is a BAM record of a possibly split mapped read, sets tags "il" to number of deleted characters from
// CIGAR string, set "cl" to number of covered bases in genome, including intron.

void intronLength(BamAlignmentRecord & record, Stats & stats,
		BamTagsDict & tagDict) {
	bool ok = true;
	String<CigarElement<> > cigar = record.cigar;
	if (length(cigar) > 0) {
		int anzN = 0;
		unsigned coverageLength = 0;
		for (unsigned long i = 0; i < length(cigar); i++) {
			unsigned cLen = length(stats.cigarHisto);
			char op = ((CigarElement<> ) cigar[i]).operation;
			resize(stats.cigarHisto, std::max(cLen, (unsigned) op + 1), 0);
			stats.cigarHisto[((CigarElement<> ) cigar[i]).operation]++;
			coverageLength += ((CigarElement<> ) cigar[i]).count;
			if (((CigarElement<> ) cigar[i]).operation != 'M') {
				unsigned len = length(stats.intronLHisto);
				resize(stats.intronLHisto,
						std::max(len, ((CigarElement<> ) cigar[i]).count + 1),
						0);
				anzN += ((CigarElement<> ) cigar[i]).count;
				stats.intronLHisto[((CigarElement<> ) cigar[i]).count]++;
			} else {
				unsigned len = length(stats.exonLHisto);
				resize(stats.exonLHisto,
						std::max(len, ((CigarElement<> ) cigar[i]).count + 1),
						0);
				stats.exonLHisto[((CigarElement<> ) cigar[i]).count]++;
			}
		}
		ok = setTagValue(tagDict, "il", anzN); //intron length
		if (! ok)
		{
			std::cerr << "error writing tag (1)" << std::endl;
		}
		//setTagValue(tagDict, "nl", record.tLen - anzN, 'i'); //inner length target length - introns for this read
		ok = setTagValue(tagDict, "cl", coverageLength);
		if (! ok)
		{
			std::cerr << "error writing tag (2)" << std::endl;
		}

	}
	
}

// Set tag "ic" to the number of introns in a record.

void intronCount(BamAlignmentRecord & record, Stats & stats,
		BamTagsDict & tagDict) {

	String<CigarElement<> > cigar = record.cigar;
	if (length(cigar) > 0) {
		unsigned anzN = 0;
		for (unsigned long i = 0; i < length(cigar); i++) {
			if (((CigarElement<> ) cigar[i]).operation != 'M')
				anzN++;
		}
		bool ok=setTagValue(tagDict, "ic", anzN);
		if (! ok)
		{
			std::cerr << "error writing tag (3)" << std::endl;
		}
		unsigned len = length(stats.intronCHisto);
		resize(stats.intronCHisto, std::max(len, anzN + 1), 0);
		stats.intronCHisto[anzN]++;
	}
}

// Compute average quality value of quality string.

void averageQV(BamAlignmentRecord & record, Stats & stats) {

	CharString scores = record.qual;

	if (length(scores) > 0) {
		unsigned len = length(stats.avrgQuality);
		resize(stats.avrgQuality,
				std::max(len, unsigned(length(record.seq)) + 1), 0);
		unsigned long sum = 0;
		for (unsigned long i = 0; i < length(scores); i++) {
			sum += scores[i];
			stats.avrgQuality[i + 1] += scores[i];
		}
//		setTagValue(tagDict, "aq", float(int(sum)) / float(length(scores)),
//				'f');
	}
}

//TODO rename is not working for me....
int myrename(char * oldfp, char * newfp) {
	char cmd[2000];
#if defined(_WIN32)
	sprintf(cmd, "move %s %s", oldfp, newfp);
#else
	sprintf(cmd, "mv %s %s", oldfp, newfp);
#endif
	return system(cmd);
}

// Update the length histogram in stats given the record.

void lengthHist(Stats& stats, BamAlignmentRecord& record) {
	//Lengths histogram
	unsigned len = length(stats.lengthsHisto);
	unsigned seqLen = length(record.seq);
	resize(stats.lengthsHisto, std::max(len, seqLen + 1), 0);
	for (unsigned i = 1; i <= seqLen; i++) {
		stats.lengthsHisto[i] += 1;
	}
}

// Create temporary BAM file.
//	createTmpFile(tmpFPnam, options, bFile);

void createTmpFile(char * tmpFPnam, Options const & options,
		BamFileOut & outFile) {
	char * tmpDir = NULL;
	char defaultTMPDIR[4095] = "/tmp";
	tmpDir = getenv("TMPDIR");
	if (tmpDir == NULL) {
		tmpDir = defaultTMPDIR;
	}
	if (options.tmpDir != "") {
		tmpDir = toCString(options.tmpDir);
		std::cerr << "tmpdir  " << tmpDir << " is set" << std::endl;
	}
	std::cerr << "this is options.tempDir " << options.tmpDir << " is set"
			<< std::endl;
	tmpFPnam = strcpy(tmpFPnam, tmpDir);
	tmpFPnam = strcat(tmpFPnam, "/samsBamStatXXXXXX");
#if defined(_WIN32)
	tmpFPnam = strcat(toCString(options.outFileName), ".tmp");
#else
	if (mkstemp(tmpFPnam) == 0) {
		std::cerr << "Cannot create temp file" << std::endl;
		//the tmp file has not been created yet so we can exit
		exit(-1);
	}
#endif
	//options.inFormat = FORMAT_SAM;
	if (options.inFormat == FORMAT_SAM) {
		tmpFPnam = strcat(tmpFPnam, ".sam");
	} else {
		tmpFPnam = strcat(tmpFPnam, ".bam");
	}
	std::cerr << "Using " << tmpFPnam << " as temp file" << std::endl;
	if (options.verbosity >= 2)
		printf("Tempname #1: %s\n", tmpFPnam);
	bool success = open(outFile, toCString(tmpFPnam));
    // Open BamFileOut for writing. Give inFile to share its BamIoContext
	if (!success)
	{	
		std::cerr << "Could not create temp file" << std::endl;
		//the tmp file has not been created yet so we can exit
		exit(-1);
	}
}
//	printStatsIDsorted(stats, reader);

void printStatsIDsorted(Stats stats, BamFileIn & bamFileIn) {


    TBamContext const & bamContext = context(bamFileIn);

	// Print stats
	std::cout << ">>Basic Statistics" << std::endl;
	std::cout << "#Measure\tValue" << std::endl;
	std::cout << "time\t" << sysTime() - programStartTime << "\n";
	std::cout << "num records\t" << stats.numRecords << std::endl;
	std::cout << "aligned records\t" << stats.alignedRecords << std::endl;
	std::cout << "aligned record %\t"
			<< 100.0 * stats.alignedRecords / stats.numRecords << std::endl;
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>NH Histogram" << std::endl;
	std::cout << "#NH value\tcount" << std::endl;
	for (unsigned i = 1; i < length(stats.nhHisto); ++i) {
		if (stats.nhHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.nhHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Flags Histogram" << std::endl;
	std::cout << "#Flag\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.flagsHisto); ++i) {
		if (stats.flagsHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.flagsHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Cigar Histogram" << std::endl;
	std::cout << "#Cigar character\tcount" << std::endl;
	for (unsigned i = 1; i < length(stats.cigarHisto); ++i) {
		if (stats.cigarHisto[i] > 0) {
			std::cout << (char) (i) << '\t';
			std::cout << stats.cigarHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Chromosome coverage" << std::endl;
	std::cout << "#Chromosome\tChr length\treads\tPercentage" << std::endl;
	std::cout << "*\t0\t" << stats.rIDs[0] << "\t0" << std::endl;
	for (unsigned i = 1; i < length(stats.rIDs); i++) {
		std::cout << contigNames(bamContext)[i - 1] << "\t"
				<< contigLengths(bamContext)[i - 1] << "\t" << stats.rIDs[i];
		std::cout << "\t"
				<< 100.0 * double(stats.rIDs[i])
						/ double(contigLengths(bamContext)[i - 1]) << std::endl;
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Intron count Distribution" << std::endl;
	std::cout << "#Number of introns\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.intronCHisto); ++i) {
		if (stats.intronCHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.intronCHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Intron Lengths Distribution" << std::endl;
	std::cout << "#Intron length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.intronLHisto); ++i) {
		if (stats.intronLHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.intronLHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Exon Lengths Distribution" << std::endl;
	std::cout << "#Exon length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.exonLHisto); ++i) {
		if (stats.exonLHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.exonLHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Per position stats" << std::endl;
	std::cout << "#position\t";
	std::cout << "quality" << std::endl;
	for (unsigned i = 1; i < length(stats.lengthsHisto); ++i) {
		std::cout << i << '\t';
		double e = stats.avrgQuality[i] / stats.lengthsHisto[i];
		std::cout << e << std::endl;
	}
	std::cout << ">>END_MODULE" << std::endl;


}

// Why did we do this???
// the only thing that comes to mind is that I believe I had once problem with writing to /dev/null
template<class T> const T& notZero(const T& a, const T& b) {
	if (a != 0)
		return a;
	return b;
}

/*
 * Analyze read that are sorted by ID
 * - check if alignment where unique
 * - can we check if there are biases on the flow cell, i.e. analyze the chip position in correlation with the number of reads
 *   or their quality value
 * - check NM/NH values
 * - analyze cigar string (intron length distribution, length as field
 * - count per chromosome unique/non-unique
 * - number of mutations
 * - quality string analysis (min, mean)
 * - quality over postion
 *
 *
 */
int analyze_idSorted(BamFileIn & reader, Options const & options, BamHeader  & header, TBamContext const & bamContext) {
	//StringSet<CharString> refNames;
	//NameStoreCache<StringSet<CharString> > refNamesCache(refNames);
	// context(refNames, refNamesCache);
	double programStartTime = sysTime();
	String<__uint64> qualSum;
	Align<Dna5String> align;
	//
	int writeFail = 0;
	Stats stats;
	resize(stats.flagsHisto, 65536, 0);
	char tmpFPnam[2000];
	BamFileOut bFile(context(reader));
	BamFileOut bamFileStdOut(context(reader), std::cout, Sam());

	StringSet<CharString> seqIds;

	// check if already run
	// Since there is a bug in IGV that doesn't allow for multiple versions of a given program to be exectued we are removing all old records.
	BamHeader newHeader;
	// newHeader.sequenceInfos = header.sequenceInfos;
	unsigned myIdx = 0;
	for (unsigned i = 0; i < length(header); i++) {
		String<char> porgID = "ID";
		if (findTagKey(myIdx, porgID, header[i])) {
			String<char> tagValue;
			getTagValue(tagValue, myIdx, header[i]);
			if (tagValue == PROGNAME_ID) {
				if (options.verbosity >= 1)
					std::cerr << "already run " << std::endl;
				if (options.forceRedo) {
					erase(header, myIdx);
				} else {
					std::cerr
							<< "Program has already been run, try using -f to force execution "
							<< "\n";
					exit(-1);
				}
			} else {
				appendValue(newHeader, header[i]);
			}
		} else {
			appendValue(newHeader, header[i]);
		}

	}

	unsigned newIdx = length(newHeader);
	resize(newHeader, length(newHeader) + 1);


	newHeader[newIdx].type = BAM_HEADER_PROGRAM;
	resize(newHeader[newIdx].tags, 3);
	newHeader[newIdx].tags[0].i1 = "ID";
	newHeader[newIdx].tags[0].i2 = PROGNAME_ID;
	newHeader[newIdx].tags[1].i1 = "VN";
	newHeader[newIdx].tags[1].i2 = VERSION;
	newHeader[newIdx].tags[2].i1 = "CL";
	newHeader[newIdx].tags[2].i2 = toCString(options.cmd);

	// Create temporary file to ensure that the file is ordered by ID and not overwrite output file
	// if file is not properly sorted the header section will not be written correctly
	createTmpFile(tmpFPnam, options, bFile);

	writeHeader(bFile, newHeader);
//	if (options.verbosity >= 2)
//		writeHeader(bamFileStdOut, newHeader);

	resize(stats.rIDs, length(contigNames(bamContext)) + 1, 0);
	StringSet<BamAlignmentRecord> idSet;
	BamAlignmentRecord oldRec;
	oldRec.qName = "";
	BamAlignmentRecord record;

	// Read alignments.
	if (options.verbosity >= 2)
		std::cerr << "Reading alignments" << std::endl;
	unsigned readCounter = 0;



	while (!atEnd(reader)) {
		readCounter++;
		if (!(readCounter % 100000)) {
			if (readCounter == 100000) {
				std::cerr
						<< "N(seq)\ttime[sec]\tRSS[kB]\tShared Memory[kB]\tPrivate Memory[kB]\n";
			}
			std::cerr << readCounter << "\t" << sysTime() - programStartTime
					<< "\t";
#if defined(_WIN32)
#else
			int tSize = 0, resident = 0, share = 0;
			std::ifstream buffer("/proc/self/statm");
			buffer >> tSize >> resident >> share;
			buffer.close();
			long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
			double rss = resident * page_size_kb;
			std::cerr << rss << "\t";
			double shared_mem = share * page_size_kb;
			std::cerr << shared_mem << "\t";
			std::cerr << rss - shared_mem;
#endif
			std::cerr << "\n";
		}
		// Read alignment record.
		readRecord(record, reader );
		stats.numRecords++;
		stats.alignedRecords += !hasFlagUnmapped(record);
		stats.flagsHisto[record.flag]++;

		BamTagsDict tagDict(record.tags);

		//Lengths histogram
		lengthHist(stats, record);

		// count reference ids
		stats.rIDs[(record.rID == BamAlignmentRecord::INVALID_REFID ? 0 : record.rID + 1)]++;


		/*
		 * we basically can only handle one group at a time any other property has to be calculated for each record
		 */
		// MD,NM distance to reference (should already be calculated... We would need the reference genome for this and compare the it
		// therefore we assume that this has already been correctly calculated.
		// H0, H1, H2 this would be cool to have, but it is too expensive to calculate
		// ni number of introns
		intronCount(record, stats, tagDict);

		// il total intron length
		intronLength(record, stats, tagDict);

		// av average quality value
		averageQV(record, stats);

		if (options.verbosity >= 2)
			writeRecord(bamFileStdOut, record);
		/*
		 * IH tag: Query hit index, indicating the alignment record is the i-th one stored
		 * NH tag: Number of reported alignments that contains the query in the current record
		 */
		//Only occurs the first time
		if (oldRec.qName == "") {
			oldRec = record;
			appendValue(idSet, record);
		} else {
			// same query Name => collect in set
			if (oldRec.qName == record.qName) {
				appendValue(idSet, record);
			} else {
				// different Name => write out old...
				if (oldRec.rID == BamAlignmentRecord::INVALID_REFID) {
					BamTagsDict tagDictid(oldRec.tags);
					setTagValue(tagDictid, "NH", 0);
					writeRecord(bFile, oldRec);
				} else {
					unsigned nID = length(idSet);
					for (unsigned i = 0; i < nID; i++) {
						BamAlignmentRecord bar;
						clear(bar);
						bar = idSet[i];
						BamTagsDict tagDictid(bar.tags);
						bool ok = setTagValue(tagDictid, "NH", nID);
						ok = setTagValue(tagDictid, "IH", i);
						// CC reference name of the next hit; "=" for the same chromosome
						// TODO once Strings can be stored, store them as Strings.
						// i.e. currently SEQAN is not able to handle strings...!!!
						// XA  iMDZ49 NM iiciilicli1NHiIHi
						if (! ok)
						{
							std::cerr << "error writing tag (4)" << std::endl;
						}
						if (i < nID - 1) {
							if (idSet[i].rID == idSet[i + 1].rID) {
								ok = setTagValue(tagDictid, "CC", '=');
							} else {
								ok = setTagValue(tagDictid, "CC", idSet[i + 1].rID);
							}
						}
						writeRecord(bFile, bar);
					}
					unsigned len = length(stats.nhHisto);
					resize(stats.nhHisto, std::max(len, nID + 1), 0);
					stats.nhHisto[nID] += 1;
				}
				clear(idSet);
				appendValue(idSet, record);
			}
		}
		oldRec = record;

	}
	// NH tag
	//finally write out last group
	int nID = length(idSet);
	for (int i = 0; i < nID; i++) {
		BamAlignmentRecord bar = idSet[i];
		BamTagsDict tagDictid(bar.tags);
		bool ok = setTagValue(tagDictid, "NH", nID);
		if (! ok)
		{
			std::cerr << "error writing tag (5)" << std::endl;
		}

		writeRecord(bFile, bar);

	}

	close(bFile);
	if (options.outFileName != "" && writeFail == 0) {
		if (myrename(tmpFPnam, toCString(options.outFileName)) != 0) {
			perror("there was a problem renaming the output file...\n");
		}
	} else {
		if (writeFail != 0) {
			std::cerr
					<< "WARNINIG: there was a problem writing to tmp file ; removing tmp file\n";
		} else {
			std::cerr << "WARNINIG: no output file given; removing tmp file\n";
		}
		remove(tmpFPnam);
	}

	// Print stats
	printStatsIDsorted(stats, reader);
	return 0;
}

bool errorMessageSent = false;


/*
 * Here we deal with the sorting and organizing different order of events...
 */
ArgumentParser::ParseResult
doWork(BamFileIn & reader, 
		Options const & options) {
	String<__uint64> qualSum;

// Read header.
	BamHeader header ;
	if (options.verbosity >= 2)
		std::cerr << "Reading header" << std::endl;
	readHeader(header, reader);

	std::cerr << "old header" << std::endl;
	BamFileOut bamFileStdErr(std::cerr, Sam());
	if (options.verbosity >= 2){
		writeHeader(bamFileStdErr, header);
	}
	TBamContext const & bamContext = context(reader);

	if (length(contigNames(bamContext)) == 0) {
		std::cerr
				<< "There was no information on sequences in the header or header is missing. Please correct\n"
				<< std::endl;
		return ArgumentParser::PARSE_ERROR;

	}
	int ret = 0;
	ret = analyze_idSorted(reader, options, header, bamContext);
	if (ret != 0) {
		//TODO remove temp file if not -vv
	}

	return ArgumentParser::PARSE_OK;
}

ArgumentParser::ParseResult
mainWithOptions(Options & options) {
	//typedef Iterator<String<CharString> >::Type TIterator;
	programStartTime = sysTime();
	std::cout << ">>Parameters\n";
	std::cout << "#Name\tvalue" << std::endl;
	std::cout << "Version\t" << VERSION << "\n";
	std::cout << "input file: \t\"" << options.inFileName << "\"" << std::endl;
	std::cout << "output file: \t\"" << options.outFileName << "\""
			<< std::endl;
	if (options.inFormat == FORMAT_SAM)
		std::cout << "file format: \tSAM" << std::endl;
	else
		std::cout << "file format: \tBAM" << std::endl;
	if (options.forceRedo)
		std::cout << "forceRedo\tyes" << std::endl;
	else
		std::cout << "forceRedo\tno" << std::endl;
	std::cout << "tmpDir\t" << options.tmpDir << "" << std::endl;
	if (options.strandSpecific)
		std::cout << "strand specific\tyes" << std::endl;
	else
		std::cout << "strand specific\tno" << std::endl;
	std::cout << ">>END_MODULE" << std::endl;

	StringSet<CharString> seqIds;
	String<char, MMap<> > seqMMapString;


	seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(options.inFileName)))
	{
		std::cerr << "ERROR: Could not open " << options.inFileName << "\n";
        return ArgumentParser::PARSE_ERROR;
	}
    return doWork(bamFileIn, options);
}

#endif  // #ifndef SANDBOX_JAGLA_APPS_SAMBAMSTAT_SAMBAMSTAT_H_
