// ==========================================================================
//     bam2roi - transform a bam/sam file into a region of interest (ROI) file
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Bernd Jagla <bernd.jagla@pasteur.fr>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#ifndef SANDBOX_MY_SANDBOX_APPS_SAMBAMSTATS_ROI_H_
#define SANDBOX_MY_SANDBOX_APPS_SAMBAMSTATS_ROI_H_


#include <fstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io/read_bam.h>
#include <seqan/bam_io/write_bam.h>
#include <seqan/bam_io/write_sam.h>

#define VERSION_ROI_H "0.2"

struct ROI {
	int startPos;          // start position
	int endPos;            // end position
	int rId;               // Reference to BAM header header.sequenceInfos[roi.rId].i1
	seqan::CharString chr; // reference sequence name
	seqan::CharString dir; // strandness of roi (+/-/.)
	int roiID;             // running number of ROI
	seqan::String<unsigned> count; // coverage vector
	seqan::String<unsigned> regionName; // name of the region
	unsigned len;          //length of ROI
	unsigned countMax;     // max coverage
	ROI(seqan::CharString refName) : startPos(seqan::maxValue<int>()), endPos(0)
	{
		dir=".";
		rId = -2;
		len = 0;
		roiID = 0;
		countMax = 0;
	}
	ROI() : startPos(seqan::maxValue<int>()), endPos(0)
		{
		dir=".";
		rId = -2;
		len = 0;
		roiID = 0;
		countMax = 0;
	}
};



// initialize global variables

int RoiAnalyzer::roiN=0;
std::fstream RoiAnalyzer::outFile;

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_SAMBAMSTATS_ROI_H_
