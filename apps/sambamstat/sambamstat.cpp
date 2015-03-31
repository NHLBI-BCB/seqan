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

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

#include "sambamstat.h"

using namespace seqan;





// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    ArgumentParser parser("sambamstat");
    Options options;
    setupCommandLineParser(parser);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    ArgumentParser::ParseResult ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != ArgumentParser::PARSE_OK)
    {
        std::cerr << "Invalid usage!" << std::endl;
        return ret;
    }
    if (options.showHelp || options.showVersion)
        return 0;
    
    if (options.inFormat == FORMAT_AUTO)
    {
        std::ifstream guessIn(toCString(options.inFileName), std::ios::binary | std::ios::in);
        if (!guessIn.good())
        {
            std::cerr << "Could not open " << options.inFileName << std::endl;
            return 1;
        }
        CharString magic;
        resize(magic, 2);
        guessIn.read(&magic[0], 2);
        if (magic != "\x1f\x8b")  // Is not gzip-compressed.
            options.inFormat = FORMAT_SAM;
        else
            options.inFormat = FORMAT_BAM;
    }




    // Finally, launch the program.
    ret = mainWithOptions(options);
    return ret;
}
