// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Code to convert between SAM and BAM format tags (textual <-> binary).
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignTagsSamToBam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TBuffer>
struct AssignTagsSamToBamOneTagHelper_
{
    TTarget     &target;
    TBuffer     buffer;
    char        typeC;
    std::string tmpBuffer;

    AssignTagsSamToBamOneTagHelper_(TTarget &target, TBuffer buffer, char typeC):
        target(target),
        buffer(buffer),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type)
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        union {
            char raw[sizeof(Type)];
            Type i;
        } tmp;

        tmp.i = lexicalCast<Type>(buffer);
        append(target, toRange(&tmp.raw[0], &tmp.raw[sizeof(Type)]));
        return true;
    }

    // we have to make this workaround until lexical_cast<float|double> can cope with Segments
    bool operator() (float)
    {
        if (BamTypeChar<float>::VALUE != typeC)
            return false;

        union {
            char raw[sizeof(float)];
            float i;
        } tmp;

        assign(tmpBuffer, buffer);
        tmp.i = lexicalCast<float>(tmpBuffer);
        append(target, toRange(&tmp.raw[0], &tmp.raw[sizeof(float)]));
        return true;
    }
};

template <typename TTarget, typename TForwardIter>
void _assignTagsSamToBamOneTag(TTarget & target, TForwardIter & iter, CharString & buffer)
{
    readUntil(target, iter, CountDownFunctor<>(2));

    clear(buffer);
    readUntil(buffer, iter, CountDownFunctor<>(3));
    SEQAN_ASSERT_EQ(buffer[0], ':');
    SEQAN_ASSERT_EQ(buffer[2], ':');

    char typeC = buffer[1];
    appendValue(target, typeC);

    switch (typeC)
    {
        case 'Z':
        case 'H':
            // BAM string
            // TODO(holtgrew): Could test on even length in case of 'H'.
            readUntil(target, iter, OrFunctor<IsTab, IsNewline>());
            appendValue(target, '\0');
            break;

        case 'B':
        {
            // BAM array

            // Read type.
            readOne(typeC, iter);
            appendValue(target, typeC);

            // Read array contents.
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());

            size_t len = length(buffer);

            // Count number of entries (== number of commas after type character).
            __uint32 nEntries = 0;
            for (size_t i = 0; i != len; ++i)
                if (buffer[i] == ',')
                    ++nEntries;

            // Write out array length.
            appendRawPod<__uint32>(target, nEntries);

            // Write out array values.
            typedef typename Infix<CharString>::Type TBufferInfix;
            size_t startPos = 0;
            for (unsigned i = 0; i < nEntries; ++i)
            {
                SEQAN_ASSERT_LT(startPos, len);
                if (buffer[startPos] == ',')
                    ++startPos;

                // search end of current entry
                size_t endPos = startPos;
                for (; endPos < len; ++endPos)
                    if (buffer[endPos] == ',' || buffer[endPos] == '\t')
                        break;

                AssignTagsSamToBamOneTagHelper_<TTarget, TBufferInfix> func(target,
                                                                            infix(buffer, startPos, endPos),
                                                                            typeC);
                if (!tagApply(func, BamTagTypes()))
                    SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);

                startPos = endPos;
            }
            break;
        }

        default:
        {
            // BAM simple value
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());

            AssignTagsSamToBamOneTagHelper_<TTarget, CharString&> func(target, buffer, typeC);
            if (!tagApply(func, BamTagTypes()))
                SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);
        }
    }
}

/*!
 * @fn assignTagsSamToBam
 * @headerfile <seqan/bam_io.h>
 * @brief Assign tags in SAM format to tags in BAM format.
 *
 * @signature void assignTagsBamToSam(bamTags, samTags);
 *
 * @param[out] bamTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the target BAM tags.
 * @param[in]  samTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the source SAM tags.
 *
 * @see assignTagsBamToSam
 */

/**
.Function.assignTagsSamToBam
..cat:BAM I/O
..summary:Assign tags in SAM format to tags in BAM format.
..signature:assignTagsSamToBam(bamTags, samTags)
..param.bamTags:Destination BAM tags.
...type:Shortcut.CharString
..param.samTags:Source SAM tags.
...type:Shortcut.CharString
..returns:$void$
..include:seqan/bam_io.h
*/

template <typename TTarget, typename TSource>
void assignTagsSamToBam(TTarget & target, TSource const & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        return;

    typedef typename Iterator<TSource const, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());

    CharString buffer;

    while (!atEnd(it))
    {
        if (value(it) == '\t')
            skipOne(it);

        _assignTagsSamToBamOneTag(target, it, buffer);
    }
}

// ----------------------------------------------------------------------------
// Function assignTagsBamToSam()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSourceIter>
struct AssignTagsBamToSamOneTagHelper_
{
    TTarget     &target;
    TSourceIter &it;
    char        typeC;

    AssignTagsBamToSamOneTagHelper_(TTarget &target, TSourceIter &it, char typeC):
        target(target),
        it(it),
        typeC(typeC)
    {}

    template <typename Type>
    bool operator() (Type)
    {
        if (BamTypeChar<Type>::VALUE != typeC)
            return false;

        appendNumber(target, reinterpret_cast<Type const &>(*it));
        it += sizeof(Type);
        return true;
    }

    bool operator() (char)
    {
        if (BamTypeChar<char>::VALUE != typeC)
            return false;

        appendValue(target, getValue(it));
        ++it;
        return true;
    }
};

template <typename TTarget, typename TSourceIter>
void _assignTagsBamToSamOneTag(TTarget & target, TSourceIter & it)
{   
    // Copy tag name.
    SEQAN_ASSERT_NOT(atEnd(it));
    appendValue(target, *it++);
    SEQAN_ASSERT_NOT(atEnd(it));
    appendValue(target, *it++);

    // Add ':'.
    appendValue(target, ':');
    
    char typeC = *it++;
    char c = FunctorLowcase<char>()(typeC);

    // The only integer type supported is a 32bit signed int (SAM Format Spec, 28 Feb 2014, Section 1.5)
    // This sucks as this projection is not identically reversible
    if (c == 'c' || c == 's' || c == 'i')
        appendValue(target, 'i');
    else
        appendValue(target, typeC);

    // Add ':'.
    appendValue(target, ':');

    switch (typeC)
    {
        case 'Z':
        case 'H':
            // BAM string
            SEQAN_ASSERT_NOT(atEnd(it));
            while (*it != '\0')
            {
                appendValue(target, *it);
                ++it;
                SEQAN_ASSERT_NOT(atEnd(it));
            }
            ++it;
            break;

        case 'B':
        {
            // BAM array
            typeC = *it++;
            appendValue(target, typeC);
            AssignTagsBamToSamOneTagHelper_<TTarget, TSourceIter> func(target, it, typeC);
            
            // Read array length.
            union {
                char raw[4];
                unsigned len;
            } tmp;
            for (unsigned i = 0; i < 4; ++i)
            {
                SEQAN_ASSERT_NOT(atEnd(it));
                tmp.raw[i] = *it++;
            }
            for (unsigned i = 0; i < tmp.len; ++i)
            {
                appendValue(target, ',');
                if (!tagApply(func, BamTagTypes()))
                    SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);
            }
            break;
        }

        default:
        {
            // BAM simple value
            AssignTagsBamToSamOneTagHelper_<TTarget, TSourceIter> func(target, it, typeC);
            if (!tagApply(func, BamTagTypes()))
                SEQAN_ASSERT_FAIL("Invalid tag type: %c!", typeC);
        }
    }
}

/*!
 * @fn assignTagsBamToSam
 * @headerfile <seqan/bam_io.h>
 * @brief Assign tags in BAM format to tags in SAM format.
 *
 * @signature void assignTagsBamToSam(samTags, bamTags);
 *
 * @param[out] samTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the target SAM tags.
 * @param[in]  bamTags A sequence of <tt>char</tt> (e.g. @link CharString @endlink) for the source BAM tags.
 *
 * @see assignTagsSamToBam
 */

/**
.Function.assignTagsBamToSam
..cat:BAM I/O
..summary:Assign tags in BAM format to tags in SAM format.
..signature:assignTagsSamToBam(bamTags, samTags)
..param.samTags:Destination SAM tags.
...type:Shortcut.CharString
..param.bamTags:Source BAM tags.
...type:Shortcut.CharString
..returns:$void$
..include:seqan/bam_io.h
..see:Function.assignTagsSamToBam
*/

template <typename TTarget, typename TSource>
void assignTagsBamToSam(TTarget & target, TSource const & source)
{
    // Handle case of empty source sequence.
    if (empty(source))
        clear(target);

    clear(target);

    typedef typename Iterator<TSource const, Rooted>::Type TSourceIter;
    TSourceIter it = begin(source, Rooted());

    bool first = true;
    while (!atEnd(it))
    {
        if (!first)
            appendValue(target, '\t');
        first = false;
        _assignTagsBamToSamOneTag(target, it);
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_SAM_CONVERSION_H_
