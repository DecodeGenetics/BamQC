#ifndef QUALITY_CHECK_H_
#define QUALITY_CHECK_H_

// ==========================================================================
// CLASS QualityCheck
// ==========================================================================

class QualityCheck
{
public:
    // Count per read position
    seqan::String<double> avgqualcount;
    seqan::String<unsigned> scposcount_5prime;
    seqan::String<unsigned> scposcount_3prime;
    seqan::String<seqan::String<uint64_t> > dnacount;

    // Count per read -> histogram
    seqan::String<unsigned> averageQual;
    seqan::String<unsigned> Ncount;
    seqan::String<uint64_t> GCcount;
    seqan::String<unsigned> insertSize;
    seqan::String<unsigned> mapQ;
    seqan::String<unsigned> readLength;
    seqan::String<unsigned> mismatch;

    // Constructors
    QualityCheck();
    QualityCheck(int isize);

    // Functions
    void avgQualPerPos();
    int check_read_len(seqan::CharString & seq, seqan::CharString & qual);
    void get_count(seqan::CharString & rec, seqan::CharString & qual);
    void cigar_count(seqan::BamAlignmentRecord & record);
    void map_Q(__uint8 & mapq);
    void insert_size(int & tlen);
    void mis_match(seqan::BamTagsDict & tagsDict);

private:
    // Count per read position
    seqan::String<uint64_t> qualcount;

    // Count per read -> histogram
    unsigned DIcount;
    unsigned qualcount_readnr;

    // Functions
    void resize_strings(seqan::CharString & dnaseq);
    void read_counts(seqan::CharString & record, seqan::CharString & qual);
    void read_length(seqan::CharString & record);
};

QualityCheck::QualityCheck(): qualcount_readnr(0)
{
    resize(dnacount, 5);
    resize(insertSize, 1001, 0);
}

QualityCheck::QualityCheck(int isize): qualcount_readnr(0)
{
    resize(dnacount, 5);
    resize(insertSize, isize + 1, 0);
}

// -----------------------------------------------------------------------------
// FUNCTION check_read_len()
// -----------------------------------------------------------------------------

int QualityCheck::check_read_len(seqan::CharString & seq, seqan::CharString & qual)
{
    // COUNTS PER READ POSITION
    if (length(seq) != length(qual))
    {
        std::cerr << "ERROR: length of sequence and quality is not the same" << "\n";
        return 1;
    }
    return 0;
}

// -----------------------------------------------------------------------------
// FUNCTION resize_strings()
// -----------------------------------------------------------------------------

void QualityCheck::resize_strings(seqan::CharString & dnaseq)
{
    // count for nucleotides and qualtype for each read position:
    // Increase number of counters if dnaseq is longer than the previous reads.
    if (length(qualcount) < length(dnaseq)) //This should only be true once...
    {
        unsigned oldSize = length(qualcount);

        resize(qualcount, length(dnaseq), 0);
        resize(avgqualcount, length(dnaseq), 0);
        resize(Ncount, length(dnaseq)+1,0);
        resize(GCcount, length(dnaseq)+1,0);
        resize(scposcount_5prime, length(dnaseq), 0); //check if n is only assigned to added size.
        resize(scposcount_3prime, length(dnaseq), 0); //check if n is only assigned to added size.

        for (unsigned j = 0; j < 5; ++j)
        {
            resize(dnacount[j], length(dnaseq), 0);
        }
    }
}

// -----------------------------------------------------------------------------
// FUNCTION get_count()
// -----------------------------------------------------------------------------

void QualityCheck::get_count(seqan::CharString & seq, seqan::CharString & qual)
{
    resize_strings(seq);
    read_counts(seq, qual);
    read_length(seq);
}

// -----------------------------------------------------------------------------
// FUNCTION read_counts()
// -----------------------------------------------------------------------------

void QualityCheck::read_counts(seqan::CharString & seq, seqan::CharString & qual)
{
    unsigned cntN = 0;
    unsigned cntGC = 0;
    unsigned avgQual = 0;

    qualcount_readnr += 1;

    unsigned j = 0;
    seqan::Iterator<seqan::CharString, seqan::Rooted>::Type itseq = begin(seq);
    seqan::Iterator<seqan::CharString, seqan::Rooted>::Type itEndseq = end(seq);
    for (; itseq != itEndseq; goNext(itseq))
    {
        dnacount[(int) seqan::ordValue((seqan::Dna5)*itseq)][j] += 1;

        if (*itseq == 'N') // count number of Ns per read
        {
            cntN += 1;
        }
        if ((*itseq == 'C') || (*itseq == 'G')) // count number of GCs per read
        {
            cntGC += 1;
        }
        ++j;
    }
    j = 0;
    seqan::Iterator<seqan::CharString, seqan::Rooted>::Type itqual = begin(qual);
    seqan::Iterator<seqan::CharString, seqan::Rooted>::Type itEndqual = end(qual);
    for (; itqual != itEndqual; goNext(itqual))
    {
        qualcount[j] += (int) seqan::ordValue(*itqual)-33; // sums up the quality value for all reads per postion
        avgQual += (int) seqan::ordValue(*itqual)-33; // sums up all quality values for the read
        ++j;
    }

    Ncount[cntN] += 1; // count number of Ns per read

    GCcount[cntGC] += 1;

    if (length(averageQual) <= ceil(double(avgQual) / length(seq)))
    {
        resize(averageQual, ceil(double(avgQual) / length(seq)) + 1, 0);
    }
    averageQual[(int)round(double(avgQual) / length(seq))] += 1; // histogram: frequency of reads with x-average quality value.
}

void QualityCheck::read_length(seqan::CharString & seq)
{
    unsigned lseq = length(seq);
    if (length(readLength) <= lseq)
    {
        resize(readLength, lseq + 1, 0);
    }
    readLength[lseq] += 1; // histogram of read lengths
}

void QualityCheck::map_Q(__uint8 & mapq)
{
    if (length(mapQ) <= mapq)
    {
        resize(mapQ, mapq + 1, 0);
    }
    mapQ[mapq] += 1; // histogram of mapping Qualities
}

void QualityCheck::insert_size(int & tlen)
{
    unsigned index = abs(tlen);

    if (index >= length(insertSize))
    {
        index = length(insertSize);
    }
    insertSize[index] += 1; // histogram of insert size
}

void QualityCheck::mis_match(seqan::BamTagsDict & tagsDict)
{

    for (unsigned tagid = 0; tagid < length(tagsDict); ++tagid)
    {
        if (getTagKey(tagsDict, tagid) == "NM")
        {
            char tagType = getTagType(tagsDict, tagid);
            if (tagType == 'c' || tagType == 'C' || tagType == 'i' || tagType == 'I' || tagType == 's' || tagType == 'S')
            {
                unsigned x = 0;
                extractTagValue(x, tagsDict, tagid); // x gives number of mismatches + deletions + insertions
                unsigned mmcount = x - DIcount; // x - deletions -insertions gives number of mismatches

                if (length(mismatch) <= mmcount)
                {
                    resize(mismatch, mmcount + 1, 0);
                }
                mismatch[mmcount] += 1; // histogram of number of mismatches
            }
        }
    }
}

void QualityCheck::cigar_count(seqan::BamAlignmentRecord & record)
{
    // function
    int c_begin = 0;                        // begin position of cigar in seq
    int c_end = 0;                          // end position of cigar in seq
    unsigned Scount_begin = 0;                   // counting soft-clipping at the beginning
    unsigned Scount_end = 0;                     // counting soft-clipping at the end
    int cigarlength = length(record.cigar); // cigar is a paired string with 'operation' (ie. S D I M) and 'count'
    DIcount = 0;                            // DIcount is used later in the function mis_match

    if (record.cigar[0].operation == 'S')
    {
        Scount_begin += record.cigar[0].count;
        for (unsigned j = 0; j < record.cigar[0].count; ++j)
        {
            scposcount_5prime[j] +=1; // if Soft-clipping at read position
        }
    }
    else if (record.cigar[cigarlength-1].operation == 'S')
    {
        Scount_end += record.cigar[cigarlength-1].count;
        for (unsigned j = length(record.seq) -record.cigar[cigarlength-1].count; j < length(record.seq); ++j)
        {
            scposcount_3prime[j] +=1; // if Soft-clipping at read position
        }
    }

    for (int i = 0; i < cigarlength; ++i)
    {
        if ((record.cigar[i].operation == 'D') || (record.cigar[i].operation == 'I')) //TODO: check if N,P,=,X matter here
        {
            DIcount += record.cigar[i].count; // counts deletions and insertions for read. Later used in function mis_match
        }

    }
}

void QualityCheck::avgQualPerPos()
{
    for (unsigned i = 0; i < length(qualcount); ++i)
    {
        avgqualcount[i] = (qualcount[i]/(double)qualcount_readnr); // sum of all read qualities per position / number of reads TODO: check if correct
    }
}

#endif  // QUALITY_CHECK_H_
