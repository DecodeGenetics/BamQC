#ifndef TRIPLET_COUNTING_H_
#define TRIPLET_COUNTING_H_

#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// STRUCT TripletCountingOptions
// --------------------------------------------------------------------------

struct TripletCountingOptions
{
    // Quality criteria
    seqan::CharString flags;
    int minBaseQ;
    char minBaseQAscii;
    __uint32 minMapQ;
    __uint32 maxClipped;
    int minAlignmentScore;

    TripletCountingOptions() :
        flags("1100xxxx000x"), minBaseQ(20), minBaseQAscii(53), minMapQ(60), maxClipped(0), minAlignmentScore(50)
    {}
};

// ==========================================================================
// STRUCT Counts
// ==========================================================================

struct TripletCounts
{
    size_t forwardFirst [4];
    size_t forwardSecond [4];
    size_t reverseFirst [4];
    size_t reverseSecond [4];

    TripletCounts()
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            forwardFirst[i] = 0;
            forwardSecond[i] = 0;
            reverseFirst[i] = 0;
            reverseSecond[i] = 0;
        }
    }
};

// --------------------------------------------------------------------------

inline int
contextToIndex(seqan::DnaString & context)
{
    return ((int)context[0] << 4) + ((int)context[1] << 2) + (int)context[2];
}

// --------------------------------------------------------------------------
// STRUCT Genome
// --------------------------------------------------------------------------

struct Genome
{
    seqan::CharString filename;
    seqan::SequenceStream stream;
    seqan::CharString chromName;
    seqan::Dna5String chrom;
    unsigned count = 0;
};

// --------------------------------------------------------------------------

bool
openFastaFile(Genome & genome)
{
    open(genome.stream, toCString(genome.filename), seqan::SequenceStream::READ, seqan::SequenceStream::FASTA);
    if (!isGood(genome.stream))
    {
        std::cerr << "ERROR: Could not open fasta file " << genome.filename << std::endl;
        return 1;
    }
    return 0;
}

// --------------------------------------------------------------------------

bool
readFastaRecord(Genome & genome, seqan::CharString recordChrom)
{
    seqan::CharString chrom;
    seqan::CharString chromName;
    if (readRecord(chromName, chrom, genome.stream) != 0)
    {
        std::cout << "ERROR: Could not read fasta record at line: " << genome.count+1 << " from " << genome.filename << std::endl;
        std::cout << "recordChrom: " << recordChrom << std::endl;
        std::cout << "Last chromName was: " << genome.chromName << std::endl;
        return 1;
    }
    genome.count +=1;
    genome.chrom = chrom;
    std::string s = toCString(chromName);
    seqan::CharString chromName_tmp = s.substr(0, s.find(' '));
    std::string ss = toCString(chromName_tmp);
    genome.chromName = ss.substr(0, ss.find('\t'));
    return 0;
}

// ==========================================================================
// FUNCTION alignmentScore()
// ==========================================================================

int
alignmentScore(seqan::BamAlignmentRecord & read)
{
    seqan::BamTagsDict tagsDict(read.tags);

    unsigned tagIdx = 0;
    if (!findTagKey(tagIdx, tagsDict, "AS"))
    {
        std::cerr << "ERROR: Read " << read.qName << " has no AS tag." << std::endl;
        return -1;
    }

    int alignmentScore = 0;
    if (!extractTagValue(alignmentScore, tagsDict, tagIdx))
    {
        std::cerr << "ERROR: Could not read AS tag for read " << read.qName << std::endl;
        return -1;
    }

    return alignmentScore;
}

// ==========================================================================
// FUNCTION checkFlagsAndQuality()
// ==========================================================================

int
checkFlagsAndQuality(seqan::BamAlignmentRecord & read, TripletCountingOptions & options)
{
    typedef seqan::Iterator<seqan::String<seqan::CigarElement<> > >::Type TIter;

    // Check flags.
    if ((hasFlagMultiple(read) && options.flags[0] == '0') || (!hasFlagMultiple(read) && options.flags[0] == '1')) return 0;
    if ((hasFlagAllProper(read) && options.flags[1] == '0') || (!hasFlagAllProper(read) && options.flags[1] == '1')) return 0;
    if (hasFlagUnmapped(read)) return 0;
    if ((hasFlagNextUnmapped(read) && options.flags[3] == '0') || (!hasFlagNextUnmapped(read) && options.flags[3] == '1')) return 0;
    if ((hasFlagRC(read) && options.flags[4] == '0') || (!hasFlagRC(read) && options.flags[4] == '1')) return 0;
    if ((hasFlagNextRC(read) && options.flags[5] == '0') || (!hasFlagNextRC(read) && options.flags[5] == '1')) return 0;
    if ((hasFlagFirst(read) && options.flags[6] == '0') || (!hasFlagFirst(read) && options.flags[6] == '1')) return 0;
    if ((hasFlagLast(read) && options.flags[7] == '0') || (!hasFlagLast(read) && options.flags[7] == '1')) return 0;
    if ((hasFlagSecondary(read) && options.flags[8] == '0') || (!hasFlagSecondary(read) && options.flags[8] == '1')) return 0;

    // Check mapping quality and alignment score of the read.
    if (read.mapQ < options.minMapQ) return 0;
    int as = alignmentScore(read);
    if (as < 0) return -1;
    if (as < options.minAlignmentScore) return 0;

    // Check clipping of read.
    __uint32 clipped = 0;
    for (TIter it = begin(read.cigar); it != end(read.cigar); ++it)
    {
        if ((*it).operation == 'S' || (*it).operation == 'H')
            clipped += (*it).count;
    }
    if (clipped > options.maxClipped) return 0;

    return 1;
}

// ==========================================================================
// FUNCTION countPosition()
// ==========================================================================

void
countPosition(TripletCounts & counts, seqan::Dna5 base, seqan::BamAlignmentRecord & read)
{
    // Count the position.
    if (hasFlagRC(read))
    {
        if (hasFlagFirst(read)) counts.reverseFirst[base] += 1;
        else counts.reverseSecond[base] += 1;
    }
    else
    {
        if (hasFlagFirst(read)) counts.forwardFirst[base] += 1;
        else counts.forwardSecond[base] += 1;
    }
    
}

// ==========================================================================
// FUNCTION countBasesInTriplets()
// ==========================================================================

void
countBasesInTriplets(seqan::String<TripletCounts> & counts, seqan::BamAlignmentRecord & read, seqan::Dna5String & chrom, char minBaseQAscii)
{
    typedef seqan::Iterator<seqan::String<seqan::CigarElement<> > >::Type TIter;

    // Initialize iterator over CIGAR string.
    TIter it = begin(read.cigar);
    size_t cigarCount = (*it).count - 1;
    
    // Iterate over the read.
    size_t chromPos = read.beginPos + 1;
    for (size_t readPos = 1; readPos < length(read.seq) - 1; ++readPos, ++chromPos, --cigarCount)
    {
        // Shift positions according to CIGAR string.
        while (cigarCount == 0)
        {
            ++it;
            SEQAN_ASSERT_NEQ(it, end(read.cigar));
            if ((*it).operation == 'D' || (*it).operation == 'N' || (*it).operation == 'H' || (*it).operation == 'P')
                chromPos += (*it).count;
            else if ((*it).operation == 'S' || (*it).operation == 'H' || (*it).operation == 'I')
                readPos += (*it).count;
            else
                cigarCount = (*it).count;
        }
        if (readPos >= length(read.seq) - 1) break;

        // Check base calling quality.
        if (read.qual[readPos] < minBaseQAscii) continue;

        // Determine context and base.
        seqan::Dna5 base = read.seq[readPos];
        if (base == 'N' || read.seq[readPos-1] == 'N' || read.seq[readPos+1] == 'N') continue;
        seqan::DnaString context = infix(chrom, chromPos-1, chromPos + 2);

        // Check prefix and suffix of context for an exact match.
        if (read.seq[readPos - 1] != context[0]) continue;
        if (read.seq[readPos + 1] != context[2]) continue;

        countPosition(counts[contextToIndex(context)], base, read);
    }
}

// ==========================================================================
// FUNCTION tripletCounting()
// ==========================================================================

template<typename TNameStore>
bool
tripletCounting(seqan::String<TripletCounts> & counts,
                seqan::BamAlignmentRecord & record,
                TNameStore & nameStore,
                Genome & genome,
                TripletCountingOptions & options)
{
    int res = checkFlagsAndQuality(record, options);
    if (res == -1) return 1;
    if (res == 0) return 0;

    seqan::CharString recordChrom = nameStore[record.rID];
    while (recordChrom != genome.chromName)
    {
        // Read the next chromosome of reference genome.
        if (readFastaRecord(genome, recordChrom) != 0) return 1;
    }

    if (recordChrom == genome.chromName)
        countBasesInTriplets(counts, record, genome.chrom, options.minBaseQAscii);

    return 0;
}

// --------------------------------------------------------------------------
// FUNCTION writeTripletCounts()
// --------------------------------------------------------------------------

void
writeTripletCounts(std::ofstream & outfile, seqan::String<TripletCounts> & counts, char base, size_t index)
{
    outfile << "triplet_counts_" << base << "_1st_FW";
    for (unsigned i = 0; i < length(counts); ++i)
        outfile << " " << counts[i].forwardFirst[index];
    outfile << std::endl;
    outfile << "triplet_counts_" << base << "_1st_RC";
    for (unsigned i = 0; i < length(counts); ++i)
        outfile << " " << counts[i].reverseFirst[index];
    outfile << std::endl;
    outfile << "triplet_counts_" << base << "_2nd_FW";
    for (unsigned i = 0; i < length(counts); ++i)
        outfile << " " << counts[i].forwardSecond[index];
    outfile << std::endl;
    outfile << "triplet_counts_" << base << "_2nd_RC";
    for (unsigned i = 0; i < length(counts); ++i)
        outfile << " " << counts[i].reverseSecond[index];
    outfile << std::endl;
}

// --------------------------------------------------------------------------

void
writeTripletCounts(std::ofstream & outfile, seqan::String<TripletCounts> & counts)
{
    char bases[] = {'A', 'C', 'G', 'T'};
    
    for (size_t i = 0; i < 4; ++i)
        writeTripletCounts(outfile, counts, bases[i], i);
}

// --------------------------------------------------------------------------

#endif  // TRIPLET_COUNTING_H_
