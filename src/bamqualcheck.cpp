#include <seqan/bam_io.h>
#include <seqan/index.h>

#include "CommandLineParser.hpp"
#include "QualityCheck.hpp"
#include "OverallNumbers.hpp"
#include "ReadQualityHasher.hpp"
#include "TripletCounting.hpp"

// -----------------------------------------------------------------------------
// STRUCT Counts
// -----------------------------------------------------------------------------

struct Counts
{
    OverallNumbers all;
    QualityCheck r1, r2;
    seqan::String<TripletCounts> tripletCounts;
    std::vector<ReadQualityHasher> sps;

    Counts(ProgramOptions & opt)
    {
        r1 = QualityCheck(opt.isize);
        r2 = QualityCheck(opt.isize);
        resize(tripletCounts, 64);

        // Initialize ReadQualityHasher.
        size_t qsize = opt.q_cutoff.size();
        size_t ksize = opt.klist.size();
        sps = std::vector<ReadQualityHasher>(qsize*ksize, ReadQualityHasher(opt));
        for (size_t i = 0; i < qsize; i++)
        {
            for (size_t j = 0; j < ksize; j++)
            {
                sps[i*ksize+j].setQualityCutoff(opt.q_cutoff[i]);
                sps[i*ksize+j].setK(opt.klist[j]);
            }
        }
    }
};

// -----------------------------------------------------------------------------
// FUNCTION getSampleIdAndLaneNames()
// -----------------------------------------------------------------------------

void getSampleIdAndLaneNames(seqan::CharString & id, std::map<seqan::CharString, unsigned> & laneNames, seqan::BamHeader & header)
{
    seqan::CharString lanename;
    for (unsigned i = 0; i < length(header.records); ++i)
    {
        if (header.records[i].type == seqan::BAM_HEADER_READ_GROUP)
        {
            for (unsigned j = 0; j < length(header.records[i].tags); ++j)
            {
                if (seqan::getValueI1(header.records[i].tags[j]) == "ID")
                {
                    lanename = seqan::getValueI2(header.records[i].tags[j]);
                    unsigned l = laneNames.size();
                    laneNames[lanename] = l;
                }
                if (seqan::getValueI1(header.records[i].tags[j]) == "SM")
                {
                    id = seqan::getValueI2(header.records[i].tags[j]);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// FUNCTION getLane()
// -----------------------------------------------------------------------------

unsigned getLane(seqan::BamAlignmentRecord & record,
                 seqan::BamTagsDict & tagsDict,
                 std::map<seqan::CharString, unsigned> & laneNames,
                 seqan::BamIOContext<seqan::StringSet<seqan::CharString> > & context)
{
    for (unsigned tagid = 0; tagid < length(tagsDict); ++tagid)
    {
        if (getTagKey(tagsDict, tagid) == "RG")
        {
            if (getTagType(tagsDict, tagid) == 'Z')
            {
                seqan::CharString readlane;
                readlane = getTagValue(tagsDict, tagid);
                readlane = infix(readlane, 1, length(readlane)-1);
                return laneNames[readlane];
                break;
            }
            else
            {
                std::cout << "Read does not have Z" << "\n";
                if (write2(std::cout, record, context, seqan::Sam()) != 0)
                {
                    std::cerr << "ERROR: Could not write record to stdout \n";
                }
                return -1;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// FUNCTION initChroms()
// -----------------------------------------------------------------------------

std::set<int> initChroms(ProgramOptions & opt, seqan::StringSet<seqan::CharString> & nameStore)
{
    seqan::StringSet<seqan::CharString> chromSet;
    resize(chromSet, 22);
    seqan::strSplit(chromSet, opt.chroms, ',');

    std::set<int> chrIdset;
    int chrId;
    typedef seqan::Iterator<seqan::StringSet<seqan::String<char> >, seqan::Standard>::Type TIterator;
    for(TIterator it = begin(chromSet, seqan::Standard()); it != end(chromSet, seqan::Standard()); ++it)
    {
        if (seqan::getIdByName(nameStore,*it, chrId))
        {
            chrIdset.insert(chrId);
        }
    }
    return chrIdset;
}


// -----------------------------------------------------------------------------
// FUNCTION openFilesAndInit()
// -----------------------------------------------------------------------------

int openFilesAndInit(seqan::Stream<seqan::Bgzf> & inStream,
                      seqan::StringSet<seqan::CharString> & nameStore,
                      seqan::BamIOContext<seqan::StringSet<seqan::CharString> > & context,
                      seqan::CharString & sampleId,
                      std::map<seqan::CharString, unsigned> & laneNames,
                      Genome & genome,
                      std::set<int> & chrIdset,
                      std::ofstream & outFile,
                      ProgramOptions & opt)
{
    // Typedefs for name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

    // Open BGZF Stream for reading.
    if (!open(inStream, toCString(opt.bamFile), "r")) //check if toCstring is needed
    {
        std::cerr << "ERROR: Could not open " << opt.bamFile << " for reading.\n";
        return 1;
    }

    outFile.open(toCString(opt.outputFile), std::ios::out | std::ios::binary);
    if (!outFile.good())
    {
        std::cerr << "ERROR: Could not open output file " << opt.outputFile << '\n';
        return 1;
    }

    TNameStoreCache nameStoreCache(nameStore);
    context = TBamIOContext(nameStore, nameStoreCache);

    seqan::BamHeader header;
    if (readRecord(header, context, inStream, seqan::Bam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from BAM file " << opt.bamFile << "\n";
        return 1;
    }

    // Initialize sample ID and lane names (read groups).
    getSampleIdAndLaneNames(sampleId, laneNames, header);
    seqan::clear(header);

    // Initialize reference genome and set of main chromosomes.
    genome.filename = opt.referenceFile;
    openFastaFile(genome);
    chrIdset = initChroms(opt, nameStore);

    return 0;
}


// -----------------------------------------------------------------------------
// FUNCTION printString()
// -----------------------------------------------------------------------------

template <typename TString>
void printString(TString & str, std::ofstream & outFile)
{
   typedef typename seqan::Iterator<TString>::Type TIterator;
   for(TIterator it = begin(str); it != end(str); ++it)
   {
       outFile << " " << *it;
   }
   outFile << std::endl;
}

// =============================================================================
// FUNCITON writeOutput()
// =============================================================================

void writeOutput(std::ofstream & outFile,
                 seqan::CharString & sampleId,
                 std::map<seqan::CharString, unsigned> & laneNames,
                 seqan::String<Counts> & counts,
                 std::vector<size_t> q_cutoff,
                 std::vector<int> klist)
{
    for (std::map<seqan::CharString, unsigned>::iterator it = laneNames.begin(); it != laneNames.end(); ++it)
    {
        outFile << "sample_id " << sampleId << std::endl;
        outFile << "lane " << it->first << std::endl;
        int lid = it->second;

        counts[lid].r1.avgQualPerPos();
        counts[lid].r2.avgQualPerPos();

        outFile << "total_read_pairs " << counts[lid].all.readcount/2 << std::endl;
        outFile << "total_bps " << counts[lid].all.totalbps << std::endl;
        outFile << "supplementary_alignments " << counts[lid].all.supplementary << std::endl;
        outFile << "marked_duplicate " << counts[lid].all.duplicates << std::endl;
        outFile << "QC_failed " << counts[lid].all.QCfailed << std::endl;
        outFile << "not_primary_alignment " << counts[lid].all.not_primary_alignment << std::endl;
        outFile << "both_reads_unmapped " << counts[lid].all.bothunmapped << std::endl;
        outFile << "first_read_unmapped " << counts[lid].all.firstunmapped << std::endl;
        outFile << "second_read_unmapped " << counts[lid].all.secondunmapped << std::endl;
        outFile << "discordant_read_pairs " << counts[lid].all.discordant << std::endl;
        outFile << "genome_coverage_histogram"; printString(counts[lid].all.poscov, outFile);
        outFile << "insert_size_histogram"; printString(counts[lid].r1.insertSize, outFile);
        outFile << "read_length_histogram_first"; printString(counts[lid].r1.readLength, outFile);
        outFile << "read_length_histogram_second"; printString(counts[lid].r2.readLength, outFile);
        outFile << "N_count_histogram_first"; printString(counts[lid].r1.Ncount, outFile);
        outFile << "N_count_histogram_second"; printString(counts[lid].r2.Ncount, outFile);
        outFile << "GC_content_histogram_first"; printString(counts[lid].r1.GCcount, outFile);
        outFile << "GC_content_histogram_second"; printString(counts[lid].r2.GCcount, outFile);
        outFile << "average_base_qual_histogram_first"; printString(counts[lid].r1.averageQual, outFile);
        outFile << "average_base_qual_histogram_second"; printString(counts[lid].r2.averageQual, outFile);
        outFile << "mapping_qual_histogram_first"; printString(counts[lid].r1.mapQ, outFile);
        outFile << "mapping_qual_histogram_second"; printString(counts[lid].r2.mapQ, outFile);
        outFile << "mismatch_count_histogram_first"; printString(counts[lid].r1.mismatch, outFile);
        outFile << "mismatch_count_histogram_second"; printString(counts[lid].r2.mismatch, outFile);
        outFile << "Ns_by_position_first"; counts[lid].r1.get_DNA_by_position(4, outFile);
        outFile << "Ns_by_position_second"; counts[lid].r2.get_DNA_by_position(4, outFile);
        outFile << "As_by_position_first"; counts[lid].r1.get_DNA_by_position(0, outFile);
        outFile << "As_by_position_second"; counts[lid].r2.get_DNA_by_position(0, outFile);
        outFile << "Cs_by_position_first"; counts[lid].r1.get_DNA_by_position(1, outFile);
        outFile << "Cs_by_position_second"; counts[lid].r2.get_DNA_by_position(1, outFile);
        outFile << "Gs_by_position_first"; counts[lid].r1.get_DNA_by_position(2, outFile);
        outFile << "Gs_by_position_second"; counts[lid].r2.get_DNA_by_position(2, outFile);
        outFile << "Ts_by_position_first"; counts[lid].r1.get_DNA_by_position(3, outFile);
        outFile << "Ts_by_position_second"; counts[lid].r2.get_DNA_by_position(3, outFile);
        outFile << "average_base_qual_by_position_first"; printString(counts[lid].r1.avgqualcount, outFile);
        outFile << "average_base_qual_by_position_second"; printString(counts[lid].r2.avgqualcount, outFile);
        outFile << "soft_clipping_5_prime_by_position_first"; printString(counts[lid].r1.scposcount_5prime, outFile);
        outFile << "soft_clipping_3_prime_by_position_first"; printString(counts[lid].r1.scposcount_3prime, outFile);
        outFile << "soft_clipping_5_prime_by_position_second"; printString(counts[lid].r2.scposcount_5prime, outFile);
        outFile << "soft_clipping_3_prime_by_position_second"; printString(counts[lid].r2.scposcount_3prime, outFile);
        counts[lid].all.ten_most_abundant_kmers(outFile);
        outFile << "8mer_count"; printString(counts[lid].all.eightmercount, outFile);
        for (size_t i = 0; i < q_cutoff.size(); i++)
        {
            for (size_t j = 0; j < klist.size(); j++)
            {
                outFile << klist[j] << "mer_count_after_qual_clipping_"<< q_cutoff[i] << " " << counts[lid].sps[i*klist.size()+j].get_sumCount() << std::endl;
                outFile << "distinct_"<< klist[j] << "mer_count_after_qual_clipping_" << q_cutoff[i] << " " << counts[lid].sps[i*klist.size()+j].F0() << std::endl;
                outFile << "unique_"<< klist[j] << "mer_count_after_qual_clipping_" << q_cutoff[i] << " " << counts[lid].sps[i*klist.size()+j].f1() << std::endl;
                outFile << klist[j] << "mer_F2_after_qual_clipping_" << q_cutoff[i] << " " << counts[lid].sps[i*klist.size()+j].F2() << std::endl;
            }
        }
        writeTripletCounts(outFile, counts[lid].tripletCounts);
    }
}

// =============================================================================
// *** MAIN ***
// =============================================================================

int main(int argc, char const ** argv)
{
    // Typedefs for name store and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

    ProgramOptions opt;
    TripletCountingOptions tripletCountingOptions;

    seqan::ArgumentParser::ParseResult res = parseCommandLine(opt, argc, argv); // Parse the command line.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Open files and initialize context.

    seqan::Stream<seqan::Bgzf> inStream;
    TNameStore nameStore;
    TBamIOContext context;

    seqan::CharString sampleId;
    std::map<seqan::CharString, unsigned> laneNames;

    Genome genome;
    std::set<int> chrIdset;

    std::ofstream outFile;

    if (openFilesAndInit(inStream, nameStore, context, sampleId, laneNames, genome, chrIdset, outFile, opt) == 1)
        return 1;

    // Initialize all counts for each lane.
    seqan::String<Counts> counts;
    unsigned lanecount = laneNames.size();
    resize(counts, lanecount, Counts(opt));

    // Iterate the input BAM file.
    seqan::BamAlignmentRecord record;
    while (!atEnd(inStream))
    {
        // Read record from BAM file.
        if (readRecord(record, context, inStream, seqan::Bam()) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM File " << opt.bamFile << "\n";
            return 1;
        }

        // Get lane of the record.
        seqan::BamTagsDict tagsDict(record.tags);
        unsigned l = getLane(record, tagsDict, laneNames, context);
        if (l == -1) return 1;

        // Check the four highest flags and continue in some cases.
        if (hasFlagSupplementary(record))
        {
            counts[l].all.supplementary += 1;
            continue;
        }
        if (hasFlagSecondary(record))
        {
            counts[l].all.not_primary_alignment += 1;
            continue;
        }
        if (hasFlagDuplicate(record))
        {
            counts[l].all.duplicates += 1;
        }
        if (hasFlagQCNoPass(record))
        {
            counts[l].all.QCfailed += 1;
        }

        // Triplet counting.
        if (tripletCounting(counts[l].tripletCounts, record, nameStore, genome, tripletCountingOptions) != 0)
            return 1;

        // Check if read is in reverse complement.
        if (hasFlagRC(record))
        {
            reverseComplement(record.seq);
            reverse(record.qual);
            reverse(record.cigar);
        }

        // Counts over all chromosomes.
        counts[l].all.readcount +=1;
        counts[l].all.totalbps += length(record.seq);
        if (seqan::hasFlagFirst(record))
        {
            counts[l].r1.check_read_len(record.seq, record.qual); // stores sequence as String of Dna5 and checks if length of seq and qual is the same
            counts[l].r1.get_count(record.seq, record.qual); // resize_string, read_count, read_length
            if (seqan::hasFlagUnmapped(record))
            {
                counts[l].all.firstunmapped += 1; // counts number of unmapped reads
                if (seqan::hasFlagNextUnmapped(record))
                {
                    counts[l].all.bothunmapped += 1; // counts number of reads where both first and second are unmapped
                }
            }
        }
        else if (seqan::hasFlagLast(record))
        {
            counts[l].r2.check_read_len(record.seq, record.qual);
            counts[l].r2.get_count(record.seq, record.qual); // resize_string, read_count, read_length
            if (seqan::hasFlagUnmapped(record))
            {
                counts[l].all.secondunmapped += 1;
            }
        }
        else
        {
            std::cerr << "ERROR: No first or second flag in read in:  " << opt.bamFile << "\n";
            return 1;
        }

        // Counts only over the specified main chromosomes.
        if (chrIdset.count(record.rID) != 0)
        {
            if (seqan::hasFlagFirst(record))
            {
                if (!seqan::hasFlagUnmapped(record))
                {
                    counts[l].r1.cigar_count(record); //soft-clipping per read positions and histogram: soft-clipping at beginning and end seperately
                    counts[l].r1.map_Q(record.mapQ); // histogram: mapping quality
                    counts[l].r1.mis_match(tagsDict); // histogram: mismatches
                    if (!seqan::hasFlagNextUnmapped(record))
                    {
                        if (chrIdset.count(record.rNextId) != 0)
                        {
                            counts[l].r1.insert_size(record.tLen);
                        }
                    }
                }
            }
            else if (seqan::hasFlagLast(record))
            {
                if (!seqan::hasFlagAllProper(record) && (!seqan::hasFlagUnmapped(record) || !seqan::hasFlagNextUnmapped(record)))
                {
                    counts[l].all.discordant += 1;
                }
                if (!seqan::hasFlagUnmapped(record))
                {
                    counts[l].r2.cigar_count(record);
                    counts[l].r2.map_Q(record.mapQ);
                    counts[l].r2.mis_match(tagsDict);
                }
            }
            if (!seqan::hasFlagUnmapped(record))
            {
                counts[l].all.coverage(record);
            }
        }

        // Count adapter 8-mers.
        counts[l].all.count8mers(record.seq);

        if (!seqan::hasFlagQCNoPass(record) && !seqan::hasFlagDuplicate(record))
        {
            RunBamStream(counts[l].sps, record.seq, record.qual, opt.q_cutoff.size(), opt.klist.size());
        }
    }

    writeOutput(outFile, sampleId, laneNames, counts, opt.q_cutoff, opt.klist);
    return 0;
}

