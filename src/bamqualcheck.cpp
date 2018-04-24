#include <seqan/bam_io.h>
#include <seqan/index.h>

#include "CommandLineParser.hpp"
#include "QualityCheck.hpp"
#include "OverallNumbers.hpp"
#include "ReadQualityHasher.hpp"
#include "TripletCounting.hpp"

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
                 seqan::String<OverallNumbers> & rall,
                 seqan::String<QualityCheck> & r1,
                 seqan::String<QualityCheck> & r2,
                 ProgramOptions & opt,
                 std::vector<std::vector<ReadQualityHasher> > & sps,
                 size_t ksize, size_t qsize,
                 seqan::String<TripletCounts> & counts)
{
    for (std::map<seqan::CharString, unsigned>::iterator it=laneNames.begin(); it!=laneNames.end(); ++it)
    {
        outFile << "sample_id " << sampleId << std::endl;
        outFile << "lane " << it->first << std::endl;
        int lid = it->second;

        r1[lid].avgQualPerPos();
        r2[lid].avgQualPerPos();

        outFile << "total_read_pairs " << rall[lid].readcount/2 << std::endl;
        outFile << "total_bps " << rall[lid].totalbps << std::endl;
        outFile << "supplementary_alignments " << rall[lid].supplementary << std::endl;
        outFile << "marked_duplicate " << rall[lid].duplicates << std::endl;
        outFile << "QC_failed " << rall[lid].QCfailed << std::endl;
        outFile << "not_primary_alignment " << rall[lid].not_primary_alignment << std::endl;
        outFile << "both_reads_unmapped " << rall[lid].bothunmapped << std::endl;
        outFile << "first_read_unmapped " << rall[lid].firstunmapped << std::endl;
        outFile << "second_read_unmapped " << rall[lid].secondunmapped << std::endl;
        outFile << "discordant_read_pairs " << rall[lid].discordant << std::endl;
        outFile << "genome_coverage_histogram"; printString(rall[lid].poscov, outFile);
        outFile << "insert_size_histogram"; printString(r1[lid].insertSize, outFile);
        outFile << "read_length_histogram_first"; printString(r1[lid].readLength, outFile);
        outFile << "read_length_histogram_second"; printString(r2[lid].readLength, outFile);
        outFile << "N_count_histogram_first"; printString(r1[lid].Ncount, outFile);
        outFile << "N_count_histogram_second"; printString(r2[lid].Ncount, outFile);
        outFile << "GC_content_histogram_first"; printString(r1[lid].GCcount, outFile);
        outFile << "GC_content_histogram_second"; printString(r2[lid].GCcount, outFile);
        outFile << "average_base_qual_histogram_first"; printString(r1[lid].averageQual, outFile);
        outFile << "average_base_qual_histogram_second"; printString(r2[lid].averageQual, outFile);
        outFile << "mapping_qual_histogram_first"; printString(r1[lid].mapQ, outFile);
        outFile << "mapping_qual_histogram_second"; printString(r2[lid].mapQ, outFile);
        outFile << "mismatch_count_histogram_first"; printString(r1[lid].mismatch, outFile);
        outFile << "mismatch_count_histogram_second"; printString(r2[lid].mismatch, outFile);
        outFile << "Ns_by_position_first"; r1[lid].get_DNA_by_position(4, outFile);
        outFile << "Ns_by_position_second"; r2[lid].get_DNA_by_position(4, outFile);
        outFile << "As_by_position_first"; r1[lid].get_DNA_by_position(0, outFile);
        outFile << "As_by_position_second"; r2[lid].get_DNA_by_position(0, outFile);
        outFile << "Cs_by_position_first"; r1[lid].get_DNA_by_position(1, outFile);
        outFile << "Cs_by_position_second"; r2[lid].get_DNA_by_position(1, outFile);
        outFile << "Gs_by_position_first"; r1[lid].get_DNA_by_position(2, outFile);
        outFile << "Gs_by_position_second"; r2[lid].get_DNA_by_position(2, outFile);
        outFile << "Ts_by_position_first"; r1[lid].get_DNA_by_position(3, outFile);
        outFile << "Ts_by_position_second"; r2[lid].get_DNA_by_position(3, outFile);
        outFile << "average_base_qual_by_position_first"; printString(r1[lid].avgqualcount, outFile);
        outFile << "average_base_qual_by_position_second"; printString(r2[lid].avgqualcount, outFile);
        outFile << "soft_clipping_5_prime_by_position_first"; printString(r1[lid].scposcount_5prime, outFile);
        outFile << "soft_clipping_3_prime_by_position_first"; printString(r1[lid].scposcount_3prime, outFile);
        outFile << "soft_clipping_5_prime_by_position_second"; printString(r2[lid].scposcount_5prime, outFile);
        outFile << "soft_clipping_3_prime_by_position_second"; printString(r2[lid].scposcount_3prime, outFile);
        rall[lid].ten_most_abundant_kmers(outFile);
        outFile << "8mer_count"; printString(rall[lid].adapterKmercount, outFile);
        for (size_t i = 0; i < qsize; i++) {
            for (size_t j = 0; j < ksize; j++) {
                outFile << opt.klist[j] << "mer_count_after_qual_clipping_"<< opt.q_cutoff[i] << " " << sps[lid][i*ksize+j].get_sumCount() << std::endl;
                outFile << "distinct_"<< opt.klist[j] << "mer_count_after_qual_clipping_" << opt.q_cutoff[i] << " " << sps[lid][i*ksize+j].F0() << std::endl;
                outFile << "unique_"<< opt.klist[j] << "mer_count_after_qual_clipping_" << opt.q_cutoff[i] << " " << sps[lid][i*ksize+j].f1() << std::endl;
                outFile << opt.klist[j] << "mer_F2_after_qual_clipping_" << opt.q_cutoff[i] << " " << sps[lid][i*ksize+j].F2() << std::endl;
            }
        }
    }
    writeTripletCounts(outFile, counts);
}

// =============================================================================
// *** MAIN ***
// =============================================================================

int main(int argc, char const ** argv)
{
    // Typedefs for name store, cache, and BAM I/O context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

    ProgramOptions opt;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(opt, argc, argv); // Parse the command line.
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }

    seqan::Stream<seqan::Bgzf> inStream; // Open BGZF Stream for reading.
    if (!open(inStream, toCString(opt.bamFile), "r")) //check if toCstring is needed
    {
        std::cerr << "ERROR: Could not open " << opt.bamFile << " for reading.\n";
        return 1;
    }
    std::ofstream outFile(toCString(opt.outputFile), std::ios::out | std::ios::binary);
    if (!outFile.good())
    {
        std::cerr << "ERROR: Could not open output file " << opt.outputFile << '\n';
        return 1;
    }

    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);

    seqan::BamHeader header;
    if (readRecord(header, context, inStream, seqan::Bam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from BAM file " << opt.bamFile << "\n";
        return 1;
    }

    // set up chromset
    seqan::StringSet<seqan::String<char> > chromSet;
    seqan::String<char> seq = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22";
    resize(chromSet, 22);
    seqan::strSplit(chromSet, seq);
    int chrId;
    std::set<int> chrIdset;
    typedef seqan::Iterator<seqan::StringSet<seqan::String<char> >, seqan::Standard>::Type TIterator;
    for(TIterator it = begin(chromSet, seqan::Standard()); it != end(chromSet, seqan::Standard()); ++it)
    {
        if (seqan::getIdByName(nameStore,*it, chrId))
        {
            chrIdset.insert(chrId);
        }
    }

    // Initialize lane names (read groups).
    seqan::CharString sampleId;
    std::map<seqan::CharString, unsigned> laneNames;
    getSampleIdAndLaneNames(sampleId, laneNames, header);
    unsigned lanecount = laneNames.size();
    seqan::clear(header);

    // Initialize counts over all reads.
    seqan::String<OverallNumbers> rall;
    resize(rall, lanecount);

    // Initialize counts for first and second reads in pair.
    seqan::String<QualityCheck> r1, r2;
    resize(r1, lanecount, QualityCheck(opt.isize));
    resize(r2, lanecount, QualityCheck(opt.isize));

    // Init ReadQualityHasher
    size_t qsize = opt.q_cutoff.size();
    size_t ksize = opt.klist.size();
    std::vector<std::vector<ReadQualityHasher> > sps(lanecount, std::vector<ReadQualityHasher> (qsize*ksize, ReadQualityHasher(opt)));
    for (unsigned k = 0;  k < lanecount; k++)
    {
        for (size_t i = 0; i < qsize; i++)
        {
            for (size_t j = 0; j < ksize; j++)
            {
              sps[k][i*ksize+j].setQualityCutoff(opt.q_cutoff[i]);
              sps[k][i*ksize+j].setK(opt.klist[j]);
            }
        }
    }

    // Init triplet counting
    seqan::String<TripletCounts> counts;
    resize(counts, 64);
    TripletCountingOptions tripletCountingOptions;

    // Reference genome.
    Genome genome;
    genome.filename = opt.referenceFile;
    openFastaFile(genome);

    // Read record
    seqan::BamAlignmentRecord record;
    while (!atEnd(inStream))
    {
        //READ RECORD
        if (readRecord(record, context, inStream, seqan::Bam()) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM File " << opt.bamFile << "\n";
            return 1;
        }
        seqan::BamTagsDict tagsDict(record.tags);
        unsigned l = getLane(record, tagsDict, laneNames, context);
        if (l == -1) return 1;

        // Check the four highest flags and continue in some cases.
        if (hasFlagSupplementary(record)) {
            rall[l].supplementary +=1;
            continue;
        }
        if (hasFlagSecondary(record)){
            rall[l].not_primary_alignment +=1;
            continue;
        }
        if (hasFlagDuplicate(record)){
            rall[l].duplicates +=1;
        }
        if (hasFlagQCNoPass(record)){
            rall[l].QCfailed +=1;
        }

        // Triplet counting.
        if (tripletCounting(counts, record, nameStore, genome, tripletCountingOptions) != 0) {
            return 1;
        }

        if (hasFlagRC(record)) // check if read is in reverse complement
        {
            reverseComplement(record.seq);
            reverse(record.qual);
            reverse(record.cigar);
        }
        // All chromosomes
        rall[l].readcount +=1;
        rall[l].totalbps += length(record.seq);
        if (seqan::hasFlagFirst(record))
        {
            r1[l].check_read_len(record.seq, record.qual); // stores dequence as String of Dna5 and checks if length of seq and qual is the same
            r1[l].get_count(record.seq, record.qual); // resize_string, read_count, read_length
            if (seqan::hasFlagUnmapped(record))
            {
                rall[l].firstunmapped +=1; // counts number of unmapped reads
                if (seqan::hasFlagNextUnmapped(record)) {
                    rall[l].bothunmapped +=1; // counts number of reads where both first and second are unmapped
                }
            }
        }
        else if (seqan::hasFlagLast(record))
        {
            r2[l].check_read_len(record.seq, record.qual);
            r2[l].get_count(record.seq, record.qual); // resize_string, read_count, read_length
            if (seqan::hasFlagUnmapped(record)){
                rall[l].secondunmapped +=1; //
            }
        }
        else
        {
            std::cerr << "ERROR: No first or second flag in read in:  " << opt.bamFile << "\n";
            return 1;
        }

        // only chromosome chr1, chr2, ... chr22
        if (chrIdset.count(record.rID)!=0)
        {
            if (seqan::hasFlagFirst(record))
            {
                if (!seqan::hasFlagUnmapped(record))
                {
                    r1[l].cigar_count(record); //soft-clipping per read positions and histogram: soft-clipping at beginning and end seperately
                    r1[l].map_Q(record.mapQ); // histogram: mapping quality
                    r1[l].mis_match(tagsDict); // histogram: mismatches
                    if (!seqan::hasFlagNextUnmapped(record))
                    {
                        if (chrIdset.count(record.rNextId)!=0)
                        {
                            r1[l].insert_size(record.tLen);
                        }
                    }
                }
            }
            else if (seqan::hasFlagLast(record))
            {
                if (!seqan::hasFlagAllProper(record) && (!seqan::hasFlagUnmapped(record) || !seqan::hasFlagNextUnmapped(record))) {
                    rall[l].discordant +=1;
                }
                if (!seqan::hasFlagUnmapped(record)){
                    r2[l].cigar_count(record);
                    r2[l].map_Q(record.mapQ);
                    r2[l].mis_match(tagsDict);
                }
            }
            if (!seqan::hasFlagUnmapped(record))
            {
                rall[l].coverage(record);
            }
        }

        rall[l].countAdapterKmers(record.seq); //count adapter 8-mers

        if (!seqan::hasFlagQCNoPass(record) && !seqan::hasFlagDuplicate(record)) {
            RunBamStream(sps[l], record.seq, record.qual, qsize, ksize);
        }
    }

    writeOutput(outFile, sampleId, laneNames, rall, r1, r2, opt, sps, ksize, qsize, counts);
    return 0;
}

