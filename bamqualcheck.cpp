#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/index.h>

#include "StreamCounter.hpp"
#include "RepHash.hpp"
#include "TripletCounting.hpp"
#include "version.h"

struct ProgramOptions
{
    // Input and output files.
    seqan::CharString bamFile;
    seqan::CharString referenceFile;
    seqan::CharString outputFile;

    // Parameters of kmerstream.
    std::vector<int> klist;
    double e;
    std::vector<size_t> q_cutoff;
    size_t q_base;
    int seed;

    // Maximum insert size to be counted in histogram.
    int isize;

    ProgramOptions() :
        referenceFile("genome.fa"), e(0.01), q_base(33), seed(0)
    {}
};


seqan::ArgumentParser::ParseResult
parseCommandLine(ProgramOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bamqualcheck");

    // Set short description, version, and date.
    setShortDescription(parser, "String Modifier");
    setDate(parser, "April 2018");
    setVersion(parser, version);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAMFILE\\fP");
    addDescription(parser, "Program for bam quality checks. ");

    // Add the required argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE,"bamfile", "False", 1));
    setValidValues(parser, 0,"bam sam");

    // Add general options.
    addSection(parser, "General options");
    addOption(parser, seqan::ArgParseOption("r", "reference", "Reference genome filename.", seqan::ArgParseArgument::STRING, "FILENAME"));
    setDefaultValue(parser, "r", "genome.fa");
    addOption(parser, seqan::ArgParseOption("i", "insert-size", "Upper bound for the insert size in insert size histogram.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "i", "1000");
    addOption(parser, seqan::ArgParseOption("o", "output-file", "Output filename.", seqan::ArgParseArgument::OUTPUTFILE, "OUT"));

    // Add kmerstream options.
    addSection(parser, "Kmerstream options");
    addOption(parser, seqan::ArgParseOption("k", "kmer-size", "Size of k-mers, single value.", seqan::ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "k", "32");
    addOption(parser, seqan::ArgParseOption("q", "quality-cutoff", "Comma separated list, keep k-mers with bases above quality threshold in PHRED.", seqan::ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "q", "17");
    addOption(parser, seqan::ArgParseOption("e", "error-rate", "Error rate guaranteed.", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "e", "0.01");
    setMinValue(parser, "e", "0");
    addOption(parser, seqan::ArgParseOption("s", "seed", "Seed value for the randomness (use time based randomness).", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "s", "1");

    // Set reference and output files as a required options.
    setRequired(parser, "r");
    setRequired(parser, "o");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Get argument and option values.
    getArgumentValue(options.bamFile, parser, 0);
    std::string qcut, kmer;
    getOptionValue(options.referenceFile, parser, "r");
    getOptionValue(qcut, parser, "q");
    getOptionValue(kmer, parser, "k");
    getOptionValue(options.e, parser, "e");
    getOptionValue(options.seed, parser, "s");
    getOptionValue(options.isize, parser, "i");
    getOptionValue(options.outputFile, parser, "o");

    // Parse the list of quality cutoff values.
    std::stringstream sq(qcut);
    size_t j;
    while (sq >> j)
    {
        options.q_cutoff.push_back(j);

        if (sq.peek() == ',')
            sq.ignore();
    }

    // Parse the list of k-mer sizes.
    std::stringstream sk(kmer);
    int i;
    while (sk >> i)
    {
        options.klist.push_back(i);

        if (sk.peek() == ',')
            sk.ignore();
    }

    return seqan::ArgumentParser::PARSE_OK;
}

class ReadQualityHasher {
    public:
        ReadQualityHasher(const ProgramOptions &opt) : hf(), k(0), sc(opt.e, opt.seed), q_cutoff(0), q_base(opt.q_base) {
            if (opt.seed != 0) {
                hf.seed(opt.seed);
            }
        }

        void setQualityCutoff(size_t q) {
            q_cutoff = q;
        }

        void setK(size_t _k) {
            k = _k;
            hf.init(k);
        }

        void operator()(const char* s, size_t l, const char* q, size_t ql) {
            // create hashes for all k-mers
            // operate on hashes
            size_t i=0, j=0;
            bool last_valid = false;
            if (l < k) {
                return;
            }
            while (j < l)
            {
                char c = s[j];
                if (c != 'N' && c != 'n' && (q[j] >= (char) (q_base+q_cutoff))) {
                    if (last_valid) {
                        hf.update(s[i], s[j]);
                        i++;
                        j++;
                    }
                    else {
                        if (i + k -1 == j) {
                            hf.init(s+i); // start the k-mer at position i
                            last_valid = true;
                            j++;
                        }
                        else {
                            j++; // move along
                        }
                    }
                }
                else {
                    // invalid character, restart
                    j++;
                    i = j;
                    last_valid = false;
                }
                if (last_valid) {
                    handle(hf.hash());
                }
            }
        }

        void handle(uint64_t val) {
            sc(val);
        }

        std::string humanReport() {
            return sc.humanReport();
        }

        std::string report() {
            return sc.report();
        }

        size_t F0() {
            return sc.F0();
        }

        size_t f1() {
            return sc.f1();
        }

        size_t get_sumCount() {
            return sc.get_sumCount(); // This function is added to StreamCounter.hpp
        }

        size_t F2() {
            return sc.F2();
        }

        bool join(const ReadQualityHasher& o) {
            return sc.join(o.sc);
        }

    private:
        size_t q_cutoff, q_base;
        RepHash hf;
        size_t k;
        StreamCounter sc;
};


class QualityCheck
{
public:
    QualityCheck();
    QualityCheck(int isize);
    unsigned maxtlen;
    seqan::String<seqan::Dna5> dnaseq;

    //functions
    void avgQualPerPos();
    int check_read_len(seqan::CharString & seq, seqan::CharString & qual);
    void get_count(seqan::CharString & rec, seqan::CharString & qual);
    void cigar_count(seqan::BamAlignmentRecord & record);
    void map_Q(__uint8 & mapq);
    void insert_size(int & tlen);
    void mis_match(seqan::BamTagsDict & tagsDict);
    void clear_var();

    //print functions
    void get_DNA_by_position(int c, std::ofstream & outFile);
    void get_avg_base_qual_by_pos(std::ofstream & outFile);
    void get_soft_clip_by_pos_5prime(std::ofstream & outFile);
    void get_soft_clip_by_pos_3prime(std::ofstream & outFile);
    void get_N_count_histogram(std::ofstream & outFile);
    void get_GC_count_histogram(std::ofstream & outFile);
    void get_soft_clipping_begin_histogram(std::ofstream & outFile);
    void get_soft_clipping_end_histogram(std::ofstream & outFile);
    void get_average_base_qual_histogram(std::ofstream & outFile);
    void get_insert_size_histogram(std::ofstream & outFile);
    void get_mapping_qual_histogram(std::ofstream & outFile);
    void get_mismatch_count_histogram(std::ofstream & outFile);
    void get_read_length_histogram(std::ofstream & outFile);

private:
    seqan::CharString softclipping;

    //Count per read position
    seqan::String<seqan::String<uint64_t> > dnacount;
    seqan::String<uint64_t> qualcount;
    seqan::String<double> avgqualcount;
    seqan::String<unsigned> scposcount_5prime;
    seqan::String<unsigned> scposcount_3prime;

    //Count per read -> histogram
    unsigned DIcount;
    unsigned qualcount_readnr;
    seqan::String<unsigned> Ncount;
    seqan::String<uint64_t> GCcount;
    seqan::String<unsigned> sccount_begin;
    seqan::String<unsigned> sccount_end;
    seqan::String<unsigned> averageQual;
    seqan::String<unsigned> mapQ;
    seqan::String<unsigned> readLength;
    seqan::String<unsigned> insertSize;
    seqan::String<unsigned> mismatch;

    //functions
    void resize_strings();
    void read_counts(seqan::CharString & record, seqan::CharString & qual);
    void read_length(seqan::CharString & record);
};

class OverallNumbers
{
public:
    OverallNumbers();
    unsigned supplementary;
    unsigned duplicates;
    unsigned QCfailed;
    unsigned not_primary_alignment;
    unsigned readcount;
    uint64_t totalbps;
    unsigned bothunmapped;
    unsigned firstunmapped;
    unsigned secondunmapped;
    unsigned discordant;

    void coverage(seqan::BamAlignmentRecord & record);
    void countAdapterKmers(seqan::CharString & seq);
    void get_poscov(std::ofstream & outFile);
    void get_kmer_histogram(std::ofstream & outFile);
    void get_adapterKmercount(std::ofstream & outFile);
    void ten_most_abundant_kmers(std::ofstream & outFile);

private:
    bool first;
    unsigned vsize;
    unsigned covsize;
    int shift;
    int id;
    seqan::Shape<seqan::Dna, seqan::UngappedShape<8> > myShape;
    seqan::String<unsigned> v1;
    seqan::String<unsigned> v2;
    seqan::String<unsigned> poscov;
    seqan::String<uint64_t> adapterKmercount;
    seqan::String<unsigned> kmer_histogram;
    seqan::String<uint64_t> ten_max_kmers; //todo move to class and make seperate print function
    seqan::StringSet<seqan::DnaString> adapterSet;
    void update_coverage();
    void update_vectors();

};

template <typename TString>
void print(TString & str, std::ofstream & outFile){
   typedef typename seqan::Iterator< TString, seqan::Standard>::Type TIterator;
   for(TIterator it = begin(str, seqan::Standard()); it != end(str, seqan::Standard()); ++it)
    {
      outFile << " " << *it;
    }
    outFile << std::endl;
}

void RunBamStream(std::vector<ReadQualityHasher> & sps, seqan::CharString & seq, seqan::CharString & qual, size_t & qsize, size_t & ksize);

int main(int argc, char const ** argv)
{

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

    typedef seqan::StringSet<seqan::CharString> TNameStore; // Setup name store, cache, and BAM I/O context.
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
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

    //Check if more then one lane.
    unsigned lanecount = 0;
    seqan::CharString lanename;
    seqan::CharString pn;
    std::map<seqan::CharString, int> laneNames;
    for (unsigned i = 0; i < length(header.records); ++i)
    {
        if (header.records[i].type == seqan::BAM_HEADER_READ_GROUP)
        {
            for (unsigned j = 0; j < length(header.records[i].tags); ++j)
            {
                if (seqan::getValueI1(header.records[i].tags[j]) == "ID")
                {
                    lanename = seqan::getValueI2(header.records[i].tags[j]);
                    laneNames[lanename] = lanecount;
                    lanecount +=1;
                }
                if (seqan::getValueI1(header.records[i].tags[j]) == "SM")
                {
                    pn = seqan::getValueI2(header.records[i].tags[j]);
                }
            }
        }
    }

    seqan::clear(header);

    seqan::String<QualityCheck> r1;
    seqan::String<QualityCheck> r2;
    resize(r1,lanecount, QualityCheck(opt.isize));
    resize(r2,lanecount, QualityCheck(opt.isize));


    seqan::String<OverallNumbers> rall;
    resize(rall,lanecount);

    size_t qsize = opt.q_cutoff.size();
    size_t ksize = opt.klist.size();
    // Init ReadQualityHasher
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
        unsigned l;
        for (unsigned tagid = 0; tagid < length(tagsDict); ++tagid)
        {
            if (getTagKey(tagsDict, tagid) == "RG")
            {
                if (getTagType(tagsDict, tagid) == 'Z')
                {
                    seqan::CharString readlane;
                    readlane = getTagValue(tagsDict, tagid);
                    readlane = infix(readlane, 1, length(readlane)-1);
                    l = laneNames[readlane];
                    break;
                }
                else
                {
                    std::cout << "Read does not have Z" << "\n";
                    if (write2(std::cout, record, context, seqan::Sam()) != 0)
                    {
                        std::cerr << "ERROR: Could not write record to stdout \n";
                    }
                    return 1;
                }
            }
        }

        // Check the four highest flags and continue if any is set.
        if (!hasFlagSupplementary(record)) {
            rall[l].supplementary +=1;
            continue;
        }
        if (hasFlagDuplicate(record)){
            rall[l].duplicates +=1;
            continue;
        }
        if (hasFlagQCNoPass(record)){
            rall[l].QCfailed +=1;
            continue;
        }
        if (hasFlagSecondary(record)){
            rall[l].not_primary_alignment +=1;
            continue;
        }

        // Triplet counting.
        if (tripletCounting(counts, record, nameStore, genome, tripletCountingOptions) != 0) {
            return 1;
        }

        if (hasFlagRC(record)) // check if read is reversed complemented
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

    for (std::map<seqan::CharString,int>::iterator it=laneNames.begin(); it!=laneNames.end(); ++it)
    {
            outFile << "pn " << pn << std::endl;
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
            outFile << "genome_coverage_histogram"; rall[lid].get_poscov(outFile);
            outFile << "insert_size_histogram"; r1[lid].get_insert_size_histogram(outFile);
            outFile << "read_length_histogram_first"; r1[lid].get_read_length_histogram(outFile);
            outFile << "read_length_histogram_second"; r2[lid].get_read_length_histogram(outFile);
            outFile << "N_count_histogram_first"; r1[lid].get_N_count_histogram(outFile);
            outFile << "N_count_histogram_second"; r2[lid].get_N_count_histogram(outFile);
            outFile << "GC_content_histogram_first"; r1[lid].get_GC_count_histogram(outFile);
            outFile << "GC_content_histogram_second"; r2[lid].get_GC_count_histogram(outFile);
            outFile << "average_base_qual_histogram_first"; r1[lid].get_average_base_qual_histogram(outFile);
            outFile << "average_base_qual_histogram_second"; r2[lid].get_average_base_qual_histogram(outFile);
            outFile << "mapping_qual_histogram_first"; r1[lid].get_mapping_qual_histogram(outFile);
            outFile << "mapping_qual_histogram_second"; r2[lid].get_mapping_qual_histogram(outFile);
            outFile << "mismatch_count_histogram_first"; r1[lid].get_mismatch_count_histogram(outFile);
            outFile << "mismatch_count_histogram_second"; r2[lid].get_mismatch_count_histogram(outFile);
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
            outFile << "average_base_qual_by_position_first"; r1[lid].get_avg_base_qual_by_pos(outFile);
            outFile << "average_base_qual_by_position_second"; r2[lid].get_avg_base_qual_by_pos(outFile);
            outFile << "soft_clipping_5_prime_by_position_first"; r1[lid].get_soft_clip_by_pos_5prime(outFile);
            outFile << "soft_clipping_3_prime_by_position_first"; r1[lid].get_soft_clip_by_pos_3prime(outFile);
            outFile << "soft_clipping_5_prime_by_position_second"; r2[lid].get_soft_clip_by_pos_5prime(outFile);
            outFile << "soft_clipping_3_prime_by_position_second"; r2[lid].get_soft_clip_by_pos_3prime(outFile);
            rall[lid].ten_most_abundant_kmers(outFile);
            outFile << "8mer_count"; rall[lid].get_adapterKmercount(outFile);
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

    return 0;
}

QualityCheck::QualityCheck(): maxtlen(1000), softclipping('S'), qualcount_readnr(0)
{
    resize(insertSize, maxtlen +1, 0);
}

QualityCheck::QualityCheck(int isize): maxtlen(isize), softclipping('S'), qualcount_readnr(0)
{
    resize(insertSize, maxtlen +1, 0);
}

int QualityCheck::check_read_len(seqan::CharString & seq, seqan::CharString & qual)
{
    // COUNTS PER READ POSITION
    dnaseq = seq;
    if (length(dnaseq) != length(qual))
    {
        std::cerr << "ERROR: length of sequence and quality is not the same" << "\n";
        return 1;
    }
    return 0;
}

void QualityCheck::resize_strings()
{
    // count for nucleotides and qualtype for each read position:
    // Increase number of counters if dnaseq is longer than the previous reads.
    if (length(dnacount) < length(dnaseq)) //This should only be true once...
    {
        unsigned oldSize = length(dnacount);
        resize(dnacount, length(dnaseq), 0);
        resize(qualcount, length(dnacount), 0);
        resize(avgqualcount, length(dnacount), 0);
        resize(Ncount, length(dnacount)+1,0);
        resize(GCcount, length(dnacount)+1,0);
        resize(scposcount_5prime, length(dnacount), 0); //check if n is only assigned to added size.
        resize(scposcount_3prime, length(dnacount), 0); //check if n is only assigned to added size.

        for (unsigned j = oldSize; j < length(dnacount); ++j)
        {
            resize(dnacount[j], 5, 0);
        }
    }
}

void QualityCheck::get_count(seqan::CharString & seq, seqan::CharString & qual)
{
    resize_strings();
    read_counts(seq, qual);
    read_length(seq);
}

void QualityCheck::read_counts(seqan::CharString & seq, seqan::CharString & qual)
{
    unsigned cntN = 0;
    unsigned cntGC = 0;
    unsigned avgQual = 0;

    qualcount_readnr+=1;

    unsigned j = 0;
    seqan::Iterator<seqan::Dna5String, seqan::Rooted>::Type itseq = begin(dnaseq);
    seqan::Iterator<seqan::Dna5String, seqan::Rooted>::Type itEndseq = end(dnaseq);
    for (; itseq != itEndseq; goNext(itseq))
    {
        dnacount[j][(int) seqan::ordValue(*itseq)] += 1;

        if (*itseq == 'N') // count number of Ns per read
        {
            cntN +=1;
        }
        if ((*itseq == 'C') || (*itseq == 'G')) // count number of GCs per read
        {
            cntGC +=1;
        }
        ++j;
    }
    j = 0;
    seqan::Iterator<seqan::CharString, seqan::Rooted>::Type itqual = begin(qual);
    seqan::Iterator<seqan::CharString, seqan::Rooted>::Type itEndqual = end(qual);
    for (; itqual != itEndqual; goNext(itqual))
    {
        qualcount[j]+= (int) seqan::ordValue(*itqual)-33; // sums up the quality value for all reads per postion
        avgQual += (int) seqan::ordValue(*itqual)-33; // sums up all quality values for the read
        ++j;
    }

    Ncount[cntN] += 1; //count number of Ns per read

    GCcount[cntGC] +=1;

    if (length(averageQual) <= ceil(double(avgQual)/length(seq))){
        resize(averageQual, ceil(double(avgQual)/length(seq))+1, 0);
        }
    averageQual[(int)round(double(avgQual)/length(seq))] +=1; // histogram: frequency of reads with x-average quality value.
}

void QualityCheck::read_length(seqan::CharString & seq)
{
    unsigned lseq = length(seq);
    if (length(readLength) <= lseq){
        resize(readLength, lseq +1, 0);
        }
    readLength[lseq] += 1; //histogram of read lengths
}

void QualityCheck::map_Q(__uint8 & mapq)
{
    if (length(mapQ) <= mapq){
        resize(mapQ, mapq +1, 0);
        }
    mapQ[mapq] +=1; //histogram of mapping Qualities
}

void QualityCheck::insert_size(int & tlen)
{
    unsigned index = abs(tlen);

    if (index >  maxtlen){
        index = maxtlen;
        }
    insertSize[index] += 1; //histogram of insert size
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
                extractTagValue(x, tagsDict, tagid); //x gives number of mismatches + deletions + insertions
                unsigned mmcount = x - DIcount; //x - deletions -insertions gives number of mismatches

                if (length(mismatch) <= mmcount){
                    resize(mismatch, mmcount +1, 0);
                    }
                mismatch[mmcount] +=1; // histogram of number of mismatches
            }
        }
    }
}

void QualityCheck::cigar_count(seqan::BamAlignmentRecord & record)
{
    // function
    int c_begin = 0;                        // begin position of cigar in seq
    int c_end = 0;                          // end position of cigar in seq
    unsigned Scount_begin = 0;                   //counting soft-clipping at the beginning
    unsigned Scount_end = 0;                     //counting soft-clipping at the end
    int cigarlength = length(record.cigar); // cigar is a paired string with 'operation' (ie. S D I M) and 'count'
    DIcount = 0;                            // DIcount is used later in the function mis_match

    if (record.cigar[0].operation == 'S')
    {
        Scount_begin += record.cigar[0].count;
        for (unsigned j = 0; j < record.cigar[0].count; ++j)
        {
            scposcount_5prime[j] +=1; //if Soft-clipping at read position
        }
    }
    else if (record.cigar[cigarlength-1].operation == 'S')
    {
        Scount_end += record.cigar[cigarlength-1].count;
        for (unsigned j = length(record.seq) -record.cigar[cigarlength-1].count; j < length(record.seq); ++j)
        {
            scposcount_3prime[j] +=1; //if Soft-clipping at read position
        }
    }

    for (int i = 0; i < cigarlength; ++i)
    {
        if ((record.cigar[i].operation == 'D') || (record.cigar[i].operation == 'I')) //TODO: check if N,P,=,X matter here
        {
            DIcount += record.cigar[i].count; // counts deletions and insertions for read. Later used in function mis_match
        }

    }
    if (length(sccount_begin) <= Scount_begin){
        resize(sccount_begin, Scount_begin+1, 0);
        }
    sccount_begin[Scount_begin] += 1; //histogram: frequency of reads with x-number of soft-clipping at the beginning

    if (length(sccount_end) <= Scount_end){
        resize(sccount_end, Scount_end+1, 0);
        }
    sccount_end[Scount_end] += 1; //histogram: frequency of reads with x-number of soft-clipping at the end
}

void QualityCheck::avgQualPerPos()
{
    for (unsigned i = 0; i < length(qualcount); ++i)
    {
        avgqualcount[i] = (qualcount[i]/(double)qualcount_readnr); // sum of all read qualities per position / number of reads TODO: check if correct
    }
}

void QualityCheck::get_DNA_by_position(int c, std::ofstream & outFile)
{
    for (unsigned j = 0; j < length(dnacount); ++j)
    {
        outFile << " " << dnacount[j][c];
    }
    outFile << std::endl;
}

void QualityCheck::get_avg_base_qual_by_pos(std::ofstream & outFile)
{
    print(avgqualcount, outFile);
}

void QualityCheck::get_soft_clip_by_pos_5prime(std::ofstream & outFile)
{
    print(scposcount_5prime, outFile);
}
void QualityCheck::get_soft_clip_by_pos_3prime(std::ofstream & outFile)
{
    print(scposcount_3prime, outFile);
}
void QualityCheck::get_N_count_histogram(std::ofstream & outFile)
{
    print(Ncount, outFile);
}
void QualityCheck::get_GC_count_histogram(std::ofstream & outFile)
{
    print(GCcount, outFile);
}
void QualityCheck::get_soft_clipping_begin_histogram(std::ofstream & outFile)
{
    print(sccount_begin, outFile);
}
void QualityCheck::get_soft_clipping_end_histogram(std::ofstream & outFile)
{
    print(sccount_end, outFile);
}
void QualityCheck::get_average_base_qual_histogram(std::ofstream & outFile)
{
    print(averageQual, outFile);
}
void QualityCheck::get_insert_size_histogram(std::ofstream & outFile)
{
    print(insertSize, outFile);
}
void QualityCheck::get_mapping_qual_histogram(std::ofstream & outFile)
{
    print(mapQ, outFile);
}
void QualityCheck::get_mismatch_count_histogram(std::ofstream & outFile)
{
    print(mismatch, outFile);
}
void QualityCheck::get_read_length_histogram(std::ofstream & outFile)
{
    print(readLength, outFile);
}
void QualityCheck::clear_var()
{
    seqan::clear(dnacount);
    seqan::clear(qualcount);
}

OverallNumbers::OverallNumbers() : supplementary(0), duplicates(0), QCfailed(0), not_primary_alignment(0), readcount(0), totalbps(0), bothunmapped(0), firstunmapped(0), secondunmapped(0), discordant(0),
first(true), vsize(1000), covsize(100), shift(0), id(0)
{
    resize(adapterKmercount, 65536, 0); //TTTTTTTT = 65535
}

void OverallNumbers::update_vectors()
{
    seqan::clear(v1);
    resize(v1,vsize, 0);
    seqan::swap(v1, v2);
}

void OverallNumbers::update_coverage()
{
    for (unsigned i = 0; i < vsize; i++)
    {
        if (v1[i] > covsize){
            poscov[covsize]+=1;
        }
        else {
            poscov[v1[i]]+=1;
        }
    }
}

void OverallNumbers::coverage(seqan::BamAlignmentRecord & record)
{
    unsigned beginpos = record.beginPos;

    if (first) // move to initializer
    {
        first = false;
        id = record.rID;
        shift = beginpos;
        resize(poscov, covsize +1, 0);
        resize(v1,vsize, 0);
        resize(v2,vsize, 0);
    }
    if (seqan::isNotEqual(id, record.rID) || ((beginpos - shift) > 2*vsize))
    {
        id = record.rID;
        update_coverage();
        update_vectors();
        update_coverage();
        seqan::clear(v1);
        resize(v1,vsize, 0);
        shift = beginpos;
    }

    unsigned pos = beginpos - shift;

    if ((pos > vsize) && (pos < 2*vsize)) // then we update_coverage and set v1 = v2 and v2 = {0}
    {
        update_coverage();
        update_vectors();
        shift += vsize;
    }
    int c = 0;
    for (unsigned i = 0; i < length(record.cigar); i++)
    {
        if (record.cigar[i].operation == 'S')
        {
            c += record.cigar[i].count;
        }
        if (record.cigar[i].operation == 'M' || record.cigar[i].operation == 'D')
        {
            for (unsigned j = c; j < record.cigar[i].count; j++) //replace to for (c; c < record.cigar[i].count; c++) and skip c +=
            {
                if (pos + j < vsize) {
                    v1[pos + j] +=1;
                }
                else {
                    v2[pos -vsize + j] +=1;
                }
            }
            c += record.cigar[i].count;
        }
    }
}

void OverallNumbers::countAdapterKmers(seqan::CharString & seq)
{
    typedef seqan::Iterator<seqan::DnaString>::Type TIterator;
    typedef seqan::Iterator<seqan::CharString>::Type TCharIterator;
    seqan::DnaString text = seq;
    int index;
    unsigned skip = 0;
    hashInit(myShape, begin(text));

    TCharIterator itSeq = seqan::begin(seq);
    for (unsigned i = 0; i < 7; ++i, ++itSeq)
    {
        if (*itSeq == 'N') skip = 8;
        if (skip > 0) --skip;
    }

    TIterator itEnd = seqan::end(text) - length(myShape) + 1;
    for (TIterator it = seqan::begin(text); it < itEnd; ++it, ++itSeq)
    {
        index = hashNext(myShape, it);
        if (*itSeq == 'N') skip = 8;
        if (skip > 0) --skip;
        else
        {
            if (index < 65536){
                ++adapterKmercount[index];
            }
            else
                std::cout << "hash_index too big\n";
        }
    }
}

void OverallNumbers::get_poscov(std::ofstream & outFile)
{
    print(poscov, outFile);
}

void OverallNumbers::get_adapterKmercount(std::ofstream & outFile)
{
    print(adapterKmercount, outFile);
}

void RunBamStream(std::vector<ReadQualityHasher> & sps, seqan::CharString & seq, seqan::CharString & qual, size_t & qsize, size_t & ksize) {

    for (size_t i = 0; i < qsize*ksize; i++) {
        sps[i](seqan::toCString(seq), seqan::length(seq),
        seqan::toCString(qual), seqan::length(qual));
    }
}

void OverallNumbers::ten_most_abundant_kmers(std::ofstream & outFile)
{
    seqan::String<uint64_t> v_copy(adapterKmercount);
    resize(ten_max_kmers,10,0);

    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type it = begin(v_copy);
    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itEnd = end(v_copy);
    std::nth_element(it, it+10, itEnd, std::greater<int>());

    for (int i=0; i<10; i++)
    {
        ten_max_kmers[i] = v_copy[i];
    }
    clear(v_copy);

    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itmax = begin(ten_max_kmers);
    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itmaxEnd = end(ten_max_kmers);
    std::sort(itmax, itmaxEnd, std::greater<uint64_t>());


    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itv = begin(adapterKmercount);
    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itvEnd = end(adapterKmercount);

    seqan::DnaString result; //todo move to class and make seperate print function
    unsigned q = 8;
    int previos_pos = 0;

    std::set<int> posSet;

    for (int i=0; i<10; i++)
    {
        int pos = std::find(itv, itvEnd, ten_max_kmers[i]) - itv;
        if (posSet.count(pos)==0)
        {
            posSet.insert(pos);
        }
        else
        {
            while (posSet.count(pos)!=0)
            {
                seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itvnext = begin(adapterKmercount) +pos+1;
                pos = std::find(itvnext, itvEnd, ten_max_kmers[i]) - itv;
            }
            posSet.insert(pos);
        }
        seqan::unhash(result, pos, q);
        outFile << "nr_" << i+1 << "_most_abundant_8mer "<< result << " " << ten_max_kmers[i] << std::endl;
    }
}













