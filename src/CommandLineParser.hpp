#ifndef COMMAND_LINE_PARSER_H_
#define COMMAND_LINE_PARSER_H_

#include <seqan/arg_parse.h>

#include "version.h"


// --------------------------------------------------------------------------
// CLASS ProgramOptions
// --------------------------------------------------------------------------

struct ProgramOptions
{
    // Input and output files.
    seqan::CharString bamFile;
    seqan::CharString referenceFile;
    seqan::CharString outputFile;

    // Main chromosomes.
    seqan::CharString chroms;

    // Parameters of kmerstream.
    std::vector<int> klist;
    double e;
    std::vector<size_t> q_cutoff;
    size_t q_base;
    int seed;

    // Maximum insert size to be counted in histogram.
    int isize;

    ProgramOptions() :
        chroms("chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"),
        referenceFile("genome.fa"), e(0.01), q_base(33), seed(0)
    {}
};

// --------------------------------------------------------------------------
// FUNCTION parseCommandLine()
// --------------------------------------------------------------------------

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
    addDescription(parser, "Program for bam quality checks. The program can read from bamfile or stdin(sam format)");

    // Add the required argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "bamfile", "False", 1));
    setValidValues(parser, 0, "- bam sam");

    // Add general options.
    addSection(parser, "General options");
    addOption(parser, seqan::ArgParseOption("r", "reference", "Reference genome filename.", seqan::ArgParseArgument::STRING, "FILENAME"));
    setDefaultValue(parser, "r", "genome.fa");
    addOption(parser, seqan::ArgParseOption("i", "insert-size", "Upper bound for the insert size in insert size histogram.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "i", "1000");
    addOption(parser, seqan::ArgParseOption("c", "chromosomes", "Comma separated list of the main chromosome names.", seqan::ArgParseArgument::STRING, "STRING"));
    addOption(parser, seqan::ArgParseOption("o", "output-file", "Output filename.", seqan::ArgParseArgument::OUTPUTFILE, "OUT"));

    // Add kmerstream options.
    addSection(parser, "Kmerstream options");
    addOption(parser, seqan::ArgParseOption("k", "kmer-size", "Comma-separated list of k-mer sizes.", seqan::ArgParseArgument::STRING, "STRING"));
    setDefaultValue(parser, "k", "32");
    addOption(parser, seqan::ArgParseOption("q", "quality-cutoff", "Comma-separated list of PHRED quality thresholds. K-mers with bases above quality threshold are kept.", seqan::ArgParseArgument::STRING, "STRING"));
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
    std::ostringstream errorStream;
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv, std::cout, errorStream);
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        if (seqan::isEqual(errorStream.str(), "bamqualcheck: illegal option -- -\n"))
        {
            std::cerr << "Reading from stdin" << std::endl;
        }
        else
        {
            std::cerr << errorStream.str() << std::endl;
            return res;
        }
    }
    // check if argument is -
    std::string lastarg = argv[argc-1];
    if (lastarg == "-") {
        argc -=1;
        options.bamFile = '-';
    }

    // Get argument and option values.

    std::string qcut, kmer;
    getOptionValue(options.referenceFile, parser, "r");
    getOptionValue(qcut, parser, "q");
    getOptionValue(kmer, parser, "k");
    getOptionValue(options.e, parser, "e");
    getOptionValue(options.seed, parser, "s");
    getOptionValue(options.isize, parser, "i");
    getOptionValue(options.chroms, parser, "c");
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

    if (length(options.bamFile) == 0) {
        getArgumentValue(options.bamFile, parser, 0);

     }

    return seqan::ArgumentParser::PARSE_OK;
}

#endif  // COMMAND_LINE_PARSER_H_
