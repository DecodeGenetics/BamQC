#ifndef OVERALL_NUMBERS_H_
#define OVERALL_NUMBERS_H_

// ==========================================================================
// CLASS OverallNumbers
// ==========================================================================

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
    unsigned FF_RR_orientation;
    unsigned properpair_count;

    seqan::String<unsigned> poscov;
    seqan::String<uint64_t> eightmercount;

    void coverage(seqan::BamAlignmentRecord & record);
    void count8mers(seqan::CharString & seq);
    void get_kmer_histogram(std::ofstream & outFile);
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
    seqan::String<unsigned> kmer_histogram;
    seqan::String<uint64_t> ten_max_kmers;
    void update_coverage();
    void update_vectors();

};

OverallNumbers::OverallNumbers() : supplementary(0), duplicates(0), QCfailed(0), not_primary_alignment(0), readcount(0), totalbps(0), bothunmapped(0), firstunmapped(0), secondunmapped(0), discordant(0),
FF_RR_orientation(0), properpair_count(0), first(true), vsize(1000), covsize(100), shift(0), id(0)
{
    resize(eightmercount, 65536, 0); //TTTTTTTT = 65535
    resize(poscov, covsize + 1, 0);
    resize(v1, vsize, 0);
    resize(v2, vsize, 0);
}

void OverallNumbers::update_vectors()
{
    seqan::clear(v1);
    resize(v1, vsize, 0);
    seqan::swap(v1, v2);
}

void OverallNumbers::update_coverage()
{
    for (unsigned i = 0; i < vsize; ++i)
    {
        if (v1[i] > covsize){
            poscov[covsize] += 1;
        }
        else {
            poscov[v1[i]] += 1;
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
    }

    if (seqan::isNotEqual(id, record.rID) || ((beginpos - shift) > 2*vsize))
    {
        id = record.rID;
        update_coverage();
        update_vectors();
        update_coverage();
        seqan::clear(v1);
        resize(v1, vsize, 0);
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
    for (unsigned i = 0; i < length(record.cigar); ++i)
    {
        if (record.cigar[i].operation == 'S')
        {
            c += record.cigar[i].count;
        }
        if (record.cigar[i].operation == 'M' || record.cigar[i].operation == 'D')
        {
            for (unsigned j = c; j < record.cigar[i].count; ++j) // replace to for (c; c < record.cigar[i].count; c++) and skip c +=
            {
                if (pos + j < vsize)
                {
                    v1[pos + j] += 1;
                }
                else
                {
                    v2[pos - vsize + j] += 1;
                }
            }
            c += record.cigar[i].count;
        }
    }
}

void OverallNumbers::count8mers(seqan::CharString & seq)
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
                ++eightmercount[index];
            }
            else
                std::cout << "hash_index too big\n";
        }
    }
}

void OverallNumbers::ten_most_abundant_kmers(std::ofstream & outFile)
{
    seqan::String<uint64_t> v_copy(eightmercount);
    resize(ten_max_kmers, 10, 0);

    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type it = begin(v_copy);
    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itEnd = end(v_copy);
    std::nth_element(it, it+10, itEnd, std::greater<int>());

    for (int i = 0; i < 10; ++i)
    {
        ten_max_kmers[i] = v_copy[i];
    }
    clear(v_copy);

    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itmax = begin(ten_max_kmers);
    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itmaxEnd = end(ten_max_kmers);
    std::sort(itmax, itmaxEnd, std::greater<uint64_t>());


    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itv = begin(eightmercount);
    seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itvEnd = end(eightmercount);

    seqan::DnaString result;
    unsigned q = 8;
    int previos_pos = 0;

    std::set<int> posSet;

    for (int i = 0; i < 10; ++i)
    {
        int pos = std::find(itv, itvEnd, ten_max_kmers[i]) - itv;
        if (posSet.count(pos) == 0)
        {
            posSet.insert(pos);
        }
        else
        {
            while (posSet.count(pos) != 0)
            {
                seqan::Iterator<seqan::String<uint64_t>, seqan::Standard>::Type itvnext = begin(eightmercount) + pos + 1;
                pos = std::find(itvnext, itvEnd, ten_max_kmers[i]) - itv;
            }
            posSet.insert(pos);
        }
        seqan::unhash(result, pos, q);
        outFile << "nr_" << i+1 << "_most_abundant_8mer "<< result << " " << ten_max_kmers[i] << std::endl;
    }
}

#endif  // OVERALL_NUMBERS_H_
