#ifndef READ_QUAL_HASHER_H_
#define READ_QUAL_HASHER_H_

#include "kmerstream/StreamCounter.hpp"
#include "kmerstream/RepHash.hpp"

#include "CommandLineParser.hpp"

// ==========================================================================
// CLASS ReadQualityHasher
// ==========================================================================

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

// ==========================================================================
// FUNCTION RunBamStream
// ==========================================================================

void RunBamStream(std::vector<ReadQualityHasher> & sps, seqan::CharString & seq, seqan::CharString & qual, size_t & qsize, size_t & ksize) {

    for (size_t i = 0; i < qsize*ksize; i++) {
        sps[i](seqan::toCString(seq), seqan::length(seq),
        seqan::toCString(qual), seqan::length(qual));
    }
}

#endif  // READ_QUAL_HASHER_H_
