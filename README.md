# BamQC
Quality control for BAM files


## Dependencies

BamQC depends on the SeqAn core library, version 1.4.2 (https://github.com/seqan/seqan).
Please use the 1.4.2 release of SeqAn as no other version will work.

## Installation

If you put the SeqAn library in a default include directory:

    $ git clone https://github.com/DecodeGenetics/BamQC.git
    $ cd BamQC
    $ make
    $ ./bamqualcheck

If you put the SeqAn library in another directory, you will need to edit line 14 of the Makefile and specify the path, e.g.

    CXXFLAGS+=-I../libraries/seqan-1.4.2/include

## Usage

The ``bamqualcheck`` program takes as input a BAM file. In addition, the reference genome and an output file are required as options:

    $ ./bamqualcheck -r genome.fa -o sample.bamqc sample.bam

If you are using another reference than the human genome or if your chromosomes are not named ``chr1, chr2, ...``, you will need to specify the list of main chromosomes using the option ``-c``, e.g.:

    $ ./bamqualcheck -r genome.fa -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, 17,18,19,20,21,22,X,Y -o sample.bamqc sample.bam

If your library has insert sizes exceeding 1000 bps, you need to increase the maximum in the insert size histogram using the option ``-i``, e.g.

    $ ./bamqualcheck -r genome.fa -i 3000 -o sample.bamqc sample.bam

## Output

The program output is a text file where each line corresponds to one count or hisotgram of counts.
If you are interested in only one specific measure, it can easily be extracted on the command line using ``grep``, e.g.

    $ grep insert_size_histogram sample.bamqc
