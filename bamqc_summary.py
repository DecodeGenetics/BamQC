#!/nfs/prog/bioinfo/apps-x86_64/python/3.5.2/bin/python3 -E

import os, sys
import argparse
from math import sqrt, exp
from scipy.optimize import brentq

# ------------------------------------------------------------------------------

 # thresholds stores three min and three max values for each field (for marking fields with ***, ** or *)
thresholds = {}

thresholds['mean_base_qual_per_read'] = {min:[25, 30, 32], max:[100, 100, 100]}
thresholds['std_base_qual_per_read'] = {min:[-1, -1, -1], max:[100, 10, 7]}
thresholds['mean_N_per_read'] =  {min:[-1, -1, -1], max:[20, 10, 5]}
thresholds['std_N_per_read'] = {min:[-1, -1, -1], max:[1000, 30, 20]}
thresholds['mean_GC_per_read'] = {min:[30, 39, 40], max:[50, 45, 43]}
thresholds['std_GC_per_read'] = {min:[-1, -1, -1], max:[18, 15, 12]}
thresholds['mean_base_qual_per_position'] = {min:[25, 30, 33], max:[100, 100, 100]}
thresholds['std_base_qual_per_position'] = {min:[-1, -1, -1], max:[8, 6, 4.5]}
thresholds['mean_N_per_position'] = {min:[-1, -1, -1], max:[30, 10, 5]}
thresholds['std_N_per_position'] = {min:[-1, -1, -1], max:[20, 10, 5]}
thresholds['mean_A_per_position'] = {min:[-1, 25, 27.25], max:[1000, 35, 32.5]}  # Not setting A,C,T,G params, should follow from GC params
thresholds['std_A_per_position'] = {min:[-1, -1, -1], max:[1000, 10, 5]}
thresholds['mean_C_per_position'] = {min:[-1, 15.5, 17.5], max:[1000, 25, 23]}  
thresholds['std_C_per_position'] = {min:[-1, -1, -1], max:[1000, 10, 5]}
thresholds['mean_G_per_position'] = {min:[-1, 17, 18], max:[1000, 24, 23]}
thresholds['std_G_per_position'] = {min:[-1, -1, -1], max:[1000, 10, 5]}
thresholds['mean_T_per_position'] = {min:[-1, 25, 27], max:[1000, 33, 31.5]}
thresholds['std_T_per_position'] = {min:[-1, -1, -1], max:[1000, 10, 5]}
thresholds['adapter_8_mers'] = {min:[-1, -1, -1], max:[10, 5, 3]}
thresholds['marked_duplicate'] = {min:[-1, -1, -1], max:[80, 60, 30]}
thresholds['unmapped'] = {min:[-1, -1, -1], max:[30, 20, 10]}
thresholds['both_unmapped'] = {min:[-1, -1, -1], max:[40, 30, 17]}
thresholds['first_unmapped'] = {min:[-1, -1, -1], max:[40, 30, 18]}
thresholds['second_unmapped'] = {min:[-1, -1, -1], max:[40, 30, 23]}
thresholds['FF_RR_oriented_pairs'] = {min:[-1, -1, -1], max:[1, 0.1, 0.01]}               ## <--- TODO
thresholds['proper_pairs'] = {min:[55, 70, 85], max:[1000, 1000, 1000]}
thresholds['proper_pairs_autosome'] = {min:[95, 95, 95], max:[1000, 1000, 1000]}
thresholds['mean_coverage'] = {min:[0.1, 0.1, 1], max:[10000, 100000, 100000]}
thresholds['total_bps'] = {min:[300000000, 300000000, 10000000000], max:[10000000000000, 100000000000000, 100000000000000]}
thresholds['std_coverage'] = {min:[-1, -1, -1], max:[100000, 100000, 100000]}
thresholds['mean_insert_size'] = {min:[-1, -1, 200], max:[10000, 10000, 10000]}
thresholds['adapter_insert_size'] = {min:[-1, -1, -1], max:[100, 20, 20]}
thresholds['std_times_mean_insert_size'] = {min:[-1, -1, -1], max:[10000000, 80000, 70000]}
thresholds['mapping_qual_60'] = {min:[-1, -1, -1], max:[1000, 1000, 15]}
thresholds['mapping_qual_40'] = {min:[-1, -1, -1], max:[1000, 1000, 12]}
thresholds['mapping_qual_20'] = {min:[-1, -1, -1], max:[1000, 1000, 9.5]}
thresholds['mean_mismatches'] = {min:[-1, -1, -1], max:[10, 5, 3]}                 ## <--- TODO
thresholds['mean_deletions'] = {min:[-1, -1, -1], max:[1000, 1000, 1000]}                  ## <--- TODO
thresholds['mean_insertions'] = {min:[-1, -1, -1], max:[1000, 1000, 1000]}                 ## <--- TODO
thresholds['nz_deletions'] = {min:[-1, -1, -1], max:[0.3, 0.1, 0.05]}                  ## <--- TODO
thresholds['nz_insertions'] = {min:[-1, -1, -1], max:[0.3, 0.1, 0.05]}                 ## <--- TODO
thresholds['clipped_5_prime'] = {min:[-1, -1, -1], max:[100, 6, 4]}
thresholds['clipped_3_prime'] = {min:[-1, -1, -1], max:[100, 30, 20]}
thresholds['CtoA'] = {min:[-1, 0.3, 0.4], max:[2,0.7,0.6]}
thresholds['GtoA'] = {min:[-1, 0.4, 0.45], max:[2,0.6,0.55]}
thresholds['TtoA'] = {min:[-1, 0.3, 0.4], max:[2,0.7,0.6]}
thresholds['AtoC'] = {min:[-1, 0.3, 0.4], max:[2,0.7,0.6]}
thresholds['GtoC'] = {min:[-1, 0.3, 0.4], max:[2,0.7,0.6]}
thresholds['TtoC'] = {min:[-1, 0.3, 0.425], max:[2,0.7,0.575]}

# ------------------------------------------------------------------------------

def parse_arguments():

    parser = argparse.ArgumentParser(description='Write a summary line of bam quality checks per lane.')

    # Positional arguments    
    parser.add_argument('bamqc_file', help='output file of the bam_qual_check program')
    
    # Optional arguments
    parser.add_argument('-l', dest='long', action='store_true', help='output in long format')
#    parser.add_argument('-d', dest='strictness', default='2', help='strictness 3, 2 or 1 for flags in dense output format (default 3)')
    parser.add_argument('-t', dest='title', action='store_true', help='print a header line in dense output format')

    return parser.parse_args()

# ------------------------------------------------------------------------------

def printUsage(prog_name):
    sys.stderr.write('Usage: ' + prog_name + '<bamqc file>' + '\n')

# ------------------------------------------------------------------------------

def init_lane(sample_id, lane_name, i):
    lane = {}
    lane['sample_id'] = sample_id
    lane['lane'] = lane_name
    lane['lane_id'] = i
    lane['read_length'] = 'variable'
    return lane

# ------------------------------------------------------------------------------

def read_bamqc_output(data, filename):

    lane = None
    i = 1
    with open(filename, 'r') as inputfile:
        for line in inputfile:

            line = line.split()

            if line[0] == 'sample_id':
                sample_id = line[1]
            elif line[0] == 'lane':
                if lane != None:
                    if 'read_length_first' in lane and 'read_length_second' in lane and lane['read_length_first'] == lane['read_length_second']:
                        lane['read_length'] = lane['read_length_first']   
                    data[len(data):] = [lane]                     
                lane = init_lane(sample_id, line[1], i)
                i += 1
            elif line[0][0:21] == 'read_length_histogram':
                if len(line) == 1:
                    lane['read_length' + line[0][21:]] = 0
                elif int(line[len(line)-1]) == lane['total_read_pairs']:
                    lane['read_length' + line[0][21:]] = len(line) - 2
                else:
                    lane[line[0]] = [float(x) for x in line[1:]]
            elif line[0][:2] == 'nr':
                lane[line[0]] = [line[1], int(line[2])]
            elif len(line) == 2 and ("histogram" not in line[0]):
                lane[line[0]] = int(line[1])
            else:
                lane[line[0]] = [float(x) for x in line[1:]]

    if lane != None:
        if 'read_length_first' in lane and 'read_length_second' in lane and lane['read_length_first'] == lane['read_length_second']:
            lane['read_length'] = lane['read_length_first']
        data[len(data):] = [lane]

# ------------------------------------------------------------------------------

def hist_avg(lane, name):
    first = lane[name + '_first']
    second = lane[name + '_second']
    avg = [0] * max(len(first), len(second))
    
    for i in range(0, len(first)):
        avg[i] += first[i]
    for i in range(0, len(second)):
        avg[i] += second[i]
    for i in range(0, len(avg)):
        avg[i] /= 2

    return avg

# ------------------------------------------------------------------------------

def hist_sum(lane, name):
    first = lane[name + '_first']
    second = lane[name + '_second']
    both = [0] * max(len(first), len(second))

    for i in range(0, len(first)):
        both[i] += first[i]
    for i in range(0, len(second)):
        both[i] += second[i]

    return both

# ------------------------------------------------------------------------------

def mean(liste):
    total = 0
    for i in liste:
        total += i
        
    return total/len(liste) if len(liste) != 0 else 'NaN'

# ------------------------------------------------------------------------------

def std(liste):
    sum_dev = 0
    m = mean(liste)
    for i in liste:
        sum_dev += (i-m)**2
    
    return sqrt(sum_dev/len(liste)) if len(liste) != 0 else 'NaN'

# ------------------------------------------------------------------------------

def stats_by_pos(first, second, total_bps):
    try:
        total_bps_first = float(len(first)) / (len(first)+len(second)) * total_bps
        total_bps_second = float(len(second)) / (len(first)+len(second)) * total_bps
    
        total_first = 0
        for i in first:
            total_first += i
    
        total_second = 0
        for i in second:
            total_second += i

        mean_first = float(total_first) / total_bps_first
        mean_second = float(total_second) / total_bps_second
        mean = (len(first) * mean_first + len(second) * mean_second) / (len(first) + len(second))
    
        sum_dev = 0
        for i in first:
            sum_dev += (float(i) / total_bps_first * len(first) - mean_first)**2
        for i in second:
            sum_dev += (float(i) / total_bps_second * len(second) - mean_second)**2

        std = sqrt(sum_dev / (len(first) + len(second)))
    
        return 100*mean, 100*std, total_first + total_second

    except ZeroDivisionError:
        return 'NaN', 'NaN', 0


# ------------------------------------------------------------------------------

def hist_mean(hist):
    length = 0
    total = 0
    for i in range(0, len(hist)):
        length += hist[i]
        total += i*hist[i]

    return total/length if length != 0 else 'NaN'

# ------------------------------------------------------------------------------

def hist_nz(hist):
    length = 0
    total = 0
    for i in range(0, len(hist)):
        length += hist[i]
        if i != 0:
            total += hist[i]

    return total/length if length != 0 else 'NaN'

# ------------------------------------------------------------------------------

def hist_std(hist):
    length = 0
    sum_dev = 0
    m = hist_mean(hist)
    if m == 'NaN':
        return 'NaN'

    for i in range(0, len(hist)):
        length += hist[i]
        sum_dev += hist[i] * (i-m)**2
        
    return sqrt(sum_dev/length) if length != 0 else 'NaN'

# ------------------------------------------------------------------------------

def count_lt(hist, maximum):
    total = 0
    for i in range(0, min(maximum, len(hist))):
        total += hist[i]
    return total

# ------------------------------------------------------------------------------

def stats_GC(lane):
    first = lane['GC_content_histogram_first']
    second = lane['GC_content_histogram_second']
    N_first = lane['Ns_by_position_first']
    N_second = lane['Ns_by_position_second']
    read_len_first = lane['read_length_first']
    read_len_second = lane['read_length_second']

    try:
        total_N_first = 0
        total_N_second = 0
        for i in range(0, len(N_first)):
            total_N_first += N_first[i]
        for i in range(0, len(N_second)):
            total_N_second += N_second[i]

        avg_N_first = total_N_first / float(lane['total_read_pairs']) / len(N_first)
        avg_N_second = total_N_first / float(lane['total_read_pairs']) / len(N_second)

        length = 0
        total_first = 0
        total_second = 0
        for i in range(0, len(first)):
            length += first[i]
            total_first += i*first[i]
        for i in range(0, len(second)):
            length += second[i]
            total_second += i*second[i]

        mean = 100.0 * (total_first + total_second) / (lane['total_bps'] - total_N_first - total_N_second)
    
        sum_dev = 0
        for i in range(0, len(first)):
            sum_dev += first[i] * (100.0*i/(read_len_first-avg_N_first) - mean)**2
        for i in range(0, len(second)):
            sum_dev += second[i] * (100.0*i/(read_len_second-avg_N_second) - mean)**2

        std = sqrt(sum_dev/length)

        return mean, std

    except ZeroDivisionError:
        return 'NaN', 'NaN'

# ------------------------------------------------------------------------------

def k_mer_error_rate(F0, f1, F1, k):

    if F1 == 0:
        return 0, 'NA'

    epstol = 1e-8
    x = float(F1) / F0
    e = 0
    e0 = 1
    x0 = 0
    maxit = 1000
    it = 0

    while abs(x - x0) / x + abs(e - e0) > epstol and it < maxit:
        it = it + 1
        #print x, e, abs(x-x0)/x
        x0 = x
        e0 = e
    
        x1 = -100
        while abs(x - x1) / x > epstol / 2:
            x1 = x
            x  = float(F1) / F0 * (3 * k * (1 - exp(-x * e / (3 * k))) + 1 - exp(-x * (1 - e)))

        def func(t):
            return float(f1) / F1 - (t * exp(-x * t / (3 * k)) + (1 - t) * exp(-x * (1 - t)))

        try:
            e = brentq(func, 0, 1)
        except ValueError:
            return x, 'NA'

    return x, e

# ------------------------------------------------------------------------------

def adapter_8_mers(kmer_counts):
    # hash values of all 8-mers in HiSeq and HiSeqX Universal adapter sequence
    indices = [558,2203,2235,8813,8940,8942,9078,9947,10754,11976,12005,13954,14187,15282,23348,24014,26035,26658,27565,27859,28091,28362,29627,30394,30523,32907,33318,35254,35762,35769,36314,38605,38771,39790,40994,41097,43016,44438,44726,45623,45728,45903,46683,46828,47479,47833,47907,47914,48023,51421,51840,52972,55660,55816,56043,56242,56558,56750,58844,60261,60557,60584,61021,61128]

    adapter_8mers = 0
    for i in indices:
        adapter_8mers += kmer_counts[i]

    return adapter_8mers * 100 / sum(kmer_counts)

# ------------------------------------------------------------------------------

def adapter_insert_size(insert_size_histogram, read_length):
    if read_length == 'variable' or sum(insert_size_histogram) - insert_size_histogram[0] == 0:
        return 'NA'

    adapter_present = 0
    for i in range(1, read_length):
        adapter_present += insert_size_histogram[i]

    return adapter_present * 100 / (sum(insert_size_histogram) - insert_size_histogram[0])

# ------------------------------------------------------------------------------

def triplet_seq(index):
    seq = ['A','A','A']
    for i in [2, 1, 0]:
        rest = index % 4
        if rest == 1:
            seq[i] = 'C'
        elif rest == 2:
            seq[i] = 'G'
        elif rest == 3:
            seq[i] = 'T'
        index = index >> 2;
    return ''.join(seq)

# ------------------------------------------------------------------------------

def reverse_triplet(index):

    # compute triplet sequence from triplet index
    seq = triplet_seq(index)

    # reverse complement the triplet
    seq_rc = ['A','A','A']
    for i in [0, 1, 2]:
        if seq[i] == 'A':
            seq_rc[2-i] = 'T'
        elif seq[i] == 'C':
            seq_rc[2-i] = 'G'
        elif seq[i] == 'G':
            seq_rc[2-i] = 'C'
    
    # compute index of rc triplet
    index_rc = 0
    bit = [16, 4, 1]
    for i in [0, 1, 2]:
        if seq_rc[i] == 'C':
            index_rc += bit[i]*1
        elif seq_rc[i] == 'G':
            index_rc += bit[i]*2
        elif seq_rc[i] == 'T':
            index_rc += bit[i]*3
    
    return index_rc

# ------------------------------------------------------------------------------

def compute_triplet_scores(summary, lane):

    lane['triplet_scores_A'] = []
    lane['triplet_scores_C'] = []
    
    counts = {}
    counts['1st_FW'] = {'A':lane['triplet_counts_A_1st_FW'], 'C':lane['triplet_counts_C_1st_FW'],
                        'G':lane['triplet_counts_G_1st_FW'], 'T':lane['triplet_counts_T_1st_FW']}
    counts['2nd_FW'] = {'A':lane['triplet_counts_A_2nd_FW'], 'C':lane['triplet_counts_C_2nd_FW'],
                        'G':lane['triplet_counts_G_2nd_FW'], 'T':lane['triplet_counts_T_2nd_FW']}
    counts['1st_RC'] = {'A':lane['triplet_counts_A_1st_RC'], 'C':lane['triplet_counts_C_1st_RC'],
                        'G':lane['triplet_counts_G_1st_RC'], 'T':lane['triplet_counts_T_1st_RC']}
    counts['2nd_RC'] = {'A':lane['triplet_counts_A_2nd_RC'], 'C':lane['triplet_counts_C_2nd_RC'],
                        'G':lane['triplet_counts_G_2nd_RC'], 'T':lane['triplet_counts_T_2nd_RC']}
                        
    a = {'C':0, 'G':0, 'T':0}
    c = {'A':0, 'G':0, 'T':0}
    a_all = {'C':0, 'G':0, 'T':0}
    c_all = {'A':0, 'G':0, 'T':0}

    for triplet in range(0,len(counts['1st_FW']['A'])):
        triplet_rc = reverse_triplet(triplet)
        middle = triplet_seq(triplet)[1]

        if middle != 'A':
            x = counts['2nd_FW']['A'][triplet] + counts['1st_RC']['A'][triplet] + counts['1st_FW']['T'][triplet_rc] + counts['2nd_RC']['T'][triplet_rc]
            y = counts['1st_FW']['A'][triplet] + counts['2nd_RC']['A'][triplet] + counts['2nd_FW']['T'][triplet_rc] + counts['1st_RC']['T'][triplet_rc]
            a[middle] += x
            a_all[middle] += x + y

        if middle != 'C':
            x = counts['2nd_FW']['C'][triplet] + counts['1st_RC']['C'][triplet] + counts['1st_FW']['G'][triplet_rc] + counts['2nd_RC']['G'][triplet_rc]
            y = counts['1st_FW']['C'][triplet] + counts['2nd_RC']['C'][triplet] + counts['2nd_FW']['G'][triplet_rc] + counts['1st_RC']['G'][triplet_rc]
            c[middle] += x
            c_all[middle] += x + y

    summary['CtoA'] = a['C'] / a_all['C'] if a_all['C'] != 0 else 'NaN'
    summary['GtoA'] = a['G'] / a_all['G'] if a_all['G'] != 0 else 'NaN'
    summary['TtoA'] = a['T'] / a_all['T'] if a_all['T'] != 0 else 'NaN'
    summary['AtoC'] = c['A'] / c_all['A'] if c_all['A'] != 0 else 'NaN'
    summary['GtoC'] = c['G'] / c_all['G'] if c_all['G'] != 0 else 'NaN'
    summary['TtoC'] = c['T'] / c_all['T'] if c_all['T'] != 0 else 'NaN'

# ------------------------------------------------------------------------------

def summarize(summary, lane):
    summary['sample_id'] = lane['sample_id']
    summary['lane'] = lane['lane']

    # total counts
    summary['total_read_pairs'] = lane['total_read_pairs']
    summary['total_bps'] = lane['total_bps']
    summary['read_length'] = lane['read_length']
    
    # fastq stats, average per read
    summary['mean_base_qual_per_read'] = hist_mean(hist_sum(lane, 'average_base_qual_histogram'))
    summary['std_base_qual_per_read'] = hist_std(hist_sum(lane, 'average_base_qual_histogram'))
    summary['mean_N_per_read'] = hist_mean(hist_sum(lane, 'N_count_histogram'))
    summary['std_N_per_read'] = hist_std(hist_sum(lane, 'N_count_histogram'))
    summary['mean_GC_per_read'], summary['std_GC_per_read'] = stats_GC(lane)

    # fastq stats, average per position
    summary['mean_base_qual_per_position'] = mean(hist_avg(lane, 'average_base_qual_by_position'))
    summary['std_base_qual_per_position'] = std(hist_avg(lane, 'average_base_qual_by_position'))
    summary['mean_N_per_position'], summary['std_N_per_position'], total_N = stats_by_pos(lane['Ns_by_position_first'], lane['Ns_by_position_second'], lane['total_bps'])
    summary['mean_A_per_position'], summary['std_A_per_position'], total_A = stats_by_pos(lane['As_by_position_first'], lane['As_by_position_second'], lane['total_bps']-total_N)
    summary['mean_C_per_position'], summary['std_C_per_position'], total_C = stats_by_pos(lane['Cs_by_position_first'], lane['Cs_by_position_second'], lane['total_bps']-total_N)
    summary['mean_G_per_position'], summary['std_G_per_position'], total_G = stats_by_pos(lane['Gs_by_position_first'], lane['Gs_by_position_second'], lane['total_bps']-total_N)
    summary['mean_T_per_position'], summary['std_T_per_position'], total_T = stats_by_pos(lane['Ts_by_position_first'], lane['Ts_by_position_second'], lane['total_bps']-total_N)
    
    # k-mer stats
    x, summary['32_mer_error_rate'] = k_mer_error_rate(lane['distinct_32mer_count_after_qual_clipping_17'],
                                                       lane['unique_32mer_count_after_qual_clipping_17'],
                                                       lane['32mer_count_after_qual_clipping_17'],
                                                       32)
    summary['adapter_8_mers'] = adapter_8_mers(lane['8mer_count'])

    if lane['total_read_pairs'] > 0:
        # marked duplicate
        summary['marked_duplicate'] = lane['marked_duplicate'] / 2.0 / lane['total_read_pairs'] * 100 

        # flag stats
        summary['unmapped'] = (lane['first_read_unmapped'] + lane['second_read_unmapped']) / 2.0 / lane['total_read_pairs'] * 100
        summary['both_unmapped'] = float(lane['both_reads_unmapped']) / lane['total_read_pairs'] * 100
        summary['first_unmapped'] = float(lane['first_read_unmapped']) / lane['total_read_pairs'] * 100 - summary['both_unmapped']
        summary['second_unmapped'] = float(lane['second_read_unmapped']) / lane['total_read_pairs'] * 100 - summary['both_unmapped']
        summary['FF_RR_oriented_pairs'] = float(lane['FF_RR_oriented_pairs']) / lane['total_read_pairs'] * 100
        summary['proper_pairs'] = float(lane['total_proper_pairs']) / lane['total_read_pairs'] * 100
        summary['proper_pairs_autosome'] = float(lane['total_proper_pairs_autosome']) / lane['first_and_or_second_read_mapped'] * 100
    else:
        summary['marked_duplicate'], summary['unmapped'], summary['both_unmapped'] = [0, 0, 0]
        summary['first_unmapped'], summary['second_unmapped'] = [0, 0]
        summary['FF_RR_oriented_pairs'], summary['proper_pairs'], summary['proper_pairs_autosome'] = [0, 0,  0]

    # histogram stats
    summary['mean_coverage'] = hist_mean(lane['genome_coverage_histogram'])
    summary['std_coverage'] = hist_std(lane['genome_coverage_histogram'])
    summary['mean_insert_size'] = hist_mean(lane['insert_size_histogram'])
    summary['std_insert_size'] = hist_std(lane['insert_size_histogram'])
    if summary['std_insert_size'] != 'NaN' and summary['mean_insert_size'] != 'NaN':
        summary['std_times_mean_insert_size'] = summary['std_insert_size'] * summary['mean_insert_size']
    else:
        summary['std_times_mean_insert_size'] = 'NaN'
    summary['adapter_insert_size'] = adapter_insert_size(lane['insert_size_histogram'], lane['read_length'])

    # aligner behaviour
    if lane['total_read_pairs'] > 0:
        summary['mapping_qual_60'] = count_lt(hist_sum(lane, 'mapping_qual_histogram'), 60) / 2.0 / lane['total_read_pairs'] * 100
        summary['mapping_qual_40'] = count_lt(hist_sum(lane, 'mapping_qual_histogram'), 40) / 2.0 / lane['total_read_pairs'] * 100
        summary['mapping_qual_20'] = count_lt(hist_sum(lane, 'mapping_qual_histogram'), 20) / 2.0 / lane['total_read_pairs'] * 100
        summary['clipped_5_prime'] = (lane['soft_clipping_5_prime_by_position_first'][0] + lane['soft_clipping_5_prime_by_position_second'][0]) / 2.0 / lane['total_read_pairs'] * 100
        summary['clipped_3_prime'] = (lane['soft_clipping_3_prime_by_position_first'][-1] + lane['soft_clipping_3_prime_by_position_second'][-1]) / 2.0 / lane['total_read_pairs'] * 100
    else:
        summary['mapping_qual_60'], summary['mapping_qual_40'], summary['mapping_qual_20'] = [0, 0, 0]
        summary['clipped_5_prime'], summary['clipped_3_prime'] = [0, 0]
    if lane['total_bps'] > 0:
        summary['mean_mismatches'] = hist_mean(hist_sum(lane, 'mismatch_count_histogram'))
        summary['mean_deletions'] = hist_mean(hist_sum(lane, 'deletion_count_histogram'))
        summary['mean_insertions'] = hist_mean(hist_sum(lane, 'insertion_count_histogram'))
        summary['nz_deletions'] = hist_nz(hist_sum(lane, 'deletion_count_histogram'))
        summary['nz_insertions'] = hist_nz(hist_sum(lane, 'insertion_count_histogram'))
    else:
        summary['mean_mismatches'] = 0
        summary['mean_deletions'] = 0
        summary['mean_insertions'] = 0
        summary['nz_deletions'] = 0
        summary['nz_insertions'] = 0

    # triplet score
    compute_triplet_scores(summary, lane)

# ------------------------------------------------------------------------------

def threshold(summary, field):
    if field not in thresholds:
        return ''
        
    if field == '32_mer_error_rate':
        field += '_low_cov' if summary['mean_coverage'] < 15 else '_high_cov'

    if summary[field] == 'NA':
        return ''
    elif summary[field] == 'NaN':
        return ' ****'
    elif summary[field] < thresholds[field][min][0] or summary[field] > thresholds[field][max][0]:
        return ' ***'
    elif summary[field] < thresholds[field][min][1] or summary[field] > thresholds[field][max][1]:
        return ' **'
    elif summary[field] < thresholds[field][min][2] or summary[field] > thresholds[field][max][2]:
        return ' *'
    else:
        return ''

# ------------------------------------------------------------------------------

def get_flags(summary, strictness):
    flags = set()

    if (len(threshold(summary, 'mean_base_qual_per_read')) > strictness or
        len(threshold(summary, 'std_base_qual_per_read')) > strictness or
        len(threshold(summary, 'mean_base_qual_per_position')) > strictness or
        len(threshold(summary, 'std_base_qual_per_position')) > strictness):
        flags.add('Q')
    if (len(threshold(summary, 'mean_N_per_read')) > strictness or
        len(threshold(summary, 'std_N_per_read')) > strictness or
        len(threshold(summary, 'mean_N_per_position')) > strictness or
        len(threshold(summary, 'std_N_per_position')) > strictness):
        flags.add('N')
    if (len(threshold(summary, 'mean_GC_per_read')) > strictness or
        len(threshold(summary, 'std_GC_per_read')) > strictness):
        flags.add('G')
    if (len(threshold(summary, 'mean_A_per_position')) > strictness or
        len(threshold(summary, 'std_A_per_position')) > strictness or
        len(threshold(summary, 'mean_C_per_position')) > strictness or
        len(threshold(summary, 'std_C_per_position')) > strictness or
        len(threshold(summary, 'mean_G_per_position')) > strictness or
        len(threshold(summary, 'std_G_per_position')) > strictness or
        len(threshold(summary, 'mean_T_per_position')) > strictness or
        len(threshold(summary, 'std_T_per_position')) > strictness):
        flags.add('B')
    if len(threshold(summary, '32_mer_error_rate')) > strictness:
        flags.add('K')
    if len(threshold(summary, 'marked_duplicate')) > strictness:
        flags.add('D')
    if (len(threshold(summary, 'unmapped')) > strictness or
        len(threshold(summary, 'both_unmapped')) > strictness or
        len(threshold(summary, 'first_unmapped')) > strictness or
        len(threshold(summary, 'second_unmapped')) > strictness):
        flags.add('U')
    if len(threshold(summary, 'proper_pairs_autosome')) > strictness:
        flags.add('P')
    if len(threshold(summary, 'FF_RR_oriented_pairs')) > strictness:
        flags.add('o')
    if (len(threshold(summary, 'mean_coverage')) > strictness or
        len(threshold(summary, 'total_bps')) > strictness or 
        len(threshold(summary, 'std_coverage')) > strictness):
        flags.add('C')
    if (len(threshold(summary, 'mean_insert_size')) > strictness or
        len(threshold(summary, 'std_insert_size')) > strictness):
        flags.add('I')
    if (len(threshold(summary, 'adapter_8mers')) > strictness or
        len(threshold(summary, 'adapter_insert_size')) > strictness):
        flags.add('A')
    if (len(threshold(summary, 'mapping_qual_60')) > strictness or
        len(threshold(summary, 'mapping_qual_40')) > strictness or
        len(threshold(summary, 'mapping_qual_20')) > strictness):
        flags.add('M')
    if len(threshold(summary, 'mean_mismatches')) > strictness:
        flags.add('m')
    if len(threshold(summary, 'nz_deletions')) > strictness:
        flags.add('d')
    if len(threshold(summary, 'nz_insertions')) > strictness:
        flags.add('i')
    if (len(threshold(summary, 'clipped_5_prime')) > strictness or
        len(threshold(summary, 'clipped_3_prime')) > strictness):
        flags.add('c')
    if (len(threshold(summary, 'CtoA')) > strictness or
        len(threshold(summary, 'GtoA')) > strictness or
        len(threshold(summary, 'TtoA')) > strictness or
        len(threshold(summary, 'AtoC')) > strictness or
        len(threshold(summary, 'GtoC')) > strictness or
        len(threshold(summary, 'TtoC')) > strictness):
        flags.add('O')

    if len(flags) == 0:
        flags.add('-')

    return flags

# ------------------------------------------------------------------------------

def write_header():
    print(('SAMPLE_ID' + '\t' + 'LANE' + '\t' + 'FAILURE_FLAGS'+ '\t' + 'JOINT_CALLING_FLAGS'+ '\t' + 'STRICT_FLAGS'), end=' ')
    print(('\t' + 'TOTAL_BPS' + '\t' + 'TOTAL_READ_PAIRS' + '\t' + 'READ_LENGTH'), end=' ')

    # fastq stats, average per read
    print(('\t' + 'MEAN_BASE_QUAL_PER_READ' + '\t' + 'STD_BASE_QUAL_PER_READ'), end=' ')
    print(('\t' + 'MEAN_N_COUNT_PER_READ' + '\t' + 'STD_N_COUNT_PER_READ'), end=' ')
    print(('\t' + 'MEAN_GC_CONTENT_PER_READ' + '\t' + 'STD_GC_CONTENT_PER_READ'), end=' ')

    # fastq stats, average per position
    print(('\t' + 'MEAN_BASE_QUAL_PER_POSITION' + '\t' + 'STD_BASE_QUAL_PER_POSITION'), end=' ')
    print(('\t' + 'MEAN_N_PER_POSITION' + '\t' + 'STD_N_PER_POSITION'), end=' ')
    print(('\t' + 'MEAN_A_PER_POSITION' + '\t' + 'STD_A_PER_POSITION'), end=' ')
    print(('\t' + 'MEAN_C_PER_POSITION' + '\t' + 'STD_C_PER_POSITION'), end=' ')
    print(('\t' + 'MEAN_G_PER_POSITION' + '\t' + 'STD_G_PER_POSITION'), end=' ')
    print(('\t' + 'MEAN_T_PER_POSITION' + '\t' + 'STD_T_PER_POSITION'), end=' ')
    
    # misc stats
    print(('\t' + '32_MER_ERROR_RATE'), end=' ')
    print(('\t' + 'ADAPTER_8_MERS'), end=' ')
    print(('\t' + 'MARKED_DUPLICATE'), end=' ')

    # alignment stats
    print(('\t' + 'UNMAPPED'), end=' ')
    print(('\t' + 'BOTH_UNMAPPED'), end=' ')
    print(('\t' + 'FIRST_UNMAPPED'), end=' ')
    print(('\t' + 'SECOND_UNMAPPED'), end=' ')
    print(('\t' + 'PROPER_PAIRS'), end=' ')
    print(('\t' + 'PROPER_PAIRS_AUTOSOME'), end=' ')
    print(('\t' + 'FF_RR_PAIRS'), end=' ')
    print(('\t' + 'MEAN_COVERAGE' + '\t' + 'STD_COVERAGE'), end=' ')
    print(('\t' + 'MEAN_INSERT_SIZE' + '\t' + 'STD_INSERT_SIZE' + '\t' + 'ADAPTER_INSERT_SIZE'), end=' ')
    print(('\t' + 'MAPPING_QUAL_60' + '\t' + 'MAPPING_QUAL_40' + '\t' + 'MAPPING_QUAL_20'), end=' ')
    print(('\t' + 'MEAN_MISMATCHES' + '\t' + 'MEAN_DELETIONS' + '\t' + 'MEAN_INSERTIONS'), end=' ')
    print(('\t' + 'NZ_DELETIONS' + '\t' + 'NZ_INSERTIONS'), end=' ')
    print(('\t' + 'CLIPPED_5_PRIME' + '\t' + 'CLIPPED_3_PRIME'), end=' ')

    # N>A and N>C scores
    print(('\t' + 'C>A' + '\t' + 'G>A' + '\t' + 'T>A'), end=' ')
    print(('\t' + 'A>C' + '\t' + 'G>C' + '\t' + 'T>C'), end=' ')

    print()

# ------------------------------------------------------------------------------

def write_line(summary):
    print((summary['sample_id'] + '\t' + summary['lane'] + '\t' + ''.join(summary['flags3']) +  '\t' + ''.join(summary['flags2'])+ '\t' + ''.join(summary['flags1'])), end=' ')
    print(('\t' + str(summary['total_bps']) + '\t' + str(summary['total_read_pairs']) + '\t' + str(summary['read_length'])), end=' ')

    # read stats, average per read
    print(('\t' + str(summary['mean_base_qual_per_read']) + '\t' + str(summary['std_base_qual_per_read'])), end=' ')
    print(('\t' + str(summary['mean_N_per_read']) + '\t' + str(summary['std_N_per_read'])), end=' ')
    print(('\t' + str(summary['mean_GC_per_read']) + '\t' + str(summary['std_GC_per_read'])), end=' ')

    # read stats, average per position
    print(('\t' + str(summary['mean_base_qual_per_position']) + '\t' + str(summary['std_base_qual_per_position'])), end=' ')
    print(('\t' + str(summary['mean_N_per_position']) + '\t' + str(summary['std_N_per_position'])), end=' ')
    print(('\t' + str(summary['mean_A_per_position']) + '\t' + str(summary['std_A_per_position'])), end=' ')
    print(('\t' + str(summary['mean_C_per_position']) + '\t' + str(summary['std_C_per_position'])), end=' ')
    print(('\t' + str(summary['mean_G_per_position']) + '\t' + str(summary['std_G_per_position'])), end=' ')
    print(('\t' + str(summary['mean_T_per_position']) + '\t' + str(summary['std_T_per_position'])), end=' ')

    # misc stats    
    print(('\t' + str(summary['32_mer_error_rate'])), end=' ')
    print(('\t' + str(summary['adapter_8_mers'])), end=' ')
    print(('\t' + str(summary['marked_duplicate'])), end=' ')

    # alignment stats
    print(('\t' + str(summary['unmapped'])), end=' ')
    print(('\t' + str(summary['both_unmapped']) + '\t' + str(summary['first_unmapped']) + '\t' + str(summary['second_unmapped'])), end=' ')
    print(('\t' + str(summary['proper_pairs']) + '\t' + str(summary['proper_pairs_autosome']) + '\t' + str(summary['FF_RR_oriented_pairs'])), end=' ')
    print(('\t' + str(summary['mean_coverage']) + '\t' + str(summary['std_coverage'])), end=' ')
    print(('\t' + str(summary['mean_insert_size']) + '\t' + str(summary['std_insert_size'])), end=' ')
    print(('\t' + str(summary['adapter_insert_size'])), end=' ')
    print(('\t' + str(summary['mapping_qual_60']) + '\t' + str(summary['mapping_qual_40']) + '\t' + str(summary['mapping_qual_20'])), end=' ')
    print(('\t' + str(summary['mean_mismatches']) + '\t' + str(summary['mean_deletions']) + '\t' + str(summary['mean_insertions'])), end=' ')
    print(('\t' + str(summary['nz_deletions']) + '\t' + str(summary['nz_insertions'])), end=' ')
    print(('\t' + str(summary['clipped_5_prime']) + '\t' + str(summary['clipped_3_prime'])), end=' ')

    # N>A and N>C scores
    print(('\t' + str(summary['CtoA']) + '\t' + str(summary['GtoA']) + '\t' + str(summary['TtoA'])), end=' ')
    print(('\t' + str(summary['AtoC']) + '\t' + str(summary['GtoC']) + '\t' + str(summary['TtoC'])), end=' ')

    print()

def write_txt(summary):
    print('# BamQC summary for')
    print(('Sample ID: ' + summary['sample_id']))
    print(('Lane: ' + summary['lane']))
    print()
    print('# Total counts')
    print(('Total basepairs: ' + '{:,}'.format(summary['total_bps'])))
    print(('Total read pairs: ' + '{:,}'.format(summary['total_read_pairs'])))
    print(('Read length: ' + str(summary['read_length'])))
    print()
    print('Labels: *** fails quality check, ** fails strict quality check, * outlier')
    print()
    print('# Read stats (per read)')
    print(('Mean of mean base calling quality:      ' + '{:8.4f}'.format(summary['mean_base_qual_per_read']) + threshold(summary, 'mean_base_qual_per_read')))
    print(('Std dev of mean base calling quality:   ' + '{:8.4f}'.format(summary['std_base_qual_per_read']) + threshold(summary, 'std_base_qual_per_read')))
    print(('Mean percent N:                         ' + '{:8.4f}'.format(summary['mean_N_per_read']) + threshold(summary, 'mean_N_per_read')))
    print(('Std dev of percent N:                   ' + '{:8.4f}'.format(summary['std_N_per_read']) + threshold(summary, 'std_N_per_read')))
    print(('Mean percent of GC bases:               ' + '{:8.4f}'.format(summary['mean_GC_per_read']) + threshold(summary, 'mean_GC_per_read')))
    print(('Std dev of percent GC:                  ' + '{:8.4f}'.format(summary['std_GC_per_read']) + threshold(summary, 'std_GC_per_read')))
    print()
    print('# Read stats (per read position)')
    print(('Mean of mean base calling quality:      ' + '{:8.4f}'.format(summary['mean_base_qual_per_position']) + threshold(summary, 'mean_base_qual_per_position')))
    print(('Std dev of mean base calling quality:   ' + '{:8.4f}'.format(summary['std_base_qual_per_position']) + threshold(summary, 'std_base_qual_per_position')))
    print(('Mean percent N:                         ' + '{:8.4f}'.format(summary['mean_N_per_position']) + threshold(summary, 'mean_N_per_position')))
    print(('Standard deviation of percent N:        ' + '{:8.4f}'.format(summary['std_N_per_position']) + threshold(summary, 'std_N_per_position')))
    print(('Mean percent A:                         ' + '{:8.4f}'.format(summary['mean_A_per_position']) + threshold(summary, 'mean_A_per_position')))
    print(('Standard deviation of percent A:        ' + '{:8.4f}'.format(summary['std_A_per_position']) + threshold(summary, 'std_A_per_position')))
    print(('Mean percent C:                         ' + '{:8.4f}'.format(summary['mean_C_per_position']) + threshold(summary, 'mean_C_per_position')))
    print(('Standard deviation of percent C:        ' + '{:8.4f}'.format(summary['std_C_per_position']) + threshold(summary, 'std_C_per_position')))
    print(('Mean percent G:                         ' + '{:8.4f}'.format(summary['mean_G_per_position']) + threshold(summary, 'mean_G_per_position')))
    print(('Standard deviation of percent G:        ' + '{:8.4f}'.format(summary['std_G_per_position']) + threshold(summary, 'std_G_per_position')))
    print(('Mean percent T:                         ' + '{:8.4f}'.format(summary['mean_T_per_position']) + threshold(summary, 'mean_T_per_position')))
    print(('Standard deviation of percent T:        ' + '{:8.4f}'.format(summary['std_T_per_position']) + threshold(summary, 'std_T_per_position')))
    print()
    print('# K-mer statistics')
    print(('Estimated 32-mer error rate:            ' + '{:8.4f}'.format(summary['32_mer_error_rate']) + threshold(summary, '32_mer_error_rate')))
    print(('Percent of Universal adapter 8-mers:    ' + '{:8.4f}'.format(summary['adapter_8_mers']) + threshold(summary, 'adapter_8_mers')))
    print()
    print('# Duplicates, proper pairs and orientation')
    print(('Percent marked as duplicate:            ' + '{:8.4f}'.format(summary['marked_duplicate']) + threshold(summary, 'marked_duplicate')))
    print(('Percent proper pairs:                   ' + '{:8.4f}'.format(summary['proper_pairs']) + threshold(summary, 'proper_pairs')))
    print(('Percent proper pairs autosome:                   ' + '{:8.4f}'.format(summary['proper_pairs_autosome']) + threshold(summary, 'proper_pairs_autosome')))
    print(('Percent FF/RR oriented pairs:           ' + '{:8.4f}'.format(summary['FF_RR_oriented_pairs']) + threshold(summary, 'FF_RR_oriented_pairs')))
    print()
    print('# Flag stats')
    print(('Percent unmapped reads:                 ' + '{:8.4f}'.format(summary['unmapped']) + threshold(summary, 'unmapped')))
    print(('Percent both reads in pair unmapped:    ' + '{:8.4f}'.format(summary['both_unmapped']) + threshold(summary, 'both_unmapped')))
    print(('Percent only first unmapped in pair:    ' + '{:8.4f}'.format(summary['first_unmapped']) + threshold(summary, 'first_unmapped')))
    print(('Percent only second unmapped in pair:   ' + '{:8.4f}'.format(summary['second_unmapped']) + threshold(summary, 'second_unmapped')))
    print()
    print('# Coverage and insert size histogram stats')
    print(('Mean coverage:                          ' + '{:8.4f}'.format(summary['mean_coverage']) + threshold(summary, 'mean_coverage')))
    print(('Std dev of coverage:                    ' + '{:8.4f}'.format(summary['std_coverage']) + threshold(summary, 'std_coverage')))
    print(('Mean insert size:                       ' + '{:8.4f}'.format(summary['mean_insert_size']) + threshold(summary, 'mean_insert_size')))
    print(('Std dev of insert size:                 ' + '{:8.4f}'.format(summary['std_insert_size']) + threshold(summary, 'std_times_mean_insert_size')))
    print(('Percent insert size < read length:      ' + '{:8.4f}'.format(summary['adapter_insert_size']) + threshold(summary, 'adapter_insert_size')))
    print()
    print('# Alignment stats')
    print(('Percent reads with mapping quality < 60:' + '{:8.4f}'.format(summary['mapping_qual_60']) + threshold(summary, 'mapping_qual_60')))
    print(('Percent reads with mapping quality < 40:' + '{:8.4f}'.format(summary['mapping_qual_40']) + threshold(summary, 'mapping_qual_40')))
    print(('Percent reads with mapping quality < 20:' + '{:8.4f}'.format(summary['mapping_qual_20']) + threshold(summary, 'mapping_qual_20')))
    print(('Mean mismatches per read pair:          ' + '{:8.4f}'.format(summary['mean_mismatches']) + threshold(summary, 'mean_mismatches')))
    print(('Mean deletions per read pair:           ' + '{:8.4f}'.format(summary['mean_deletions']) + threshold(summary, 'mean_deletions')))
    print(('Mean insertions per read pair:          ' + '{:8.4f}'.format(summary['mean_insertions']) + threshold(summary, 'mean_insertions')))
    print(('Fraction of read pairs with deletions:           ' + '{:8.4f}'.format(summary['nz_deletions']) + threshold(summary, 'mean_deletions')))
    print(('Fraction of read pairs with insertions:          ' + '{:8.4f}'.format(summary['nz_insertions']) + threshold(summary, 'mean_insertions')))
    print(('Percent of reads clipped at 5\'-end:     ' + '{:8.4f}'.format(summary['clipped_5_prime']) + threshold(summary, 'clipped_5_prime')))
    print(('Percent of reads clipped at 3\'-end:     ' + '{:8.4f}'.format(summary['clipped_3_prime']) + threshold(summary, 'clipped_3_prime')))
    print()
    print('# Fraction of 1st-FW + 2nd-RC reads for mismatches (expect 0.5)')
    print(('C>A fraction:                           ' + '{:8.4f}'.format(summary['CtoA']) + threshold(summary, 'CtoA')))
    print(('G>A fraction:                           ' + '{:8.4f}'.format(summary['GtoA']) + threshold(summary, 'GtoA')))
    print(('T>A fraction:                           ' + '{:8.4f}'.format(summary['TtoA']) + threshold(summary, 'TtoA')))
    print(('A>C fraction:                           ' + '{:8.4f}'.format(summary['AtoC']) + threshold(summary, 'AtoC')))
    print(('G>C fraction:                           ' + '{:8.4f}'.format(summary['GtoC']) + threshold(summary, 'CtoC')))
    print(('T>C fraction:                           ' + '{:8.4f}'.format(summary['TtoC']) + threshold(summary, 'TtoC')))

# ------------------------------------------------------------------------------

if __name__ == '__main__':

    # parse command line
    args = parse_arguments()

    # read the input data
    data = []
    read_bamqc_output(data, args.bamqc_file)
    
    if args.title == True and args.long != True:
        write_header()

    for lane in data:
        # compute summary statistics
        summary = {}
        summarize(summary, lane)

        # print summary statistics
        if args.long != True:
            summary['flags3'] = get_flags(summary, 3)
            summary['flags2'] = get_flags(summary, 2)
            summary['flags1'] = get_flags(summary, 1)
            write_line(summary)
        else:
            write_txt(summary)

