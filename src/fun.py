import re
from Bio import SeqIO
import pysam
import parasail
import math
from itertools import product
from array import array

# Define a function to parse the motifs from the population allele file
def parse_pop_alleles(file_path, read_length):
    # Initialize an empty list to store the parsed motifs
    motifs_db = []
    
    # Initialize an empty set to store unique motif names
    motifs_names = set()
    
    # Initialize an empty set to store unique chromosomes
    chromosomes = set()
    
    # Open the file at the specified file path in read mode
    with open(file_path, 'r') as file:

        # Initialize variables to store the parsed data
        old_name = ""
        name = None
        old_chromosome = None
        chromosome = None
        old_start = None
        old_end = None
        start = None
        end = None
        lengths = array('h') 
        frequencies = []
        patterns = []

        # Iterate over each line in the file
        for line in file:
            # Check if the line is not empty or just whitespace
            if line.strip():
                # Split the line into frequency, name, and location using regex (split by whitespace)
                frequency, name, location = re.split(r'\s+', line.strip(), maxsplit=2)
                
                # Use regex to match and capture chromosome, start, end, and repeat patterns from the location
                match = re.match(r'chr([\dXY]+):g\.(\d+)_(\d+)(.*)', location)
                if match:
                    # Extract the captured groups from the regex match
                    chromosome, start, end, repeats = match.groups()
                    
                    # Find all repeat patterns and their counts within the repeats string
                    repeat_patterns = re.findall(r'([A-Z]+)\[([A-Z\d]+)\]', repeats)

                    # Convert repeat patterns to integers, and replace 'E' with 151 and remove rows containing 'B'
                    b = False
                    e = False
                    e_idx = -1
                    
                    # Initialize an empty list to store the parsed repeat patterns
                    repeats = []

                    length = 0

                    # Iterate over each repeat pattern
                    i = 0
                    for repeat in repeat_patterns:
                        if repeat[1] == 'E':
                            e = True
                            e_idx = i
                            repeats.append([repeat[0], -1])
                        elif repeat[1] == 'B':
                            b = True
                        else:
                            repeats.append([repeat[0], int(repeat[1])])
                            length += len(repeat[0]) * int(repeat[1])
                        i += 1
                    if b:
                        continue
                    if e:
                        count = int(math.ceil(read_length - length) / len(repeats[e_idx][0])) + 1
                        repeats[e_idx][1] = count
                        length += len(repeats[e_idx][0]) * count
                    
                    # Handle the first line
                    if old_name == "":
                        old_name = name
                        old_chromosome = chromosome
                        old_start = start
                        old_end = end
                    
                    if old_name != name:
                        # Append a dictionary with the parsed data to the motifs list
                        motifs_db_append(motifs_db, motifs_names, chromosomes, frequencies, old_name, old_chromosome, old_start, old_end, patterns, lengths)

                        # Update the variables for the next iteration
                        old_name = name
                        old_chromosome = chromosome
                        old_start = start
                        old_end = end
                        patterns = []
                        frequencies = []
                        lengths = array('h') 
                        patterns.append(repeats)
                        frequencies.append(float(frequency))
                        lengths.append(length)
                    else:
                        patterns.append(repeats)
                        frequencies.append(float(frequency))
                        lengths.append(length)
        # Append the last motif to the motifs list
        motifs_db_append(motifs_db, motifs_names, chromosomes, frequencies, name, chromosome, start, end, patterns, lengths)

    return motifs_db, chromosomes

# Define a function to append a motif to the motifs list
def motifs_db_append(motifs_db, motifs_names, chromosomes, frequencies, name, chromosome, start, end, patterns, lengths):
    
    pattern, repeats = process_patterns(patterns)

    motifs_db.append({
                        'name': name,
                        'chromosome': chromosome,
                        'start': int(start)-1,
                        'end': int(end),
                        'patterns': pattern,
                        'repeats': repeats,
                        'frequencies': frequencies,
                        'left_flang': None,
                        'right_flang': None,
                        'lengths': lengths,
                    })
    motifs_names.add(name)
    chromosomes.add(chromosome)

def process_patterns(patterns):

    number_of_motifs = len(patterns)
    number_of_patterns = len(patterns[0])

    pattern = []
    repeat = [array('h') for _ in range(number_of_patterns)]

    for i in range(number_of_patterns):
        pattern.append(patterns[0][i][0])

        for j in range(number_of_motifs):
            repeat[i].append(patterns[j][i][1])

    return pattern, repeat

def generate_tuples(n, bounds):
    # Create ranges for each xi
    ranges = [range(1, x + 1) for x in bounds]
    
    # Generate all possible tuples using itertools.product
    tuples = list(product(*ranges))
    
    return tuples

def pseudo_counts_alleles(m, read_length):
    
    patterns = m['patterns']
    repeats = m['repeats']
    frequencies = m['frequencies']

    population = []
    for r in repeats:
        population.append(set(r))
    
    pattern_bound = []
    pattern_length = []
    for p in patterns:
        l = len(p)
        pattern_length.append(l)
        pattern_bound.append(math.ceil(read_length / l))

    number_of_patterns = len(patterns)
    tuples = generate_tuples(number_of_patterns, pattern_bound)
    
    ix = 0
    for t in tuples:

        ix += 1
        if ix % 20: continue

        variant = array('h')
        same = 0
        length = 0
        for i in range(number_of_patterns):
            length += t[i] * pattern_length[i]
            if t[i] in population[i]:
                same += 1
            variant.append(t[i])
        if length > read_length or same == number_of_patterns:
            continue

        for i in range(number_of_patterns):
            m['repeats'][i].append(variant[i])
        m['frequencies'].append(1)
        m['lengths'].append(length)

# Define a function to parse the reference genome
def parse_genome(file_path, valid_chromosomes):
    # Initialize an empty dictionary to store the reference genome
    reference_genome = {}
    
    for record in SeqIO.parse(file_path, "fasta"):
        chromosome = record.id.split()[0][3:] 
        if chromosome in valid_chromosomes:
            reference_genome[chromosome] = str(record.seq)

    return reference_genome

# Define a function to fetch the flangs for a motif
def fetch_flangs(motifs_db, flang_length, chromosomes, file_path):
    # Parse the reference genome
    reference_genome = parse_genome(file_path, chromosomes)

    for motif in motifs_db:
        # Extract the motif data
        chromosome = motif['chromosome']
        start = motif['start']
        end = motif['end']
        
        # Get the sequence from the reference genome
        chromosome_sequence = reference_genome[chromosome]
        
        # Extract 1000 bases before the start position
        motif['start'] = max(0, start - flang_length)
        motif['left_flang'] = chromosome_sequence[motif['start']:start]
        
        # Extract 3 bases after the end position
        motif['end'] = min(len(chromosome_sequence), end + flang_length)
        motif['right_flang'] = chromosome_sequence[end:motif['end']]
    
# Define a function to convert Phred quality scores to probabilities
def phred_to_prob(quality_scores):
    # Apply the transformation formula to each score in the list
    return [10**(-q / 10) for q in quality_scores]

# Fetching relevant reads for the motifs
def filter_relevant_reads(motif, bam_file_path, read_length):
    relevant_reads = []

    bam_file = pysam.AlignmentFile(bam_file_path, "rb")

    chr = "chr" + motif['chromosome']
    start = motif['start'] + read_length
    end = motif['end'] - read_length

    for read in bam_file.fetch(chr, start, end):
        read_info = {
            'sequence': read.query_sequence,
            'quality_scores': phred_to_prob(read.query_qualities)
        }
        relevant_reads.append(read_info)
    
    return relevant_reads

# Create the motif sequence by concatenating the left flang, patterns, and right flang
def create_motif_sequence(motif, pattern_index):

    seq = motif['left_flang']

    number_of_patterns = len(motif['patterns'])
    for i in range(number_of_patterns):
        seq += motif['patterns'][i] * motif['repeats'][i][pattern_index]
    
    seq += motif['right_flang']

    return seq

# function to mapp the reads to the motifs
def mapping_reads_to_motifs(motif, motifs_num, relevant_reads):
    mapping = [[[] for _ in range(len(relevant_reads))] for _ in range(motifs_num)]
    mapped_indices = [None for _ in range(motifs_num)]

    k = 0
    for i in range(motifs_num):
        indices = set()
        ref_seq = create_motif_sequence(motif, i)
        j = 0
        for r in relevant_reads:
            query_seq = r['sequence']
            res = parasail.sg_dx_trace(query_seq, ref_seq, 2, 1, parasail.blosum62)
            result, flags, s, e = extraction(res, indices)
            mapping[i][j] = (result, flags, s, e)
            j += 1

        mapped_indices[i] = sorted(indices)

    return mapping

# Get the mapping indeces
def extraction(mapping, mapped_indices):
    indices = []
    flags = array('h')

    comp = mapping.traceback.comp
    ref = mapping.traceback.ref
    query = mapping.traceback.query

    start = 0
    while query[start] == '-':
        start += 1

    end = len(query) - 1
    while query[end] == '-':
        end -= 1
    end += 1

    r_idx = -1
    m_idx = start - 1
    begin = -1
    finish = -1
    for i in range(start, end):
        cmp = 0
        if comp[i] != ' ':
            cmp = 1

        if query[i] != '-':
            r_idx += 1
        
        if ref[i] != '-':
            m_idx += 1
            if query[i] != '-':
                mapped_indices.add(m_idx)
                flags.append(r_idx)
            else:
                flags.append(-1)
            if begin == -1: begin = m_idx
            finish = m_idx

        indices.append(array('h', [r_idx, m_idx, cmp]))

    return indices, flags, begin, finish


# Probabilities calculation
# E step
def E_step(s_given_r_matrix, motif, read, mapping, pr_base_given_position_motif, motif_count, reads_count, b_g_b): 

    r_given_s_matrix = [[0 for _ in range(motif_count)] for _ in range(reads_count)]
    denominators = [0 for _ in range(reads_count)]

    for i in range(reads_count):
        sum = 0
        for j in range(motif_count):
            tmp = pr_r_given_s(motif, read, j, i, mapping, pr_base_given_position_motif, b_g_b) 
            r_given_s_matrix[i][j] = tmp
            sum += tmp * motif['frequencies'][j]
        denominators[i] = sum

    for j in range(reads_count):
        denominator = denominators[j]
        for i in range(motif_count): 
            nominator = r_given_s_matrix[j][i] * motif['frequencies'][i]

            if denominator == 0:
                s_given_r_matrix[i][j] = 0
            else:
                s_given_r_matrix[i][j] = nominator / denominator    

# Calculate the probability of a read given a motif by multiplying the probabilities of the all mapped bases
def pr_r_given_s(motif, read, motif_index, read_index, mapping, pr_base_given_position_motif, b_g_b):
    prod = 1

    interval = mapping[motif_index][read_index][0]
    for i in interval:
        if i[2] == 0:
            prod *= 0.05
        else:
            seq = create_motif_sequence(motif, motif_index)
            if seq == read[read_index]['sequence'][i[0]]:
                prod *= pr_base_given_motif(motif_index, read_index, i[0], i[1], pr_base_given_position_motif, b_g_b, read[read_index]['sequence'][i[0]])
            else:
                prod *= pr_base_given_motif(motif_index, read_index, i[0], i[1], pr_base_given_position_motif, b_g_b, "Q")

    return prod

# Calculate the probability of a base given a motif by summing the probabilities of the base given each of the four bases
def pr_base_given_motif(motif_index, read_index, base_index, m_base_index, pr_base_given_position_motif, b_g_b, flag): 
    sum = 0

    for i in range(4):
        sum += b_g_b[read_index][base_index][i] * pr_base_given_position_motif[motif_index][m_base_index][i]

    return sum

# Calculate the probability of a base given another base by checking if they are the same or not
def pr_base_given_base(reads, read_index, base_index, idx):
    bases = ["A", "C", "G", "T"]

    base = reads[read_index]['sequence'][base_index]
    quality = reads[read_index]['quality_scores'][base_index]

    if base == bases[idx]:
        return 1 - quality
    else:
        return quality / 3
    
def pre_pr_base_given_base(reads, read_length, reads_count):

    pr = [[[0 for _ in range(4)] for _ in range(read_length)] for _ in range(reads_count)]

    for i in range(reads_count):
        for j in range(read_length):
            for k in range(4):
                pr[i][j][k] = pr_base_given_base(reads, i, j, k)

    return pr

def pre_pr_base_gvien_motif(pr_base_given_position_motif, motif_count, motif, reads_count, mapping, relevant_reads):
    for i in range(motif_count):
            for j in range(motif['lengths'][i]):

                counts = [0, 0, 0, 0]
                sum = 0
                for l in range(reads_count):
                    start = mapping[i][l][2]
                    finish = mapping[i][l][3]
                    if j >= start and j <= finish:
                        pos = mapping[i][l][1][j - start]
                        if pos != -1:
                            sum += 1
                            if relevant_reads[l]['sequence'][pos] == "A":
                                counts[0] += 1
                            if relevant_reads[l]['sequence'][pos] == "C":
                                counts[1] += 1
                            if relevant_reads[l]['sequence'][pos] == "G":
                                counts[2] += 1
                            if relevant_reads[l]['sequence'][pos] == "T":
                                counts[3] += 1

                if sum == 0: counts = [0, 0, 0, 0] 
                else:
                    for k in range(4):
                        counts[k] /= sum

                pr_base_given_position_motif[i][j] = counts


                                
        














