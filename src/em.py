from fun import *
import os
import time
import shutil

def EM(output_folder, pseudo_counts, flank_length, read_length, alleles_file, reference_genome_file, reference_alignment_file):

    # Global variables
    output_folder = output_folder  + "frequencies/"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)

    init_start = time.time()

    # Parse the population alleles
    # Create a list of dictionaries containing the parsed motifs
    # with the moif structure: {'name': str, 'chromosome': int, 'start': int, 'end': int, 'patterns': dict, 'frequencies': array, 'left_flang': str, 'right_flang': str}
    # and set of motifs names and chromosomes
    # [start, end) - 0-based indexing to reference genome
    motifs_db, chromosomes = parse_pop_alleles(alleles_file, read_length)

    # Get the flangs
    fetch_flangs(motifs_db, flank_length, chromosomes, reference_genome_file)

    init_end = time.time()

    print(f"Init Time: {init_end - init_start}\n")

    # For each motif in the motifs_db
    m_idx = 0
    for motif in motifs_db:
        
        # if motif['name'] != "HD":                           # Testing purposes
        #     continue

        print(m_idx, motif['name'])
        
        parse_start = time.time()

        if pseudo_counts: pseudo_counts_alleles(motif, read_length)

        # print(f"Number of variants {len(motif['repeats'][0])}")
        # print(f"Pseudo counts done {time.time() - parse_start}")

        # Normalize the frequencies
        sum = 0
        for f in motif['frequencies']:
            sum += f

        motif_count = len(motif['frequencies'])
        for i in range(motif_count):
            motif['frequencies'][i] = motif['frequencies'][i] / sum
            motif['lengths'][i] += 2 * flank_length

        # print(f"Motif frequencies normalized {time.time() - parse_start}")

        # Get all the relevant reads
        # list of relevant reads for each motif
        relevant_reads = filter_relevant_reads(motif, reference_alignment_file, read_length)
        reads_count = len(relevant_reads)

        # print(f"Relevant reads count: {reads_count}")
        # print(f"Relevant reads filtered {time.time() - parse_start}")

        # s_given_r_matrix - probability of a motif given a read initialized with the motif frequencies
        s_given_r_matrix = [[0 for _ in range(reads_count)] for _ in range(motif_count)]
        for i in range(motif_count):
            for j in range(reads_count):
                s_given_r_matrix[i][j] = motif['frequencies'][i]

        # print(f"SR Matrix initialized {time.time() - parse_start}")

        # Precompution of the probability of a base given a base
        b_g_b = pre_pr_base_given_base(relevant_reads, read_length, reads_count)

        # print(f"Precomputation of the base given base done {time.time() - parse_start}")

        # Mapping for each motif and all its relevant reads
        mapping = mapping_reads_to_motifs(motif, motif_count, relevant_reads)

        # print(f"Mapping done {time.time() - parse_start}")
        
       # Calculate the probability of a base given a position and a motif
        pr_base_given_position_motif = [[None for _ in range(motif['lengths'][i])] for i in range(motif_count)]
        pre_pr_base_gvien_motif(pr_base_given_position_motif, motif_count, motif, reads_count, mapping, relevant_reads)

        # print(f"Precomputation of the base given position and motif done {time.time() - parse_start}")

        # Create the output folder for the motif
        name = motif['name']
        folder = output_folder + name
        # if os.path.exists(folder):
        #     shutil.rmtree(folder)
        os.makedirs(folder)

        parse_end = time.time()

        print(f"Parse Time: {parse_end - parse_start}")

        # EM Algorithm
        number_of_iterations = 0
        while True:
            # print(f"Iteration: {number_of_iterations}")

            E_step(s_given_r_matrix, motif, relevant_reads, mapping, pr_base_given_position_motif, motif_count, reads_count, b_g_b)

            # M Step
            # Chceck if the motifs frequencies have changed, if not break the loop
            pr_s = []
            for i in range(motif_count):
                sum = 0
                for j in range(reads_count):
                    sum += s_given_r_matrix[i][j]
                pr_s.append(sum / reads_count)


            diff = 0
            for i in range(motif_count):
                diff = max(diff, abs(motif['frequencies'][i] - pr_s[i]))
                motif['frequencies'][i] = pr_s[i]
                
            if diff < 0.001:
                break

            # Save the results of the iteration
            number_of_iterations += 1

        print(f"Number of iterations: {number_of_iterations}\n")

        THRESHOLD = 1e-4
        motif['frequencies'] = [0 if freq < THRESHOLD else round(freq, 4) for freq in motif['frequencies']]

        # Create the output file with the results for the motif with the name
        file = folder + "/" + name + ".txt"

        with open(file, "w") as f:
            for i in range(motif_count):
                p = ""
                number_of_patterns = len(motif['patterns'])
                for j in range(number_of_patterns):
                    p += f"{motif['patterns'][j]}[{motif['repeats'][j][i]}]"

                f.write(f"Frequency: {motif['frequencies'][i]}, Patterns {p}\n")
        f.close()

        motifs_db[m_idx] = None
        m_idx += 1

# EM()