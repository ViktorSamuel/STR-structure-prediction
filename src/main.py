import os
import shutil

from em import EM


def alleles_output(output_folder):
   
    # Loop through each subfolder
    main_folder = output_folder + "frequencies/"
    res_folder = output_folder + "alleles/"
    if os.path.exists(res_folder):
        shutil.rmtree(res_folder)
    os.makedirs(res_folder)

    for subfolder in os.listdir(main_folder):
        subfolder_path = os.path.join(main_folder, subfolder)
        file = subfolder + ".txt"
        file_path = os.path.join(subfolder_path, file)
        
        with open(file_path, 'r') as f:
            # Read file content
            lines = f.readlines()

        # Sort lines based on frequency in descending order
        sorted_lines = sorted(
            lines,
            key=lambda line: float(line.split("Frequency: ")[1].split(",")[0]),
            reverse=True
        )

        # Get the lines with the highest and second highest frequencies
        filtered_lines = sorted_lines[:2]

        if len(filtered_lines) != 2:
            filtered_lines.append(filtered_lines[0])

        if float(filtered_lines[1].split("Frequency: ")[1].split(",")[0]) < 0.17:
            filtered_lines[1] = filtered_lines[0]

        # Write filtered content to a new file in the results folder
        output_path = os.path.join(res_folder, file)  # Keep the same file name
        with open(output_path, 'w') as output_file:
            output_file.writelines(filtered_lines)


pseudo_counts = False
flank_length = 296
read_length = 148
alleles_file = "../data/pop_alleles.txt"
reference_genome_file = "../data/grch38_decoy.fa"
reference_alignment_file = "../data/HG002.GRCh38.selected_w_pairs.bam"
output_folder = "../output/"
    
EM(output_folder, pseudo_counts, flank_length, read_length, alleles_file, reference_genome_file, reference_alignment_file)

alleles_output(output_folder)




