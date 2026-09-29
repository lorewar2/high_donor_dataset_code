import random
import os
from collections import defaultdict
import csv
import pysam
import sys
from pathlib import Path

UNIQUE_CELLS = 500 # number of cells required
BAR_CODE_MIN_READ = 10_000 # Min number of reads corrosponding to cell
SEED = 255

def main():
    #analysis_code()
    extract_code()
    return

# extract code
def extract_code():
    random.seed(SEED)
    cb_dict = defaultdict(list)
    folder_path = sys.argv[1]
    bam_files = list(Path(folder_path).glob("*.bam"))
    bam_file_path = bam_files[0]
    cell_file_path = folder_path + "/barcodes_high_reads.tsv"
    # from the cell file randomly select UNIQUE CELLS
    cells = []
    with open(cell_file_path, "r", encoding="utf-8") as f:
        for cell in f:
            cells.append(cell.split("\t")[0])
    selected_cells = random.sample(cells, min(UNIQUE_CELLS, len(cells)))
    print(selected_cells)
    # load the selected reads
    for selected_cell in selected_cells:
        cb_dict[selected_cell] = []
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
            # Iterate through all reads in the BAM file
            break_index = 0
            for read in bam_file.fetch():
                if read.has_tag("CB"):  # Check if read has "CB" tag
                    cb_tag = read.get_tag("CB")
                    if cb_tag not in cb_dict:
                        continue
                    else:
                        cb_dict[cb_tag].append(read)
                if break_index > 500_000_000:
                    break
                break_index += 1
                if break_index % 1_000_000 == 0:
                    print(break_index)
    # modify the CB tag
    modified_cb_dict = modify_cb_tags(cb_dict, folder_path)
    # save the reads with modified CB tag
    save_modified_reads(modified_cb_dict , folder_path + "/final.bam", folder_path + "/final.tsv", bam_file_path)
    return

def save_modified_reads(modified_reads, output_bam_path, output_barcodes_path, template_path):
    # Get the list of unique CB tags
    unique_cb_tags = list(modified_reads.keys())
    # Open a new BAM file for writing
    bamfile = pysam.AlignmentFile(template_path, "rb")
    with pysam.AlignmentFile(output_bam_path, "wb", template=bamfile) as out_bam:
        # Write each read to the output BAM file
        for cb_tag, reads in modified_reads.items():
            for read in reads:
                out_bam.write(read)
    # Write the unique CB tags to a barcodes.tsv file
    with open(output_barcodes_path, "w") as barcode_file:
        for cb_tag in unique_cb_tags:
            barcode_file.write(f"{cb_tag}\n")
        
def modify_cb_tags(sampled_reads, modify_with):
    # Dictionary to store the modified sampled reads with updated CB tags
    modified_reads = {}
    for cb_tag, reads in sampled_reads.items():
        # Modify the CB tag from ending with "-1" to "-2"
        if cb_tag.endswith("-1"):
            new_cb_tag = cb_tag[:-2] + modify_with
        else:
            print(f"Warning: CB tag {cb_tag} does not end with '-1'. Skipping.")
            new_cb_tag = cb_tag + modify_with
        # Update the CB tag in each read
        modified_reads[new_cb_tag] = []
        for read in reads:
            # Create a copy of the read to avoid modifying the original object (if needed)
            read_copy = read
            # Update the "CB" tag within the read to match the new CB tag
            read_copy.set_tag("CB", new_cb_tag, value_type="Z")
            # Add the modified read to the new CB tag list
            modified_reads[new_cb_tag].append(read_copy)
    return modified_reads

# analysis code
def analysis_code():
    cb_dict = defaultdict(list)
    folder_path = sys.argv[1]
    bam_files = list(Path(folder_path).glob("*.bam"))
    bam_file_path = bam_files[0]
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        # Iterate through all reads in the BAM file
        break_index = 0
        for read in bam_file.fetch():
            if read.has_tag("CB"):  # Check if read has "CB" tag
                cb_tag = read.get_tag("CB")
                if cb_tag not in cb_dict:
                    cb_dict[cb_tag] = 1
                else:
                    cb_dict[cb_tag] += 1
            if break_index > 500_000_000:
                break
            break_index += 1
            if break_index % 1_000_000 == 0:
                print(break_index)
    cb_dict = dict(sorted(cb_dict.items(), key=lambda item: item[1]))
    num_greater_than_min = 0
    with open(folder_path + '/barcodes_high_reads.tsv', 'w') as file:
        for key, value in cb_dict.items():
            if value > BAR_CODE_MIN_READ:
                num_greater_than_min += 1
                file.write(f"{key}\t{value}\n")
                print(f"{key}\t{value}")
    print(num_greater_than_min)
    return

if __name__ == "__main__":
    main()