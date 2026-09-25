import random
import os
from collections import defaultdict
import csv
import pysam

SEED = 10
INPUT_FILE_PATH = "data.bam"
INPUT_BAR_PATH = "data.tsv"
OUTPUT_BAM_PATH = "with_doubs.bam5"
OUTPUT_BAR_PATH = "with_doubs.tsv5"
REQUIRED_CELL_COUNT = 2500
START_INDEX = 0
END_INDEX = 1_000_000_000
DOUBLET_VALUE = 5 # for 5% percent == 10 # for 10% percent == 5
FIRST_CELL_COUNT = 50 # minority cluster cell count (10% of total?)

def main():
    # Random initialize with seed
    random.seed(SEED)
    normal_list, doublet_list = open_bar_code_file_get_doublet_cells(INPUT_BAR_PATH)
    paired_doublet_list = generate_doublet_list(doublet_list)
    normal_by_donor, doublet_by_donor = read_all_bam_files(doublet_list, normal_list)
    doublet_list = attach_reads_using_doublet_list(paired_doublet_list, doublet_by_donor)
    save_modified_reads(normal_by_donor, doublet_list, OUTPUT_BAM_PATH, OUTPUT_BAR_PATH, INPUT_FILE_PATH)
    return

def generate_doublet_list(doublet_list):
    paired_doublet_list = []
    doublet_list_donor = [[] for _ in range(71)]
    # put doublet list by donor
    for entry in doublet_list:
        bar_code = entry.split("-")[0]
        donor = int(entry.split("-")[1])
        doublet_list_donor[donor].append(entry)
    # remove empty
    doublet_list_donor = [arr for arr in doublet_list_donor if arr]
    # Select pairs until list empty
    doublet_list_donor_empty = False
    while (doublet_list_donor_empty == False):
        two_selected_donors = random.sample(range(len(doublet_list_donor)), 2)
        cell_index_1 = random.sample(range(len(doublet_list_donor[two_selected_donors[0]])), 1)[0]
        cell_index_2 = random.sample(range(len(doublet_list_donor[two_selected_donors[1]])), 1)[0]
        # new required data
        cb_tag_1 = doublet_list_donor[two_selected_donors[0]][cell_index_1]
        cb_tag_2 = doublet_list_donor[two_selected_donors[1]][cell_index_2]
        new_cb_tag = "{}-{}-{}".format(doublet_list_donor[two_selected_donors[0]][cell_index_1].split("-")[0], doublet_list_donor[two_selected_donors[0]][cell_index_1].split("-")[1], doublet_list_donor[two_selected_donors[1]][cell_index_2].split("-")[1])
        #print(cb_tag_1, cb_tag_2, new_cb_tag)
        paired_doublet_list.append((cb_tag_1, cb_tag_2, new_cb_tag))
        # delete the entry
        del(doublet_list_donor[two_selected_donors[0]][cell_index_1])
        del(doublet_list_donor[two_selected_donors[1]][cell_index_2])
        # remove empty
        doublet_list_donor = [arr for arr in doublet_list_donor if arr]
        if len(doublet_list_donor) < 2:
            doublet_list_donor_empty = True
    print(doublet_list_donor)
    return paired_doublet_list

def attach_reads_using_doublet_list(paired_doublet_list, doublet_dics_by_donor):
    list_with_doublets = []
    # go thorugh the paired doublet list get reads and append them together
    for (cb_1, cb_2, new_cb) in paired_doublet_list:
        # get donor 1
        donor_1 = int(cb_1.split("-")[1])
        # get donor 2
        donor_2 = int(cb_2.split("-")[1])
        appended_reads = doublet_dics_by_donor[donor_1][cb_1] + doublet_dics_by_donor[donor_2][cb_2]
        modified_reads = []
        for read in appended_reads:
            # Create a copy of the read to avoid modifying the original object (if needed)
            read_copy = read
            # Update the "CB" tag within the read to match the new CB tag
            read_copy.set_tag("CB", new_cb, value_type="Z")
            modified_reads.append(read_copy)
        list_with_doublets.append((new_cb, modified_reads))
    return list_with_doublets

def open_bar_code_file_get_doublet_cells(barcode_file):
    curr_donor = 0
    prev_donor = 0
    final_list_of_doublets = []
    final_list_of_cells = []
    temp_list = []
    cell_index = 0
    donor_index = 0
    first_cell = True
    with open(barcode_file, 'r') as file:
        for line in file:
            split_line = line.split("\t")
            cell_id = split_line[0].split("-")[0]
            curr_donor = int(split_line[0].split("-")[1])
            if curr_donor == prev_donor:
                if first_cell and (cell_index < FIRST_CELL_COUNT):
                    temp_list.append(cell_id.strip())
                elif not first_cell and (cell_index < REQUIRED_CELL_COUNT):
                    temp_list.append(cell_id.strip())
            else:
                first_cell = False
                # process stuff, get 10 percent of cells and save it in final_list of doublets
                percent_10_size = int(len(temp_list) / DOUBLET_VALUE)
                doublet_selected = random.sample(temp_list, percent_10_size)
                for value in temp_list:
                    if value in doublet_selected:
                        final_list_of_doublets.append("{}-{}".format(value, prev_donor))
                    else:
                        final_list_of_cells.append("{}-{}".format(value, prev_donor))
                # clear the list
                temp_list.clear()
                temp_list.append(cell_id.strip())
                prev_donor = curr_donor
                donor_index += 1
                cell_index = 0
            cell_index += 1
        percent_10_size = int(len(temp_list) / DOUBLET_VALUE)
        doublet_selected = random.sample(temp_list, percent_10_size)
        for value in temp_list:
            if value in doublet_selected:
                final_list_of_doublets.append("{}-{}".format(value, prev_donor))
            else:
                final_list_of_cells.append("{}-{}".format(value, prev_donor))
    return final_list_of_cells, final_list_of_doublets

def read_all_bam_files(doublet_list, normal_list):
    doublet_by_donor = [defaultdict(list) for _ in range(71)]
    normal_by_donor = [defaultdict(list) for _ in range(71)]
    bam_file_path = "{}".format(INPUT_FILE_PATH)
    index = 0
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        # Iterate through all reads in the BAM file
        for read in bam_file.fetch():
            index += 1
            if index % 1_000_000 == 0:
                print(index)
            if index < START_INDEX:
                continue
            if index > END_INDEX:
                break
            if read.has_tag("CB"):  # Check if read has "CB" tag
                cb_tag = read.get_tag("CB")
                donor_index = int(cb_tag.split("-")[1])
                if cb_tag.strip() in doublet_list:
                    doublet_by_donor[donor_index][cb_tag].append(read)
                elif cb_tag.strip() in normal_list:
                    normal_by_donor[donor_index][cb_tag].append(read)

    return (normal_by_donor, doublet_by_donor)

def save_modified_reads(all_list, doublet_list, output_bam_path, output_barcodes_path, template_path):
    # Get the list of unique CB tags
    unique_cb_tags = []
    # Open a new BAM file for writing
    bamfile = pysam.AlignmentFile(template_path, "rb")
    with pysam.AlignmentFile(output_bam_path, "wb", template=bamfile) as out_bam:
        # Write each read to the bam file
        for donor_index in range(71):
            unique_cb_tags.extend(list(all_list[donor_index].keys()))
            # all file
            for cb_tag, reads in all_list[donor_index].items():
                for read in reads:
                    out_bam.write(read)
        # doublet file
        for entry in doublet_list:
            unique_cb_tags.append(entry[0])
            for read in entry[1]:
                out_bam.write(read)
            
    # Write the unique CB tags to a barcodes.tsv file
    with open(output_barcodes_path, "w") as barcode_file:
        for cb_tag in unique_cb_tags:
            barcode_file.write(f"{cb_tag}\n")

if __name__ == "__main__":
    main()