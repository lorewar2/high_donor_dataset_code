from collections import defaultdict

CELLS_PER_DONOR = 500
NUMBER_OF_DONORS = 250

def main():
    selects_cells_and_donors()
    return

def selects_cells_and_donors():
    donor_barcode_dict = defaultdict(list)
    filtered_barcodes_path = "merged.tsv"
    # get the available barcodes 
    with open(filtered_barcodes_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            donor = line.split("-")[1]
            if donor not in donor_barcode_dict:
                if len(donor_barcode_dict) < NUMBER_OF_DONORS:
                    donor_barcode_dict[donor] = []
            else:
                if len(donor_barcode_dict[donor]) < CELLS_PER_DONOR:
                    donor_barcode_dict[donor].append(line)
    with open("{}_{}.tsv".format(NUMBER_OF_DONORS, CELLS_PER_DONOR), "w") as barcode_file:
        for (donor, line) in donor_barcode_dict.items():
            for entry in line:
                barcode_file.write(f"{entry}\n")
if __name__ == "__main__":
    main()