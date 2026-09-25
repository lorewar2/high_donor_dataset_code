import pandas as pd
import os
import shutil

# Replace 'your_file.tsv' with the path to your TSV file
file_path = 'HCA2.tsv'

# Read the TSV file
df = pd.read_csv(file_path, sep='\t')

pd.set_option('display.max_rows', None)
# Display the contents of the DataFrame
print(df.columns)
df_bam = df[df["file_format"] == "bam"]
print(df_bam[[ "file_name", "bundle_uuid"]])

for index, row in df.iterrows():
    file_name = row['file_name']
    folder_name = row['bundle_uuid']
    if folder_name == "286eb67a-7ea7-54dc-a456-07d366cb30a8":
        print("Skipping")
        continue
    #print(row['donor_organism.biomaterial_core.biomaterial_id'], row['file_name'])
    file_path = "./" + folder_name + "/" + file_name
    #print(file_path)
    
    if os.path.isfile(file_path):
        # Create a directory with the corresponding name (excluding file extension)
        print("file_available", file_name)
        newpath = row['donor_organism.biomaterial_core.biomaterial_id']
        os.makedirs(newpath, exist_ok=True)
        source = file_path
        destination = "./" + newpath + "/" + file_name
        print(destination)
        try:
            os.rename(source, destination)
        except:
            print("error")
