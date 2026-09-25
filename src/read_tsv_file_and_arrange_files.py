import pandas as pd
import os
import shutil

# Replace 'your_file.tsv' with the path to your TSV file
file_path = 'HCA.tsv'

# Read the TSV file
df = pd.read_csv(file_path, sep='\t')

pd.set_option('display.max_rows', None)
# Display the contents of the DataFrame
print(df.columns)
#df_bam = df[df["file_format"] == "bam"]
print(df[["donor_organism.biomaterial_core.biomaterial_id", "file_name"]])

for index, row in df.iterrows():
    file_name = row['file_name']
    #print(row['donor_organism.biomaterial_core.biomaterial_id'], row['file_name'])
    if os.path.isfile(file_name):
        # Create a directory with the corresponding name (excluding file extension)
        print("file_available", file_name)
        newpath = row['donor_organism.biomaterial_core.biomaterial_id']
        os.makedirs(newpath, exist_ok=True)
        source = "./" + file_name
        destination = "./" + newpath + "/" + file_name
        print(destination)
        #try:
        #    os.rename(source, destination)
        #except:
        #    print("what")