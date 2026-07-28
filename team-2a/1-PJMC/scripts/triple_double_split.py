import os
import argparse
import sys

#to use, run this in your run folder:
#python ../../scripts/triple_double_split.py "$(pwd)"

parser = argparse.ArgumentParser()
parser.add_argument("data_folder", type=str)
args = parser.parse_args()

data_folder = args.data_folder  

if not os.path.isdir(data_folder):
    print(f"Error: '{data_folder}' is not a valid directory.")
    sys.exit(1)

print(f"Received folder: {data_folder}")
#===PATHS====


#now only needs to be {data_folder}/csv/...
input_file = f'{data_folder}/csv/compiled_totalGammas.csv'

#let's not use a subdirectory
output_directory = f'{data_folder}/csv'

#make directory in the datafolder if needed
os.makedirs(output_directory, exist_ok=True)

#12 and 8 columns files
twelve_col_file = output_directory + '/triples_data.csv' 
eight_col_file = output_directory + '/doubles_data.csv'

#process the file 
with open(input_file, 'r') as input_f, \
     open(twelve_col_file, 'w') as output_twelve, \
     open(eight_col_file, 'w') as output_eight: 
    
    for line in input_f:
        strip = line.strip()
        if not strip:
            continue 

        col = strip.split(',')
        col_nums = len(col)
        
        if col_nums == 12:
            output_twelve.write(strip + "\n")
        elif col_nums == 8:
            output_eight.write(strip + "\n")

print(f"data is now separated into triples and doubles")

