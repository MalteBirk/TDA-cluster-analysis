import os
import argparse
parser = argparse.ArgumentParser(description='find folder')
parser.add_argument('-l','--list', action='append', help='<Required> Set flag', required=True)
args = parser.parse_args()
folder = args.list[0]
organism = args.list[1]

all_files = os.listdir(folder)

for file in all_files:
    gene_file = open(folder + "/" + file, "r")
    for line in gene_file:
        line = line.split(" ", 1)[1]
        if organism in line:
            print(file)
        break