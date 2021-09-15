import os
import configparser
import subprocess
from blast_analysis import BlastAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class HousekeepingAnalysis:

    def __init__(self, gene, folder_name):
        self.gene = gene
        self.folder_name = folder_name
        self.blast_functions_object = BlastAnalysis(self.gene, self.folder_name)
        self.housekeeping_list = config_parser.get(self.gene, "housekeeping_genes").split(",")
        self.extracted_path = self.gene + "_" + "extracted_housekeeping_genes"
        self.aligned_path = self.gene + "_" + "alligned_housekeeping_genes"
        self.extracted_protein_path = self.gene + "_" + "extracted_housekeeping_proteins"
        self.aligned_protein_path = self.gene + "_" + "alligned_housekeeping_proteins"

    def extract_genes(self):
        # NOTE: Writing over and inside original files, fix so that new files are made.
        self.blast_functions_object.check_directory_path(self.extracted_path)
        self.blast_functions_object.check_directory_path(self.aligned_path)
        cds_files = os.listdir("cds_fasta" + "_" + self.folder_name)

        fs = []
        try:
            housekeeping_index = 0
            housekeeping_number_list = []
            for gene in self.housekeeping_list:
                housekeeping_number_list.append(gene)
                housekeeping_number_list.append(housekeeping_index)
                fs.append(open(self.extracted_path + "/" + gene, 'w'))
                housekeeping_index += 1
            for cds_file in cds_files:
                flag = False
                infile = open("cds_fasta" + "_" + self.folder_name + "/" + cds_file)
                for line in infile:
                    if line.startswith(">"):
                        gene_name = line.split(" ")[1]
                        if gene_name.startswith("[gene="):
                            gene_name = gene_name.split("=")[1][:-1]
                            if gene_name in self.housekeeping_list:
                                file_position = housekeeping_number_list.index(gene_name) + 1
                                fs[housekeeping_number_list[file_position]].write(line)
                                flag = True
                            else:
                                flag = False
                        else:
                            flag = False
                    elif flag == True:
                        fs[housekeeping_number_list[file_position]].write(line)
                
                infile.close()
        
        finally:
            for f in fs:
                f.close()

    def extract_proteins(self):
        # NOTE: Writing over and inside original files, fix so that new files are made.
        self.blast_functions_object.check_directory_path(self.extracted_protein_path)
        self.blast_functions_object.check_directory_path(self.aligned_protein_path)
        cds_files = os.listdir("protein_cds_fasta" + "_" + self.folder_name)

        fs = []
        try:
            housekeeping_index = 0
            housekeeping_number_list = []
            for gene in self.housekeeping_list:
                housekeeping_number_list.append(gene)
                housekeeping_number_list.append(housekeeping_index)
                fs.append(open(self.extracted_protein_path + "/" + gene, 'w'))
                housekeeping_index += 1
            for cds_file in cds_files:
                flag = False
                infile = open("protein_cds_fasta" + "_" + self.folder_name + "/" + cds_file)
                for line in infile:
                    if line.startswith(">"):
                        gene_name = line.split(" ")[1]
                        if gene_name.startswith("[gene="):
                            gene_name = gene_name.split("=")[1][:-1]
                            if gene_name in self.housekeeping_list:
                                file_position = housekeeping_number_list.index(gene_name) + 1
                                fs[housekeeping_number_list[file_position]].write(line)
                                flag = True
                            else:
                                flag = False
                        else:
                            flag = False
                    elif flag == True:
                        fs[housekeeping_number_list[file_position]].write(line)
                
                infile.close()
        
        finally:
            for f in fs:
                f.close()