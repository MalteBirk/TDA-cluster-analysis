import os
import configparser
import subprocess
import shutil

from blast_analysis import BlastAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class HousekeepingAnalysis:

    def __init__(self, gene, folder_name):
        self.gene = gene
        self.folder_name = folder_name
        self.blast_functions_object = BlastAnalysis(self.gene, self.folder_name)

        folder_name = folder_name.split("/")
        self.housekeeping_list = config_parser.get(self.gene, "housekeeping_genes").split(",")
        self.extracted_path = self.blast_functions_object.directory_formatter(self.gene + "_" + "extracted_housekeeping_genes")
        self.aligned_path = self.blast_functions_object.directory_formatter(self.gene + "_" + "aligned_housekeeping_genes")
        self.extracted_protein_path = self.blast_functions_object.directory_formatter(self.gene + "_" + "extracted_housekeeping_proteins")
        self.aligned_protein_path = self.blast_functions_object.directory_formatter(self.gene + "_" + "aligned_housekeeping_proteins")
        
        self.cds_folder = self.blast_functions_object.directory_formatter("cds_fasta" + "_" + folder_name[-1])
        self.protein_cds_folder = self.blast_functions_object.directory_formatter("protein_cds_fasta" + "_" + folder_name[-1])

    def extract_genes(self):
        # NOTE: Writing over and inside original files, fix so that new files are made or remove old folders.
        # Remove old files
        if os.path.exists(self.extracted_path):
            shutil.rmtree(self.extracted_path)
        if os.path.exists(self.aligned_path):
            shutil.rmtree(self.aligned_path)

        self.blast_functions_object.check_directory_path(self.extracted_path)
        self.blast_functions_object.check_directory_path(self.aligned_path)
        cds_files = os.listdir(self.cds_folder)

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
                infile = open(self.cds_folder + "/" + cds_file)
                for line in infile:
                    if line.startswith(">"):
                        # Some genes don't have 
                        try:
                            gene_name = line.split(" ")[1]
                        except:
                            print(line)
                            print(cds_file)
                            continue
                            #exit(1)
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
        finally:
            for f in fs:
                f.close()
        
        housekeeping_files = os.listdir(self.extracted_path)
        for file in housekeeping_files:
            self.alignment(self.extracted_path + "/" + file, self.aligned_path, file)


    def extract_proteins(self):
        # NOTE: Writing over and inside original files, fix so that new files are made.
        # Remove old files
        if os.path.exists(self.extracted_protein_path):
            shutil.rmtree(self.extracted_protein_path)
        if os.path.exists(self.aligned_protein_path):
            shutil.rmtree(self.aligned_protein_path)
        self.blast_functions_object.check_directory_path(self.extracted_protein_path)
        self.blast_functions_object.check_directory_path(self.aligned_protein_path)
        cds_files = os.listdir(self.protein_cds_folder)

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
                infile = open(self.protein_cds_folder + "/" + cds_file)
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
        
        housekeeping_files = os.listdir(self.extracted_protein_path)
        for file in housekeeping_files:
            self.alignment(self.extracted_protein_path + "/" + file, self.aligned_protein_path, file)
    
    def alignment(self, path, second_location, gene_name):
        alignment = subprocess.run(["mafft", "--quiet", "--adjustdirection", path], capture_output=True)
        alignment_file = open(second_location + "/" + gene_name + ".aln", "wb")
        alignment_file.write(alignment.stdout)
        alignment_file.close()