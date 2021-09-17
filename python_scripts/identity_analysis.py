import os
import configparser
import subprocess

from blast_analysis import BlastAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class IdentityAnalysis:

    def __init__(self, gene, folder_name):
        self.gene = gene
        self.folder_name = folder_name
        self.blast_functions_object = BlastAnalysis(self.gene, self.folder_name)
        
        folder_name = folder_name.split("/")
        self.cds_blast_folder = self.blast_functions_object.directory_formatter(self.gene + "_Nucleotide_Blast")
        self.reference_folder = self.blast_functions_object.directory_formatter("reference" + self.gene)
        self.identity_analysis_folder = self.blast_functions_object.directory_formatter("reference" + self.gene)

    def blast_to_length_comparison(self):
        self.blast_functions_object.check_directory_path(self.identity_analysis_folder)
        url = config_parser.get(self.gene, "mibig_url")
        filename = url.split("/")[-1]
        # NOTE: Many more options here. When more than TDA gene is tested investigate this!!
        if filename.endswith(".gbk"):
            filename = filename.split(".")[0]
        ref_gene = open(self.reference_folder + "/" + filename + ".ffn", "rb")
        seq = ""
        for line in ref_gene:
            if line.startswith(b">"):
                name = line.rstrip()
            else:
                seq = line.rstrip()
                self.identity_and_length(name, len(seq), self.cds_blast_folder)
    
    def identity_and_length(self, name, seq_length, folder):
        self.identity_list = []
        self.length_list = []
        print(folder)
        blast_files = os.listdir(folder)
        for file in blast_files:
            blast_file = open(folder + "/" + file,"r")
            for blast_result in blast_file:
                if blast_result.startswith(str(name[1:], "UTF-8")):
                    # TODO: Fix this.
                    print(blast_result.split("\t")[2])
                    exit(1)
                    try:
                        identity = blast_result.split("\t")[2]
                        length = (blast_result.split("\t")[3]/seq_length)*100
                        self.identity_list.append(identity)
                        self.length_list.append(length)
                    except:
                        print("Error")
                    break
        print(self.identity_list)
        print(self.length_list)
        # TODO: Plot with: https://matplotlib.org/stable/tutorials/introductory/pyplot.html