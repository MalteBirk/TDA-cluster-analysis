import os
import subprocess
import configparser

from file_extractor import FileExtractor
from blast_analysis import BlastAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class MainScript:

    def __init__(self, test):
        self.test = test
    
    def gbk_splitter_script(self):
        script_directory = os.listdir()
        if "seqIO_extract.py" not in script_directory:
            subprocess.run(["wget", 
            "https://raw.githubusercontent.com/davised/seqIO_extract/master/seqIO_extract.py"])
        
    def file_extractor(self, gene):
        extractor_object = FileExtractor(gene)
        extractor_object.extract_reference_genes()
    
    def blast_analysis(self):
        blast_object = BlastAnalysis()
        blast_object.gene_blast()


test_object = MainScript(".")

test_object.gbk_splitter_script()
# Extract files and make directories
# NOTE: Make it possible to have multiple config files and parse them on the command line.
for sect in config_parser.sections():
    test_object.file_extractor(sect)
    test_object.blast_analysis(sect)
