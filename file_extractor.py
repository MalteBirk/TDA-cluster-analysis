import os
import subprocess

class FileExtractor:

    def __init__(self,test):
        self.info = test
    
    def extract_reference_genes(self):
        os.mkdir("referenceTDA")
        subprocess.run(["wget", "https://mibig.secondarymetabolites.org/repository/BGC0000932/BGC0000932.gbk", "-P", "referenceTDA/"])
