import os
import subprocess
import configparser

from file_extractor import FileExtractor
from blast_analysis import BlastAnalysis
from housekeeping_analysis import HousekeepingAnalysis
from identity_analysis import IdentityAnalysis
from organisation_analysis import OrganisationAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class MainScript:

    def __init__(self, name):
    # Have some argparses here?
        self.folder_name = name

    def initator(self):
        if os.path.exists("../Results"):
            pass
        else:
            os.mkdir("../Results")
        self._gbk_splitter_script()
        # Extract files and make directories
        # NOTE: Make it possible to have multiple config files and parse them on the command line.
        for gene in config_parser.sections():
            if os.path.exists("../Results/" + config_parser.get(gene, "gene_name")):
                pass
            else:
                os.mkdir("../Results/" + config_parser.get(gene, "gene_name"))
            self.folder_name = "../Results/" + config_parser.get(gene, "gene_name") + "/" + \
                               config_parser.get(gene, "gene_name") + "_" \
                               + config_parser.get(gene, "organism_list").replace(",","_")
            #self._file_extractor(gene)
            #print("starting_blast_analysis")
            #self._blast_analysis(gene)
            #print("extracting_housekeeping_genes")
            #self._housekeeping_analysis(gene)
            #print("finding gene organisation")
            #self._organisation_analysis(gene)
            print("Doing identity analysis")
            self._identity_analysis(gene)

    def _gbk_splitter_script(self):
        script_directory = os.listdir()
        if "seqIO_extract.py" not in script_directory:
            subprocess.run(["wget", 
            "https://raw.githubusercontent.com/davised/seqIO_extract/master/seqIO_extract.py"])
        
    def _file_extractor(self, gene):
        extractor_object = FileExtractor(gene, self.folder_name)
        extractor_object.extract_reference_genes()
    
    def _blast_analysis(self, gene):
        blast_object = BlastAnalysis(gene, self.folder_name)
        blast_object.gene_blast()
        blast_object.protein_blast()
        blast_object.extract_genes()
        blast_object.extract_protein()
    
    def _housekeeping_analysis(self, gene):
        housekeeping_object = HousekeepingAnalysis(gene, self.folder_name)
        housekeeping_object.extract_genes()
        housekeeping_object.extract_proteins()
    
    def _organisation_analysis(self, gene):
        organisation_object = OrganisationAnalysis(gene, self.folder_name)
        #organisation_object.extract_gene_organisation()
        organisation_object.find_organisation()

    def _identity_analysis(self, gene):
        identity_object = IdentityAnalysis(gene, self.folder_name)
        #identity_object.blast_to_length_comparison()
        identity_object.tree_maker()


if __name__ == "__main__":
    test_object = MainScript(" ")
    test_object.initator()

