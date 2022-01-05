import os
import subprocess
import configparser
import argparse
import shutil

from file_extractor import FileExtractor
from blast_analysis import BlastAnalysis
from housekeeping_analysis import HousekeepingAnalysis
from identity_analysis import IdentityAnalysis
from organisation_analysis import OrganisationAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class MainScript:

    def __init__(self, name, args):
        self.folder_name = name
        self.args = args

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
            self._file_extractor(gene)
            self._local_file_handler()
            print("starting_blast_analysis")
            self._blast_analysis(gene)
            print("extracting_housekeeping_genes")
            self._housekeeping_analysis(gene)
            print("finding gene organisation")
            self._organisation_analysis(gene)
            print("Doing identity analysis")
            self._identity_analysis(gene)

    def _local_file_handler(self):
        if self.args.local_genome is not None:
            genomes_folder_name = self.folder_name.split("/")
            genomes_folder_name[-1] = "genomes" + "_" + genomes_folder_name[-1]
            self.genomes_folder_name = "/".join(genomes_folder_name)
            local_genome_files = os.listdir(self.args.local_genome)
            for file in local_genome_files:
                shutil.copyfile(self.args.local_genome + "/" + file, self.genomes_folder_name + "/" + file)
        
        if self.args.local_pgb is not None:
            local_pgb_files = os.listdir(self.args.local_pgb)
            cds_folder_name = self.folder_name.split("/")
            cds_folder_name[-1] = "protein_cds_fasta" + "_" + cds_folder_name[-1]
            self.cds_folder_name = "/".join(cds_folder_name)

            for file in local_pgb_files:
                filename = file.split(".")[0] + "_cds_from_genomic.fna"
                outfile = open(self.cds_folder_name + "/" + filename, "w")
                genbank_file = open(self.args.local_pgb + "/" + file, "r")
                flag = False
                translation = ""
                locus_tag = ""
                location = ""
                gene = ""
                for line in genbank_file:
                    if "ORIGIN" in line:
                        flag = False
                        if translation != "":
                            outfile.write(">Jyllinge_" + locus_tag + " [gene=" + gene + "]" + " [location=" + location + "]" + "\n")
                            for i in range(0,len(translation),60):
                                outfile.write(translation[i:i+60] + "\n")
                        translation = ""
                    if " CDS " in line:
                        if translation != "":
                            outfile.write(">Jyllinge_" + locus_tag + " [gene=" + gene + "]" + " [location=" + location + "]" + "\n")
                            for i in range(0,len(translation),60):
                                outfile.write(translation[i:i+60] + "\n")
                        location = line.split(" ")[-1].rstrip()
                        flag = False
                        gene = ""
                    if " /gene=" in line:
                        gene = line.split("=")[-1].replace('"',"").rstrip()
                    if " /locus_tag=" in line:
                        locus_tag = line.split("=")[-1].replace('"',"").rstrip()
                        name = locus_tag.rsplit("_", 1)[0] + "_lcl|"
                        number = "WP_" + locus_tag.rsplit("_", 1)[1] + ".1_1_1"
                        locus_tag = name + number
                    if flag == True:
                        translation += line.split(" ")[-1].replace('"',"").rstrip()
                    if " /translation=" in line:
                        translation = line.split("=")[-1].replace('"',"").rstrip()
                        flag = True
                if translation != "":
                    outfile.write(">Jyllinge_" + locus_tag + " [gene=" + gene + "]" + " [location=" + location + "]" + "\n")
                    for i in range(0,len(translation),60):
                        outfile.write(translation[i:i+60] + "\n")

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
        housekeeping_object.protein_blast()
        housekeeping_object.extract_housekeeping_protein()
    
    def _organisation_analysis(self, gene):
        organisation_object = OrganisationAnalysis(gene, self.folder_name)
        organisation_object.extract_gene_organisation()
        organisation_object.find_organisation()

    def _identity_analysis(self, gene):
        identity_object = IdentityAnalysis(gene, self.folder_name)
        identity_object.blast_to_length_comparison()
        identity_object.tree_maker()
        identity_object.blast_identity_to_operons()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='find local folder')
    parser.add_argument('-local_pgb', '--local_pgb', action='store', required = False, dest = "local_pgb")
    parser.add_argument('-local_genome' '--local_genome', action='store', required = False, dest = "local_genome")
    args = parser.parse_args()
    test_object = MainScript(" ", args)
    test_object.initator()
