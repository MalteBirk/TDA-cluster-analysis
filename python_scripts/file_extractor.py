import os
import subprocess
import configparser
import fileinput

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class FileExtractor:

    def __init__(self, gene):
        self.gene = gene
        self.folder_name = config_parser.get(self.gene, "gene_name") + "_" \
            + config_parser.get(self.gene, "organism_list").replace(",","_")
        self.cds_list = os.listdir("cds_fasta" + "_" + self.folder_name)

    def extract_reference_genes(self):
        self.reference_gene = "../"+"reference"+self.gene + "/"
        url = config_parser.get(self.gene, "mibig_url")
        filename = url.split("/")[-1]
        gene_name = filename.split(".")[0]
        if os.path.exists(self.reference_gene):
            pass
        else:
            os.mkdir(self.reference_gene)
        
        # Get reference gene.
        subprocess.run(["wget", url, "-P", self.reference_gene])

        # NOTE: Make it remove the additional reference gene after script has run through.
        file = subprocess.run(["python3", "seqIO_extract.py", 
                                "--fna", self.reference_gene + filename, "all"], 
                                capture_output = True)
        self._fna_file_writer(gene_name, file)

        file = subprocess.run(["faidx", self.reference_gene + gene_name + ".fna", 
                        gene_name + config_parser.get(self.gene, "fasta_extraction")], 
                        capture_output = True)
        self._fna_file_writer(config_parser.get(self.gene, "gene_names"), file)
        self._gene_list_writer()
        
        # ffn file
        file = subprocess.run(["python3", "seqIO_extract.py", 
                                "--outname", "gene", "--ffn", 
                                "--matchtype", "exact", "--searchname", 
                                "gene", self.reference_gene + filename, 
                                self.reference_gene + "geneList.txt"],
                                capture_output = True)
        outfile = open(self.reference_gene + gene_name + ".ffn", "wb")
        outfile.write(file.stdout)
        outfile.close()
        
        # ffa file
        file = subprocess.run(["python3", "seqIO_extract.py", 
                        "--outname", "gene", 
                        "--matchtype", "exact", "--searchname", 
                        "gene", self.reference_gene + filename, 
                        self.reference_gene + "geneList.txt"],
                        capture_output = True)
        outfile = open(self.reference_gene + gene_name + ".faa", "wb")
        outfile.write(file.stdout)
        outfile.close()

        os.remove(self.reference_gene + filename)
        os.remove(self.reference_gene + gene_name + ".fna")
        os.remove(self.reference_gene + gene_name + ".fna.fai")
        os.remove(self.reference_gene + "geneList.txt")

        #self._ncbi_download()
        #self._coding_sequence_processing()
        self._translate_coding_sequences()

    def _fna_file_writer(self, gene_name, file):
        outfile = open(self.reference_gene + gene_name + ".fna", "wb")
        outfile.write(file.stdout)
        outfile.close()
    
    def _gene_list_writer(self):
        gene_list = config_parser.get(self.gene, "gene_list")
        gene_list = gene_list.split(",")
        outfile = open(self.reference_gene + "geneList.txt", "w")
        for gene in gene_list:
            outfile.write(gene + "\n")
        outfile.close()

    def _ncbi_download(self):      
        print("downloading cds from", \
            config_parser.get(self.gene, "organism_list").replace(","," "))
        subprocess.run(["ncbi-genome-download", "-p",
                        config_parser.get(self.gene, "ncbi_parallel"), 
                        "-F", "cds-fasta", "-l", "all", 
                        "--flat-output", "-g", config_parser.get(self.gene, "organism_list"), 
                        "-o", "cds_fasta" + "_" + self.folder_name, "bacteria"])
        subprocess.run(["gunzip", "-r", "cds_fasta" + "_" + self.folder_name])

        print("downloading whole genomes from", \
            config_parser.get(self.gene, "organism_list").replace(","," "))
        subprocess.run(["ncbi-genome-download", "-p",
                config_parser.get(self.gene, "ncbi_parallel"), 
                "-F", "fasta", "-l", "all", 
                "--flat-output", "-g", config_parser.get(self.gene, "organism_list"), 
                "-o", "genomes" + "_" + self.folder_name, "bacteria"])
        subprocess.run(["gunzip", "-r", "genomes" + "_" + self.folder_name])

    def _coding_sequence_processing(self):
        genome_list = os.listdir("genomes" + "_" + self.folder_name)
        for i in range(0, len(self.cds_list)):
            genome_name = self.cds_list[i].split("_", 3)
            genome_name = "_".join(genome_name[0:3])
            # Needs to seek here if the plasmid should be included.
            if genome_list[i].startswith(genome_name):
                genome_file = open("genomes" + "_" + self.folder_name + "/" + genome_list[i])
                # Some highly weird names sometime. Not a good overall pattern anywhere.
                strain_name = genome_file.readline().split(" ", 1)[1].split(",")[0]
                genome_file.close()
            cds_filename = "cds_fasta" + "_" + self.folder_name + "/" + self.cds_list[i]
            # Uses print to change the lines
            for line in fileinput.input(cds_filename, inplace = True):
                line = line.rstrip('\r\n')
                if line.startswith(">"):
                    print(">" + strain_name + " " + line[1:])
                else:
                    print(line)

    def _translate_coding_sequences(self):
        if os.path.exists("protein_cds_fasta" + "_" + self.folder_name):
            pass
        else:
            os.mkdir("protein_cds_fasta" + "_" + self.folder_name)
        for cds in self.cds_list:
            subprocess.run(["transeq", "-sequence",
                "cds_fasta" + "_" + self.folder_name + "/" + cds, 
                "-trim", "boolean", 
                "-outseq", "protein_cds_fasta" + "_" + self.folder_name + "/" + cds], 
                capture_output=True)

# ncbi-genome-download -p 40 -F s-fasta -l all  --flat-output -g Phaeobacter,Epibacterium -o cds_fasta bacteria

#ncbi-genome-download -p 40 -F cds-fasta -l all  --flat-output -g Phaeobacter,Epibacterium,Pseudovibrio,Paracoccus,Stappia -o cds_fasta bacteria
#gunzip cds_fasta/*

#get all genomes, only for naming purposes
#ncbi-genome-download -p 40 -F fasta -l all  --flat-output -g Phaeobacter,Epibacterium,Pseudovibrio,Paracoccus,Stappia -o genomes bacteria
#gunzip genomes/*