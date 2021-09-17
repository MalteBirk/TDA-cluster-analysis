import os
import configparser
import subprocess

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class BlastAnalysis:

    def __init__(self, gene, folder_name):
        self.gene = gene
        self.folder_name = folder_name
        self.reference_gene = "../"+"reference"+self.gene + "/"
        self.url = config_parser.get(self.gene, "mibig_url")
        
        folder_name = folder_name.split("/")
        self.cds_blast_folder = self.directory_formatter(gene + "_Nucleotide_Blast")
        self.cds_folder = self.directory_formatter("cds_fasta" + "_" + folder_name[-1])
        self.extracted_gene_folder = self.directory_formatter(gene + "_" + "Extracted_cds")
        self.gene_align_folder = self.directory_formatter(gene + "_" + "Aligned_cds")
        
        self.protein_blast_folder = self.directory_formatter(gene + "_" + "Protein_Blast")
        self.protein_cds_folder = self.directory_formatter("protein_cds_fasta" + "_" + folder_name[-1])
        self.extracted_protein_folder = self.directory_formatter(gene + "_" + "Extracted_protein_cds")
        self.protein_align_folder = self.directory_formatter(gene + "_" + "Aligned_protein_cds")

    def gene_blast(self):
        #TDA_Nucleotide_Blast
        self.check_directory_path(self.cds_blast_folder)
        
        filename = self.url.split("/")[-1]
        gene_name = filename.split(".")[0]
        # Make user able to change identity here later?
        blast_statement = 'blastn  -perc_identity 50 -qcov_hsp_perc 50 -subject \
                        '+ self.cds_folder +'/{} -query '+ self.reference_gene + gene_name + ".ffn"\
                         +' -outfmt 6 -word_size 5 > '+ self.cds_blast_folder +'/{.}.blast'
        cds_list = subprocess.Popen(["ls", self.cds_folder], 
                                stdout=subprocess.PIPE)
        # NOTE: Change the 12 to be altered in the config file.
        parallel_blast = subprocess.Popen(["parallel", "-j", "12",
                        blast_statement],
                        stdin = cds_list.stdout, stderr=subprocess.DEVNULL)
        parallel_blast.wait()
        print("finished blastn")

        count = 0
        # Removes empties
        for file in os.listdir(self.cds_blast_folder):
            if os.path.getsize(self.cds_blast_folder + "/" + file) == 0:
                count += 1
                os.remove(self.cds_blast_folder + "/" + file)
        print(count)
        exit(1)

    def protein_blast(self):

        self.check_directory_path(self.protein_blast_folder)

        filename = self.url.split("/")[-1]
        gene_name = filename.split(".")[0]
        # Make user able to change identity here later?
        blast_statement = 'blastp -subject \
                        '+ self.protein_cds_folder +'/{} -query '+ self.reference_gene + gene_name + ".faa"\
                         +' -outfmt 6 -word_size 5 -max_target_seqs 1 > '+ self.protein_blast_folder +'/{.}.blast'
        cds_list = subprocess.Popen(["ls", self.protein_cds_folder], 
                                stdout=subprocess.PIPE)
        
        parallel_blast = subprocess.Popen(["parallel", "-j", "12",
                        blast_statement],
                        stdin = cds_list.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        parallel_blast.wait()
        print("finished blastp")

        count = 0
        # Removes empties
        for file in os.listdir(self.protein_blast_folder):
            if os.path.getsize(self.protein_blast_folder + "/" + file) == 0:
                count += 1
                os.remove(self.protein_blast_folder + "/" + file)
        print(count)

    def extract_genes(self):
        # NOTE: Consider running the extraction parallel if there are far too many genes.

        self.check_directory_path(self.extracted_gene_folder)
        self.check_directory_path(self.gene_align_folder)

        gene_list = config_parser.get(self.gene, "gene_list")
        gene_list = gene_list.split(",")
        blast_results = os.listdir(self.cds_blast_folder)
        cds_files = os.listdir(self.cds_folder)
        
        for gene in gene_list:
            fasta_gene = b""
            for blast in blast_results:
                blast_list = subprocess.run(["grep", gene, self.cds_blast_folder + "/" + blast], 
                capture_output = True)
                # If the gene could not be found by BLAST
                if blast_list.stdout == b'':
                    continue
                # Saves the name of the gene and the name of the cds file.
                gene_name = blast_list.stdout.split(b'\t')[1]
                gene_name = ">" + str(gene_name, "UTF-8")
                common_file_name = blast.split("_", 3)
                common_file_name = "_".join(common_file_name[0:3])
                for file in cds_files:
                    if file.startswith(common_file_name):
                        cds_file_name = file
                        break        
                faidx_extract = subprocess.run(["faidx", "--regex", gene_name, self.cds_folder + 
                "/" + cds_file_name],
                capture_output=True)
                fasta_gene += faidx_extract.stdout
            self.file_writing_and_alignment(self.extracted_gene_folder + "/" + gene + ".fasta", fasta_gene, self.gene_align_folder + "/" + gene + ".aln")
    
    def extract_protein(self):
        # NOTE: Consider running the extraction parallel if there are far too many genes.
        # NOTE: Very similar to the extract genes.

        self.check_directory_path(self.extracted_protein_folder)
        self.check_directory_path(self.protein_align_folder)

        gene_list = config_parser.get(self.gene, "gene_list")
        gene_list = gene_list.split(",")
        blast_results = os.listdir(self.protein_blast_folder)
        cds_files = os.listdir(self.protein_cds_folder)
        
        for gene in gene_list:
            fasta_gene = b""
            for blast in blast_results:
                blast_list = subprocess.run(["grep", gene, self.protein_blast_folder + "/" + blast], 
                capture_output = True)
                # If the gene could not be found by BLAST
                if blast_list.stdout == b'':
                    continue
                # Saves the name of the gene and the name of the cds file.
                gene_name = blast_list.stdout.split(b'\t')[1]
                gene_name = ">" + str(gene_name, "UTF-8")
                common_file_name = blast.split("_", 3)
                common_file_name = "_".join(common_file_name[0:3])
                for file in cds_files:
                    if file.startswith(common_file_name):
                        cds_file_name = file
                        break        
                faidx_extract = subprocess.run(["faidx", "--regex", gene_name, self.protein_cds_folder
                + "/" + cds_file_name],
                capture_output=True)
                fasta_gene += faidx_extract.stdout
            self.file_writing_and_alignment(self.extracted_protein_folder + "/" + gene + ".fasta", fasta_gene, self.protein_align_folder + "/" + gene + ".aln")

    def file_writing_and_alignment(self, first_location, fasta_gene, second_location):
        gene_file = open(first_location, "wb")
        gene_file.write(fasta_gene)
        gene_file.close()
        alignment = subprocess.run(["mafft", "--quiet", "--adjustdirection", first_location], capture_output=True)
        alignment_file = open(second_location, "wb")
        alignment_file.write(alignment.stdout)
        alignment_file.close()
    
    def check_directory_path(self, path):
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)
    
    def directory_formatter(self, name):
        folder = self.folder_name.split("/")
        folder[-1] = name
        folder = "/".join(folder)
        return folder