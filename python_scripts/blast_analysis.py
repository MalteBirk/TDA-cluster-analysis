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
        self.gene_folder = self.gene + "_" + "Nucleotide_Blast"
        self.protein_folder = self.gene + "_" + "Protein_Blast"
        self.extracted_gene_folder = self.gene +"_" + "Extracted_cds"
        self.gene_align_folder = self.gene + "_" + "Aligned_cds"
        self.extracted_protein_folder = self.gene +"_" + "Extracted_protein_cds"
        self.protein_align_folder = self.gene + "_" + "Aligned_protein_cds"

    def gene_blast(self):
        self.check_directory_path(self.gene_folder)
        
        cds_folder = "cds_fasta" + "_" + self.folder_name
        filename = self.url.split("/")[-1]
        gene_name = filename.split(".")[0]
        # Make user able to change identity here later?
        blast_statement = 'blastn  -perc_identity 50 -qcov_hsp_perc 50 -subject \
                        '+ cds_folder +'/{} -query '+ self.reference_gene + gene_name + ".ffn"\
                         +' -outfmt 6 -word_size 5 > '+ self.gene_folder +'/{.}.blast'
        cds_list = subprocess.Popen(["ls", cds_folder], 
                                stdout=subprocess.PIPE)
        # NOTE: Change the 12 to be altered in the config file.
        parallel_blast = subprocess.Popen(["parallel", "-j", "12",
                        blast_statement],
                        stdin = cds_list.stdout)
        parallel_blast.wait()
        print("finished blastn")

        count = 0
        # Removes empties
        for file in os.listdir(self.gene_folder):
            if os.path.getsize(self.gene_folder + "/" + file) == 0:
                count += 1
                os.remove(self.gene_folder + "/" + file)
        print(count)

    def protein_blast(self):
        self.check_directory_path(self.protein_folder)
        
        protein_folder = "protein_cds_fasta" + "_" + self.folder_name
        filename = self.url.split("/")[-1]
        gene_name = filename.split(".")[0]
        # Make user able to change identity here later?
        blast_statement = 'blastp -subject \
                        '+ protein_folder +'/{} -query '+ self.reference_gene + gene_name + ".faa"\
                         +' -outfmt 6 -word_size 5 -max_target_seqs 1 > '+ self.protein_folder +'/{.}.blast'
        cds_list = subprocess.Popen(["ls", protein_folder], 
                                stdout=subprocess.PIPE)
        
        parallel_blast = subprocess.Popen(["parallel", "-j", "12",
                        blast_statement],
                        stdin = cds_list.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        parallel_blast.wait()
        print("finished blastp")

        count = 0
        # Removes empties
        for file in os.listdir(self.protein_folder):
            if os.path.getsize(self.protein_folder + "/" + file) == 0:
                count += 1
                os.remove(self.protein_folder + "/" + file)
        print(count)

    def extract_genes(self):
        # NOTE: Consider running the extraction parallel if there are far too many genes.
        self.check_directory_path(self.extracted_gene_folder)
        self.check_directory_path(self.gene_align_folder)

        gene_list = config_parser.get(self.gene, "gene_list")
        gene_list = gene_list.split(",")
        blast_results = os.listdir(self.gene_folder)
        cds_files = os.listdir("cds_fasta" + "_" + self.folder_name)
        
        for gene in gene_list:
            fasta_gene = b""
            for blast in blast_results:
                blast_list = subprocess.run(["grep", gene, self.gene_folder + "/" + blast], 
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
                faidx_extract = subprocess.run(["faidx", "--regex", gene_name, "cds_fasta" + "_" + 
                self.folder_name + "/" + cds_file_name],
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
        blast_results = os.listdir(self.protein_folder)
        cds_files = os.listdir("protein_cds_fasta" + "_" + self.folder_name)
        
        for gene in gene_list:
            fasta_gene = b""
            for blast in blast_results:
                blast_list = subprocess.run(["grep", gene, self.protein_folder + "/" + blast], 
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
                faidx_extract = subprocess.run(["faidx", "--regex", gene_name, "protein_cds_fasta" + "_" + 
                self.folder_name + "/" + cds_file_name],
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
