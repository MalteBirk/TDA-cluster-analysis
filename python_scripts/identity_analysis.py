import os
import configparser
import subprocess
import plotly_express as px

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
        self.protein_blast_folder = self.blast_functions_object.directory_formatter(self.gene + "_Protein_Blast")
        self.reference_folder = self.blast_functions_object.directory_formatter("reference" + self.gene)
        self.identity_analysis_folder = self.blast_functions_object.directory_formatter("reference" + self.gene)
        self.blast_identities_figures_folder = self.blast_functions_object.directory_formatter(self.gene + "_figures/blast_identities")
        self.tree_figures_folder = self.blast_functions_object.directory_formatter(self.gene + "_figures/trees")
        self.operon_identity_analysis_folder = self.blast_functions_object.directory_formatter(self.gene + "_figures/operon_identities")

        self.gene_align_folder = self.blast_functions_object.directory_formatter(gene + "_" + "aligned_cds")
        self.protein_align_folder = self.blast_functions_object.directory_formatter(gene + "_" + "aligned_protein_cds")
        self.aligned_housekeeping_path = self.blast_functions_object.directory_formatter(self.gene + "_" + "aligned_housekeeping_genes")
        self.aligned_housekeeping_protein_path = self.blast_functions_object.directory_formatter(self.gene + "_" + "aligned_housekeeping_proteins")

    def blast_to_length_comparison(self):
        self.blast_functions_object.check_directory_path(self.blast_identities_figures_folder)
        self.blast_functions_object.check_directory_path(self.identity_analysis_folder)
        self.blast_functions_object.check_directory_path(self.operon_identity_analysis_folder)
        url = config_parser.get(self.gene, "mibig_url")
        filename = url.split("/")[-1]
        # NOTE: Many more options here. When more than TDA gene is tested investigate this!!
        if filename.endswith(".gbk"):
            self.filename = filename.split(".")[0]
        self._ref_gene_analyser(".ffn", self.cds_blast_folder, "_nucleotide")
        self._ref_gene_analyser(".faa", self.protein_blast_folder, "_protein")


    def tree_maker(self):
        self.blast_functions_object.check_directory_path(self.tree_figures_folder)
        # MANUEL for now. Rscript messing up
        #self._tree_figure_writer(self.gene_align_folder, "-nt")
        #self._tree_figure_writer(self.protein_align_folder, None)
        #self._tree_figure_writer(self.aligned_housekeeping_path, "-nt")
        self._tree_figure_writer(self.aligned_housekeeping_protein_path, None)

    def blast_identity_to_operons(self):
        self._operon_plot("protein")

    def _tree_figure_writer(self, folder, data_type):
        alignments = os.listdir(folder)
        
        for file in alignments:
            if data_type == "-nt":
                tree_outfile = open(file.split(".")[0] + "_" + "tree" + ".tree", "wb")
                tree_file = subprocess.run(["fasttree", "-nt", "-quiet", folder + "/" + file], capture_output = True)
            else:
                tree_outfile = open(file.split(".")[0] + "Protein_" + "tree" + ".tree", "wb")
                tree_file = subprocess.run(["fasttree", "-quiet", folder + "/" + file], capture_output = True)
            tree_outfile.write(tree_file.stdout)
            tree_outfile.close()
        
        
        self._tree_data_wrangler()
        self._R_tree_analysis()
    
    def _tree_data_wrangler(self):
        current_directory = os.listdir()
        for file in current_directory:
            if file.endswith(".tree"):
                # Genus name only
                tree_file = open(file, "r")
                for line in tree_file:
                    all_data = line.split(",")
                df_name = file.split("_")[0]
                filename = "".join(file.split(".")[:-1])
                wrangled_tree_file = open(filename + ".wrangled", "w")
                data_frame = open(df_name + ".tsv", "w")
                count = 0
                for tree_info in all_data:
                    
                    name = "_".join(tree_info.split("_", 2)[:2])
                    if "lcl" in tree_info:
                        other_name = "".join(tree_info.split("_lcl", 1)[:1])
                        last_of_name = "_".join(other_name.split("_", 2)[2:])
                    else:
                        other_name = "".join(tree_info.split(":", 1)[0])
                        last_of_name = "_".join(other_name.split("_", 2)[2:])
                    tree_value = tree_info.split(":", 1)[1]
                    genus_name = name.replace(")","").replace("(","")
                    combined_genus_name = genus_name + "_" + last_of_name
                    # For the weird case
                    second_flag = False
                    if (name + "_" + last_of_name).startswith("Phaeobacter_italicus_strain_CECT_7645"):
                        last_of_name = last_of_name + "_" + str(count)
                        combined_genus_name = genus_name + "_" + last_of_name
                        count += 1
                        operon = "Incomplete"
                        second_flag = True

                    if tree_value.endswith(";") or tree_value.endswith("\n"):
                        wrangled_tree_file.write(name + "_" + last_of_name + ":" + tree_value)
                    else:
                        wrangled_tree_file.write(name + "_" + last_of_name + ":" + tree_value + ",")
                    
                    flag = False

                    if flag == False:
                        cluster_file = open("complete_cluster_list", "r")
                        for line in cluster_file:
                            if line.rstrip() == combined_genus_name:
                                operon = "Complete"
                                flag = True

                    if flag == False:
                        cluster_file = open("semi_complete_cluster_list", "r")
                        for line in cluster_file:
                            if line.rstrip() == combined_genus_name:
                                operon = "SemiComplete"
                                flag = True
                    
                    if flag == False:
                        cluster_file = open("incomplete_cluster_list", "r")
                        for line in cluster_file:
                            if line.rstrip() == combined_genus_name:
                                operon = "Incomplete"
                                flag = True
                    
                    if flag == False and second_flag == False:
                        print(combined_genus_name)
                        operon = "NA"
                    data_frame.write(combined_genus_name + "\t" + genus_name.split("_")[0] + "\t" + operon + "\n")
            
                # Inserting for known producers
                known_producer_list = open("known_producer_file.txt", "r")
                producer_list = []
                for i in known_producer_list:
                    i = i.rstrip()
                    producer_list.append(i)
                data_frame.close()
                data_frame = open(df_name + ".tsv", "r")
                known_producer_data_frame = open(df_name + "KnownProducer" + ".tsv", "w")
                #known_producer_data_frame.write("Name" + "\t" + "Genus" + "\t" + "Cluster" + "\t" + "KnownProducer" + "\n")
                count = 0
                for line in data_frame:
                    known_flag = False
                    non_producer_flag = False
                    known_producer = ""
                    written = False
                    for producer in producer_list:
                        if "_" + producer + "_" in line.split("\t")[0] or "_" + producer + " " in line.split("\t")[0] or " " + producer + "_" in line.split("\t")[0]:
                            known_producer_data_frame.write(line.rstrip() + known_producer + "\n")
                            written = True
                        if known_flag == True and known_producer in line.split("\t")[0]:
                            known_producer = "\t" + "Producer"
                        if non_producer_flag == True and known_producer in line.split("\t")[0]:
                            known_producer = "\t" + "Non-Producer"
                        if producer.startswith(">producers"):
                            known_flag = True
                        if producer.startswith(">non_producers"):
                            known_flag = False
                            non_producer_flag = True
                    if written == False:
                        known_producer_data_frame.write(line.rstrip() + "\t" + "Unknown" + "\n")

    def _ref_gene_analyser(self, file_ending, folder, info_type):
        # NOTE: Something off with tdaA for proteins!!!!
        ref_gene = open(self.reference_folder + "/" + self.filename + file_ending, "rb")
        seq = ""
        for line in ref_gene:
            if line.startswith(b">"):
                name = line.rstrip()
            else:
                seq = line.rstrip()
                self._identity_and_length(name, len(seq), folder, info_type)
        #self._R_analysis()

    def _identity_and_length(self, name, seq_length, folder, info_type):
        blast_files = os.listdir(folder)
        if os.path.exists(str(name[1:], "UTF-8") + "_data_frame"):
            os.remove(str(name[1:], "UTF-8") + "_data_frame")
        data_frame_file = open(str(name[1:], "UTF-8") + info_type + "_data_frame", "w")
        for file in blast_files:
            blast_file = open(folder + "/" + file, "r")
            for blast_result in blast_file:
                if blast_result.startswith(str(name[1:], "UTF-8")):
                    identity = blast_result.split("\t")[2]
                    length = (int(blast_result.split("\t")[3])/seq_length)*100
                    # Expects it as the first name
                    genus_name = blast_result.split("\t")[1].split("_")[0]
                    # TODO: Look more in depth at this once more species start to emerge
                    species_name = blast_result.split("\t")[1].split("|")[0]
                    data_frame_file.write(file + "\t" 
                                          + str(identity) + "\t" 
                                          + str(length) + "\t" 
                                          + genus_name + "\t" 
                                          + species_name + "\n")
        data_frame_file.close()

    def _R_analysis(self):
        if os.path.exists("info_data_frame"):
            os.remove("info_data_frame")
        self.blast_functions_object.check_directory_path(self.blast_identities_figures_folder)
        # Information might be redundant.
        information_file = open("R_information_file", "w")
        information_file.write(self.blast_identities_figures_folder + "/")
        # NOTE: Remove capture_output if you want to get messages from R.
        subprocess.run(["Rscript", "../R_scripts/point_plot_and_histogram.R"])
        
        current_working_directory = os.listdir(".")
        for file in current_working_directory:
            if file.endswith("_data_frame") or file == "R_information_file":
                os.remove(file)
            if file.endswith("_Rfig.png"):
                os.rename(file, self.blast_identities_figures_folder + "/" + file)
    
    def _R_tree_analysis(self):
        information_file = open("R_information_file", "w")
        information_file.write(self.tree_figures_folder + "/")
        # NOTE:!!! Try to remove | in R to get more informative trees.
        # NOTE: !!! This will only work, since the metadata file keeps that part of the name.
        #subprocess.run(["Rscript", "../R_scripts/tree_maker.R"])
        exit(1)
        current_working_directory = os.listdir(".")
        for file in current_working_directory:
            # NOTE: Mute here if you want to look at tree files. Might be necessary to delete some information.
            if file.endswith("_tree.tree") or file.endswith("_tree_protein.tree") or file == "R_information_file" or file.endswith(".tsv") or file.endswith(".wrangled"):
                os.remove(file)
            if file.endswith("_Rfig.png"):
                os.rename(file, self.tree_figures_folder + "/" + file)

    def _operon_plot(self, type):
        all_files = os.listdir()
        gene_list = config_parser.get(self.gene, "gene_list").split(",")
        for gene in gene_list:
            new_data_frame = open(gene + "_identity_operon_data_frame", "w")
            for file in all_files:
                if file.startswith(gene + "_" + type) and file.endswith("data_frame"):
                    data_frame = open(file, "r")
                if file.startswith(gene) and file.endswith("KnownProducer.tsv"):
                    operon_file = open(file, "r")
                    operon_frame = []
                    for line in operon_file:
                        operon_frame.append(line)
            for line in data_frame:
                #print(line)
                identity = line.split("\t")[1]
                name = line.split("\t")[4].split("_lcl")[0]
                genus = line.split("\t")[3]
                    
                for entry in operon_frame:
                   # print(entry.split("\t")[0])
                    #if entry.split("\t")[0] == "Phaeobacter_inhibens_2.10" and name == "Phaeobacter_inhibens_2.10":
                    #    print("my name", name)
                    #    print("entry", entry.split("\t")[0])
                    if entry.split("\t")[0] == name:
                        operon = entry.split("\t")[2]
                        #print(entry.split("\t")[0])
                        #print(operon)
                        
                new_data_frame.write(name + "\t" + genus + "\t" + identity + "\t" + operon + "\n")        

            data_frame.close()
            operon_file.close()
            new_data_frame.close()
