import os
import configparser
import subprocess
#import matplotlib.pyplot as plt

from blast_analysis import BlastAnalysis

config_parser = configparser.ConfigParser()
config_parser.read("../config.ini")

class OrganisationAnalysis:

    def __init__(self, gene, folder_name):
        self.gene = gene
        self.folder_name = folder_name
        folder_name = folder_name.split("/")
        self.gene_list = config_parser.get(self.gene, "gene_list").split(",")
        self.length = int(config_parser.get(self.gene, "approximate_gene_length"))

        self.blast_functions_object = BlastAnalysis(self.gene, self.folder_name)
        self.protein_cds_folder = self.blast_functions_object.directory_formatter("protein_cds_fasta" + "_" + folder_name[-1])
        self.extracted_protein_folder = self.blast_functions_object.directory_formatter(gene + "_" + "Extracted_protein_cds")

        self.organisation_analysis_folder = self.blast_functions_object.directory_formatter(self.gene + "_gene_organisation")

    def extract_gene_organisation(self):
        self.blast_functions_object.check_directory_path(self.organisation_analysis_folder)
        all_protein_cds = os.listdir(self.protein_cds_folder)
        all_extracted_proteins = os.listdir(self.extracted_protein_folder)
        for file in all_protein_cds:
            genome_file = open(self.protein_cds_folder + "/" + file, "r")
            first_line = genome_file.readline()
            organism_name = first_line.split("|")[0]
            gene_name_list = []
            # Extracts the correct gene name from the extracted proteins
            for gene_file in all_extracted_proteins:
                single_gene_file = open(self.extracted_protein_folder + "/" + gene_file, "r")
                for line in single_gene_file:
                    if line.startswith(organism_name):
                        single_gene_name = gene_file.rsplit(".", 1)[0]
                        gene_name = line.rstrip()
                        gene_name_list.append(gene_name)
                        gene_name_list.append(single_gene_name)
            # Goes back to beginning of file, and extract the lines with locations for the genes.
            genome_file.seek(0)
            #organism_name = organism_name.replace("/","")
            if os.path.exists(self.organisation_analysis_folder + "/" + file):
                os.remove(self.organisation_analysis_folder + "/" + file)
            outfile = open(self.organisation_analysis_folder + "/" + file, "w")
            #outfile = open(self.organisation_analysis_folder + "/" + organism_name[1:].rsplit("_", 1)[0] + ".txt", "w")

            for line in genome_file:
                if line.startswith(">"):
                    potential_gene_name = line.split(" ")[0]
                    if potential_gene_name in gene_name_list:
                        gene_index = gene_name_list.index(potential_gene_name)
                        outfile.write(line.rstrip() + " " + str(gene_name_list[gene_index + 1]) + "\n")
                        gene_name_list.remove(gene_name_list[gene_index + 1])
                        gene_name_list.remove(potential_gene_name)
                elif gene_name_list == []:
                    break
                else:
                    pass

    def find_organisation(self):
        all_organisations = os.listdir(self.organisation_analysis_folder)
        complete_cluster_list = []
        semi_complete_cluster_list = []
        incomplete_cluster_list = []
        for file in all_organisations:
            min_num = None
            max_num = None
            self.coupled = []
            self.name_and_pos_list = []
            organisation_file = open(self.organisation_analysis_folder + "/" + file, "r")
            found_on_same_contig = False
            contig_list = []
            larger_clusters = []
            for line in organisation_file:
                line_list = line.split(" ")
                # Different wrangling of data result in specific case for Jyllinge strains
                if line_list[0].startswith(">Jyllinge"):
                    gene_number = line_list[0].rsplit("_")[-3].split(".")[0]
                # Will be second to last position in normal cases
                else:
                    gene_number = line_list[0].rsplit("_")[-2]
                contig = line_list[0].split("|")[1].split(".")[0]
                species_name = line_list[0]
                contig_list.append(contig)
                for information in line_list:
                    if information.startswith("[location"):
                        location = information.replace("[","").replace("]","").replace(")","").replace("<", "").replace(">", "")
                        if "complement" in location:
                            location_of_gene = location.rsplit("(")[-1].split("..")
                        else:
                            location_of_gene = location.rsplit("=")[-1].split("..")
                        if min_num is None:
                            min_num = min(location_of_gene)
                        else:
                            if min(location_of_gene) < min_num:
                                min_num = min(location_of_gene)
                        if max_num is None:
                            max_num = max(location_of_gene)
                        else:
                            if max(location_of_gene) > max_num:
                                max_num = max(location_of_gene)
                
                
                # NOTE: Figure out the ordering and get it right.
                self.name_and_pos_list.append(line_list[-1].rstrip())
                self.name_and_pos_list.append(gene_number)
            
            if len(contig_list) > 0 :
                # Expects them to be on same contig from previous findings.
                if line_list[0].startswith(">Jyllinge"):
                    found_on_same_contig = True
                else:
                    found_on_same_contig = contig_list.count(contig_list[0]) == len(contig_list)

            opereon_length = int(max_num) - int(min_num)

            for i in range(0, len(self.name_and_pos_list), 2):
                self._couple_finder(1)
                self._couple_finder(-1)
            if len(self.coupled) != 0:
                larger_clusters = [self.coupled[0]]
            
            for i in range(0, len(self.coupled) - 1):
                overlap = self.overlap((self.coupled[i]), (self.coupled[i + 1]))
                if overlap != []:
                    pair_list = self.coupled[i].copy()
                    overlapping_list = self.coupled[i + 1].copy()
                    overlapping_list.remove(overlap[0])
                    combined_list = pair_list
                    combined_list.extend(overlapping_list)
                    flag = False
                    for gene in combined_list:
                        for cluster in larger_clusters:
                            if gene in cluster:
                                cluster.extend(combined_list)
                                flag = True
                    if flag == False:
                        larger_clusters.append(combined_list)
            only_uniques_clusters = []  
            for i in range(len(larger_clusters)):
                reduced_cluster = list(set(larger_clusters[i]))
                only_uniques_clusters.append(reduced_cluster)

            gene_list = []
            flag = False
            for cluster in only_uniques_clusters:
                for i in cluster:
                    if i not in gene_list:
                        gene_list.append(i)
                    else:
                        print("Warning, duplicate found in", file)
                        flag = True
                    if flag == True:
                        break
                if flag == True:
                    break
            
            if found_on_same_contig == False:
                incomplete_cluster_list.append([species_name, file, "NA"])
            elif opereon_length > int(self.length) * 5:
                incomplete_cluster_list.append([species_name, file, opereon_length])
            elif len(only_uniques_clusters) > 1:
                semi_complete_cluster_list.append([species_name, file, opereon_length])
            else:
                complete_cluster_list.append([species_name, file, opereon_length])
        
        # For now, just species name.
        if os.path.exists("complete_cluster_list"):
            os.remove("complete_cluster_list")
        outfile = open("complete_cluster_list", "w")
        outfile.write("species\n")
        for i in complete_cluster_list:
            species_name = self._cluster_name_wrangler(i[0])
            outfile.write(species_name + "\n")
        
        if os.path.exists("semi_complete_cluster_list"):
            os.remove("semi_complete_cluster_list")
        outfile = open("semi_complete_cluster_list", "w")
        outfile.write("species\n")
        for i in semi_complete_cluster_list:
            species_name = self._cluster_name_wrangler(i[0])
            outfile.write(species_name + "\n")
        
        if os.path.exists("incomplete_cluster_list"):
            os.remove("incomplete_cluster_list")
        outfile = open("incomplete_cluster_list", "w")
        outfile.write("species\n")
        for i in incomplete_cluster_list:
            species_name = self._cluster_name_wrangler(i[0])
            outfile.write(species_name + "\n")
    
    def _couple_finder(self, orientation):
        for i in range(0, len(self.name_and_pos_list), 2):
            if str(int(self.name_and_pos_list[i + 1]) + orientation) in self.name_and_pos_list:
                index = self.name_and_pos_list.index(str(int(self.name_and_pos_list[i + 1]) + orientation))
                # Checks if the couple is already there
                flag = False
                for j in self.coupled:
                    if self.name_and_pos_list[i] in j and self.name_and_pos_list[index - 1] in j:
                        flag = True
                if flag == True:
                    continue
                self.coupled.append([self.name_and_pos_list[i],self.name_and_pos_list[index - 1]])
    
    def _cluster_name_wrangler(self, cluster):
        species_name = cluster.split("_lcl")[0].replace(">","")
        if "(" in species_name:
            species_name = species_name.split("(")[0]
        return species_name

    def overlap(self, list1, list2):
        return list(set(list1) & set(list2))
