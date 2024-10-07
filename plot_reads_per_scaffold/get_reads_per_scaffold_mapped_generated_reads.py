import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import sys
import os
import random
#GCA_to_name_dict = {"GCA_000149925_1_ASM14992v1_genomic":"Puccinia graminis", "GCA_002873125_1_ASM287312v1_genomic":"Puccinia coronata", "GCA_001624995_1_ASM162499v1_genomic":"Puccinia horiana", "GCA_007896445_1_ASM789644v1_genomic":"Puccinia hordei","GCA_021901695_1_Pst134E36_v1_pri_genomic":"Puccinia striiformis", "GCA_026914185_1_ASM2691418v1_genomic":"Puccinia triticina", "GCA_004348175_1_PNOVO_Noble_KM_genomic":"Puccinia novopanici", "GCA_019395275_1_ASM1939527v1_genomic":"Puccinia brachypodii", "GCA_001263375_1_ASM126337v1_genomic":"Puccinia sorghi"}
GCA_to_name_dict = {"GCA_000149925_1_ASM14992v1_genomic":"Puccinia graminis", "GCA_002873125_1_ASM287312v1_genomic":"Puccinia coronata", "GCA_001624995_1_ASM162499v1_genomic":"Puccinia horiana", "GCA_007896445_1_ASM789644v1_genomic":"Puccinia hordei","GCA_021901695_1_Pst134E36_v1_pri_genomic":"Puccinia striiformis", "GCA_026914185_1_ASM2691418v1_genomic":"Puccinia triticina", "GCA_004348175_1_PNOVO_Noble_KM_genomic":"Puccinia novopanici", "GCA_019395275_1_ASM1939527v1_genomic":"Puccinia brachypodii", "GCA_001263375_1_ASM126337v1_genomic":"Puccinia sorghi", "GCA_921294245_1_PGI_MELIAE_v1_genomic":"Brassicogethes aeneus", "GCA_019925205_1_ASM1992520v1_genomic":"Metopolophium dirhodum","GCA_000142945_1_ASM14294v1_genomic":"Phytophthora infestans", "GCA_016880985_1_ASM1688098v1_genomic":"Phytophthora idaei", "GCA_033557915_1_HumboldtARI_Pten_1_0_genomic":"Phytophthora tentaculata", "GCA_000247585_2_PP_INRA-310_V2_genomic":"Phytophthora nicotianae", "GCA_020882245_1_ASM2088224v1_genomic":"Rhopalosiphum padi", "GCA_003676215_3_ASM367621v3_genomic":"Rhopalosiphum maidis", "GCA_036289425_1_JGU_Rn_01_genomic":"Rhopalosiphum nymphaeae", "GCA_019425605_1_ASM1942560v1_genomic":"Sitobion avenae", "GCA_008086715_1_ASM808671v1_genomic":"Sitobion miscanthi"}
#File is generated using calculate_nr_reads_to_subsample.py
number_to_subsample_file = "nr_reads_to_subsample_from_generated.txt"
#Create a dictionary which contain all taxid for species for which reads have been mapped to
#These taxid function as the key to access a dictionary which contains information 
nr_to_subsample_dict = {}
with open(number_to_subsample_file, "r") as sample_file:
    for line in sample_file:
        if line[0] == "R":
            current_taxid = line.split()[-1].strip("\n")[:-1]
            nr_to_subsample_dict[current_taxid] = {}
        else:
            nr = line.strip("\n").split()[-1]
            GCA = line.split()[-2].strip(":")
            nr_to_subsample_dict[current_taxid][GCA] = int(float(nr))

generated_GCA_dict = {"GCA_019925205_1_ASM1992520v1_genomic":["GCA_019925205_1_ASM1992520v1_genomic"], "GCA_921294245_1_PGI_MELIAE_v1_genomic":["GCA_921294245_1_PGI_MELIAE_v1_genomic"], "GCA_000142945_1_ASM14294v1_genomic": ["GCA_000142945_1_ASM14294v1_genomic", "GCA_000247585_2_PP_INRA-310_V2_genomic", "GCA_016880985_1_ASM1688098v1_genomic", "GCA_033557915_1_HumboldtARI_Pten_1_0_genomic"], "GCA_020882245_1_ASM2088224v1_genomic": ["GCA_003676215_3_ASM367621v3_genomic", "GCA_020882245_1_ASM2088224v1_genomic", "GCA_036289425_1_JGU_Rn_01_genomic"], "GCA_019425605_1_ASM1942560v1_genomic": ["GCA_008086715_1_ASM808671v1_genomic", "GCA_019425605_1_ASM1942560v1_genomic"], "GCA_007896445_1_ASM789644v1_genomic":["GCA_000149925_1_ASM14992v1_genomic", "GCA_001263375_1_ASM126337v1_genomic", "GCA_001624995_1_ASM162499v1_genomic", "GCA_002873125_1_ASM287312v1_genomic", "GCA_004348175_1_PNOVO_Noble_KM_genomic", "GCA_007896445_1_ASM789644v1_genomic", "GCA_019395275_1_ASM1939527v1_genomic", "GCA_021901695_1_Pst134E36_v1_pri_genomic", "GCA_026914185_1_ASM2691418v1_genomic"],"GCA_002873125_1_ASM287312v1_genomic":["GCA_000149925_1_ASM14992v1_genomic", "GCA_001263375_1_ASM126337v1_genomic", "GCA_001624995_1_ASM162499v1_genomic", "GCA_002873125_1_ASM287312v1_genomic", "GCA_004348175_1_PNOVO_Noble_KM_genomic", "GCA_007896445_1_ASM789644v1_genomic", "GCA_019395275_1_ASM1939527v1_genomic", "GCA_021901695_1_Pst134E36_v1_pri_genomic", "GCA_026914185_1_ASM2691418v1_genomic"]}

pair_list = [["44670", "GCA_019925205_1_ASM1992520v1_genomic", "Metopolophium dirhodum", "insert_176_insstd_88"]]#, ["4787", "GCA_000142945_1_ASM14294v1_genomic", "Phytophthora infestans", "insert_196_insstd_96"], ["40932", "GCA_020882245_1_ASM2088224v1_genomic", "Rhopalosiphum padi", "insert_180_insstd_86"], ["44664", "GCA_019425605_1_ASM1942560v1_genomic", "Sitobion avenae", "insert_215_insstd_99"],["27344", "GCA_002873125_1_ASM287312v1_genomic", "Puccinia coronata", "insert_245_insstd_115"], ["27345", "GCA_007896445_1_ASM789644v1_genomic", "Puccinia hordei", "insert_203_insstd_95"], ["1431903", "GCA_921294245_1_PGI_MELIAE_v1_genomic", "Brassicogethes aeneus", "insert_170_insstd_87"]]

for pair in pair_list:
    tax_id = pair[0]
    name = pair[2]
    GCA = pair[1]
    #The following file is genrated using get_contig_lengths_from_fasta.py and contains information about scaffold names and lengths of a genome
    contig_lengths_file = "scaffold_summaries/{}.txt".format(GCA) 
    generated_info = pair[3]
    minid = "97"
    amb = "toss"
    #Get 
    generated_GCA_list = generated_GCA_dict[GCA]
    for generated_GCA in generated_GCA_list:
        bam_file = "map_files/bam/{}_generated_perfect_NEAT_10x_{}_mapped_to_{}_amb_{}_minid_0_{}.bam".format(generated_GCA, generated_info, tax_id, amb, minid)
        out_bam = "map_files/bam/{}_generated_perfect_NEAT_10x_{}_mapped_to_{}_amb_{}_minid_0_{}_subsample_{}.bam".format(generated_GCA, generated_info, tax_id, amb, minid, str(nr_to_subsample_dict[tax_id][generated_GCA]))
        #Count total number of mapped pair end reads
        map_count = 0
        with open(bam_file, "r") as bam_f:
            for line in bam_f:
                if line[0] != "@":
                    nextline = next(bam_f)
                    map_count += 1
        print(map_count)
        if map_count >= nr_to_subsample_dict[tax_id][generated_GCA]:
            #Randomly decide what pair end reads to subsample
            read_line_to_subsample = random.sample(range(1,map_count+1), nr_to_subsample_dict[tax_id][generated_GCA])
            line_count = 0
            contig_count_dict = {}
            with open(bam_file, "r") as bam_f:
                with open(out_bam, "w") as out_bam:
                    for line in bam_f:
                        if line[0] != "@":
                            nextline = next(bam_f)
                            line_count += 1
                            #If this is one of the random pair end reads to subsample, increase count of number of mapped reads for the scaffold the read pair mapped
                            if line_count in read_line_to_subsample:
                                contig1 = line.strip("\n").split("\t")[2].split()[0]
                                contig2 = nextline.strip("\n").split("\t")[2].split()[0]
                                if contig1 == contig2:
                                    if contig1 in contig_count_dict.keys():
                                        contig_count_dict[contig1] += 1
                                    else:
                                        contig_count_dict[contig1] = 1
            #Go through all scaffolds and get their, names, lengths, and number of mapped reads
            contig_names = []
            contig_lengths = []
            contig_reads = []
            
            with open(contig_lengths_file, "r") as cl_file:
                for line in cl_file:
                    split_line = line.strip("\n").split(",")
                    contig = split_line[0]
                    length = int(split_line[1])
                    try:
                        reads = contig_count_dict[contig]
                    except:
                        reads = 0
                    contig_lengths.append(length)
                    contig_reads.append(reads)
                    contig_names.append(contig)
            #Save files if needed for future plotting
            np_scaff_length = np.array(contig_lengths)
            np_nr_mapped = np.array(contig_reads)
            np.savetxt("numpy_files/reads_per_contig_{}_subsampled_generated_from_{}_mapped_to_{}_minid{}_amb_{}.txt".format(str(nr_to_subsample_dict[tax_id][generated_GCA]),generated_GCA,tax_id, minid, amb), np_nr_mapped)
            np.savetxt("numpy_files/contig_length_{}_subsampled_generated_from_{}_mapped_to_{}.txt".format(str(nr_to_subsample_dict[tax_id][generated_GCA]),generated_GCA,tax_id), np_scaff_length)
            #Used to place label when ploting
            min_reads = np_nr_mapped.min()
            max_contig = np_scaff_length.max()
            #Calculate pearson
            res = stats.pearsonr(contig_lengths, contig_reads)
            pearson = res.statistic
            #Plot number of mapped reads per scaffold (y-axsis) vs scaffold length (x-axsis)
            plt.text(max_contig*0.78, min_reads, 'PCC: %.2f' % pearson, fontsize=16)
            plt.scatter(contig_lengths, contig_reads, color="black")
            plt.xlabel("Scaffold length")
            plt.ylabel("Number of mapped reads")
            plt.title("{}".format(str(nr_to_subsample_dict[tax_id][generated_GCA])), fontsize=18)
            plt.savefig("coverage_plots/reads_vs_scaffold_length/reads_vs_scaffold_length_{}_subsampled_reads_from_{}_generated_NEAT_error{}_mapped_to_{}_minid0{}_ambigious_{}.png".format(str(nr_to_subsample_dict[tax_id][generated_GCA]),generated_GCA,generated_info ,tax_id, minid, amb))
            plt.clf()