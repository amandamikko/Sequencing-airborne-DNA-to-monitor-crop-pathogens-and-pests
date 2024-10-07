import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import sys
import os

#Specify species to evaluate   
pair_list = [["44670", "GCA_019925205_1_ASM1992520v1_genomic", "Metopolophium dirhodum"]] #, ["1431903", "GCA_921294245_1_PGI_MELIAE_v1_genomic", "Brassicogethes aeneus"], ["4787", "GCA_000142945_1_ASM14294v1_genomic", "Phytophthora infestans"], ["40932", "GCA_020882245_1_ASM2088224v1_genomic", "Rhopalosiphum padi"], ["44664", "GCA_019425605_1_ASM1942560v1_genomic", "Sitobion avenae"], ["27344", "GCA_002873125_1_ASM287312v1_genomic", "Puccinia coronata"], ["27345", "GCA_007896445_1_ASM789644v1_genomic", "Puccinia hordei"]]
#Go through all species and evaluate coverage depth of the mapped air filter reads across the genome
#results will be written to a output file specified by the variable cov_depth_file
for pair in pair_list:
    tax_id = pair[0]
    #files should be generated using snakefile
    input = "map_files/basecov/basecov_merged_mapped_to_{t}_amb_toss_minid_0_97.txt".format(t=tax_id)  
    #output file                                             
    cov_depth_file = "coverage_coord/coverage_coordinates_merged_mapped_to_{t}_amb_toss_minid_0_97_summarize_stack_hight.txt".format(t=tax_id)                                             
    scaffold = ""                                                            
    found_start = False                                                      
    start_coord = ""                                                         
    coords = []             
    #Got through a basecov-file to find all areas with coverage and determine the maximum coverage of each covered  area
    prev_scaffold_name = ""
    with open(input, "r") as basecov:                                        
        for line in basecov:                                                 
            if line[0] != "#":                                                                         
                split_line = line.strip("\n").split()                        
                scaffold_name = split_line[0].split()[0]   
                #If we are no longer looking at the same scaffold and previous coordinate have coverage, add last region of previous scaffold to list of regions with coverage
                if prev_scaffold_name != scaffold_name and found_start == True:
                    coords.append([curr_scaff, str(start_coord), str(prev_coord), str(max_stack_hight)])
                    found_start = False               
                #If we have not found a start but there is coverage on the current coordinate, start a new region and document maximum coverage
                if int(split_line[-1]) != 0 and found_start == False:        
                    found_start = True                                       
                    start_coord = int(split_line[-2])
                    curr_scaff = scaffold_name    
                    max_stack_hight = int(split_line[-1])
                #If we have found a start and there is coverage on the current coordinate, then look at weither the surrent coverage is higher than the current maximum and update which coordinate was the previous one (for next iteration)                
                if int(split_line[-1]) != 0 and found_start == True:
                    prev_coord = int(split_line[-2]) 
                    if  int(split_line[-1]) > max_stack_hight:
                        max_stack_hight = int(split_line[-1])  
                # If we have found a start and there is no coverage on the current coordinate, append start coordinate and previous coordinate to coordinate list                          
                if int(split_line[-1]) == 0 and found_start == True:         
                    found_start = False                                      
                    stop_coord = prev_coord
                    if scaffold_name == curr_scaff:                
                        coords.append([split_line[0], str(start_coord), str(stop_coord), str(max_stack_hight)])
                    else:
                        coords.append([curr_scaff, str(start_coord), str(stop_coord), str(max_stack_hight)])
                prev_scaffold_name = scaffold_name

    #When done, if applicable, add regions which were potentially at the end of the last scaffold
    if found_start == True:                     
        found_start = False                                                  
        stop_coord = int(split_line[-2])                                     
        coords.append([split_line[0], str(start_coord), str(stop_coord), str(max_stack_hight)])    
    #Go through scaffolds and regions and write them to output                          
    with open(cov_depth_file, "w") as out_file:                                
        out_file.write("scaffold_name,start_coordinate,stop_coordinate,max_stack_hight\n") 
        for coord in coords:                                                                                                   
            out_file.write(coord[0] + "," + coord[1] + "," + coord[2] + "," + coord[3] + "\n")

#NEW
#Get number of mapped air filter reads per scaffold (Each counted read can not map to the same location as more than one other read) and plot this value against scafold length
#output file generated so that it can be used by calculate_nr_reads_to_subsample.py
summary_out_file = "Summary_total_amount_of_mapped_reads_per_species.txt"

for pair in pair_list:
    #Specify parameters and file names
    tax_id = pair[0]
    name = pair[2]
    GCA = pair[1]
    minid = "97"
    amb = "toss"
    out_bam = "map_files/bam/merged_all_weeks_mapped_to_{}_minid0{}_ambigious_{}_only_reads_from_stacks_of_max_2_hight.sam".format(tax_id, minid, amb)
    coverage_info_file = "coverage_coord/coverage_coordinates_merged_mapped_to_{}_amb_toss_minid_0_97_summarize_stack_hight.txt".format(tax_id) #For all weeks mapping
    scaffold_lengths_file = "scaffold_summaries/{}.txt".format(GCA) 
    bam_file = "map_files/bam/merged_mapped_to_{}_amb_{}_minid_0_{}.sorted.sam".format(tax_id, amb, minid)
    cov_dict = {}
    min_overlap = 0
    #Get information about where reads were mapped from file produced by analyze_stacks_from_mapping.py
    #Create a dict and store scaffold names paired with another dictionary containing information about maximum coverage at different coordinates of the scaffold
    with open(coverage_info_file, "r") as coverage_info_f:
        coverage_info_f.readline()
        for line in coverage_info_f:
            scaff = line.split(",")[0]
            start = int(line.split(",")[1])
            stop = int(line.split(",")[2])
            cov = int(line.split(",")[3].strip("\n"))
            if scaff not in cov_dict.keys():
                cov_dict[scaff] = {(start, stop): cov}
            else:
                cov_dict[scaff][(start, stop)] = cov
    #Go through bam file and count reads per scaffold and make sure that each counted read pair is not from a region with more than two reads coverage
    reads_per_scaffold = {}
    total_reads = 0
    with open(out_bam, "w") as out_f:
        with open(bam_file, "r") as f:
                for line in f:
                    if line[0] != "@":
                        nextline = next(f)
                        #Get mapping for forward and reverese read
                        scaff_name_map = line.strip("\n").split("\t")[2].split()[0]
                        scaff_name_map2 = nextline.strip("\n").split("\t")[2].split()[0]
                        #Make sure forward and revere map to same scaffold
                        if scaff_name_map == scaff_name_map2:
                            start_coord1 = int(line.strip("\n").split("\t")[3])
                            read_length1 = line.strip("\n").split("\t")[8]
                            stop_coord1 = start_coord1 + int(read_length1)
                            #Establish witch coordinate is the start coordinate of the first read
                            if stop_coord1 < start_coord1:
                                temp_start = start_coord1
                                temp_stop = stop_coord1
                                start_coord1 = temp_stop
                                stop_coord1 = temp_start
                            start_coord2 = int(nextline.strip("\n").split("\t")[3])
                            read_length2 = nextline.strip("\n").split("\t")[8]
                            stop_coord2 = start_coord2 + int(read_length2)
                            #Establish witch coordinate is the start coordinate of the secound read
                            if stop_coord2 < start_coord2:
                                temp_start = start_coord2
                                temp_stop = stop_coord2
                                start_coord2 = temp_stop
                                stop_coord2 = temp_start
                            not_found1 = True
                            not_found2 = True
                            #Go through regions to which reads have mapped and find to which region each read has mapped
                            for mappable_region in cov_dict[scaff_name_map].keys():
                                start = mappable_region[0]
                                stop = mappable_region[1]
                                if start_coord1 >= start + min_overlap and start_coord1 <= stop - min_overlap:
                                    cov1 = cov_dict[scaff_name_map][mappable_region] 
                                    map_region1 = mappable_region
                                    not_found1 = False
                                elif stop_coord1 >= start + min_overlap and stop_coord1 <= stop - min_overlap:
                                    cov1 = cov_dict[scaff_name_map][mappable_region] 
                                    map_region1 = mappable_region
                                    not_found1 = False
                                if start_coord2 >= start + min_overlap and start_coord2 <= stop - min_overlap:
                                    cov2 = cov_dict[scaff_name_map][mappable_region] 
                                    map_region2 = mappable_region
                                    not_found2 = False
                                elif stop_coord2 >= start + min_overlap and stop_coord2 <= stop - min_overlap:
                                    cov2 = cov_dict[scaff_name_map][mappable_region] 
                                    map_region2 = mappable_region
                                    not_found2 = False
                            try:
                                #If the reads are from different mappable regions, check the maximum coverage of either region
                                if map_region1 != map_region2 and not not_found1 and not not_found2:
                                    cov = max(cov1, cov2)
                                    #If maximum coverage is less or equal to two increase count for number of mapped reads of the scaffold where the reads map
                                    if cov <= 2:
                                        total_reads += 1
                                        if scaff_name_map in reads_per_scaffold.keys():
                                            reads_per_scaffold[scaff_name_map] += 1
                                        else:
                                            reads_per_scaffold[scaff_name_map] = 1
                                        #Write bam enries to output file
                                        out_f.write(line)
                                        out_f.write(nextline)
                                #If reads are from the same region, check if maximum coverage of the region is no more than two, if so increase count for number of mapped reads of the scaffold where the reads map
                                elif map_region1 == map_region2 and not not_found1 and not not_found2:
                                    cov = cov1
                                    if cov <= 2:
                                        total_reads += 1
                                        if scaff_name_map in reads_per_scaffold.keys():
                                            reads_per_scaffold[scaff_name_map] += 1
                                        else:
                                            reads_per_scaffold[scaff_name_map] = 1
                                        #Write bam enries to output file
                                        out_f.write(line)
                                        out_f.write(nextline)
                            except: 
                                pass
                    else:
                        #Wtite header to output file
                        out_f.write(line)

    #Document total number of reads mapped to species
    with open(summary_out_file, "a") as out_f:
        out_f.write("{},{},{},{}\n".format(tax_id,minid,amb,str(total_reads)))

    #Go through all scaffolds and get their, names, lengths, and number of mapped reads
    scaffold_names = []
    scaffold_lengths = []
    scaffold_reads = []
    with open(scaffold_lengths_file, "r") as cl_file:
        for line in cl_file:
            split_line = line.strip("\n").split(",")
            scaffold = split_line[0]
            length = int(split_line[1])
            try:
                reads = reads_per_scaffold[scaffold]
            except:
                reads = 0
            scaffold_lengths.append(length)
            scaffold_reads.append(reads)
            scaffold_names.append(scaffold)
    #Save files if needed for future plotting
    np_scaff_length = np.array(scaffold_lengths)
    np_nr_mapped = np.array(scaffold_reads)
    np.savetxt("numpy_files/reads_per_scaffold_{}_minid{}_amb_{}.txt".format(tax_id, minid, amb), np_nr_mapped)
    np.savetxt("numpy_files/scaffold_length_{}.txt".format(tax_id), np_scaff_length)
    #Calculate pearson
    res = stats.pearsonr(np_scaff_length, np_nr_mapped)
    pearson = res.statistic
    #Used to place label when ploting
    min_reads = np_nr_mapped.min()
    max_scaffold = np_scaff_length.max()
    print(max_scaffold)
    #Plot number of mapped reads per scaffold (y-axsis) vs scaffold length (x-axsis)
    plt.scatter(np_scaff_length, np_nr_mapped, color="black")
    plt.text(max_scaffold*0.78, min_reads, 'PCC: %.2f' % pearson, fontsize=16)
    plt.xlabel("Scaffold length")
    plt.ylabel("Number of mapped reads")
    plt.title("{}".format(name))
    plt.savefig("plots/reads_vs_scaffold_length/reads_vs_scaffold_length_all_weeks_mapped_to_{}_minid0{}_ambigious_{}_v2.png".format(tax_id, minid, amb))
    plt.clf()