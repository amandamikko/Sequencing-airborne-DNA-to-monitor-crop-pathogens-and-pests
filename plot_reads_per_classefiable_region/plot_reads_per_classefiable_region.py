import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import sys
import os

species_name_dict = {"40559" : "Botrytis cinerea", "984962" : "Heterobasidion irregulare", "53485" : "Pyrenophora teres",
                           "161013" : "Thrips palmi", "34373" : "Blumeria graminis", "13684" : "Parastagonospora nodorum",
                           "5599" : "Alternaria alternata", "1047171" : "Zymoseptoria tritici", "45130" : "Bipolaris sorokiniana",
                           "138532" : "Curtobacterium flaccumfaciens pv. flaccumfaciens", "5518" : "Fusarium graminearum", "36044" : "Erysiphe pisi",
                           "5113" : "Epichloe typhina", "112498" : "Ramularia collo-cygni", "120017" : "Ustilago hordei",
                           "45151" : "Pyrenophora tritici-repentis", "263140" : "Sitodiplosis mosellana", "27350" : "Puccinia striiformis",
                           "208348" : "Puccinia triticina", "5297" : "Puccinia graminis", "220672" : "Plenodomus biglobosus",
                           "13164" : "Myzus persicae", "5180" : "Sclerotinia sclerotiorum", "48100" : "Alternaria solani",
                           "27337" : "Verticillium dahliae", "294": "Pseudomonas fluorescens",
                           "282267" : "Fusarium asiaticum",
                           "1634478" : "Aquanectria penicillioides", "64609" : "Ilyonectria destructans",
                           "182845" : "Calonectria ilicicola", "1079257" : "Ilyonectria robusta", "78403" : "Neonectria neomacrospora",
                           "1028729" : "Fusarium pseudograminearum", "1715230": "Fusarium solani", "660027": "Fusarium oxysporum",
                           "28447" : "Clavibacter michiganensis", "294":"Pseudomonas fluorescens", "7029":"Acyrthosiphon pisum ",
                           "1431903": "Brassicogethes aeneus", "4787": "Phytophthora infestans", "27344": "Puccinia coronata",
                           "27345": "Puccinia hordei", "40932": "Rhopalosiphum padi", "44670": "Metopolophium dirhodum",
                           "44664": "Sitobion avenae",}

category_coulour_dict = {"Arthropoda" : "#a42255", "Bacteria" : "#88ccee", "Fungi": "#d8af39",
                         "Nematodes":  "#e15b64", "Oomycetes": "#117733", "Viruses": "#44aa99"}

tax_id_to_type = {"53485":"Fungi", "34373":"Fungi", "1047171":"Fungi", "45151":"Fungi", "7029":"Arthropoda", "27350":"Fungi"}
minid = "97"
amb = "toss"

min_overlap = 1
taxid_GCA_list = [["45151", "GCA_000149985_1_ASM14998v1_genomic", "_insert_237_insstd_98_conf01_min_hit_10"]]#, ["27350", "GCA_021901695_1_Pst134E36_v1_pri_genomic", "_insert_259_insstd_98_conf01_min_hit_10"], ["53485", "GCA_014334815_1_ASM1433481v1_genomic", "_insert_224_insstd_97_conf01_min_hit_10"], ["34373", "GCA_905067625_1_Bgtriticale_THUN12_genome_v1_2_genomic", "_insert_263_insstd_102_conf01_min_hit_10"], ["1047171", "GCA_000219625_1_MYCGR_v2_0_genomic", "_insert_280_insstd_95_conf01_min_hit_10"], ["7029", "GCA_005508785_2_pea_aphid_22Mar2018_4r6ur_v2_genomic", "_insert_247_insstd_86_conf01_min_hit_10"]]

for pair in taxid_GCA_list:
    taxid = pair[0]
    GCA = pair[1]
    error = pair[2]
    name = species_name_dict[taxid]
    mapped_file = "map_files/bam/no_PCR_merged_conf01_min_hit_10_subseq_{t}_mapped_to_{t}_amb_toss_minid_0_97.sorted.sam".format(t=taxid)
    scaff_file = "coverage_coord/coverage_coordinates_{g}_generated_perfect_NEAT_10x{e}_subseq_{t}_mapped_to_{t}_amb_toss_minid_0_97.txt".format(g=GCA,t=taxid, e=error)
    scaff_dict = {}
    # uses the covcoord file created from get_covcoord.py (scaff_coord)
    with open(scaff_file, "r") as scaff_coord:
        for line in scaff_coord:
            scaff = line.split(",")[0]
            start = int(line.split(",")[1])
            stop = int(line.split(",")[2])
            if scaff not in scaff_dict.keys():
                scaff_dict[scaff] = {(start, stop): 0}
            else:
                scaff_dict[scaff][(start, stop)] = 0


    with open(mapped_file, "r") as mapped_stats:
        for mapped_line in mapped_stats:
            if mapped_line[0] != "@":
                nextline = next(mapped_stats)
                #Check that read is only mapping to 1 contig
                scaff_name_map = mapped_line.strip("\n").split()[3]
                scaff_name_map2 = nextline.strip("\n").split()[3]
                if scaff_name_map == scaff_name_map2:
                    if scaff_name_map in scaff_dict.keys():
                        start_coord = int(mapped_line.strip("\n").split()[4])
                        read_length = mapped_line.strip("\n").split()[9]
                        stop_coord = start_coord + int(read_length)
                        #Find which coordinate is smallest and assing that to be the start coordinate
                        if stop_coord < start_coord:
                            temp_start = start_coord
                            temp_stop = stop_coord
                            start_coord = temp_stop
                            stop_coord = temp_start
                        #Go through all mappable regions and see if start or stop coordinate falls within the region, if so, add one mapped read to that region
                        for mappable_region in scaff_dict[scaff_name_map].keys():
                            start = mappable_region[0]
                            stop = mappable_region[1]
                            if start_coord >= start + min_overlap and start_coord <= stop - min_overlap:
                                scaff_dict[scaff_name_map][mappable_region] += 1
                            elif stop_coord >= start + min_overlap and stop_coord <= stop - min_overlap:
                                scaff_dict[scaff_name_map][mappable_region] += 1
                #If it maps to multiple contigs, go though classefiable regions on both contigs 
                else:
                    if scaff_name_map in scaff_dict.keys():
                        start_coord = int(mapped_line.strip("\n").split()[4])
                        read_length = mapped_line.strip("\n").split()[9]
                        stop_coord = start_coord + int(read_length)
                        #Find which coordinate is smallest and assing that to be the start coordinate
                        if stop_coord < start_coord:
                            temp_start = start_coord
                            temp_stop = stop_coord
                            start_coord = temp_stop
                            stop_coord = temp_start
                        #Go through all mappable regions and see if start or stop coordinate falls within the region, if so, add one mapped read to that region
                        for mappable_region in scaff_dict[scaff_name_map].keys():
                            start = mappable_region[0]
                            stop = mappable_region[1]
                            if start_coord >= start + min_overlap and start_coord <= stop - min_overlap:
                                scaff_dict[scaff_name_map][mappable_region] += 1
                            elif stop_coord >= start + min_overlap and stop_coord <= stop - min_overlap:
                                scaff_dict[scaff_name_map][mappable_region] += 1
                    if scaff_name_map2 in scaff_dict.keys(): #Changed from elif to avaluate if we get same result
                        start_coord = int(nextline.strip("\n").split()[4])
                        read_length = nextline.strip("\n").split()[9]
                        stop_coord = start_coord + int(read_length)
                        #Find which coordinate is smallest and assing that to be the start coordinate
                        if stop_coord < start_coord:
                            temp_start = start_coord
                            temp_stop = stop_coord
                            start_coord = temp_stop
                            stop_coord = temp_start
                        #Go through all mappable regions and see if start or stop coordinate falls within the region, if so, add one mapped read to that region
                        for mappable_region in scaff_dict[scaff_name_map2].keys():
                            start = mappable_region[0]
                            stop = mappable_region[1]
                            if start_coord >= start + min_overlap and start_coord <= stop - min_overlap:
                                scaff_dict[scaff_name_map2][mappable_region] += 1
                            elif stop_coord >= start + min_overlap and stop_coord <= stop - min_overlap:
                                scaff_dict[scaff_name_map2][mappable_region] += 1
   #Add all classefiable region lengths and number of mapped reads on those regions to lists for plotting
    nr_mapped = []
    scaff_length = []
    for scaff_name_map in scaff_dict.keys():
        for region in scaff_dict[scaff_name_map].keys():
            count = scaff_dict[scaff_name_map][region]
            nr_mapped.append(count)
            start = region[0]
            stop = region[1]
            length = stop - start + 1
            scaff_length.append(length)
    #Convert to numpy arrays
    np_scaff_length = np.array(scaff_length)
    np_nr_mapped = np.array(nr_mapped)

    #Get latin name and concatinate it so that it can be used to save files 
    name = species_name_dict[taxid]
    split_name = name.split()
    lat_name = ""
    for part in split_name:
        lat_name += part + "_"
    lat_name = lat_name[:-1]

    np.savetxt("numpy_files/counts_{}_classefiable_regions.txt".format(lat_name), np_nr_mapped)
    np.savetxt("numpy_files/length_{}_classefiable_regions.txt".format(lat_name), np_scaff_length)

    #Plot number of mapped reads per classefiable region (y-axis) against classefiable region length (x-axis)
    colour_to_plot = category_coulour_dict[tax_id_to_type[taxid]]
    res = stats.pearsonr(np_scaff_length,np_nr_mapped)
    pearson = res.statistic
   
    min_reads = np_nr_mapped.min()
    max_contig = np_scaff_length.max()
    plt.text(max_contig*0.78, min_reads, 'PCC: %.2f' % pearson, fontsize=16)
    plt.scatter(np_scaff_length, np_nr_mapped, color="black")
    plt.xlabel("Scaffold length")
    plt.ylabel("Number of mapped reads")
    plt.title("{}".format(name))
    plt.savefig("plots/reads_vs_classefiable_region/reads_vs_contig_length_all_weeks_mapped_to_{}_minid0{}_ambigious_{}.png".format(taxid, minid, amb))
    plt.clf()

