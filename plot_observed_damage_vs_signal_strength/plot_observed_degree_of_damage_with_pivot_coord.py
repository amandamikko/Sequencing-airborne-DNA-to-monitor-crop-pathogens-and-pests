import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Extracts pivot coordinates for observed species specified in the file observed_species_info/observed_species_name_and_taxid.tsv and plot these values with average degree of observed damage

obs_date_to_week = {"2007-04-16": "2007-16",
                    "2007-04-23": "2007-17",
                    "2007-04-30": "2007-18",
                    "2007-05-07": "2007-19",
                    "2007-05-14": "2007-20",
                    "2007-05-21": "2007-21",
                    "2007-05-28": "2007-22",
                    "2007-06-04": "2007-23",
                    "2007-06-11": "2007-24",
                    "2007-06-18": "2007-25",
                    "2007-06-25": "2007-26",
                    "2007-07-02": "2007-27",
                    "2007-07-09": "2007-28",
                    "2007-07-16": "2007-29",
                    "2007-07-23": "2007-30",
                    "2007-07-30": "2007-31",
                    "2007-08-06": "2007-32",
                    "2007-08-13": "2007-33",
                    "2007-08-20": "2007-34",
                    "2007-08-27": "2007-35"}
#Create list of all week number (used when creating plots)
all_obs_dates = obs_date_to_week.keys()
week_dates = []
for d in all_obs_dates:
    week_dates.append(obs_date_to_week[d])

#Path to information about observed species
obs_pat_info = "../observed_species_info/observed_species_name_and_taxid.tsv"
#Path to output_file
output_file = "../kraken_reports/pivot_coord_obs_species.csv"

#Load file with calculated pivotcoordinates
pivot_coord = pd.read_csv('../kraken_reports/all_weeks_count_values_removed_zero_inflated_gbm_pivot_coord.csv', sep=',')
#Get taxids for which pivot coordinates have been calculated
tax_ids = list(pivot_coord.columns)
pivot_coord = pivot_coord.transpose()
#Assign correct labels
week_dates_air_filer = ["2007-15", "2007-19", "2007-23", "2007-27", "2007-29", "2007-35", "2007-37", "2007-41", "2007-45"]
pivot_coord.columns = week_dates_air_filer
pivot_coord = pivot_coord.transpose()

#Create dictionaries which describe connection between tax id and common name of species as well as common name and latin name
tax_id_to_common = {}
common_to_lat_names = {}
with open(obs_pat_info, "r") as obs_pat:
    for line in obs_pat:
        split_line = line.strip("\n").split("\t")
        if len(split_line[1].split(",")) == 1:
            taxid = split_line[1].strip()
            #Only stor info for species for which pivot coordinates could be calculated
            if taxid in tax_ids:
                common_name = split_line[2]
                tax_id_to_common[taxid] = common_name
                lat_name = split_line[0]
                common_to_lat_names[common_name] = lat_name

#Extract pivot coordinates for all observed species for which pivot coordinates could be calculated
columns_to_extract = list(tax_id_to_common.keys())
relevant_pivot = pivot_coord[columns_to_extract].copy()
#Rename columns to species names
common_list = []
for tax in columns_to_extract:
    common_list.append(tax_id_to_common[tax])
relevant_pivot.columns = common_list

#Save file if needed only for plotting later
out_pd = relevant_pivot.transpose()
out_pd.to_csv(output_file,sep=";")

#Load pivot coordinates
df_filter_summary = pd.read_csv("../kraken_reports/pivot_coord_obs_species.csv", sep=";", index_col=0)

kraken_indices = list(df_filter_summary.index.values)

#Create list of distance intervals in which you want to calculate average obseved damage
length_steps = []
length_steps.append((0,10))
for step in range(10, 100, 10):
    length_steps.append((step, step+10))

#Dictionary to stor all aberage observation values for all distance intervals
coord_to_obs_dict = {}

with open("../observed_species_info/Graderingslista_2007_Skåne.csv", "r") as csv:
    headers = csv.readline().strip('\n').split("\t")
    #Load observation values
    pat_header = headers[4:]
    nr_pat = len(pat_header)
    nr_to_obs_pat = {}
    #Create index number for each pathigen
    nr = 0
    for pat in pat_header:
        nr_to_obs_pat[pat] = nr
        nr += 1
    #Go through observation values
    for line in csv:
        #Makesure line is not blank
        if len(line.strip('\n').split("\t")) >  3 \
                and line.strip('\n').split("\t")[2] != "":
            #Replace comma so that values can be converted to float
            obs_values = line.strip('\n').replace(",", ".").split("\t")[4:]
            #Load date, length from filter station
            date = line.strip('\n').split("\t")[3]
            curr_interval = (int(line.strip('\n').split("\t")[2].split("-")[0]), int(line.strip('\n').split("\t")[2].split("-")[1]))
            obs_interval = coord_to_obs_dict.keys()
            #See if damage within interval has been reported for any pathogen
            if curr_interval in obs_interval:
                #Load current count values (Number of pathogen damage was scored) and current average observation value
                coord_obs_df = coord_to_obs_dict[curr_interval][0]
                count_obs_df = coord_to_obs_dict[curr_interval][1]
                #Go through all pathogens and update count and average observation values
                for pat in nr_to_obs_pat.keys():
                    nr = nr_to_obs_pat[pat]
                    curr_val_to_insert = obs_values[nr]
                    if curr_val_to_insert != "":
                        float_curr_val_to_insert = float(obs_values[nr])
                        curr_count = count_obs_df[pat][obs_date_to_week[date]]
                        #If pathogen damage was not yet reported insert pathogen into count_obs_df and coord_obs_df
                        if curr_count == 0:
                            count_obs_df[pat][obs_date_to_week[date]] += 1
                            coord_obs_df[pat][obs_date_to_week[date]] = float_curr_val_to_insert
                        #Else update current values
                        else:
                            curr_sum = curr_count*count_obs_df[pat][obs_date_to_week[date]]
                            coord_obs_df[pat][obs_date_to_week[date]] = (curr_count*coord_obs_df[pat][obs_date_to_week[date]]
                                                                          + float_curr_val_to_insert)\
                                                                         / \
                                                                         (count_obs_df[pat][obs_date_to_week[date]] + 1)
                            count_obs_df[pat][obs_date_to_week[date]] += 1

                #Update dict which stores values for all pathogens in the current interval
                coord_to_obs_dict[curr_interval] = [coord_obs_df, count_obs_df]
            #If damage from any pathogen was not yet reported in the current interval insert interval into coord_to_obs_dict
            else:
                #Initiate dataframes used to store number of times damage from pathogen was reported and dataframe used to store average observation values
                count_obs_df = pd.DataFrame(0, index=week_dates, columns=pat_header)
                coord_obs_df = pd.DataFrame(index=week_dates, columns=pat_header)
                #Update average values and counts
                for pat in nr_to_obs_pat.keys():
                    nr = nr_to_obs_pat[pat]
                    curr_val_to_insert = obs_values[nr]
                    if curr_val_to_insert != "":
                        float_curr_val_to_insert = float(obs_values[nr])
                        coord_obs_df[pat][obs_date_to_week[date]] = float_curr_val_to_insert
                        count_obs_df[pat][obs_date_to_week[date]] += 1
                #Save count values and average obervation values in dictionary which stores all dfs for all intervals
                coord_to_obs_dict[curr_interval] = [coord_obs_df, count_obs_df]

all_coords = coord_to_obs_dict.keys()
for pat in kraken_indices:
    nr_reads_pat = df_filter_summary.loc[pat]
    filter = nr_reads_pat.values

    #min-max normalize
    filter_min_val = np.nanmin(filter)
    filter_max_val = np.nanmax(filter)
    filter_normalized_arr = (filter - filter_min_val) / (filter_max_val - filter_min_val)
    filter_normalized_arr = np.nan_to_num(filter_normalized_arr)
    kraken_dates = nr_reads_pat.index.tolist()

    #Plot filter signal
    week_list = []
    for date in kraken_dates:
        week_list.append(int(date[-2:]))
    plt.bar(week_list, filter_normalized_arr, width=0.2, color="black", label="Airfilter reads") #, zorder=0)
    plt.ylabel("Normalized observation values")
    plt.xlabel("Week number")

    #Create string which display the species common name
    concat_reg_name = ""
    for p in pat.split(" "):
        concat_reg_name += p + "_"
    concat_reg_name = concat_reg_name[:-1]
    lat_name = common_to_lat_names[pat]


    #Get values for pathogen from all intervals
    all_coord_lists = []
    all_coord_values = []
    for coord in all_coords:
        coord_info = coord_to_obs_dict[coord]
        obs_values = coord_info[0]
        obs_count = coord_info[1]
        pat_values = obs_values[pat].tolist()
        for val in pat_values:
            all_coord_values.append(val)
        obs_dates = coord_info[1].index.values
        week_list = []
        for date in obs_dates:
            week_list.append(int(date[-2:]))
        label = "{}-{}km".format(coord[0], coord[1])
        all_coord_lists.append([pat_values, week_list, label])

    #Min max normalize observed degre of damage
    obs_min_val = np.nanmin(all_coord_values)
    obs_max_val = np.nanmax(all_coord_values)
    all_normalized_lists = []
    for list in all_coord_lists:
        obs_normalized_arr = (list[0] - obs_min_val) / (obs_max_val - obs_min_val)
        all_normalized_lists.append([obs_normalized_arr, list[1], list[2]])

    #Create string which display species clatine name
    concat_lat_name = ""
    for p in lat_name.split(" "):
        concat_lat_name += p + "_"
    concat_lat_name = concat_lat_name[:-1]
    plt.title("{}".format(lat_name))
    
    #Plot observed degree of damage and add lable
    colour_list = ["#9C4511", "#BD651A", "#DD8629", "#F5AD52", "#FED693", "#BBE4D1", "#76C7BE", "#3EA8A6", "#208288"]
    colour_num=0
    for list in all_normalized_lists:
        if (np.isnan(list[0])).sum() > 0:
            plt.scatter(list[1], list[0], label=list[2], c=colour_list[colour_num])
        colour_num += 1
    plt.legend(bbox_to_anchor=(1, 0.5), loc="center left" , fontsize=8)
    plt.tight_layout()
    plt.savefig("../plots/scatter_and_bar/scatter_{}_{}_obs_damage_with_pivot_coord_test_code.png".format(concat_lat_name, concat_reg_name.replace("/","per"))) #HERE
    plt.clf()