import pandas as pd

dict_sample_week = {"Sample_Air-20161114-3": "2007-19",
        "Sample_Air-20161115-21": "2007-29",
        "Sample_Air-20161115-23": "2007-35",
        "Sample_Air-20161115-24": "2007-27",
        "Sample_Air-20161121-40": "2007-37",
        "Sample_Air-20161122-56": "2007-15",
        "Sample_Air-20161122-64": "2007-45",
        "Sample_Air-20161124-84": "2007-41",
        "Sample_Air-20161124-86": "2007-23"}

path_to_kraken_report = "../kraken_reports/"
path_file_w_report_names = "../kraken_reports/all_report_names.txt"
output_df = path_to_kraken_report + "all_weeks_count_values_removed_zero_inflated.csv"
text_after_sample_name_in_file_name = "_conf01_min_hit_10.report"
write_read_count_to_file = True
write_count_values_to_file_zero = True

read_counts = pd.DataFrame()


#Create list of repport names
repport_names = []
with open(path_file_w_report_names, "r") as rep_names:
    for line in rep_names:
        repport_names.append(line.strip("\n"))

#Go through all reports and add count values to pandas dataframe that uses taxid as index
for name in repport_names:
    print(name)
    with open(path_to_kraken_report + name, "r") as repport:
        report_name = name.split(text_after_sample_name_in_file_name)[0]
        week_name = dict_sample_week[report_name]
        read_counts[week_name] = 0
        total_count = 0
        read_counts.loc["0", week_name] = int(0)
        for line in repport:
            split_list = line.split("\t")
            rank = split_list[3]
            species = split_list[-1].strip("\n").strip()
            tax_id = split_list[-2]
            nr_classified_clade = int(split_list[1])
            nr_classified_at = int(split_list[2])
            if rank == "S":
                #Insert new row if species not seen before
                if tax_id not in list(read_counts.index):
                    new_row = (len(list(read_counts.columns)) - 1) * [0] + [nr_classified_clade]
                    read_counts.loc[tax_id] = new_row
                    read_counts = read_counts.sort_index()
                #If taxid has been observed in previous repports update value for current report
                else:
                    read_counts.at[tax_id, week_name] = nr_classified_clade
            #Add reads not classified at species level to the unclassified taxon
            elif "S" not in rank:
                total_count += int(nr_classified_at)

        read_counts.at["0", week_name] = total_count


#Sort weeks so that they are in right order in df
read_counts = read_counts[sorted(read_counts.columns)]
if write_read_count_to_file:
    read_counts.to_csv(path_to_kraken_report + "all_weeks_count_values.csv")

#Find and remove zeroinflated taxa
max_allowed_zero = 2/3
T_F_remove_vector = read_counts.eq(0).sum(axis=1) >= len(list(read_counts.columns))*max_allowed_zero
remove_rows = read_counts[T_F_remove_vector]


#Add reads from removed taxa to unclassified
for name in repport_names:
    report_name = name.split(text_after_sample_name_in_file_name)[0]
    week_name = dict_sample_week[report_name]
    read_counts.at["0", week_name] += remove_rows[week_name].sum()

#Drop zero inflated taxa
read_counts.drop(index=list(remove_rows.index), inplace=True)

if write_count_values_to_file_zero:
    read_counts.to_csv(output_df)
