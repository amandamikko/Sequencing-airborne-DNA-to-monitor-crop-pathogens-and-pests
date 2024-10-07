# Sequencing-air-to-monitor-plant-pathogens

Code used to obtain results from the paper *Sequencing air to monitor plant pathogens*.

## Generating, classifying and mapping reads
Snakemake-file can be used to:
1. Generate reads from all genomes in a ganera and map the gnerated reads to a species in that same ganera. In the paper, this was done to evaliuate if DNA from the lowly abundant species were present in the air filters. This was done by comparing correlations between scaffold length and number of mapped reads for mapped air filter reads and correlations produced by generated reads. 

2. Generate and classify reads using Kraken2, to evaluate classification potential. How the generated reads are classified can then be visualized using Sankey and Krona plots.

3. Find classefiable regions of species and map air filter reads to that species. Number of mapped reads per classefiable reigion can after running the snakemake-file be visualized using plot_reads_per_classefiable_region.py. 

To run the snakemake-file, you must:
Download code to generate reads using NEAT (https://github.com/zstephens/neat-genreads). Path to gen_reads.py within NEAT code should be provided as the parameter path_to_NEAT in the Snakemake file.

Create references for the genomes that you wish to map reads to. These should be created in the ref_data folder in a subfolder named after the taxid of the genome. For example ref_data/27350

Download picard and provide a path to it specified by the variable path_to_picard in the Snakefile.

To dun the Snakemake pipelinge you can for example type: "snakemake --cores 15 --use-conda"
After running snakemake-file, you can:

## Visualize observed degree of damage and air filter signal strength (bar+scatter plot)
1. plot_observed_damage_vs_signal_strength/create_df_from_kraken_report.py Creates dataframe with number of reads classifed to species in each kraken reports.
2. plot_observed_damage_vs_signal_strength/perform_gbm.R Does zero raplacement so that pivot coordinates can be calculated
3. plot_observed_damage_vs_signal_strength/get_pivot_coor_from_gbm.py calculate pivot coordinates
4. plot_observed_damage_vs_signal_strength/plot_observed_degree_of_damage_with_pivot_coord.py Calculate average observed damage at different distance and plot the values together with air filter signal strength.

## Visualixe of number of mapped redads per classefiable region
This can be done by running plot_reads_per_classefiable_region/plot_reads_per_classefiable_region.py 

## Visualize reads per scaffold of air filter reads and generated reads
This can be done by running:
1. plot_reads_per_scaffold/analyze_stacks_from_mapping.py to plot reads per scaffold of the mapped air filter reads. This code makes sure that all counted reads do not overlap with more than one other read. This is done to disregard reads from repetetive and conserved regions.
2. plot_reads_per_scaffold/calculate_nr_reads_to_subsample.py to Calculate number of reads to subsample from the mapped generated reads. Requires a csv file describing number of reads from the air filters which map to the different scaffolds (genrated from analyze_stacks_from_mapping.py).
3. plot_reads_per_scaffold/get_reads_per_scaffold_mapped_generated_reads.py to plot number of mapped reads per scaffold vs scaffold length (of generated reads).
