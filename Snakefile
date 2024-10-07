GENOME_Pyrenophora_teres = ["GCA_014334815_1_ASM1433481v1_genomic"]
GENOME_Blumeria_graminis = ["GCA_905067625_1_Bgtriticale_THUN12_genome_v1_2_genomic"]
GENOME_Zymoseptoria_tritici = ["GCA_000219625_1_MYCGR_v2_0_genomic"]
GENOME_Pyrenophora_tritici_repentis = ["GCA_000149985_1_ASM14998v1_genomic"]
GENOME_Acyrthosiphon_pisum = ["GCA_005508785_2_pea_aphid_22Mar2018_4r6ur_v2_genomic"]
GENOME_Puccinia_striiformis = ["GCA_021901695_1_Pst134E36_v1_pri_genomic"]
Phytophthora_GENOMES = ["GCA_000247585_2_PP_INRA-310_V2_genomic", "GCA_033557915_1_HumboldtARI_Pten_1_0_genomic", "GCA_016880985_1_ASM1688098v1_genomic", "GCA_000142945_1_ASM14294v1_genomic"]
Sitobion_GENOMES = ["GCA_008086715_1_ASM808671v1_genomic", "GCA_019425605_1_ASM1942560v1_genomic"]
Rhopalosiphum_GENOMES = ["GCA_003676215_3_ASM367621v3_genomic", "GCA_036289425_1_JGU_Rn_01_genomic", "GCA_020882245_1_ASM2088224v1_genomic"]
Metopolophium_GENOMES = ["GCA_019925205_1_ASM1992520v1_genomic"]
Brassicogethes_GENOMES = ["GCA_921294245_1_PGI_MELIAE_v1_genomic"]
pyrenophora_GENOMES = ["GCA_000149985_1_ASM14998v1_genomic", "GCA_014334815_1_ASM1433481v1_genomic", "GCA_000465215_2_Pysem1_0_genomic", "GCA_012365135_1_ASM1236513v1_genomic"]
acyrthosiphon_genomes = ["GCA_005508785_2_pea_aphid_22Mar2018_4r6ur_v2_genomic"]
zymoseptoria_genomes = ["GCA_000983655_1_ASM98365v1_genomic", "GCA_000223825_2_ASM22382v2_genomic", "GCA_000223765_2_ASM22376v2_genomic", "GCA_000223685_2_ASM22368v2_genomic", "GCA_026122155_1_ASM2612215v1_genomic", "GCA_000219625_1_MYCGR_v2_0_genomic"]
blumeria_genomes = ["GCA_905067625_1_Bgtriticale_THUN12_genome_v1_2_genomic", "GCA_900239735_1_BGH_DH14_v4_genomic"]
puccinia_genomes = ["GCA_000149925_1_ASM14992v1_genomic", "GCA_001263375_1_ASM126337v1_genomic", "GCA_001624995_1_ASM162499v1_genomic", "GCA_002873125_1_ASM287312v1_genomic", "GCA_004348175_1_PNOVO_Noble_KM_genomic", "GCA_007896445_1_ASM789644v1_genomic", "GCA_019395275_1_ASM1939527v1_genomic", "GCA_021901695_1_Pst134E36_v1_pri_genomic", "GCA_026914185_1_ASM2691418v1_genomic"]
GENOMES_obs_species = ["GCA_014334815_1_ASM1433481v1_genomic", "GCA_905067625_1_Bgtriticale_THUN12_genome_v1_2_genomic","GCA_000219625_1_MYCGR_v2_0_genomic","GCA_000149985_1_ASM14998v1_genomic", "GCA_005508785_2_pea_aphid_22Mar2018_4r6ur_v2_genomic", "GCA_021901695_1_Pst134E36_v1_pri_genomic", "GCA_921294245_1_PGI_MELIAE_v1_genomic", "GCA_000142945_1_ASM14294v1_genomic", "GCA_002873125_1_ASM287312v1_genomic", "GCA_007896445_1_ASM789644v1_genomic", "GCA_020882245_1_ASM2088224v1_genomic", "GCA_019925205_1_ASM1992520v1_genomic", "GCA_019425605_1_ASM1942560v1_genomic"]
all_genrated_reads = puccinia_genomes + Brassicogethes_GENOMES + Metopolophium_GENOMES + Rhopalosiphum_GENOMES + Sitobion_GENOMES + Phytophthora_GENOMES

confidence_Puccinia_striiformis = ["_insert_259_insstd_98"]
confidence_Blumeria_graminis = ["_insert_263_insstd_102"]
confidence_Pyrenophora_tritici_repentis = ["_insert_237_insstd_98"]
confidence_Pyrenophora_teres = ["_insert_224_insstd_97"]
confidence_Zymoseptoria_tritici = ["_insert_280_insstd_95"]
confidence_Acyrthosiphon_pisum = ["_insert_247_insstd_86"]
confidence_Sclerotinia_sclerotiorum = ["_insert_232_insstd_97"]

obs_species_with_more_than_1000_reads_TAXID = ["45151", "1047171", "34373", "7029", "27350", "53485"]
obs_species_with_less_than_1000_reads_TAXID = ["44670", "44664","40932", "4787", "1431903", "27344", "27345"]
ambiguous = ["toss"]
COV=["10"]
minid = ["97"]
confidence = "_conf01_min_hit_10"
path_to_genomes="/proj/nobackup/snic2019-35-73/private/amanda/generate_classify/genomes/last_genomes/"
path_to_fq = "/proj/nobackup/snic2019-35-73/private/amanda/filtered_human/"
#Pattern that can list all samples in the directory specified by path_to_fq
sample_pattern_R1 = "*_unmapped_R1.fq.gz"
sample_pattern_R2 = "*_unmapped_R2.fq.gz"
#Common start of all sample files
sample_start = "Sample_Air-"
#Get all sequence names
samples = glob_wildcards(path_to_fq+"{sname}_unmapped_R1.fq.gz").sname
#Files which are needed before running the snakemake file
path_to_selmeout = "selmeout.py"
path_to_db = "/proj/nobackup/snic2019-35-73/private/amanda/renamed_header_nt_201222_pathogens"
path_to_names = path_to_db+"/taxonomy/names.dmp"
path_to_nodes = path_to_db+"/taxonomy/nodes.dmp"
path_to_NEAT="/proj/nobackup/snic2019-35-73/private/amanda/git/NEAT/gen_reads.py"
path_to_picard="picard.jar"

#TO DO, write rules (or appropriate targets) subsequencing and mapping subsequenced airfilter reads to reference genomes of the species they were classified to
#Impelemnt bam to sam conversion
ruleorder: bbmap_alig_subseq > sort_bam_generated_NEAT > run_picard_generated > bam_to_sam
rule target:
  input:
    #Classify air filter reads
    expand("kraken_outputs/{file_name}_default_conf_min_hit.output.gz", file_name=samples),
    expand("kraken_outputs/{file_name}_conf01_min_hit_10.output.gz", file_name=samples),
    #Rule to first generate reads using NEAT and then map the generated reads to evaluate if DNA from species could be present in air filters
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_196_insstd_96_mapped_to_4787_amb_{amb}_minid_0_{minid}.bam", genome_name = Phytophthora_GENOMES, amb = ambiguous, minid = minid),
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_215_insstd_99_mapped_to_44664_amb_{amb}_minid_0_{minid}.bam", genome_name = Sitobion_GENOMES, amb = ambiguous, minid = minid),
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_180_insstd_86_mapped_to_40932_amb_{amb}_minid_0_{minid}.bam", genome_name = Rhopalosiphum_GENOMES, amb = ambiguous, minid = minid),
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_176_insstd_88_mapped_to_44670_amb_{amb}_minid_0_{minid}.bam", genome_name = Metopolophium_GENOMES, amb = ambiguous, minid = minid),
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_170_insstd_87_mapped_to_1431903_amb_{amb}_minid_0_{minid}.bam", genome_name = Brassicogethes_GENOMES, amb = ambiguous, minid = minid),
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_245_insstd_115_mapped_to_27344_amb_{amb}_minid_0_{minid}.bam", genome_name = puccinia_genomes, amb = ambiguous, minid = minid),
    expand("map_files/bam/{genome_name}_generated_perfect_NEAT_10x_insert_203_insstd_95_mapped_to_27345_amb_{amb}_minid_0_{minid}.bam", genome_name = puccinia_genomes, amb = ambiguous, minid = minid),
    #Map all sequenced air filter reads to species and get coverage coordinates
    expand("map_files/bam/merged_mapped_to_{seqid}_amb_{amb}_minid_0_{minid}.sorted.sam", seqid = obs_species_with_less_than_1000_reads_TAXID, amb = ambiguous, minid = minid),
    #Get scafold lengths (used when plotting number of reads vs scafold length)
    expand("scaffold_summaries/{genome_name}.txt", genome_name=all_genrated_reads),
    #Generate and classify reads for evaluation of classification potenial
    expand("kraken_reports/{genome_name}_generated_perfect_NEAT_10x_insert_300_insstd_30_conf01_min_hit_10.report", genome_name = GENOMES_obs_species),
    #Get clasefiable regions
    expand("coverage_coord/coverage_coordinates_GCA_005508785_2_pea_aphid_22Mar2018_4r6ur_v2_genomic_generated_perfect_NEAT_{cov}x{conf}_conf01_min_hit_10_subseq_7029_mapped_to_7029_amb_{amb}_minid_0_{min}.txt", amb=ambiguous, conf=confidence_Acyrthosiphon_pisum, min=minid, cov=COV),
    expand("coverage_coord/coverage_coordinates_GCA_014334815_1_ASM1433481v1_genomic_generated_perfect_NEAT_{cov}x{conf}_conf01_min_hit_10_subseq_53485_mapped_to_53485_amb_{amb}_minid_0_{min}.txt", amb=ambiguous, conf=confidence_Pyrenophora_teres, min=minid, cov=COV),
    expand("coverage_coord/coverage_coordinates_GCA_905067625_1_Bgtriticale_THUN12_genome_v1_2_genomic_generated_perfect_NEAT_{cov}x{conf}_conf01_min_hit_10_subseq_34373_mapped_to_34373_amb_{amb}_minid_0_{min}.txt", amb=ambiguous, conf=confidence_Blumeria_graminis, min=minid, cov=COV),
    expand("coverage_coord/coverage_coordinates_GCA_000219625_1_MYCGR_v2_0_genomic_generated_perfect_NEAT_{cov}x{conf}_conf01_min_hit_10_subseq_1047171_mapped_to_1047171_amb_{amb}_minid_0_{min}.txt", amb=ambiguous, conf=confidence_Zymoseptoria_tritici, min=minid, cov=COV),
    expand("coverage_coord/coverage_coordinates_GCA_000149985_1_ASM14998v1_genomic_generated_perfect_NEAT_{cov}x{conf}_conf01_min_hit_10_subseq_45151_mapped_to_45151_amb_{amb}_minid_0_{min}.txt", amb=ambiguous, conf=confidence_Pyrenophora_tritici_repentis, min=minid, cov=COV),
    expand("coverage_coord/coverage_coordinates_GCA_021901695_1_Pst134E36_v1_pri_genomic_generated_perfect_NEAT_{cov}x{conf}_conf01_min_hit_10_subseq_27350_mapped_to_27350_amb_{amb}_minid_0_{min}.txt", amb=ambiguous, conf=confidence_Puccinia_striiformis, min=minid, cov=COV),
    #Remove PCR duplicates from reads classified to seqid and then mapped to seqid reference genome
    expand("map_files/bam/no_PCR_merged{conf}_subseq_{seqid}_mapped_to_{seqid}_amb_{amb}_minid_0_{min}.sorted.sam", amb=ambiguous, conf=confidence, min=minid, seqid=obs_species_with_more_than_1000_reads_TAXID),


rule generate_reads_NEAT_custom_coverage_insert_and_insertstd:
  input:
    i1 = path_to_genomes+"{genome_name}.fna"
  output:
    o1 = "generated_reads/NEAT/{genome_name}_generated_perfect_NEAT_{cov}x_insert_{ins}_insstd_{insstd}_read1.fq.gz",
    o2 = "generated_reads/NEAT/{genome_name}_generated_perfect_NEAT_{cov}x_insert_{ins}_insstd_{insstd}_read2.fq.gz"
  threads: 4
  log:
    "logs/generate_insert_and_insertstd/{genome_name}_generated_NEAT_{cov}x_insert_{ins}_insstd_{insstd}.log"
  params:
    path_to_NEAT=path_to_NEAT,
    path_to_genomes=path_to_genomes,
  wildcard_constraints:
    insstd="\d+",
    ins="\d+"
  conda:
    "envs/NEAT.yml"
  shell:
    """
    python {params.path_to_NEAT} -r {params.path_to_genomes}{wildcards.genome_name}.fna -o generated_reads/NEAT/{wildcards.genome_name}_generated_perfect_NEAT_{wildcards.cov}x_insert_{wildcards.ins}_insstd_{wildcards.insstd} -c {wildcards.cov} -R 126 -E 0 --force-coverage --pe {wildcards.ins} {wildcards.insstd} 2> {log}
    """

rule generate_reads_NEAT_custom_coverage_insert_and_insertstd_gzip_input:
  input:
    i1 = path_to_genomes+"{genome_name}.fna.gz"
  output:
    o1 = "generated_reads/NEAT/{genome_name}_generated_perfect_NEAT_{cov}x_insert_{ins}_insstd_{insstd}_read1.fq.gz",
    o2 = "generated_reads/NEAT/{genome_name}_generated_perfect_NEAT_{cov}x_insert_{ins}_insstd_{insstd}_read2.fq.gz"
  threads: 4
  log:
    "logs/generate_insert_and_insertstd/{genome_name}_generated_NEAT_{cov}x_insert_{ins}_insstd_{insstd}.log"
  params:
    path_to_NEAT=path_to_NEAT,
    path_to_genomes=path_to_genomes,
  wildcard_constraints:
    insstd="\d+",
    ins="\d+"
  conda:
    "envs/NEAT.yml"
  shell:
    """
    gzip -d {input.i1}
    python {params.path_to_NEAT} -r {params.path_to_genomes}{wildcards.genome_name}.fna -o generated_reads/NEAT/{wildcards.genome_name}_generated_perfect_NEAT_{wildcards.cov}x_insert_{wildcards.ins}_insstd_{wildcards.insstd} -c {wildcards.cov} -R 126 -E 0 --force-coverage --pe {wildcards.ins} {wildcards.insstd} 2> {log}
    gzip {params.path_to_genomes}{wildcards.genome_name}.fna
    """

rule classify_kraken_reads_generated:
  input:
    i1 = "generated_reads/NEAT/unique/{file_name}_R1.fq.gz",
    i2 = "generated_reads/NEAT/unique/{file_name}_R2.fq.gz"
  output:
    o1 = "kraken_reports/{file_name}_conf01_min_hit_10.report",
    o2 = "kraken_outputs/{file_name}_conf01_min_hit_10.output.gz"
  threads: 36
  log:
    "logs/kraken_logs/{file_name}.log"
  params:
    path_to_db=path_to_db
  conda:
    "envs/Kraken2.yml"
  shell:
    """
    kraken2 --paired --gzip-compressed --threads 36 --minimum-hit-groups 10 --confidence 0.1 --output kraken_outputs/{wildcards.file_name}_conf01_min_hit_10.output --report kraken_reports/{wildcards.file_name}_conf01_min_hit_10.report  --db {params.path_to_db} {input.i1} {input.i2} 2> {log}
    gzip kraken_outputs/{wildcards.file_name}_conf01_min_hit_10.output
    """
rule classify_kraken_reads:
  input:
    i1 = path_to_fq+"{file_name}_unmapped_R1.fq.gz",
    i2 = path_to_fq+"{file_name}_unmapped_R2.fq.gz"
  output:
    o1 = "kraken_reports/{file_name}_conf01_min_hit_10.report",
    o2 = "kraken_outputs/{file_name}_conf01_min_hit_10.output.gz"
  threads: 36
  log:
    "logs/kraken_logs/{file_name}.log"
  params:
    path_to_db=path_to_db
  conda:
    "envs/Kraken2.yml"
  shell:
    """
    kraken2 --paired --gzip-compressed --threads 36 --minimum-hit-groups 10 --confidence 0.1 --output kraken_outputs/{wildcards.file_name}_conf01_min_hit_10.output --report kraken_reports/{wildcards.file_name}_conf01_min_hit_10.report  --db {params.path_to_db} {input.i1} {input.i2} 2> {log}
    gzip kraken_outputs/{wildcards.file_name}_conf01_min_hit_10.output
    """

rule classify_kraken_reads_default_stringency:
  input:
    i1 = path_to_fq+"{file_name}_unmapped_R1.fq.gz",
    i2 = path_to_fq+"{file_name}_unmapped_R2.fq.gz"
  output:
    o1 = "kraken_reports/{file_name}_default_conf_min_hit.report",
    o2 = "kraken_outputs/{file_name}_default_conf_min_hit.output.gz"
  threads: 36
  log:
    "logs/kraken_logs/{file_name}.log"
  params:
    path_to_db=path_to_db
  conda:
    "envs/Kraken2.yml"
  shell:
    """
    kraken2 --paired --gzip-compressed --threads 36 --output kraken_outputs/{wildcards.file_name}_default_conf_min_hit.output --report kraken_reports/{wildcards.file_name}_default_conf_min_hit.report  --db {params.path_to_db} {input.i1} {input.i2} 2> {log}
    gzip kraken_outputs/{wildcards.file_name}_default_conf_min_hit.output
    """

rule give_reads_unique_id:
  input:
    i1 = "generated_reads/NEAT/{file_name}_read1.fq.gz",
    i2 = "generated_reads/NEAT/{file_name}_read2.fq.gz",
  output:
    o1 = "generated_reads/NEAT/unique/{file_name}_R1.fq.gz",
    o2 = "generated_reads/NEAT/unique/{file_name}_R2.fq.gz",
  log:
    "logs/give_reads_unique_id/{file_name}.log"
  threads: 1
  shell:
    """
    python give_reads_unique_id.py {input.i1} {input.i2} {output.o1} {output.o2} 2> {log}
    """

rule bbmap_alig_generated:
  input:
    i1 = "generated_reads/NEAT/unique/{file_name}_R1.fq.gz",
    i2 = "generated_reads/NEAT/unique/{file_name}_R2.fq.gz",
  output:
    o1 = "map_files/bam/{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.bam",
    o2 = "map_files/basecov/basecov_{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.txt",
  log:
    statsfile="logs/bbmap_align/{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.statsfile.txt",
    stderr="logs/bbmap_align/{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.stderr.log",
  threads: 15
  params:
    ref_tmp = "ref_data/{genomeid}",
  wildcard_constraints:
    min="\d+",
  conda:
    "envs/BBMap.yml"
  shell:
    """
    bbmap.sh \
    threads=15 \
    pairedonly=t \
    ambiguous={wildcards.amb} \
    strictmaxindel=t \
    overwrite=t \
    minid=0.{wildcards.min} \
    in1={input.i1} \
    in2={input.i2} \
    path={params.ref_tmp} \
    outm={output.o1} \
    statsfile={log.statsfile} \
    basecov={output.o2}\
    Xmx32g \
    2> {log.stderr}
    """

rule bbmap_alig_subseq:
  input:
    i1 = "seq/{file_name}_R1.fq.gz",
    i2 = "seq/{file_name}_R2.fq.gz",
  output:
    o1 = "map_files/bam/{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.bam",
    o2 = "map_files/basecov/basecov_{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.txt",
  log:
    statsfile="logs/bbmap_align/{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.statsfile.txt",
    stderr="logs/bbmap_align/{file_name}_mapped_to_{genomeid}_amb_{amb}_minid_0_{min}.stderr.log",
  threads: 15
  params:
    ref_tmp = "ref_data/{genomeid}",
  wildcard_constraints:
    min="\d+",
  conda:
    "envs/BBMap.yml"
  shell:
    """
    bbmap.sh \
    threads=15 \
    pairedonly=t \
    ambiguous={wildcards.amb} \
    strictmaxindel=t \
    overwrite=t \
    minid=0.{wildcards.min} \
    in1={input.i1} \
    in2={input.i2} \
    path={params.ref_tmp} \
    outm={output.o1} \
    statsfile={log.statsfile} \
    basecov={output.o2}\
    Xmx32g \
    2> {log.stderr}
    """

rule get_readID:
  input:
    i1 = "kraken_outputs/{file_name}.output.gz",
  output:
    o1 = "readID/{file_name}_taxids_{tax}.readIDs"
  threads: 1
  log:
    "logs/get_readID/{file_name}_taxids_{tax}.log"
  params:
    extract_mode = "clade",
    path_to_selmeout = path_to_selmeout,
    path_to_names = path_to_names,
    path_to_nodes = path_to_nodes,
  shell:
    """
    gzip -d {input.i1}
    python {params.path_to_selmeout} --mode {params.extract_mode} --names {params.path_to_names} --nodes {params.path_to_nodes} --output {output.o1} --input kraken_outputs/{wildcards.file_name}.output --tax_id {wildcards.tax}
    gzip kraken_outputs/{wildcards.file_name}.output
    """

rule get_sequences_generated:
  input:
    i1 = "readID/{file_name}_conf01_min_hit_10_taxids_{tax}.readIDs",
    i2 = "generated_reads/NEAT/unique/{file_name}_R1.fq.gz",
    i3 = "generated_reads/NEAT/unique/{file_name}_R2.fq.gz",
  output:
    o1 = "seq/{file_name}_conf01_min_hit_10_subseq_{tax}_R1.fq.gz",
    o2 = "seq/{file_name}_conf01_min_hit_10_subseq_{tax}_R2.fq.gz",
  threads: 1
  log:
    "logs/get_sequences/{file_name}_conf01_min_hit_10_subseq_{tax}.log"
  conda:
    "envs/seqtk.yml"
  shell:
    """
    python new_readID.py {input.i1} readID/{wildcards.file_name}_conf01_min_hit_10_taxids_{wildcards.tax}_1.readIDs readID/{wildcards.file_name}_conf01_min_hit_10_taxids_{wildcards.tax}_2.readIDs
    seqtk subseq {input.i2} readID/{wildcards.file_name}_conf01_min_hit_10_taxids_{wildcards.tax}_1.readIDs | gzip > {output.o1}
    seqtk subseq {input.i3} readID/{wildcards.file_name}_conf01_min_hit_10_taxids_{wildcards.tax}_2.readIDs | gzip > {output.o2}
    """

rule get_sequences:
  input:
    i1 = "readID/{file_name}_conf01_min_hit_10_taxids_{tax}.readIDs",
    i2 = path_to_fq+"{file_name}_unmapped_R1.fq.gz",
    i3 = path_to_fq+"{file_name}_unmapped_R2.fq.gz",
  output:
    o1 = "seq/{file_name}_conf01_min_hit_10_subseq_{tax}_R1.fq.gz",
    o2 = "seq/{file_name}_conf01_min_hit_10_subseq_{tax}_R2.fq.gz",
  threads: 1
  log:
    "logs/get_sequences/{file_name}_conf01_min_hit_10_{tax}.log"
  conda:
    "envs/seqtk.yml"
  shell:
    """
    seqtk subseq {input.i2} {input.i1} | gzip > {output.o1}
    seqtk subseq {input.i3} {input.i1} | gzip > {output.o2}
    """

rule merge_sequences_subseq:
  input:
    expand(["seq/{sample}{{conf}}_subseq_{{seqid}}_R1.fq.gz"], sample=samples)
  output:
    o1 = "seq/merged{conf}_subseq_{seqid}_R1.fq.gz",
    o2 = "seq/merged{conf}_subseq_{seqid}_R2.fq.gz",
  threads: 1
  wildcard_constraints:
    conf="_conf01_min_hit_10|_default_conf_min_hit",
  params:
    sample_start=sample_start
  shell:
    """
    cat seq/{params.sample_start}*{wildcards.conf}_subseq_{wildcards.seqid}_R1.fq.gz > {output.o1}
    cat seq/{params.sample_start}*{wildcards.conf}_subseq_{wildcards.seqid}_R2.fq.gz > {output.o2}
    """

rule merge_sequences:
  input:
    expand([path_to_fq+"{sample}_unmapped_R1.fq.gz"], sample=samples)
  output:
    o1 = "seq/merged_R1.fq.gz",
    o2 = "seq/merged_R2.fq.gz",
  threads: 1
  params:
    path_to_fq=path_to_fq,
    sample_pattern_R1=sample_pattern_R1,
    sample_pattern_R2=sample_pattern_R2,
  shell:
    """
    cat {params.path_to_fq}{params.sample_pattern_R1} > {output.o1}
    cat {params.path_to_fq}{params.sample_pattern_R2} > {output.o2}
    """

rule get_covoord:
  input:
    i1 = "map_files/basecov/basecov_{file_name}.txt",
  output:
    o1 = "coverage_coord/coverage_coordinates_{file_name}.txt",
  log:
    stderr="logs/coverage_coord/{file_name}.log",
  threads: 5
  shell:
    """
    python get_covcoord.py {input.i1} {output.o1} 2> {log.stderr}
    gzip {input.i1}
    """

rule sort_bam_generated_NEAT:
  input:
    i1 = "{f_name}.bam"
  output:
    o1 = "{f_name}.sorted.bam",
  threads: 5
  conda:
    "envs/samtools.yml"
  shell:
    """
    samtools sort {input.i1} -o {output.o1}
    """

rule bam_to_sam:
  input:
    i1 = "{f_name}.bam"
  output:
    o1 = "{f_name}.sam",
  threads: 5
  conda:
    "envs/samtools.yml"
  shell:
    """
    samtools view -h {input.i1} > {output.o1}
    """

rule run_picard_generated:
  input:
    i1 = "map_files/bam/{file_name}.sorted.bam"
  output:
    o1 = "map_files/bam/no_PCR_{file_name}.sorted.bam",
    o2 = "map_files/bam/no_PCR_{file_name}.sorted.sam"
  log:
    metric="logs/picard/metrics_{file_name}.log",
    stderr="logs/picard/{file_name}.log"
  threads: 5
  conda:
    "envs/picard.yml"
  shell:
    """
    java -jar picard.jar MarkDuplicates -I {input.i1} -O {output.o1} -M {log.metric} -REMOVE_DUPLICATES True 2> {log.stderr}
    samtools view -h {output.o1} > {output.o2} 
    """

rule gzip_basecov:
  input:
    i1 = "map_files/basecov/basecov_{f_name}.txt"
  output:
    o1 = "map_files/basecov/basecov_{f_name}.txt.gz",
  threads: 5
  shell:
    """
    gzip {input.i1}
    """

rule gest_scaff_lengths_gz:
  input:
    i1 = path_to_genomes+"{file_name}.fna.gz"
  output:
    o1 = "scaffold_summaries/{file_name}.txt"
  threads: 1
  params:
    path=path_to_genomes
  log:
    stderr="logs/calc_scaff_lenths/{file_name}.log"
  conda:
    "envs/python.yml"
  shell:
    """
    gzip -d {input.i1}
    python get_scaffold_lengths_from_fasta.py {params.path}{wildcards.file_name}.fna {output.o1}
    gzip {params.path}{wildcards.file_name}.fna
    """

rule gest_scaff_lengths:
  input:
    i1 = path_to_genomes+"{file_name}.fna"
  output:
    o1 = "scaffold_summaries/{file_name}.txt"
  threads: 1
  params:
    path=path_to_genomes
  log:
    stderr="logs/calc_scaff_lenths/{file_name}.log"
  conda:
    "envs/python.yml"
  shell:
    """
    python get_scaffold_lengths_from_fasta.py {input.i1} {output.o1}
    """