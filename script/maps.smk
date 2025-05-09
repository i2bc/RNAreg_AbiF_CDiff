# aim: raw MAPS sequence data analysis workflow 
#
# Author: Claire Toffano-Nioche, Joe Ueda
# I2BC, France
# claire.toffano-nioche @ i2bc.paris-saclay.fr
# joe.ueda @ i2bc.paris-saclay.fr
# 
# april 2025
# usage example on I2BC cluster:
# qsub smk_maps.qsub
# module load snakemake/snakemake-8.4.6
# snakemake -c3 -s maps.smk --configfile maps.yml --use-conda
# or
# module load snakemake/snakemake-8.4.6
# snakemake -c3 -s maps.smk --configfile maps.yml --use-conda
# 



SAMPLES, = glob_wildcards(config["dataDir"]+"{sample}"+config["pePrefix"]+"1"+config["peSuffix"]+config["fastqSuff"])
BWT2_IDX = ["1","2","3","4","rev.1","rev.2"]


rule all:
  input:   
    expand(config["cwdDir"]+config["fastqcDir"]+"{sample}"+config["pePrefix"]+"1"+config["peSuffix"]+"_fastqc.html", sample=SAMPLES),
    expand(config["cwdDir"]+config["fastqcDir"]+"{sample}"+config["pePrefix"]+"2"+config["peSuffix"]+"_fastqc.html", sample=SAMPLES),
    expand(config["cwdDir"]+config["trimDir"]+"{sample}"+config["pePrefix"]+"1"+config["trimSuff"], sample=SAMPLES),
    expand(config["cwdDir"]+config["trimDir"]+"{sample}"+config["pePrefix"]+"2"+config["trimSuff"], sample=SAMPLES),
    expand(config["cwdDir"]+config["trimDir"]+"{sample}"+"report.json", sample=SAMPLES),
    expand(config["cwdDir"]+config["trimDir"]+"{sample}"+"report.html", sample=SAMPLES),
    expand(config["cwdDir"]+config["bwt2Dir"]+config["genomeIdxPrefix"]+".{bwt2_idx}.bt2", bwt2_idx=BWT2_IDX),
    expand(config["cwdDir"]+config["bwt2Dir"]+"{sample}.bam", sample=SAMPLES),
    expand(config["cwdDir"]+config["bwt2Dir"]+"{sample}.flagstat", sample=SAMPLES),
    expand(config["cwdDir"]+config["bwt2Dir"]+"{sample}.bam.bai", sample=SAMPLES),
    expand(config["cwdDir"]+config["featureCountsDir"]+"{sample}_fc.txt", sample=SAMPLES),
    expand(config["cwdDir"]+config["featureCountsDir"]+"{sample}_fc4DEG.txt", sample=SAMPLES),
    expand(config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete.txt"),
    expand(config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+"_shortName_oldLocusTag.txt"),
    expand(config["cwdDir"]+"volcano_plot.pdf")





rule bowtie2_genome_indexing:
  output:
    expand(config["cwdDir"]+config["bwt2Dir"]+config["genomeIdxPrefix"]+".{ext}.bt2", ext=BWT2_IDX)
  input:
    fasta=config["cwdDir"]+config["genomeFna"]
  params:
    dir=config["cwdDir"]+config["bwt2Dir"]+config["genomeIdxPrefix"]
  resources:
    mem_mb=2000  
  threads : 1
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"bowtie2.yml"
  log:
    out=config["cwdDir"]+"Logs/genome_bwt2_index.stdout",
    err=config["cwdDir"]+"Logs/genome_bwt2_index.stderr"
  shell:
    "bowtie2-build --threads {threads} {input.fasta} {params.dir} 1>{log.out} 2>{log.err}"



rule fastqc:
  output:
    config["cwdDir"]+config["fastqcDir"]+"{sample}"+config["pePrefix"]+"1"+config["peSuffix"]+"_fastqc.zip",
    config["cwdDir"]+config["fastqcDir"]+"{sample}"+config["pePrefix"]+"2"+config["peSuffix"]+"_fastqc.zip",
    config["cwdDir"]+config["fastqcDir"]+"{sample}"+config["pePrefix"]+"1"+config["peSuffix"]+"_fastqc.html",
    config["cwdDir"]+config["fastqcDir"]+"{sample}"+config["pePrefix"]+"2"+config["peSuffix"]+"_fastqc.html"
  input: 
    r1=config["dataDir"]+"{sample}"+config["pePrefix"]+"1"+config["peSuffix"]+config["fastqSuff"],
    r2=config["dataDir"]+"{sample}"+config["pePrefix"]+"2"+config["peSuffix"]+config["fastqSuff"]
  params:
    outDir=config["fastqcDir"]
  log:
    out_r1=config["cwdDir"]+"Logs/{sample}"+config["pePrefix"]+"1_fastqc.stdout",
    err_r1=config["cwdDir"]+"Logs/{sample}"+config["pePrefix"]+"1_fastqc.stderr",
    out_r2=config["cwdDir"]+"Logs/{sample}"+config["pePrefix"]+"2_fastqc.stdout",
    err_r2=config["cwdDir"]+"Logs/{sample}"+config["pePrefix"]+"2_fastqc.stderr"
  conda: 
    config["cwdDir"]+config["condaEnvDir"]+"fastqc.yml"
  shell: 
    """
       fastqc --outdir {params.outDir} {input.r1} 1>{log.out_r1} 2>{log.err_r1} ;
       fastqc --outdir {params.outDir} {input.r2} 1>{log.out_r2} 2>{log.err_r2}
    """



rule trimming:
  output:
    r1trimmed = config["cwdDir"]+config["trimDir"]+"{sample}"+config["pePrefix"]+"1"+config["trimSuff"],
    r2trimmed = config["cwdDir"]+config["trimDir"]+"{sample}"+config["pePrefix"]+"2"+config["trimSuff"],
    json = config["cwdDir"]+config["trimDir"]+"{sample}"+"report.json",
    html = config["cwdDir"]+config["trimDir"]+"{sample}"+"report.html"
  input:
    r1=config["dataDir"]+"{sample}"+config["pePrefix"]+"1"+config["peSuffix"]+config["fastqSuff"],
    r2=config["dataDir"]+"{sample}"+config["pePrefix"]+"2"+config["peSuffix"]+config["fastqSuff"]
  log:
    out=config["cwdDir"]+"Logs/{sample}_trimming.stdout",
    err=config["cwdDir"]+"Logs/{sample}_trimming.stderr"
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"fastp.yml"
  shell:
    """
        fastp -i {input.r1} -I {input.r2} -o {output.r1trimmed} -O {output.r2trimmed} -j {output.json} -h {output.html} 1> {log.out} 2> {log.err}
    """




rule bowtie2_mapping:
  output:
    temp(config["cwdDir"]+config["bwt2Dir"]+"{sample}.sam")
  input:
    idx=rules.bowtie2_genome_indexing.output,
    r1trimmed = config["cwdDir"]+config["trimDir"]+"{sample}"+config["pePrefix"]+"1"+config["trimSuff"],
    r2trimmed = config["cwdDir"]+config["trimDir"]+"{sample}"+config["pePrefix"]+"2"+config["trimSuff"]
  params:
    index=config["cwdDir"]+config["bwt2Dir"]+config["genomeIdxPrefix"],
    orientation=config["bwtSequencingOrientation"]
  log:
    out=config["cwdDir"]+"Logs/{sample}_bwt2_mapping.stdout",
    err=config["cwdDir"]+"Logs/{sample}_bwt2_mapping.stderr"
  benchmark:
    "benchmarks/{sample}_bwt2_mapping.benchmark.txt"
  resources: 
    part="long",
    mem_mb=2000
  threads: 8
  conda: 
    config["cwdDir"]+config["condaEnvDir"]+"bowtie2.yml"
  shell:
    "bowtie2 -p {threads} --omit-sec-seq --sam-no-qname-trunc {params.orientation} -x {params.index} -1 {input.r1trimmed} -2 {input.r2trimmed} -S {output} 1> {log.out} 2> {log.err}"


rule bwt_sam_2_sorted_bam:
  output:
    config["cwdDir"]+config["bwt2Dir"]+"{sample}.bam"
  input:
    config["cwdDir"]+config["bwt2Dir"]+"{sample}.sam"
  log:
    out=config["cwdDir"]+"Logs/{sample}_bwt_sam_2_sorted_bam.stdout",
    err=config["cwdDir"]+"Logs/{sample}_bwt_sam_2_sorted_bam.stderr"
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"samtools.yml"
  shell:
    "samtools sort -O BAM -o {output} {input} 1> {log.out} 2> {log.err}"


rule bwt_bam_2_bai_and_stat:
  output:
    bai=config["cwdDir"]+config["bwt2Dir"]+"{sample}.bam.bai",
    fst=config["cwdDir"]+config["bwt2Dir"]+"{sample}.flagstat"
  input:
    config["cwdDir"]+config["bwt2Dir"]+"{sample}.bam"
  log:
    out_bai=config["cwdDir"]+"Logs/{sample}_bwt_bam_2_bai.stdout",
    err_bai=config["cwdDir"]+"Logs/{sample}_bwt_bam_2_bai.stderr",
    err_stat=config["cwdDir"]+"Logs/{sample}_bwt_bam_2_stat.stderr"
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"samtools.yml"
  shell:
    """
       samtools index {input} {output.bai} 1> {log.out_bai} 2> {log.err_bai} ;
       samtools flagstat {input} 1> {output.fst} 2> {log.err_stat}
    """



rule featureCount:
  output: config["cwdDir"]+config["featureCountsDir"]+"{sample}_fc.txt"
  input: config["cwdDir"]+config["bwt2Dir"]+"{sample}.bam"
  params:
    annotation=config["annotation"],
    feature=config["fc_feature"],
    tag=config["fc_tag"],
    mode=config["fc_mode"]
  log:
    out_fc=config["cwdDir"]+"Logs/{sample}_fcount.stdout",
    err_fc=config["cwdDir"]+"Logs/{sample}_fcount.stderr" 
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"count.yml"
  threads: 8
  shell:
    """
       featureCounts -T {threads} {params.mode} -t {params.feature} -g {params.tag} -a {params.annotation} -o {output} {input} 1> {log.out_fc} 2> {log.err_fc}
    """


rule prepareCounts:
  output: config["cwdDir"]+config["featureCountsDir"]+"{sample}_fc4DEG.txt"
  input: config["cwdDir"]+config["featureCountsDir"]+"{sample}_fc.txt"
  params:
    awkscript=config["cwdDir"]+config["awksc"]
  shell:
    """
    echo {wildcards.sample}
    awk -f {params.awkscript} {input} > {output}
    """
  

rule analyseDiff:
  output:
    RData=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/"+config["st_comparison"]+".RData", 
    html=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/"+config["st_comparison"]+"_report.html",
    table=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete.txt"
  input: 
    design=config["cwdDir"]+config["st_designDir"]+config["st_comparison"]+".design4sartools",
    files=expand(config["cwdDir"]+config["featureCountsDir"]+"{sample}_fc4DEG.txt", sample=SAMPLES)
  params:
    resultDir=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/",
    RData=config["st_comparison"]+".RData",
    html=config["st_comparison"]+"_report.html",
    script=config["cwdDir"]+config["st_script"],
    compName=config["st_comparison"],
    rawDir=config["cwdDir"]+config["featureCountsDir"],
    group=config["st_group"],
    condRef=config["st_condRef"]
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"sartools.yml"
  shell:
    """
    Rscript {params.script} --projectName="{params.compName}" --targetFile="{input.design}" --rawDir="{params.rawDir}" --varInt="{params.group}" --condRef="{params.condRef}"
    mv tables figures {params.html} {params.RData} {params.resultDir}
    """
  

rule addOldLocusTag:
  output:
    finalTable=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+"_shortName_oldLocusTag.txt"
  input:
    table=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+".complete.txt",
    conversionOLT=config["cwdDir"]+config["conversion_table4locusTag"],
    conversionSN=config["cwdDir"]+config["conversion_table4shortName"]
  params:
    old_script=config["cwdDir"]+config["olt_script"],
    shortName_script=config["cwdDir"]+config["sn_script"],
    only_oldLoucTag=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+"_with_OldLocusTag.txt"
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"python.yml"
  shell:
    """
    python {params.old_script} {input.table} {input.conversionOLT} {params.only_oldLoucTag}
    python {params.shortName_script} {params.only_oldLoucTag} {input.conversionSN} {output.finalTable}
    """

rule volcanoPlot:
  output:
    volcanoPlot=config["cwdDir"]+"volcano_plot.pdf"
  input: 
    updatedTable=config["cwdDir"]+config["st_dir"]+config["st_comparison"]+"/tables/"+config["st_condition"]+"_shortName_oldLocusTag.txt"
  params:
    vp_script=config["cwdDir"]+config["volcano_script"],
    path2config=config["cwdDir"]+config["volcano_config"]
  conda:
    config["cwdDir"]+config["condaEnvDir"]+"R-EnhancedVolcano_ce.yml"
  shell:
    """
    Rscript {params.vp_script} {params.path2config} {input.updatedTable}
    """
  
