# Deciphering the RNA-based regulation mechanism of phage-encoded AbiF system in Clostridioides difficile

This GitHub repository accompanies the preprint where we report the identification of a new AbiF-like system within a prophage of an hypervirulent *Clostridioides difficile* R20291 strain: ["Deciphering the RNA-based regulation mechanism of phage-encoded AbiF system in Clostridioides difficile"](https://www.biorxiv.org/).

## Contents

- [Dataset](#dataset)
- [Conservation](#conservation) of the AbiF-like system of the hypervirulent ribotype 027 strain
- [MAPS](#MAPS) analysis for the RCd22 ncRNA
- [Citation](#citation)

## Dataset

- genome and annotation of the *C. difficile* R20291 strain used: download both the genome (`fna`) and annotation (`gff3`) files from the [GCF_000027105.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000027105.1/) ncbi assembly. Add to this annotation file the identification of ncRNA identified in CD360 (`data/RCdxxx.gff`)
- [Similarity search of the AbiF-like system](zenodo_xxx or MatSupp) containing `20241029_distribution_abiF.csv`
- [genomic contexte of AbiF_like systems](zenodo_xxx or MatSupp) containing:
```bash
 9989 Positions_genes_autours_abi-2_ou_abiF_5neg.csv
10012 Positions_genes_autours_abi-2_ou_abiF_4neg.csv
10036 Positions_genes_autours_abi-2_ou_abiF_3neg.csv
10053 Positions_genes_autours_abi-2_ou_abiF_2neg.csv
10080 Positions_genes_autours_abi-2_ou_abiF_1neg.csv
10085 Positions_genes_autours_abi-2_ou_abiF_1pos.csv
10066 Positions_genes_autours_abi-2_ou_abiF_2pos.csv
10047 Positions_genes_autours_abi-2_ou_abiF_3pos.csv
10031 Positions_genes_autours_abi-2_ou_abiF_4pos.csv
10006 Positions_genes_autours_abi-2_ou_abiF_5pos.csv
```
- MAPS experiments data: stand in two parts: RNAseq fraction ([ENA_xxx](ENA_xxx)) and [proteic fraction](zenodo_xxx, MatSupp, or github/data `R20291_RCd22_Soutourina_120723.xlsx` = 145K ; supprimer les colonnes "long" ?)

- [color selection for the figure of AbiF-like system conservation](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/data/color_selection.tsv)


## Conservation 

### Similarity search of the AbiF-like system 

The AbiF-like system identified in the *C. difficile* hypervirulent ribotype 027 strain were searched in an in-house database of 47545 completely sequences genomes (chromosome level assembly, November 2023) both from Refseq and Genbank databases with PSI-BLAST (2.16.1 version). It was run against COG, CD and PFAM profiles from CDD database (E-value=1e-4), other parameters were default.
The resulting table `20241029_distribution_abiF.csv` (see [Dataset](#dataset) section) contains the number of valid hits by genomes.
The genomic context of the resulting hits including 5 proteins upstream and downstream from the location of the AbiF-like system where extracted and can be found in the corresponding 10 separate files (see [Dataset](#dataset) section)

### AbiF-like system representation in the bacterial species tree (figure)

#### Setup 

<!--
- [conda environement for figure generation](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/conda_environments/fig_conda_env.yaml)
- [preprocess script for figure](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/script/figure.R)
-->
- iTOL templates to download: [colors](https://itol.embl.de/help/dataset_color_strip_template.txt) and [simplebar](https://itol.embl.de/help/dataset_simplebar_template.txt)
- python3


#### species tree 

Create the tree (Newick format) from the taxonomy present in the AbiF resulting table `20241029_distribution_abiF.csv`:

```bash
awk -F "\t" '{if($4~/Bacteria/){print $4"\t"$3"\t"$1}}' 20241029_distribution_abiF.csv | sed 's/:/-/;s/\t/;/g' | sort > full_lineage_bact.tsv
python3 scripts/convert_lineage_to_nwk.py full_lineage_bact.tsv > full_lineage_bact.nwk
```
Firstly, the NCBI bacterial lineage (4th column) is completed by the strain name (3rd column) and the assembly ID (1st column), the characters ":" not allowed in iTOL are replaced by "-" and the separators ";" in the lineage are replaced by tab characters. 
Next, this one-line-per-species format is converted into a tree following the Newick format using the python code `convert_lineage_to_nwk.py` taken from the [Mark Watson stackoverflow response](https://stackoverflow.com/questions/26146623/convert-csv-to-newick-tree/26147279#26147279) and adapted for python 3.

Resulting file: `full_lineage_bact.nwk`


#### species colors 

For ease of viewing, only 23 taxonomic levels have been manually selected for color display : `data/color_selection.tsv`.
This selection was chosen according to different taxonomic levels to highlight *C. difficile* and *Staphylococcus aureus* and visualize those with counts of zero (if more than 160 genomes). 
The other selected taxonomic levels constitute a balanced grouping according to the number of genomes in which the AbiF_like system was searched.

Into the iTOL templates [colors](https://itol.embl.de/help/dataset_color_strip_template.txt), manually change:
- SEPARATOR SPACE for TAB 
- COLOR_BRANCHES from 0 to 1 
- DATASET_LABEL label1 to Philum
and add the `color_selection.tsv` at the end:

```bash
cat dataset_color_strip_template.txt data/color_selection.tsv > full_lineage_bact_colors.txt
```
Resulting file: `full_lineage_bact_colors.txt`


#### AbiF distribution

On the [simplebar](https://itol.embl.de/help/dataset_simplebar_template.txt) iTOL template file, uncomment the `WIDTH,1000` line and change its value from 1000 to 200.
Complete this `dataset_simplebar_template.txt` with the AbiF counts (2nd column) present in each assembly (1rst column) of the resulting AbiF table `20241029_distribution_abiF.csv`:
```bash
cp dataset_simplebar_template.txt full_lineage_bact_simplebar_abiF.txt
awk -F "\t" '{if($4~/Bacteria/){print $1","$2}}' 20241029_distribution_abiF.csv | sort >> full_lineage_bact_simplebar_abiF.txt
```
Resulting file: `full_lineage_bact_simplebar_abiF.txt`


#### Figure creation on the iTOL web page

On the [iTOL v7](https://itol.embl.de/upload.cgi) web page:

Import the species tree structure ("Upload" tab with `full_lineage_bact.nwk`)
As the tree contains too many leaves, iTOL by default compresses some of the data we want to display. 
To correct this, performe manually the steps below in the iTOL "Control panel" and "Basic" tab :
- place "Labels" to "Hide"
- change 350° to 359° on "Mode option" and "Arc"
- choose grey color (#5b5b5B) for "Line style"
- select "Un-collapse all" on "Advanced" tab and "Nodes options"

Add colors on the tree by "Upload annotation files" (in "Control panel" then "Datasets" tab) with `full_lineage_bact_colors.txt`
<!-- problème : variable COLOR #ff0000 unknown  ?? en trop ??? -->
<!-- ajouter la légende "automatique" sous le bouton "roue crantée" -->

Add the AbiF distribution : in "Datasets" table, "Upload" tab with `full_lineage_bact_simplebar_abiF.txt` 


Save the figure.

### Environment of AbiF-like (figure)

#### Setup 

- [conda environement for figure generation](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/conda_environments/fig_conda_env.yaml)
- [Rscript for figure generation](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/script/genomiccontext.R)

The figure showing the genomic context of AbiF-like systems was produced with the R script `genomiccontext.R` after a few steps of filtering and reorganizing the genomic context data of AbiF-like systems (see Process below).
If required, the R-base version used (4.2.3) and added with the ggplot2 package (3.4.4) can be installed with the conda environment file `conda_environments/env_Rbase4.2.3.yml`: `conda env create -f conda_environments/env_Rbase4.2.3.yml`

#### Process

10 input files: see "genomic contexte of AbiF_like systems" in Dataset

select COG/Pfam or CDD if > 100 occurencies from the 10 input files:
```bash
rm abi_COG_up100.txt ; for p in 5neg 4neg 3neg 2neg 1neg 1pos 2pos 3pos 4pos 5pos ; do 
   sort -t $'\t' -k9,9 Positions_genes_autours_abi-2_ou_abiF_$p.csv > Pos_abi_$p.tmp 
   join -t $'\t' -1 9 -2 1 -e "empty" -a 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3 Pos_abi_$p.tmp col9_up100_cdd_pfam_COG.txt > Pos_abi_$po.csv 
done
```
format for graph:
```bash
echo "position"$'\t'"abundance"$'\t'"COG"$'\t'"classe" > abi_COG_up100_graph.txt ; 
for p in 5neg 4neg 3neg 2neg 1neg abi 1pos 2pos 3pos 4pos 5pos ; do 
   awk -F "\t" '{if(length($14)>=1){print $14"\t"$15}}' Pos_abi_$p.csv | sort -t $'\t' -k1,1 | uniq -c | awk -v pos=${p} 'BEGIN{nbl=0;somme=0;OFS="\t"}{if(pos!="abi"){if($1>100){funcCl="";for(i=3;i<=NF;i++){funcCl=sprintf("%s %s",funcCl, $i)};print pos,$1,$2,funcCl;somme=somme+$1};nbl=nbl+$1}}END{if(pos=="abi"){print pos, 10000, "V", "Defense mechanisms"}else{print pos,nbl-somme,"others","Other COG, pfam or CDD < 100"}}' | sed 's/5neg/-5/;s/5pos/+5/;s/4neg/-4/;s/4pos/+4/;s/3neg/-3/;s/3pos/+3/;s/2neg/-3/;s/2pos/+2/;s/1neg/-1/;s/1pos/+1/' >> abi_COG_up100_graph.txt
done
```
create graph:
```bash
conda activate Rbase4.2.3 
# pour les tests : conda activate /DATA/miniconda3/envs/Rbase4.2.3
Rscript abi_context_graph.R
```

result: `abi_context_graph.png`
 
 
## MAPS
 
MAPS (MS2‐affinity purification coupled with RNA sequencing) analysis for the RCd22 ncRNA

### MAPS rnaseq fraction

Setup:
- [conda environement for maps analysis](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/conda_environments/maps_conda_env.yaml)
- [snakemake pipeline for maps analysis](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/script/maps.smk)

MAPS data (rnaseq part): 
4 fastq RCd22 : comparison of a *C. difficile* R20291 containing a plasmid with an MS2 tag followed (rcd22) or not (crtl) by RCd22 ncRNA (2 replicates for each condition, see Dataset section).

MAPS analysis: 
The analysis workflow can be found in the snakemake file (see `scripts/maps.smk`) and the parameters are set in the corresponding configuration file (`scripts/maps.yml`). The workflow basically contains the following steps: quality control of reads (fastqc, fastqscreen) and correction if necessary (fastp), mapping to the genome sequence (bowtie2), count of the number of reads per gene (featurecount), differential capture analysis (SARTools used with DESeq2 mode), and creation of graph (R-enhancedVolcano).
Software versions are described in the associated conda environments (`conda_environments/*.yml`).

Results: 
complete, up, and down tables of differentially captured RNA and the associated volcano-plot will be provided in the specified `st_dir/st_comparison/` (see the fixed values of `st_dir` and `st_comparison` in the `scripts/maps.yml`) repository.

### MAPS proteic fraction

Setup:
- [conda environement for maps analysis](https://github.com/i2bc/RNAreg_AbiF_CDiff/blob/main/conda_environments/env_Rbase4.2.3.yml)
- download the `template_script_DESeq2_CL.r` from the [SARTools](https://github.com/PF2-pasteur-fr/SARTools) github pages.

MAPS data (proteic part): 
The xlxs file was manually saved with a csv format (with tabulation separator).
Total spectrum counts columns are broken down as follows: E464 and E464.2 (columns 10 and 11) stand for the 2 RCd22 replicates while E608 and E608.2 denote control replicates (columns 14 and 15). 

From the csv file, count and ID columns for each sample were extracted as follow:
```bash
mkdir MAPSprot
awk -v col=10 'BEGIN{FS="\t";OFS="\t"}{if($col==""){$col=0};nbsplit=split($2, id , "|");if(nbsplit==3){print id[3],$col}}' R20291_RCd22_Soutourina_120723.csv | awk '{print $1"\t"$NF}' > MAPSprot/rcd22_1.txt
awk -v col=11 'BEGIN{FS="\t";OFS="\t"}{if($col==""){$col=0};nbsplit=split($2, id , "|");if(nbsplit==3){print id[3],$col}}' R20291_RCd22_Soutourina_120723.csv | awk '{print $1"\t"$NF}' > MAPSprot/rcd22_2.txt
awk -v col=14 'BEGIN{FS="\t";OFS="\t"}{if($col==""){$col=0};nbsplit=split($2, id , "|");if(nbsplit==3){print id[3],$col}}' R20291_RCd22_Soutourina_120723.csv | awk '{print $1"\t"$NF}' > MAPSprot/msctr_1.txt
awk -v col=15 'BEGIN{FS="\t";OFS="\t"}{if($col==""){$col=0};nbsplit=split($2, id , "|");if(nbsplit==3){print id[3],$col}}' R20291_RCd22_Soutourina_120723.csv | awk '{print $1"\t"$NF}' > MAPSprot/msctr_2.txt
```
We apply a differential expression analysis to these counts (in the activated conda environment):
```bash
conda activate Rbase4.3.1_ce
Rscript template_script_DESeq2_CL.r --projectName="RCd22-MAPS_MS" --targetFile="data/MassSpec.design4sartools" --rawDir="MAPSprot/" --varInt="group" --condRef="Crtl" 
```
Results:
complete, up, and down tables of differentially total spectrum counts and the associated volcano-plot will be provided in the `MAPSprot/tables` and `MAPSprot/figures` repositories.

## Reference

- conda
- Psi-blast (v2.16.1)
- Refseq and Genbank databases
- CDD database
- snakemake
- fastqc
- fastqScreen
- fastp
- bowtie2
- subread
- SARtools
- DESeq2
- R-enhancedVolcano


## Citation

If you find these pages useful for your research, please cite the relevant paper

```
@article{saunier_abiF_2025,
  title = {Deciphering the RNA-based regulation mechanism of phage-encoded AbiF system in Clostridioides difficile},
  author = {Saunier, Marion and Soutourina, Olga and ...},
  journal = {Bio},
  year = {2025},
  url = {https://www.biorxiv.org/}
}
```
