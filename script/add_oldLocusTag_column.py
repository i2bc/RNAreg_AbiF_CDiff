# Script for adding the column old_locus_tag based on a conversion table
# 
# use:
#   python add_oldLocusTag_column.py results/long/new_params/tables/maps_r20_newParams.txt table_locusTag_oldLocusTag.csv long_new_with_oldLocusTag.txt


import pandas as pd
from collections import defaultdict
import sys


# read arguments
original_MAPS_table = sys.argv[1] # original MAPS table
locusTag_oldLocuTag_table = sys.argv[2] # conversion table (csv file): locus_tag,old_locus_tag
output_MAPStable_name = sys.argv[3]  # name of the output file


def update_maps_with_oldLT(maps_file, tag_table, output_file):
    # Read the MAPS table
    maps_data = pd.read_csv(maps_file, sep='\t', header=0)  # Ensure header is read
    maps_data.insert(loc=1, column='old_locus_tag', value=['NA' for _ in range(maps_data.shape[0])])

    # Read locus_tag table as dict. format= locus_tag : old_locus_tag
    locus_tag_data = open(tag_table,"r")
    locus_tag_dict={}
    for i in locus_tag_data:
        if i.startswith("CDR20291"):
            locus_tag_dict[i.replace("\n","").split(",")[0]]=i.replace("\n","").split(",")[1]


    # Add corresponding old_locus_tag to the table
    for index, row in maps_data.iterrows():
        id_new = row.iloc[0] 
        id_old = locus_tag_dict.get(id_new, 'NA') # get old_locus_tag, otherwise NA
        
        # Update table with found old_locus_tag
        maps_data.at[index, 'old_locus_tag'] = id_old

    # Write new table
    maps_data.to_csv(output_file, sep='\t', header=True, index=False)
    

update_maps_with_oldLT(original_MAPS_table, locusTag_oldLocuTag_table, output_MAPStable_name)
