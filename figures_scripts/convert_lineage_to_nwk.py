#! python3
#
# code from Mark Watson post: 
# https://stackoverflow.com/questions/26146623/convert-csv-to-newick-tree/26147279#26147279
# 
# e.g. input file format:
# Archaea;;;;;;archaeon GW2011_AR20
# Archaea;Bathyarchaeota;Bathyarchaeia;Bathyarchaeales;Bathyarchaeaceae;Bathyarchaeum;Bathyarchaeum sp.
# Archaea;Bathyarchaeota;;;;;Bathyarchaeota archaeon
# Archaea;Euryarchaeota;Archaeoglobi;Archaeoglobales;Archaeoglobaceae;Archaeoglobus;Archaeoglobus fulgidus
#
# output format: newick tree (without distances)
#
# command line:
# python3 convert_lineage_to_nwk.py lineage.csv
#
import csv
from collections import defaultdict
from pprint import pprint
import sys

def tree(): return defaultdict(tree)

def tree_add(t, path):
  for node in path:
    t = t[node]

def pprint_tree(tree_instance):
    def dicts(t): return {k: dicts(t[k]) for k in t}
    pprint(dicts(tree_instance))

def csv_to_tree(input):
    t = tree()
    for row in csv.reader(input, quotechar='\''):
        tree_add(t, row)
    return t

def tree_to_newick(root):
    items = []
    for k in root:
        s = ''
        if len(root[k].keys()) > 0:
            sub_tree = tree_to_newick(root[k])
            if sub_tree != '':
                s += '(' + sub_tree + ')'
        s += k
        items.append(s)
    return ','.join(items)

def csv_to_weightless_newick(input):
    t = csv_to_tree(input)
    #pprint_tree(t)
    return tree_to_newick(t)

if __name__ == '__main__':
    
    input = [
        "'Phylum','Class','Order','Family','Genus','Species','Subspecies','unique_gi'", 
        "'Phylum','Class','Order','example'",
        "'Another','Test'",
    ]
    
    with open(sys.argv[1],"r") as f:
    	input = [",".join(list(map("'{}'".format,l.strip().split(";")))) for l in f.readlines()]
    	#print(input[:5])

    print(csv_to_weightless_newick(input))
