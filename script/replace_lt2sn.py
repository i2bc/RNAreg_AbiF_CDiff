# Script for replacing locus_tag with the short name in the DESeq2 table (RCd22)
#
# use:
#   python replace_lt2sn.py ../long/old_params/tables/maps_r20_oldParams.txt ../../ad_shortname/RCd22_shortNames.csv  maps_r20_long_oldParams_renamed.txt
#   


import sys

deseq_path = sys.argv[1]    # path to deseq2 output table
convert_path = sys.argv[2]  # path to conversion table (csv file, header: locusTag,shortName)
outputName = sys.argv[3]    # name of the output file 


# read deseq table
with open(deseq_path) as f:
    listLines = [line.split("\t") for line in f]
'''
for i in listLines:
    print(i)
'''
#read conversion table
with open(convert_path) as f:
    lines = [line.replace("\n","") for line in f]
dictNames = {}
for i in lines:
    if i.startswith("\"locusTag")==False:
        split=i.split(",")
        dictNames[split[0]]=split[1]
'''
for i in dictNames:
    print(i, dictNames[i])
'''


for k in listLines:
    if k[0].startswith("CDR20291"):
        id_lt=k[0]
        if id_lt in dictNames:
            k[0] = dictNames[id_lt]

outputRenamedFile = open(outputName, 'w')  

for i in listLines:
    outputRenamedFile.writelines("\t".join(i))