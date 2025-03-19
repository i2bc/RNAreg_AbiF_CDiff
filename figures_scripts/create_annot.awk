
# lineage cut following a number level 
# level is between 1 to 8:
# 'Phylum','Class','Order','Family','Genus','Species','Subspecies','unique_gi'
#
# -v group=8
# -v lincol=4
#
# awk -f create_annot.awk -v group=7 -v lincol=4 -v meancol=2 distribution_abiF.csv
#
function lin_cut(complete_lineage, level) { 
   nb=split(complete_lineage,lineage,";")
   lineage_cut=lineage[1]
   for(i=2;i<level;i++){
      lineage_cut=sprintf("%s;%s",lineage_cut,lineage[i])
   }
   if(lineage[level]==""){
      lineage_cut=sprintf("%s;%s",lineage_cut,"other")
   }else{
      lineage_cut=sprintf("%s;%s",lineage_cut,lineage[level])
   }
   return lineage_cut
}

BEGIN{
   FS="\t"
   mem=""
   nbl=0
   abiF=0
   if(group==0){group=8}
   if(lincol==0){lincol=4}
   if(meancol==0){lincol=2}
}
{
   if(NR==1){
      mem=lin_cut($lincol,group)
   }
   current_lin=lin_cut($lincol,group)
   if(current_lin!~mem){
            	#print("div par 0 ?", abiF,"/",nbl" ; "mem" next="$0)
      printf("%s,%.2f\n",mem,abiF/nbl)
      mem=current_lin
      nbl=1
      abiF=$meancol
   }else{
      nbl++
      abiF=abiF+$meancol
   }
}
END{
   printf("%s,%.2f\n",mem,abiF/nbl)
}
