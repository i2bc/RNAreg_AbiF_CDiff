#
# convert the S7 table saved (excluding the column 11) in tsv format 
# into a file for each position relative to Abi.
#
# Abi stands in pos_0.txt
# upstream or donwstream proteins stand in pos_-x.txt or pos_x.txt 
#
# command line:
# awk -f Abif_environment_2_pos.awk S7_table.tsv
#
function whereisabi(nbline, abivalue, abilines) {
   pos=-1
   i=1
   while((i<=nbline)&&(pos==-1)){ 
      if(match(abivalue,abilines[i])==1){
         pos=i
      }
      i++
   }
   return pos
}
function printabi(nbpos, pos, abilines) {
   for(i=1;i<=nbpos;i++){  
      posfile=i-pos
      print abilines[i]>"pos_"posfile".txt"
   }
   return 1
}
BEGIN{
   FS="\t"
   tmp=1
   abi=""
   nbpos=0
}
{
   # new abi:
   if($1=="==="){      
     # manage the current abi:
     abipos=whereisabi(nbpos, abi, memo)
     ok=printabi(nbpos, abipos, memo)
     # init for the next abi:
     abi=$2 
     nbpos=$3
     tmp=1
   }else{
      # memorize lines
      memo[tmp++]=$0
   }
}
END{ # last abi:
   abipos=whereisabi(nbpos, abi, memo)
   ok=printabi(nbpos, abipos, memo) 
}  
