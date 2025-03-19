{
if(NR==1){mem=$4}
if ($4 != mem) {
   printf("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",mem,count0,count1,count2,count3,count4,count5,count6,count7,count8)
   count0=0;count1=0;count2=0;count3=0;count4=0;count5=0;count6=0;count7=0;count8=0;
   mem=$4
} 
if ($2==0) { 
   count0++
} else { 
   if ($2==1) { 
      count1++ 
   } else { 
      if ($2==2) { 
         count2++ ;
      } else { 
         if ($2==3) { 
            count3++ ;
         } else { 
            if ($2==4) { 
               count4++ ;
            } else { 
               if ($2==5) { 
                  count5++ ;
               } else { 
                  if ($2==6) { 
                     count6++ ;
                  } else { 
                     if ($2==7) {
                        count7++ ;
                     } else { 
                        if ($2==8) { 
                           count8++ ;
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
}
END{
   printf("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",mem,count0,count1,count2,count3,count4,count5,count6,count7,count8)
}
