BEGIN{
 N=0; LastKey = "-1-1-1-"
}

{
 N++;
 if(N==1) {for(i=1; i<= NF; i++) printf "%s ",$i;}
 else {
   if($1!=LastKey) { printf "\n"; for(i=1; i<= NF; i++) printf "%s ",$i;}
   else { for(i=2; i<= NF; i++) printf "%s ",$i; }
 }
 LastKey = $1; 
}
