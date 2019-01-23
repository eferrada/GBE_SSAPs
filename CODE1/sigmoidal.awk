BEGIN{
  kT= 0.593;
  for ( i=1;i<=L;i=i+1 ){
     a[$2] = 0;
  }
}
{
  if ( FILENAME==ARGV[1] ) {
   ddg = $4;
   Pi  = exp(-ddg/kT)/(1+exp(-ddg/kT));

   a[$2] = a[$2] + Pi;
  }
  else {
    ddg = $4;
    Pi  = exp(-ddg/kT)/(1+exp(-ddg/kT));
    print $1"\t"$2"\t"$3"\t"int(1000*Pi/a[$2])/1000; #"\t"int(1000*$4)/1000;
  }
}
