BEGIN{
  N = 0;
  kT= 0.593;
  MAX = 500;

  a[$2] = 0;
}
{
  if ( FILENAME==ARGV[1] ) {
   N = N +1;

   d1[N] = $1;
   d2[N] = $2;
   d3[N] = $3;
   d4[N] = $4;
  }
  else {
    ddg = $4*$4;

    if ( ddg > MAX ) {
      a[$2]  = a[$2] + exp(-(MAX));
    }
    else {
       a[$2]  = a[$2] + exp(-(ddg));
    }
  }
}
END{
  for(i = 1;i<=N;i=i+1){
    ddg = d4[i]*d4[i];

    if ( ddg > MAX ) {
      Pi= exp(-(MAX));
    }
    else {
      Pi = exp(-(ddg));
    }

    T = d2[i];

    print d1[i]"\t"d2[i]"\t"d3[i]"\t"int(1000*Pi/a[T])/1000; #"\t"int(1000*d4[i])/1000;
  }
}
