BEGIN{
  N = 0;
  s = 0.0;
  kT= 0.593;
  Ne = 1;

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
    ddg = $4;

    a[$2]  = a[$2] + exp(-(ddg));
  }
}
END{
  for(i = 1;i<=N;i=i+1){
    ddg = d4[i];
    Pi = exp(-(ddg));

    T = d2[i];
    print d1[i]"\t"d2[i]"\t"d3[i]"\t"int(1000*Pi/a[T])/1000; #"\t"int(1000*d4[i])/1000;
  }
}
