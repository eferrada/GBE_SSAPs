

# CODE 1: Generating site-specific amino acid preference profiles from structure data.


## Setting up folders:

```
mkdir STRs
mkdir OUTs
mkdir MODELs
mkdir REPAIRs
mkdir MUTANT_FILEs
```

## Reconstructing and repairing initial structures:

 We use the following FoldX commands to reconstruct residue sidechains and to optimize the initial thermodynamic stability of the structure:

```
./foldx --command=ReconstructSideChains --pdb-dir=./STRs --pdb=d2gi9a_.pdb --output-dir=./OUTs/ 
./foldx -c RepairPDB --pdb-dir=./OUTs --pdb=Rebuilt_d2gi9a_.pdb --output-dir=./REPAIRs/
```

With foldx, the executable foldx software. The file "rotabase.txt" must be in the same folder as foldx. 
An academic license can be obtained from : [http://foldxsuite.crg.eu]

## Generating all possible mutations per site

The following files are required:
```
AA.dat  : Alphabetic list of the 20 proteinaceous amino acids.
scp.mut : Awk script to generate an exhaustive list of single mutations.
```

A list of all single mutant sequence are generated using the following script:

```
mkdir MUTANT_FILEs/d2gi9a_
awk -f scp.mut -v WT=MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE AA.dat > MUTANT_FILEs/d2gi9a_/mutant_file.txt
```

The following foldX command generate the single mutant models and measure stability of each, 5 times.

```
./foldx --command=BuildModel --pdb-dir=./REPAIRs/ --pdb=Rebuilt_d2gi9a__Repair.pdb --mutant-file=./MUTANT_FILEs/d2gi9a_/mutant_file.txt --output-dir=./MODELs/d2gi9a_/ --wildtype=MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE --numberOfRuns=5 --out-pdb=false
```

## Filtering FoldX output and generating SSAP profiles:

 The following scripts can be used to transform thermodynamic stability differences of mutations into SSAP profiles based on alternative biophysical models.


### Threshold stability model: sigmoidal function.

**sigmoidal.awk**:

```
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
    print $1"\t"$2"\t"$3"\t"int(1000*Pi/a[$2])/1000; 
  }
}

```

### Optimum stability model: Gaussian.

**gaussian.awk**:

```
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

    print d1[i]"\t"d2[i]"\t"d3[i]"\t"int(1000*Pi/a[T])/1000; 
  }
}
```



### Maximunm stability model : exponential decay.

**exponential.awk**: 

```
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
    print d1[i]"\t"d2[i]"\t"d3[i]"\t"int(1000*Pi/a[T])/1000; 
  }
}

```


 Example of pipeline to generate SSAP profile using the the stability models:

```
tail -1064 MUTANT_FILEs/d2gi9a_/mutant_file.txt > a.1
awk '{if($1~/\_/){print}}' MODELs/d2gi9a_/Average_d2gi9a_.fxout | cut -f1-3 > a.0
paste a.0 a.1 > a.2
awk -f scp.map -v WT=MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE a.2 | sort -n -k2 > a.3
cp a.3 a.4

awk -f exponential.awk a.3 a.4 > SSAPs/d2gi9a_.P1
awk -f gaussian.awk a.3 a.4 > SSAPs/d2gi9a_.P2
awk -f sigmoidal.awk a.3 a.4 > SSAPs/d2gi9a_.P3 
```

 The file **d2gi9a_.P***, summarizes the SSAP profile derived from FoldX and the maximum stability model:

```
#ResID  NRes    aaID     SSAP
M	1	A	0.064	
M	1	C	0.055	
M	1	D	0.036	
M	1	E	0.059
M	1	F	0.018
M	1	G	0.058
M	1	H	0.028
...
E	56	Q	0.054	
E	56	R	0.028
E	56	S	0.029
E	56	T	0.040
E	56	V	0.066
E	56	W	0.011
E	56	Y	0.052

```



