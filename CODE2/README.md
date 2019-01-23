
# CODE 2: Calculating the fraction of sites with significantly different amino acid preferences. 

Matlab scripts in this directory randomize amino acid preferences per site according to Eq.6 in the paper.
In addition, the scripts implement the method published by Doud et al. (2015), for the comparison of pairs of site-specific
amino acid preference (SSAP) profiles, based on Jensen-Shannon distance and an exact permutation test. For details, see Method section of the paper and Doud et al. (2015) MBE article.

## Example for running the code:
```
/{PATH_TO_MATLAB}/bin/matlab -nosplash < script.example.exe
```

 The script file above calls the main function **_rangeCorr_** (e.g. __rangeCorr('d2gi9a_', 'd1em7a_', 5)__ )

```
function [  ] = rangeCorr ( dom1, dom2, sampleSize )
% This function compares the preference profiles between the equivalent sites of
% dom1 and dom2, by generating "sampleSize" randomizations of replicates with Pearson correlation coefficients ranging
% between 0.5 and 1.0.
%
% dom1: ID domain 1.
% dom2: ID domain 2.
% sampleSize: Total randomization samples

```

 The script file can be a list of comparisons.
 The script calls a pair of precalculated files located at __./DATA/MATCHED/__ with the following format: 

```
./DATA/MATCHED/', dom1, '_', dom2, '.match.1'
./DATA/MATCHED/', dom1, '_', dom2, '.match.2'
```

 These files list the preference profiles per site, and per amino acid (in alphabetic order).
 The lists match the residue number of all equivalent residues in the structure pair. Equivalent residues are those found to match (< 3.5 Angstroms ) in a 3D alignment.

```
> head DATA/MATCHED/d2gi9a__d1em7a_.match.1
1	0.066
1	0.058
1	0.035
1	0.062
...

```
 There are two output files: **[dom1]_[dom2]_res.dat** and **[dom1]_[dom2]_pval.dat**

### [dom1]_[dom2]_res.dat:
 This output file list the sites found significantly different at varying correlations.  
 It has the following format:

```
[Pearson's r]  [Number of sites at p-value < 0.01]  [Number of sites at p-value < 0.05]  [List of residues with p-value < 0.01]  [List of residues with p-value < 0.05]
0.87794,6,24,38,45,52,54,18,9,38,45,52,54,18,9,17,23,16,30,37,47,27,13,7,3,43,56,48,55,36,6,24,19
0.7823,2,14,23,45,23,45,9,52,54,56,38,37,47,17,2,3,24,30
0.76486,2,23,45,52,45,52,55,54,47,27,9,10,37,43,18,38,50,13,19,7,11,35,17,53,28,29,31
0.75657,8,13,38,9,17,23,52,27,45,8,38,9,17,23,52,27,45,8,54,43,55,24,7
0.72463,5,10,11,45,54,56,23,11,45,54,56,23,47,52,38,37,55
0.73844,4,13,45,52,38,54,45,52,38,54,44,4,13,28,10,47,56,23,27
0.69948,1,4,45,45,17,54,56
0.69466,0,6,47,2,9,24,37,54
0.66736,1,4,45,45,54,36,52
0.65117,2,2,45,45,38
0.64766,0,4,38,43,45,54
0.61432,2,3,52,45,52,45,23
0.62366,1,1,45
0.59828,0,3,45,17,47
0.60914,0,1,45
0.57431,0,2,37,52
0.54776,0,4,9,45,52,53
[P-values]
0.14956
0.15457
0.14264
0.10006
...

```

 At the end of the file, a list of all p-values is printed.

### [dom1]_[dom2]_pval.dat:
  This file lists the Pearson correlation and the p-values per site.

 It has the following format:

```
[Pearson's r], [List of p-values per site ]
0.96886,0.03725,0.045698,0.023095,0.012,0.20182,0.023654,0.018529,0.023654,0.005,0.031111,0.018437,0.03725,0.02225,0.025455,0.035132,0.0175,0.019211,0.025455,0.025,0.025455,0.19317,0.049239,0.012,0.023654,0.1125,0.20759,0.0175,0.037439,0.032568,0.0175,0.059468,0.20182,0.062813,0.01125,0.025455,0.049,0.01,0.01,0.023654,0.025,0.20182,0.12284,0.01,0.039524,0.005,0.026176,0.0175,0.023704,0.026286,0.023654,0.048295,0.005,0.088673,0.005,0.019167,0.014091
0.96734,0.034459,0.039419,0.016818,0.01125,0.26348,0.01125,0.01125,0.017609,0.005,0.12717,0.014306,0.044091,0.039419,0.02875,0.034459,0.016818,0.016818,0.0191,0.034459,0.022963,0.15713,0.034934,0.0078571,0.016818,0.030909,0.12717,0.014306,0.076957,0.1164,0.01125,0.037436,0.11082,0.11882,0.0078571,0.032353,0.030909,0.01125,0.01125,0.018333,0.029828,0.21027,0.11082,0.006,0.039419,0.005,0.030909,0.01125,0.11082,0.030909,0.022963,0.039419,0.005,0.051556,0.005,0.01125,0.01125
...
```

 The input parameter and outputs are explained as part of the documentation of each script.


