function [  ] = rangeCorr ( dom1, dom2, sampleSize )
% This function compares the preference profiles between the equivalent sites of
% dom1 and dom2, using a sample (sampleSize) of replicates with correlations ranging
% between range1 and range2
%
% matchDir: Depends on the class. Eg: 'MATCHED_D'
% dom1: ID domain 1.
% dom2: ID domain 2.
% range1, range2 : lower and upper range of Confidence Intervals.
% sampleSize: Total samples

% Reading equivalent sites from previously calculated files:
D1 = {'./DATA/MATCHED/', dom1, '_', dom2, '.match.1'};
D2 = {'./DATA/MATCHED/', dom1, '_', dom2, '.match.2'};

% Preparing output files.
OUTPUT = {'./DATA/RESULTs/',dom1, '_', dom2, '_res.dat' };
PVAL = {'./DATA/RESULTs/',dom1, '_', dom2, '_pval.dat' };

a = load (strjoin(D1,''));
b = load (strjoin(D2,''));

[s1 s2] = size (a(:,2));
LENGTH = floor ( s1/20);

count = 0;
c = zeros(LENGTH*20,1);
for i = 1:LENGTH
    for j=1:20
        count = count +1;
        c(count, 1) = i;
    end
end

T= zeros(1,LENGTH);
tc = 0;

ofile = fopen ( strjoin(OUTPUT,''), 'a' );
    fprintf ( ofile, "#[Pearson's r]  [Number of sites at p-value < 0.01]  [Number of sites at p-value < 0.05]  [List of residues with p-value < 0.01]  [List of residues with p-value < 0.05]\n" );
fclose( ofile );

for i=1:100
    for k=1:sampleSize
        [q1 r1] = rndProf ( a(:,2), i/100);
        [q2 r2] = rndProf ( a(:,2), i/100);
    
        [p1 r3] = rndProf ( b(:,2), i/100);
        [p2 r4] = rndProf ( b(:,2), i/100);

        if ((((r1+r2+r3+r4)/4)>=0.50))    	
            [Rc Rr P S] = RMSDarray ( [c(:,1) a(:,2) q1' q2' ],[c(:,1) b(:,2) p1' p2'], LENGTH, 0.05);
  
            T = T+P;
            tc = tc +1;

            [v y] = sort (P);
    
            dlmwrite ( strjoin(OUTPUT,''),[(r1+r2+r3+r4)/4 size(find(P<=0.01),2)  size(find(P<=0.05),2) y(v<0.01) y(v<0.05)],'-append');
            dlmwrite ( strjoin(PVAL,''),[(r1+r2+r3+r4)/4 P],'-append');
        end
    end
  
   if (((r1+r2+r3+r4)/4)<=0.50)
        break;
   end
end

T = T/tc;
ofile = fopen ( strjoin(OUTPUT,''), 'a' );
    fprintf ( ofile, "# P-values\n" );
fclose( ofile );

dlmwrite(strjoin(OUTPUT,''),[T'], '-append');
end



