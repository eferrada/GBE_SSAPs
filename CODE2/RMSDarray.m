function [Rc Rr P S] = RMSDarray ( obs, pre, LENGTH, ALPHA )
% obs    : file listing Ne empirical preferences.
% pre    : file lisitng Nt theoretical preferences.
% LENGTH : sequence length of the protein.
% 
% Rc : RMSD corrected for the non-randomized data.
% Rr : RMSD from randomized replicate groups.
% P  : p-values from randomization test and FDR at ALPHA.
% S  : significantly different sites at alpha = ALPHA

%% Calculating RMSDs for non-randomized data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rw Rb Rc] = RMSDcorrected (obs, pre, LENGTH);

%% Opening files, defining variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe = obs; 
ft = pre; 

[re ce] = size(fe);
[rt ct] = size(ft);

% Effective number of replicates.
Na = ce-1;
Nb = ct-1;

% Creating matrix with all data.
M = [fe(:,2:ce) ft(:,2:ct)];
[rm cm] = size(M);

% Calculating possible combinations of replicates per dataset.
Ca = nchoosek([1:cm],Na);
Cb = nchoosek([1:cm],Nb);

%% Exact randomizations:
%%%%%%%%%%%%%%%%%%%%%%%%
E = zeros( re, ce, 1);
T = zeros( rt, ct, 1);
[ra ca] = size(Ca);
[rb cb] = size(Cb);
Rr = zeros ( LENGTH, ra*rb, 1);

count = 0;
for i=1:ra
    for j=1:rb
        count = count + 1;
        [Rw Rb Rr(:,count)] = RMSDcorrected ([fe(:,1) M(:,Ca(i,:))], [ft(:,1) M(:,Cb(j,:))], LENGTH);
    end
end

%% Calculating p-values
%%%%%%%%%%%%%%%%%%%%%%%
total = LENGTH * ra * rb;
for j=1:LENGTH
    P(j) = length(find(Rr(:,:)>=Rc(j)))/total;
end

FDR = mafdr(P,'BHFDR', 'true');

P = FDR;
S = find(P<=ALPHA);

end
