function [Rw Rb Rc] = RMSDcorrected (obs, pre, LENGTH)
% obs : file listing Ne empirical preferences.
% pre : file lisitng Nt theoretical preferences.
% LENGTH     : sequence length of the protein.
% 
% Rw : RMSD within
% Rb : RMSD between
% Rc : RMSD corrected

fe = obs; 
ft = pre; 

[re ce] = size(obs);
[rt ct] = size(pre);

% Ne = Number of replicates in A. 
% Nt = Number of replicates in B.
Na = ce-1;
Nb = ct-1;

%% Calculating RMSD-between (Rb)
Nbetween = Na*Nb;
Rb = zeros(LENGTH,1);
for i = 1:Na
    for j = 1:Nb
        for k = 1:LENGTH
            Rb(k) = Rb(k) + JSDiv(fe(fe(:,1)==k,i+1)',ft(ft(:,1)==k,j+1)'); 
        end
    end
end

for k = 1:LENGTH
    Rb(k) = sqrt(Rb(k)./Nbetween);
end

%% Calculating RMSD-within (Rw)
if ( Na > 1 )
Naa = nchoosek(Na,2);
else
    Naa = 1;
end

if ( Nb > 1 )
    Nbb = nchoosek(Nb,2);
else
    Nbb = 1;
end

Rw = zeros(LENGTH,1);
Rwa = zeros(LENGTH,1);
Rwb= zeros(LENGTH,1);

if ( Na > 1 )
    for i=1:Na
        for j=i+1:Na
            for k = 1:LENGTH
                Rwa(k) = Rwa(k) + JSDiv(fe(fe(:,1)==k,i+1)',fe(fe(:,1)==k,j+1)'); 
            end
        end
    end
    
    for k = 1:LENGTH
        Rwa(k) = sqrt(Rwa(k)./Naa);
    end
end

if ( Nb > 1 )
    for i=1:Nb
        for j=i+1:Nb
            for k = 1:LENGTH
                Rwb(k) = Rwb(k) + JSDiv(ft(ft(:,1)==k,i+1)',ft(ft(:,1)==k,j+1)');  
            end
        end
    end
    
    for k = 1:LENGTH
        Rwb(k) = sqrt(Rwb(k)./Nbb);
    end
end

Rw = 0.5*Rwa + 0.5*Rwb;

%% Calculating RMSD-corrected
Rc = Rb-Rw;

end

