function [Q r] = rndProf( p, f )
%% Receives:
% p = original profile (vector multiple of 20).
% f = noise in percentage
%% Returns:
% Q = new profile with white noise in the same dimensions than p.
% r = correlation between p and q.

%% Transform input profile into matrix
c = vec2mat ( p,20);
[Length s1] = size (c);

%% Initialize new profile:
q = zeros ( 56,20);

% For each position:
for i=1:Length
    % gets a multinomial random sample of size 100
    q(i,:) = mnrnd(100,c(i,:)./sum(c(i,:)));

    % adds noise as a fraction f of 100.
    q(i,:) = q(i,:) + abs(normrnd(q(i,:),f*100,1,20));

    % renormalized output profile per site.
    q(i,:) = q(i,:)./sum(q(i,:));
end

%% Obtain correlation between p and q.
P = reshape(c',[1,20*Length]);
Q = reshape(q',[1,20*Length]);
[R h] = corrcoef( reshape(c',[1,20*Length]), reshape (q',[1,20*Length]));
r = R(1,2); 
%scatter( reshape(c,[1,20*Length]), reshape (r,[1,20*Length]))

end

